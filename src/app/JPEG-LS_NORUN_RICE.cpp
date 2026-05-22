#include <cassert>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <vector>
#include <array>
#include <string>
#include <chrono>
#include <algorithm>
#include <sstream>
#include <filesystem>

#include "../lib/arithCoding.h"
#include "../lib/pgm.h"

namespace fs = std::filesystem;

//======================================================
//
//   C O M M A N D    L I N E   P A R A M E T E R S
//
//======================================================
// Je 6 Schwellen fuer GT1 und GTN => je 7 Klassen
struct cmdPars
{
    bool        decode;
    std::string inname;
    std::string outname;
    std::string outdir;
};


bool readCmdLine(int argc, char** argv, cmdPars& pars)
{
    std::stringstream err;
    int arg = 0;
    bool ok = true;

    ok = ok && (argc > ++arg);
    if (ok) {
        if      (std::string(argv[arg]) == "-e") pars.decode = false;
        else if (std::string(argv[arg]) == "-d") pars.decode = true;
        else {
            err << "ERROR: First parameter must be '-e' or '-d'." << std::endl;
            ok = false;
        }
    }

    ok = ok && (argc > ++arg);
    if (ok) {
        pars.inname = argv[arg];
        std::ifstream intest(pars.inname);
        if (!intest.good()) {
            err << "ERROR: Cannot open input file \"" << pars.inname << "\"." << std::endl;
            ok = false;
        }
    }

    ok = ok && (argc > ++arg);
    if (ok) pars.outname = argv[arg];

    if (argc > ++arg) pars.outdir = argv[arg];
    else              pars.outdir = "output";


    return ok;
}


//======================================================
//
//   G L O B A L E   T E S T P A R A M E T E R
//
//======================================================

//======================================================
//
//   E N T R O P Y   C O D I N G
//
//======================================================

class EntropyCoderBase
{
protected:
    static const unsigned N = 8;
    static const unsigned NUM_CTX = 256 + 256 + 256 + 256;
    unsigned GROUP_SIZE;

    std::array<std::array<uint8_t, N>, NUM_CTX> m_pmfAbsCtx;
    std::array<uint8_t, NUM_CTX> m_pmfSignCtx;
    std::array<uint8_t, NUM_CTX> m_pmfSigCtx;
    std::array<std::array<uint8_t, N>, NUM_CTX> m_pmfExpCtx;


   

public:
    EntropyCoderBase(unsigned groupSize)
        : GROUP_SIZE(groupSize)
    {
        for (auto& ctx : m_pmfAbsCtx)
            for (auto& p : ctx) p = 0;
        for (auto& p : m_pmfSignCtx) p = 0;
        for (auto& p : m_pmfSigCtx)  p = 0;
        for (auto& ctx : m_pmfExpCtx)
            for (auto& p : ctx) p = 0;
    }
};


class EntropyEncoder : protected EntropyCoderBase
{
public:
    EntropyEncoder(OBitstream& bs, unsigned groupSize)
        : EntropyCoderBase(groupSize), aenc(bs)
    { aenc.start(); }

    void encodeSample(PGMImage::Sample s, unsigned ctxIdx, unsigned sigIdx, int a_err, int b_err, int c_err, int d_err, int q1, int q2, int q3, int k);
    void finish() { aenc.finish(); }
    void encRice(unsigned n, unsigned r);

private:
    ArithmeticEncoder aenc;
};


class EntropyDecoder : protected EntropyCoderBase
{
public:
    EntropyDecoder(IBitstream& bs, unsigned groupSize)
        : EntropyCoderBase(groupSize), adec(bs)
    { adec.start(); }

    PGMImage::Sample decodeSample(unsigned ctxIdx, unsigned sigIdx, int a_err, int b_err, int c_err, int d_err, int q1, int q2, int q3, int k);
    unsigned decRice(unsigned r);


private:
    ArithmeticDecoder adec;
};

void EntropyEncoder::encRice(unsigned n, unsigned r)
{
    int pre = n >> r;
    int suf = n - (pre << r);
    
    // Unary encoding: encode pre ones followed by a zero
    for (int i = 0; i < pre; i++)
        aenc.encBit(1);
    aenc.encBit(0);
    
    // Fixed encoding: encode suf with r bits
    if (r > 0)
        aenc.encBits(suf, r);
}

unsigned EntropyDecoder::decRice(unsigned r)
{
    int pre = 0;
    // Unary decoding: count ones until zero
    while (adec.decBit()) {
        pre++;
    }
    
    unsigned suf = 0;
    // Fixed decoding: read r bits
    if (r > 0)
        suf = adec.decBits(r);
    
    unsigned n = (pre << r) + suf;
    return n;
}

void EntropyEncoder::encodeSample(PGMImage::Sample s, unsigned /*ctxIdx*/, unsigned /*sigIdx*/, int /*a_err*/, int /*b_err*/, int /*c_err*/, int /*d_err*/, int /*q1*/, int /*q2*/, int /*q3*/, int k)
{
    unsigned absValue = unsigned(s < 0 ? -s : s);
    encRice(absValue, k);

    if (absValue)
        aenc.encBit(s < 0);
}

PGMImage::Sample EntropyDecoder::decodeSample(unsigned /*ctxIdx*/, unsigned /*sigIdx*/,
                                               int /*a_err*/, int /*b_err*/, int /*c_err*/, int /*d_err*/, int /*q1*/, int /*q2*/, int /*q3*/, int k)
{
    unsigned absValue = decRice(k);
    
    if (absValue) {
        bool sign = adec.decBit();
        return sign ? -PGMImage::Sample(absValue) : PGMImage::Sample(absValue);
    }
    return 0;
}
//======================================================
//
//   LOCO-I PREDICTOR MIT DIAGNOSE
//
//======================================================
class Prediction
{
public:
    Prediction(int width, int height, std::vector<PGMImage::Sample>& img,
               bool isDecoding = false, const std::string& outputDir = "output")
        : m_width(width), m_height(height), m_data(img), m_org(m_data),
          m_isDecoding(isDecoding), m_outputDir(outputDir)
    {
        m_ctxs.resize(367);
        for (auto& ctx : m_ctxs) {
            ctx.N = 1;
            ctx.A = 4;
            ctx.B = 0;
            ctx.C = 0;
        }

        buildContextMap();
    }

void subtractPrediction(EntropyEncoder& eenc, int /* gsize*/)
{
    const PGMImage::Sample* orgData = m_org.data();
    //PGMImage::Sample*       data    = m_data.data();

    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {

            // Lambda-Funktion zum sicheren Zugriff auf Nachbarpixel (Ränder = 0)
            auto get = [&](int xx, int yy) -> int {
                if (xx < 0 || yy < 0 || xx >= m_width || yy >= m_height) return 0;
                return int(orgData[yy * m_width + xx]);
            };

            // --- Nachbarn ---
            int a = get(x-1,y), b = get(x,y-1), c = get(x-1,y-1), d = get(x+1,y-1), Ix = get(x,y);

            // Step 2. For the current sample, compute the local gradients according to Code segment A.1. 
            // A.3.1 Local gradient computation
            int g1 = d-b, g2 = b-c, g3 = c-a;

            // Step 3. Select the coding mode
            // A.3.2 Mode selection  - Code segment A.3 – Mode selection procedure for lossless coding
            if (g1 == 0 && g2 == 0 && g3 == 0 && 0) {
                // --- RUN ---
                // später

            } else {
                // --- Normal ---


                // Step 4. Quantize the local gradients according to the steps detailed in Code segment A.4. 
                // --- reguläre Prädiktion --- A.3.3 Local gradient quantization
                int q1 = quantize(g1), q2 = quantize(g2), q3 = quantize(g3);

                // Step 5. Check  and  change  if  necessary  the  signs  of  the  components  of  the  vector  representing  the  context,  
                // modifying accordingly the variable SIGN (see A.3.4).
                // --- SIGN --- A.3.4 Quantized gradient merging
                // A.3.4 Quantized gradient merging
                int SIGN = computeSign(q1, q2, q3);

                // kanonische Richtung herstellen
                q1 *= SIGN;
                q2 *= SIGN;
                q3 *= SIGN;

                // Quantized gradient merging

                int Q = contextID(q1, q2, q3);

                // Step 6. Compute Px according to Code segment A.5. 
                // A.4.1 Edge-detecting predictor
                // --- Prädiktor ---
                int Predictor_x;

                codeSegmentA5(a,b,c, Predictor_x);

                // Step 7. Correct Px using C[Q] and the variable SIGN, and clamp the corrected value to the interval 
                // [0..MAXVAL] according to the procedure in Code segment A.6. 
                // A.4.2 Prediction correction
                codeSegmentA6(SIGN, Predictor_x, Q);

                // Step 8. Compute  the  prediction  error  and,  if  necessary,  invert  its  sign  according  to  the  procedure  in  Code  
                // segment A.7.
                // A.4.3 Computation of Prediction error
                int Errval;
                codeSegmentA7(Errval, Ix, Predictor_x, SIGN);


                // A.4.4 Error quantization for near-lossless coding, and reconstructed value computation
                // Step 9 or  lossless  coding,  update  the  reconstructed  value  by  setting  Rx  equal  
                // to Ix. 

                //int Rx = Ix;  // bei Verlustlosigkeit ist Quantisierung identischt

                // codeSegmentA9(); // später, bei Verlustbehaftung

                // Step 10. Reduce the error to the relevant range according to Code segment A.9.
                // A.4.5 Modulo reduction of the prediction error
                //codeSegmentA9(Errval);

                // Step 11. Compute the context-dependent Golomb variable k according to the procedure in Code segment A.10. 
                // A.5.1 Estimation of k - Golomb coding variable computation
                int k = 0;
                codeSegmentA10(k, Q);

                // Step 12. Perform the error mapping according to the procedure in Code segment A.11.
                //int MErrval;
                //codeSegmentA11(MErrval, Errval, k, Q);

                // Debug: Erste 5 Pixel ausgeben

                // Step 13. Encode the mapped error value MErrval using the limited length Golomb code function LG(k, LIMIT), as 
                // specified in A.5.3
                // A.5.3 Mapped error encoding
                eenc.encodeSample(static_cast<PGMImage::Sample>(Errval),
                              /*model*/0, /*sig*/0,
                              0,0,0,0,
                              q1,q2,q3,k);

                // Step 14. Update the variables according to Code segment A.12. 
                codeSegmentA12(Q, Errval);
                
                // Step 15. Update the prediction correction value C[Q] according to the procedure in Code segment A.13.
                codeSegmentA13(Q);
                
                // Step 16. Next sample
            }
        }
    }
    eenc.finish();
}

void addPrediction(EntropyDecoder& edec)
{
    std::vector<PGMImage::Sample> tempData(m_data.size());

    for (int y = 0; y < m_height; ++y) {
        for (int x = 0; x < m_width; ++x) {

            // Lambda-Funktion zum sicheren Zugriff auf Nachbarpixel (Ränder = 0)
            auto get = [&](int xx, int yy) -> int {
                if (xx < 0 || yy < 0 || xx >= m_width || yy >= m_height) return 0;
                return int(tempData[yy * m_width + xx]);
            };

            // --- Nachbarn ---
            int a = get(x-1,y), b = get(x,y-1), c = get(x-1,y-1), d = get(x+1,y-1);

            int g1 = d-b, g2 = b-c, g3 = c-a;

            if (g1 == 0 && g2 == 0 && g3 == 0 && 0) {
                // --- RUN ---
                // später
            } else {
            // --- Normal ---


            // Step 4. Quantize the local gradients according to the steps detailed in Code segment A.4. 
            // --- reguläre Prädiktion --- A.3.3 Local gradient quantization
            int q1 = quantize(g1), q2 = quantize(g2), q3 = quantize(g3);

            // Step 5. Check  and  change  if  necessary  the  signs  of  the  components  of  the  vector  representing  the  context,  
            // modifying accordingly the variable SIGN (see A.3.4).
            // --- SIGN --- A.3.4 Quantized gradient merging
            // A.3.4 Quantized gradient merging
            int SIGN = computeSign(q1, q2, q3);

            // kanonische Richtung herstellen
            q1 *= SIGN;
            q2 *= SIGN;
            q3 *= SIGN;

            // Quantized gradient merging

            int Q = contextID(q1, q2, q3);

            // Step 6. Compute Px according to Code segment A.5. 
            // A.4.1 Edge-detecting predictor
            // --- Prädiktor ---
            int Predictor_x;

            codeSegmentA5(a,b,c, Predictor_x);

            // Step 7. Correct Px using C[Q] and the variable SIGN, and clamp the corrected value to the interval 
            // [0..MAXVAL] according to the procedure in Code segment A.6. 
            // A.4.2 Prediction correction
            codeSegmentA6(SIGN, Predictor_x, Q);


            // Step 8. Compute the context-dependent Golomb variable k according to the procedure in Code segment A.10. 
            // A.5.1 Estimation of k - Golomb coding variable computation
            int k = 0;
            codeSegmentA10(k, Q);


            // A.5.3 Mapped error decoding
            // Step 9. Decode the mapped error value MErrval: 
            int Errval = edec.decodeSample(/*ctxIdx*/0, /*sigIdx*/0,
                              0,0,0,0,
                              q1,q2,q3,k);


            // Step 10. Perform  the  inverse  of  the  error  mapping  indicated  in  Code  segment  A.11,  where  now  MErrval  is  given  
            // and Errval is computed.
            //int Errval;
            //codeSegmentA11D(Errval, MErrval, k, Q);
            
            // Step 11. Update the variables according to Code segment A.12. 
            codeSegmentA12(Q, Errval);

            // Step 12. Only for nearlossless

            // Step 13. Invert sign of Errval if the variable SIGN is negative
            if (SIGN < 0) Errval = -Errval;

            // Step 14. 
            int Rx = Predictor_x + Errval; 
            Rx %= 256; // Modulo-Operation, um sicherzustellen, dass Rx im Bereich [0..255] liegt

            if (Rx < 0) Rx += 256; // Falls Rx negativ ist, korrigieren
            else if (Rx > 255) Rx -= 256; // Falls Rx über 255 liegt, korrigieren

            if (Rx < 0) Rx = 0;
            else if (Rx > 255) Rx = 255;

            // Debug: Erste 5 Pixel ausgeben

            // Step 15. ap Rx  using  the  inverse  point  transform  Pt  specified  by  the  parameter  A
            int Ix = Rx; // bei Verlustlosigkeit ist die Punkttransformation identisch

            tempData[y * m_width + x] =
                static_cast<PGMImage::Sample>(Ix);

            // --- Context Update ---

            biasupdateLocoI(Q);
            }
        }
    }

    m_data = tempData;
}


private:
    struct Context { int A=0, B=0, N=1, C=0; };


    inline void updateContextLocoI(int ctxIdx, int error)
    {
        Context& ctx = m_ctxs[ctxIdx];
        ctx.B += error;
        ctx.A += std::abs(error); 
        ctx.N += 1;

        if (ctx.N == 64) {
            ctx.A >>= 1;
            if (ctx.B >= 0) {
                ctx.B >>= 1;
            } else {
                ctx.B = -((1-ctx.B) >> 1);
            }
            ctx.N >>= 1;
        }
        ctx.N += 1;
    }

    inline void biasupdateLocoI(int ctxIdx)
    {
        // MAX_C maximum allowed value of C[0..364], equal to 127 
        // MIN_C minimum allowed value of C[0..364], equal to –128 
        Context& ctx = m_ctxs[ctxIdx];
        if (ctx.B <= -ctx.N) {
            ctx.B += ctx.N;
            if (ctx.C > -128) ctx.C -= 1;
            if (ctx.B <= -ctx.N) ctx.B = -ctx.N + 1; // prevent overflow
        } else if (ctx.B > 0) {
            ctx.B -= ctx.N;
            if (ctx.C < 127) ctx.C += 1;
            if (ctx.B > 0) ctx.B = 0; // prevent overflow
        }
    }


    inline void codeSegmentA5(int a, int b, int c, int& Px) {
        // Code segment A.5
        if (c >= std::max(a,b)) {
             Px = std::min(a,b);
        } else if (c <= std::min(a,b)) {
            Px = std::max(a,b);
        } else {
             Px = a + b - c;
        }

    }

    inline void codeSegmentA6(int SIGN, int& Px, int Q) {
        // Code segment A.6
        if (SIGN == 1)
            Px += m_ctxs[Q].C; // Q ist hier 0, da Kontext-Update erst später erfolgt
        else
            Px -= m_ctxs[Q].C; // Q ist hier 0, da Kontext-Update erst später erfolgt
        if (Px < 0) Px = 0;
        else if (Px > 255) Px = 255;
    }

    inline void codeSegmentA7(int& Errval, int Ix, int Px, int SIGN) {
        // Code segment A.7
        Errval = Ix - Px;
        if (SIGN == -1)
            Errval = -Errval;
    }

    inline void codeSegmentA8(int& Rx, int Ix) {
        // Code segment A.8
        Rx = Ix; // bei Verlustlosigkeit ist Quantisierung identisch
    }
    
    inline void codeSegmentA9(int& Errval) {
        // Code segment A.9
        if (Errval < 0)
            Errval += 256;     // Errval is now in [0..255]
        if (Errval >= (256+1) / 2)
            Errval -= 256;     // Errval is now in [-128..127]
    }

    inline void codeSegmentA10(int& k, int Q) {
        // Code segment A.10
        k = 0;
        while ((m_ctxs[Q].N << k) < m_ctxs[Q].A) k++;
    }

    inline void codeSegmentA11(int& MErrval, int Errval, int k, int Q) {
        // Code segment A.11
        if (k == 0 && ((2 * m_ctxs[Q].B) <= -m_ctxs[Q].N)) {
            if (Errval >= 0) {
                MErrval = 2 * Errval + 1;
            } else {
                MErrval = -2 * (Errval + 1);
            }
        } else {
            if (Errval >= 0) {
                MErrval = (2 * Errval);
            } else {
                MErrval = (-2*(Errval)) - 1;
            }
        }
    }

    inline void codeSegmentA11D(int& Errval, int MErrval, int k, int Q) {
        // Code segment A.11 inverse
        if (k == 0 && ((2 * m_ctxs[Q].B) <= -m_ctxs[Q].N)) {
            if (MErrval % 2 == 1) {
                Errval = MErrval / 2;
            } else {
                Errval = -(MErrval / 2) - 1;
            }
        } else {
            if (MErrval % 2 == 0) {
                Errval = MErrval / 2;
            } else {
                Errval = -(MErrval + 1) / 2;
            }
        }
    }

    inline void codeSegmentA12(int Q, int Errval) {
        // Code segment A.12
        updateContextLocoI(Q, Errval);
    }

    inline void codeSegmentA13(int Q) {
        // Code segment A.13
        biasupdateLocoI(Q); // Kontext-Index ist hier 0, da Kontext-Update erst später erfolgt
    }

    int contextMap[9][9][9];

    void buildContextMap()
    {
        int id = 0;
        for (int q1 = -4; q1 <= 4; ++q1)
        {
            for (int q2 = -4; q2 <= 4; ++q2)
            {
                for (int q3 = -4; q3 <= 4; ++q3)
                {
                    bool keep =
                        (q1 > 0) ||
                        (q1 == 0 && q2 > 0) ||
                        (q1 == 0 && q2 == 0 && q3 >= 0);
                    if (keep)
                    {
                        contextMap[q1 + 4]
                                  [q2 + 4]
                                  [q3 + 4] = id;
                        ++id;
                    }
                }
            }
        }
    }

    inline int computeSign(int q1, int q2, int q3) const
    {
        if ((q1 < 0) ||
            (q1 == 0 && q2 < 0) ||
            (q1 == 0 && q2 == 0 && q3 < 0))
        {
            return -1;
        }

        return +1;
    }

    inline int contextID(int q1, int q2, int q3) const
    {
        return contextMap[q1 + 4]
                         [q2 + 4]
                         [q3 + 4];
    }

    inline int quantize(int grad) const
    {
        constexpr int T1 = 3;
        constexpr int T2 = 7;
        constexpr int T3 = 21;

        if      (grad <= -T3) return -4;
        else if (grad <= -T2) return -3;
        else if (grad <= -T1) return -2;
        else if (grad <   0 ) return -1;
        else if (grad ==  0 ) return  0;
        else if (grad <   T1) return  1;
        else if (grad <   T2) return  2;
        else if (grad <   T3) return  3;
        else                  return  4;
    }

    int contextModel(int a, int b, int c, int d) const
    {
        return std::abs(a) + std::abs(b) + std::abs(c) + std::abs(d);
    }

    
    bool          m_isDecoding;
    std::string   m_outputDir;
    int           bias_used;
    int           total_pixels;
    const int     m_alphabetSize = 256;
    int           m_width;
    int           m_height;
    int           RUNindex = 0;
    int           Nn[2] = {0, 0};
    const int     range = 256;
    const int     qbpp = 8;
    const int     bpp = 8;
    const int     limit = 32;
    const int J[32] = {
    0, 0, 0, 0,
    1, 1, 1, 1,
    2, 2, 2, 2,
    3, 3, 3, 3,
    4, 4,
    5, 5,
    6, 6,
    7, 7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15
    };

    std::vector<PGMImage::Sample>& m_data;
    std::vector<PGMImage::Sample>  m_org;
    std::vector<Context>           m_ctxs;
};


//======================================================
//
//   M A I N   E N C O D I N G + D E C O D I N G
//
//======================================================
void encode(const std::string& inname, const std::string& outname,
            unsigned groupSize, const std::string& outputDir = "output")
{
    PGMImage img;
    img.read(inname);

    std::filesystem::create_directories(outputDir);
    const std::filesystem::path outpath =
        std::filesystem::path(outputDir) / outname;

    std::ofstream stream(outpath, std::ios::binary);
    OBitstream bs(stream);

    bs.addFixed<unsigned>(img.getWidth(),  16);
    bs.addFixed<unsigned>(img.getHeight(), 16);

    EntropyEncoder eenc(bs, groupSize);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData(), false, outputDir);
    pred.subtractPrediction(eenc, groupSize);

    bs.byteAlign();
    stream.close();
}

void decode(const std::string& inname, const std::string& outname,
            unsigned groupSize, const std::string& outputDir = "output")
{
    std::ifstream stream(inname, std::ios::binary);
    IBitstream bs(stream);

    int width  = bs.getFixed<unsigned>(16);
    int height = bs.getFixed<unsigned>(16);

    PGMImage img(width, height);

    EntropyDecoder edec(bs, groupSize);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData(), true, outputDir);
    pred.addPrediction(edec);

    std::filesystem::create_directories(outputDir);
    const std::filesystem::path outpath =
        std::filesystem::path(outputDir) / outname;

    img.write(outpath.string());
}

int main(int argc, char** argv)
{
    cmdPars pars;
    if (!readCmdLine(argc, argv, pars)) return 1;

    fs::create_directories(pars.outdir);

    std::string fullOutputPath = pars.outdir + "/" + pars.outname;

    unsigned groupSize = 1;

    try {
        if (pars.decode) {
            decode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        } else {
            encode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Fehler: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}