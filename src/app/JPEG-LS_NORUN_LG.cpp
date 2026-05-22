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
    // GT1-Schwellen
    int gt1_t1, gt1_t2, gt1_t3, gt1_t4, gt1_t5, gt1_t6;
    // GTN-Schwellen
    int gtn_t1, gtn_t2, gtn_t3, gtn_t4, gtn_t5, gtn_t6;
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

    auto readOptInt = [&](int& dest, int defaultVal, const char* name) {
        if (argc > ++arg) {
            try {
                dest = std::stoi(argv[arg]);
                if (dest < 1 || dest > 99999) {
                    err << "ERROR: " << name << " must be between 1 and 99999." << std::endl;
                    ok = false;
                }
            } catch (...) {
                err << "ERROR: " << name << " must be an integer." << std::endl;
                ok = false;
            }
        } else {
            dest = defaultVal;
        }
    };

    // GT1-Schwellen (Default: binaer 1,2,4,8,16,32)
    readOptInt(pars.gt1_t1,  1, "GT1_Threshold1");
    readOptInt(pars.gt1_t2,  2, "GT1_Threshold2");
    readOptInt(pars.gt1_t3,  4, "GT1_Threshold3");
    readOptInt(pars.gt1_t4,  8, "GT1_Threshold4");
    readOptInt(pars.gt1_t5, 16, "GT1_Threshold5");
    readOptInt(pars.gt1_t6, 32, "GT1_Threshold6");

    // GTN-Schwellen (Default: binaer 1,2,4,8,16,32)
    readOptInt(pars.gtn_t1,  1, "GTN_Threshold1");
    readOptInt(pars.gtn_t2,  2, "GTN_Threshold2");
    readOptInt(pars.gtn_t3,  4, "GTN_Threshold3");
    readOptInt(pars.gtn_t4,  8, "GTN_Threshold4");
    readOptInt(pars.gtn_t5, 16, "GTN_Threshold5");
    readOptInt(pars.gtn_t6, 32, "GTN_Threshold6");

    // Validierung GT1: streng aufsteigend
    if (ok) {
        int vals[6] = {
            pars.gt1_t1, pars.gt1_t2, pars.gt1_t3,
            pars.gt1_t4, pars.gt1_t5, pars.gt1_t6
        };
        for (int i = 1; i < 6; i++) {
            if (vals[i] <= vals[i-1]) {
                err << "ERROR: GT1-Schwellen muessen streng aufsteigend sein: gt1_t"
                    << i << "=" << vals[i-1] << " >= gt1_t" << (i+1) << "=" << vals[i] << std::endl;
                ok = false;
                break;
            }
        }
    }

    // Validierung GTN: streng aufsteigend
    if (ok) {
        int vals[6] = {
            pars.gtn_t1, pars.gtn_t2, pars.gtn_t3,
            pars.gtn_t4, pars.gtn_t5, pars.gtn_t6
        };
        for (int i = 1; i < 6; i++) {
            if (vals[i] <= vals[i-1]) {
                err << "ERROR: GTN-Schwellen muessen streng aufsteigend sein: gtn_t"
                    << i << "=" << vals[i-1] << " >= gtn_t" << (i+1) << "=" << vals[i] << std::endl;
                ok = false;
                break;
            }
        }
    }

    if (!ok) {
        std::string pname = argv[0];
        std::cerr << err.str() << std::endl;
        std::cerr << "Usage: " << pname
                  << " -e|-d inFile outFile [outputDir]"
                  << " [gt1_t1..gt1_t6] [gtn_t1..gtn_t6]" << std::endl;
        std::cerr << "  GT1 t1..t6: 1-99999, streng aufsteigend, Default: 1 2 4 8 16 32" << std::endl;
        std::cerr << "  GTN t1..t6: 1-99999, streng aufsteigend, Default: 1 2 4 8 16 32" << std::endl;
        return false;
    }
    return true;
}


//======================================================
//
//   G L O B A L E   T E S T P A R A M E T E R
//
//======================================================
// GT1-Schwellen
int g_gt1_t1 =  1, g_gt1_t2 =  2, g_gt1_t3 =  4;
int g_gt1_t4 =  8, g_gt1_t5 = 16, g_gt1_t6 = 32;

// GTN-Schwellen
int g_gtn_t1 =  1, g_gtn_t2 =  2, g_gtn_t3 =  4;
int g_gtn_t4 =  8, g_gtn_t5 = 16, g_gtn_t6 = 32;

std::ofstream g_bpp_results_file;


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
    std::array<uint8_t, NUM_CTX> m_pmfGt1Ctx;
    std::array<uint8_t, NUM_CTX> m_pmfGt2Ctx;
    std::array<uint8_t, NUM_CTX> m_pmfGt3Ctx;
    std::array<uint8_t, NUM_CTX> m_pmfGt4Ctx;
    std::array<std::array<uint8_t, N>, NUM_CTX> m_pmfExpCtx;

    inline unsigned mappedCtx(unsigned ctxIdx) const {
        ctxIdx = 1;
        return ctxIdx;
    }

    static inline unsigned get_gt1_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);

        if      (sum <= 1) return 0;
        else if (sum <= 5) return 1;
        else if (sum <= 10) return 2;
        else if (sum <= 17) return 3;
        else if (sum <= 29) return 4;
        else if (sum <= 61) return 5;
        else                return 6;
    }

    static inline unsigned get_gtN_index(int a_err, int b_err, int c_err, int d_err, unsigned binIdx)
    {
        int sum = std::abs(a_err) + std::abs(b_err) + std::abs(c_err) + std::abs(d_err);
        unsigned bin = std::min<unsigned>(N-1, binIdx);
        if      (sum <= 1)  return 0  + bin;
        else if (sum <= 2)  return 1+N  + bin;  // Abstand muss >= N sein!
        else if (sum <= 4)  return 1+2*N  + bin;
        else if (sum <= 8)  return 1+3*N  + bin;
        else if (sum <= 16) return 1+4*N  + bin;  // Abstand muss >= N sein!
        else if (sum <= 32) return 1+5*N  + bin;
        else                return 1+6*N  + bin;
    }


    inline unsigned get_sign_index(int q1, int q2, int q3,
                                                   int a_err, int b_err)
    {
        unsigned ctx = 0;

        // Bit 0: Gradient-Vorzeichen
        int gradient_sign = (std::abs(q1) > 0) ? (q1 < 0 ? -1 : 1) :
                            (std::abs(q2) > 0) ? (q2 < 0 ? -1 : 1) : 
                                                 (q3 < 0 ? -1 : 1);
        if (gradient_sign < 0) ctx |= 0b001;

        // Bit 1: Fehler-Magnitude (ob a_err oder b_err größer)
        if (std::abs(a_err) > std::abs(b_err)) ctx |= 0b010;

        // Bit 2: Fehler-Vorzeichen
        int residual_sign = (a_err < 0 || b_err < 0) ? -1 : 1;
        if (residual_sign < 0) ctx |= 0b100;

        return ctx;
    }

        static inline unsigned get_gt2_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);
        if      (sum <= 2)  return 0;
        else if (sum <= 7)  return 1;
        else if (sum <= 14) return 2;
        else if (sum <= 23) return 3;
        else if (sum <= 38) return 4;
        else if (sum <= 72) return 5;
        else                return 6;
    }

    static inline unsigned get_gt3_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);
        if      (sum <= 3)  return 0;
        else if (sum <= 9)  return 1;
        else if (sum <= 18) return 2;
        else if (sum <= 30) return 3;
        else if (sum <= 50) return 4;
        else if (sum <= 90) return 5;
        else                return 6;
    }

    static inline unsigned get_gt4_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);
        if      (sum <= 4)  return 0;
        else if (sum <= 11) return 1;
        else if (sum <= 22) return 2;
        else if (sum <= 35) return 3;
        else if (sum <= 60) return 4;
        else if (sum <= 110) return 5;
        else                return 6;
    }

public:
    EntropyCoderBase(unsigned groupSize)
        : GROUP_SIZE(groupSize)
    {
        for (auto& ctx : m_pmfAbsCtx)
            for (auto& p : ctx) p = 0;
        for (auto& p : m_pmfSignCtx) p = 0;
        for (auto& p : m_pmfSigCtx)  p = 0;
        for (auto& p : m_pmfGt1Ctx)  p = 0;
        for (auto& p : m_pmfGt2Ctx)  p = 0;
        for (auto& p : m_pmfGt3Ctx)  p = 0;
        for (auto& p : m_pmfGt4Ctx)  p = 0;
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
    void encLimitedGolomb(unsigned n, unsigned k, unsigned limit, unsigned qbpp);

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
    unsigned decLimitedGolomb(unsigned k, unsigned limit, unsigned qbpp);


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

void EntropyEncoder::encLimitedGolomb(
    unsigned n,
    unsigned k,
    unsigned limit,
    unsigned qbpp)
{
    unsigned q = n >> k;

    unsigned maxq = limit - qbpp - 1;

    if (q < maxq)
    {
        // unary: q zeros + 1
        for (unsigned i = 0; i < q; ++i)
            aenc.encBit(0);

        aenc.encBit(1);

        // remainder
        if (k > 0)
            aenc.encBits(n & ((1u << k) - 1), k);
    }
    else
    {
        // escape
        for (unsigned i = 0; i < maxq; ++i)
            aenc.encBit(0);

        aenc.encBit(1);

        // explicit binary of n-1
        aenc.encBits(n - 1, qbpp);
    }
}

unsigned EntropyDecoder::decLimitedGolomb(
    unsigned k,
    unsigned limit,
    unsigned qbpp)
{
    unsigned maxq = limit - qbpp - 1;

    unsigned q = 0;

    while (q < maxq && adec.decBit() == 0)
        ++q;

    if (q < maxq)
    {
        unsigned r = 0;

        if (k > 0)
            r = adec.decBits(k);

        return (q << k) | r;
    }
    else
    {
        return adec.decBits(qbpp) + 1;
    }
}

void EntropyEncoder::encodeSample(
    PGMImage::Sample merrval,
    unsigned /*ctxIdx*/,
    unsigned /*sigIdx*/,
    int /*a_err*/,
    int /*b_err*/,
    int /*c_err*/,
    int /*d_err*/,
    int /*q1*/,
    int /*q2*/,
    int /*q3*/,
    int k)
{
    constexpr unsigned LIMIT = 32;
    constexpr unsigned QBPP  = 8;

    encLimitedGolomb(
        unsigned(merrval),
        k,
        LIMIT,
        QBPP);
}
PGMImage::Sample EntropyDecoder::decodeSample(
    unsigned /*ctxIdx*/,
    unsigned /*sigIdx*/,
    int /*a_err*/,
    int /*b_err*/,
    int /*c_err*/,
    int /*d_err*/,
    int /*q1*/,
    int /*q2*/,
    int /*q3*/,
    int k)
{
    constexpr unsigned LIMIT = 32;
    constexpr unsigned QBPP  = 8;

    return PGMImage::Sample(
        decLimitedGolomb(
            k,
            LIMIT,
            QBPP));
}

//======================================================
//
//   S I G N I F I C A N C E   I N D E X
//
//   6 Schwellen => 7 Klassen (0..6)
//   Klasse 0: sum <= t1
//   Klasse 1: t1 < sum <= t2
//   ...
//   Klasse 6: sum > t6
//
//======================================================
inline unsigned compute_significant_Idx(
    int /*q1*/, int /*q2*/, int /*q3*/,
    int a_err, int b_err, int c_err, int d_err,
    int /*a*/, int /*b*/, int /*c*/, int /*d*/)
{
    const int T1 = 1;
    const int T2 = 3;
    const int T3 = 8;
    const int T4 = 14;
    const int T5 = 29;
    const int T6 = 61;

    int sum = std::abs(a_err) +
              std::abs(b_err) +
              std::abs(c_err) +
              std::abs(d_err);

    if (sum <= T1) return 0;
    else if (sum <= T2) return 1;
    else if (sum <= T3) return 2;
    else if (sum <= T4) return 3;
    else if (sum <= T5) return 4;
    else if (sum <= T6) return 5;
    else return 6;
}


inline int moduloRange(int v, int range = 256)
{
    return (v % range + range) % range;
}

inline int toSigned(int v, int range = 256)
{
    if (v >= range / 2)
        return v - range;
    return v;
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

        fs::create_directories(outputDir);
        if (isDecoding) {
            m_decDebugFile.open(outputDir + "/dec_debug.txt", std::ios::out);
            m_dec_first5.open(outputDir + "/dec_first5_values.txt", std::ios::out);
            m_kHistFile.open(outputDir + "/k_histogram.txt", std::ios::out);
        } else {
            m_encDebugFile.open(outputDir + "/enc_debug.txt", std::ios::out);
            m_enc_first5.open(outputDir + "/enc_first5_values.txt", std::ios::out);
        }
        bias_used    = 0;
        total_pixels = 0;
        RUNindex = 0;
        Nn[0] = 0;
        Nn[1] = 0;
        m_pixel_counter = 0;
        m_dec_pixel_counter = 0;
    }


    ~Prediction()
    {
        if (m_encDebugFile.is_open()) m_encDebugFile.close();
        if (m_decDebugFile.is_open()) m_decDebugFile.close();
        if (m_enc_first5.is_open()) m_enc_first5.close();
        if (m_dec_first5.is_open()) m_dec_first5.close();
        if (m_kHistFile.is_open()) m_kHistFile.close();
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
                codeSegmentA9(Errval);

                // Step 11. Compute the context-dependent Golomb variable k according to the procedure in Code segment A.10. 
                // A.5.1 Estimation of k - Golomb coding variable computation
                int k = 0;
                codeSegmentA10(k, Q);

                // k sammeln
                m_kHistogram[k]++;
                m_totalK += k;
                m_totalPixels++;


                // Step 12. Perform the error mapping according to the procedure in Code segment A.11.
                int MErrval;
                codeSegmentA11(MErrval, Errval, k, Q);

                // Debug: Erste 5 Pixel ausgeben
                if (m_pixel_counter < 800) {
                    writeDebugFirstFive(m_pixel_counter, q1, q2, q3, SIGN, Errval, MErrval, k, Q, Ix);
                    m_pixel_counter++;
                }

                // Step 13. Encode the mapped error value MErrval using the limited length Golomb code function LG(k, LIMIT), as 
                // specified in A.5.3
                // A.5.3 Mapped error encoding
                eenc.encodeSample(static_cast<PGMImage::Sample>(MErrval),
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
    writeKHistogram();
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
            int MErrval = edec.decodeSample(/*ctxIdx*/0, /*sigIdx*/0,
                              0,0,0,0,
                              q1,q2,q3,k);


            // Step 10. Perform  the  inverse  of  the  error  mapping  indicated  in  Code  segment  A.11,  where  now  MErrval  is  given  
            // and Errval is computed.
            int Errval;
            codeSegmentA11D(Errval, MErrval, k, Q);
            
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
            if (m_dec_pixel_counter < 800) {
                writeDebugFirstFive(m_dec_pixel_counter, q1, q2, q3, SIGN, Errval, MErrval, k, Q, -1, Rx);
                m_dec_pixel_counter++;
            }

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

    void writeDebugFirstFive(int iteration, int q1, int q2, int q3, int sign, 
                              int errval, int merrval, int k, int Q, int Ix = -1, int Rx = -1)
    {
        if (!m_enc_first5.is_open() && !m_dec_first5.is_open()) return;
        std::ostringstream oss;
        oss << "Iteration " << iteration << ": "
            << "q1=" << q1 << " | q2=" << q2 << " | q3=" << q3 << " | "
            << "sign=" << sign << " | errval=" << errval << " | "
            << "merrval=" << merrval << " | k=" << k << " | Q=" << Q;
        
        if (Ix != -1) oss << " | Ix=" << Ix;
        if (Rx != -1) oss << " | Rx=" << Rx;
        
        oss << std::endl;
        
        if (m_enc_first5.is_open())
            m_enc_first5 << oss.str();
        if (m_dec_first5.is_open())
            m_dec_first5 << oss.str();
    }

void writeKHistogram()
{
    if (!m_kHistFile.is_open()) {
        std::cerr << "Histogram file not open!\n";
        return;
    }

    std::cout << "writing histogram...\n";

    double avgK = (m_totalPixels > 0)
        ? double(m_totalK) / m_totalPixels
        : 0.0;

    m_kHistFile << "Average k-value: " << avgK << "\n";
    m_kHistFile << "Total pixels processed: " << m_totalPixels << "\n\n";

    m_kHistFile << "k-value distribution:\n";
    m_kHistFile << "k\tCount\tPercentage\n";
    m_kHistFile << "---\t-----\t----------\n";

    for (auto& [k, count] : m_kHistogram)
    {
        double perc = 100.0 * count / m_totalPixels;

        m_kHistFile << k << "\t"
                    << count << "\t"
                    << perc << "%\n";
    }

    m_kHistFile.flush();   // 🔥 WICHTIG
}


    void writeDebugSubtract(int idx, int x, int y, int a, int b, int c, int d,
                            int g1, int g2, int g3, int q1, int q2, int q3,
                            int sign, int pred_med, int pred_corr, int actual,
                            int error, int ctxIdx, const Context& ctx)
    {
        if (!m_encDebugFile.is_open()) return;
        m_encDebugFile
            << "idx=" << idx << " | (x,y)=(" << x << "," << y << ")"
            << " | [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "]"
            << " | grads={" << g1 << "," << g2 << "," << g3 << "}"
            << " | q={" << q1 << "," << q2 << "," << q3 << "}"
            << " | sign=" << sign << " | pred_med=" << pred_med
            << " | pred_corr=" << pred_corr << " | actual=" << actual
            << " | error=" << error << " | ctx=" << ctxIdx
            << " | ctx(A,B,N,C)=(" << ctx.A << "," << ctx.B << "," << ctx.N << "," << ctx.C << ")"
            << std::endl;
    }

    void writeDebugAdd(int idx, int x, int y, int a, int b, int c, int d,
                       int g1, int g2, int g3, int q1, int q2, int q3,
                       int sign, int pred_med, int pred_corr, int reconstructed,
                       int error, int ctxIdx, const Context& ctx)
    {
        if (!m_decDebugFile.is_open()) return;
        m_decDebugFile
            << "idx=" << idx << " | (x,y)=(" << x << "," << y << ")"
            << " | [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "]"
            << " | grads={" << g1 << "," << g2 << "," << g3 << "}"
            << " | q={" << q1 << "," << q2 << "," << q3 << "}"
            << " | sign=" << sign << " | pred_med=" << pred_med
            << " | pred_corr=" << pred_corr << " | reconstructed=" << reconstructed
            << " | error=" << error << " | ctx=" << ctxIdx
            << " | ctx(A,B,N,C)=(" << ctx.A << "," << ctx.B << "," << ctx.N << "," << ctx.C << ")"
            << std::endl;
    }

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

    std::ofstream m_encDebugFile;
    std::ofstream m_decDebugFile;
    std::ofstream m_enc_first5;
    std::ofstream m_dec_first5;
    std::ofstream m_kHistFile;
    bool          m_isDecoding;
    std::string   m_outputDir;
    int           m_pixel_counter = 0;
    int           m_dec_pixel_counter = 0;
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
    std::map<int, int> m_kHistogram;
    long long m_totalK = 0;
    long long m_totalPixels = 0;

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
    const std::filesystem::path outpath = std::filesystem::path(outputDir) / outname;
    std::ofstream stream(outpath, std::ios::binary);
    OBitstream bs(stream);

    bs.addFixed<unsigned>(img.getWidth(),  16);
    bs.addFixed<unsigned>(img.getHeight(), 16);

    EntropyEncoder eenc(bs, groupSize);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData(), false, outputDir);
    pred.subtractPrediction(eenc, groupSize);
    bs.byteAlign();
    stream.close();

    unsigned long filesize = static_cast<unsigned long>(std::filesystem::file_size(outpath));
    unsigned long pixels   = img.getWidth() * img.getHeight();
    double bpp = (double)(filesize * 8) / pixels;

    if (g_bpp_results_file.is_open()) {
        g_bpp_results_file
            << "GT1=[" << g_gt1_t1 << "," << g_gt1_t2 << "," << g_gt1_t3 << ","
            << g_gt1_t4 << "," << g_gt1_t5 << "," << g_gt1_t6 << "]"
            << " GTN=[" << g_gtn_t1 << "," << g_gtn_t2 << "," << g_gtn_t3 << ","
            << g_gtn_t4 << "," << g_gtn_t5 << "," << g_gtn_t6 << "]"
            << " | File=" << inname
            << " | Size=" << filesize << " bytes"
            << " | Pixels=" << pixels
            << " | BPP=" << std::fixed << std::setprecision(4) << bpp
            << std::endl;
        g_bpp_results_file.flush();
    }

    // Ausgabe-Format fuer PS-Regex:
    // "  GT1=[1,2,4,8,16,32] GTN=[1,2,4,8,16,32]: 4.4401 BPP"
    std::cout << "  GT1=["
              << g_gt1_t1 << "," << g_gt1_t2 << "," << g_gt1_t3 << ","
              << g_gt1_t4 << "," << g_gt1_t5 << "," << g_gt1_t6 << "]"
              << " GTN=["
              << g_gtn_t1 << "," << g_gtn_t2 << "," << g_gtn_t3 << ","
              << g_gtn_t4 << "," << g_gtn_t5 << "," << g_gtn_t6 << "]"
              << ": " << std::fixed << std::setprecision(4) << bpp << " BPP"
              << std::endl;
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
    const std::filesystem::path outpath = std::filesystem::path(outputDir) / outname;
    img.write(outpath.string());
}


int main(int argc, char** argv)
{
    cmdPars pars;
    if (!readCmdLine(argc, argv, pars)) return 1;

    fs::create_directories(pars.outdir);

    // GT1-Schwellen setzen
    g_gt1_t1 = pars.gt1_t1;  g_gt1_t2 = pars.gt1_t2;
    g_gt1_t3 = pars.gt1_t3;  g_gt1_t4 = pars.gt1_t4;
    g_gt1_t5 = pars.gt1_t5;  g_gt1_t6 = pars.gt1_t6;

    // GTN-Schwellen setzen
    g_gtn_t1 = pars.gtn_t1;  g_gtn_t2 = pars.gtn_t2;
    g_gtn_t3 = pars.gtn_t3;  g_gtn_t4 = pars.gtn_t4;
    g_gtn_t5 = pars.gtn_t5;  g_gtn_t6 = pars.gtn_t6;

    std::string fullOutputPath = pars.outdir + "/" + pars.outname;
    std::string bpp_log_file   = pars.outdir + "/threshold_bpp_results.txt";

    if (!g_bpp_results_file.is_open()) {
        g_bpp_results_file.open(bpp_log_file, std::ios::app);
        if (!g_bpp_results_file.is_open())
            std::cerr << "WARNING: Cannot open BPP results file: " << bpp_log_file << std::endl;
    }

    unsigned groupSize = 1;

    try {
        if (pars.decode) {
            std::cout << "Dekodiere: " << pars.inname << " -> " << fullOutputPath << std::endl;
            decode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        } else {
            std::cout << "Kodiere: " << pars.inname << " -> " << fullOutputPath
                      << " (GT1=[" << g_gt1_t1 << "," << g_gt1_t2 << "," << g_gt1_t3 << ","
                      << g_gt1_t4 << "," << g_gt1_t5 << "," << g_gt1_t6 << "]"
                      << " GTN=[" << g_gtn_t1 << "," << g_gtn_t2 << "," << g_gtn_t3 << ","
                      << g_gtn_t4 << "," << g_gtn_t5 << "," << g_gtn_t6 << "])" << std::endl;
            encode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        }
        std::cout << "OK" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Fehler: " << e.what() << std::endl;
        if (g_bpp_results_file.is_open()) g_bpp_results_file.close();
        return 1;
    }

    if (g_bpp_results_file.is_open()) g_bpp_results_file.close();
    return 0;
}