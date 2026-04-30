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
        m_ctxs.resize(365);
        for (auto& ctx : m_ctxs) {
            ctx.N = 1;
            ctx.A = std::max(2, (256 + 32) / 64);
            ctx.B = 0;
            ctx.C = 0;
        }
        fs::create_directories(outputDir);
        if (isDecoding)
            m_decDebugFile.open(outputDir + "/dec_debug.txt", std::ios::out);
        else
            m_encDebugFile.open(outputDir + "/enc_debug.txt", std::ios::out);
        bias_used    = 0;
        total_pixels = 0;
    }

    ~Prediction()
    {
        if (m_encDebugFile.is_open()) m_encDebugFile.close();
        if (m_decDebugFile.is_open()) m_decDebugFile.close();
    }

    void subtractPrediction(EntropyEncoder& eenc, int gsize)
    {
        const PGMImage::Sample* orgData = m_org.data();
        PGMImage::Sample*       data    = m_data.data();
        std::vector<int16_t> residualRaw(m_width * m_height);
        int debug = 0, zero_residuals = 0, small_residuals = 0;
        (void)gsize;

        for (int y = 0; y < m_height; ++y) {
            for (int x = 0; x < m_width; ++x) {
                total_pixels++;

                auto get = [&](int xx, int yy) -> int {
                    if (xx < 0 || yy < 0 || xx >= m_width || yy >= m_height) return 0;
                    return int(orgData[yy * m_width + xx]);
                };
                auto folderrget = [&](int xx, int yy) -> int {
                    if (xx < 0 || yy < 0 || xx >= m_width || yy >= m_height) return 0;
                    return int(data[yy * m_width + xx]);
                };

                int a = get(x-1,y), b = get(x,y-1), c = get(x-1,y-1), d = get(x+1,y-1);
                int g1 = d-b, g2 = b-c, g3 = c-a;
                int q1 = quantize(g1), q2 = quantize(g2), q3 = quantize(g3);

                int sign = 1;
                if      (q1 < 0)                        sign = -1;
                else if (q1 == 0 && q2 < 0)             sign = -1;
                else if (q1 == 0 && q2 == 0 && q3 < 0)  sign = -1;

                int q1_n = q1, q2_n = q2, q3_n = q3;
                if (sign < 0) { q1_n = -q1; q2_n = -q2; q3_n = -q3; }

                int ctxIdx = q1_n*81 + q2_n*9 + q3_n;
                Context& ctx = m_ctxs[ctxIdx];

                int pred_med;
                if      (c >= std::max(a,b)) pred_med = std::min(a,b);
                else if (c <= std::min(a,b)) pred_med = std::max(a,b);
                else                         pred_med = a+b-c;

                int pred_corr = pred_med + (sign > 0 ? ctx.C : -ctx.C);
                pred_corr = std::clamp(pred_corr, 0, 255);
                if (std::abs(pred_corr - pred_med) > 0) bias_used++;

                int actual     = int(orgData[y * m_width + x]);
                int error_code = actual - pred_corr;
                residualRaw[y * m_width + x] = static_cast<int16_t>(error_code);

                int e_fold = error_code;
                if (e_fold < 0)    e_fold += 256;
                if (e_fold >= 128) e_fold -= 256;
                data[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(e_fold, -128, 127));

                int fe_a = folderrget(x-1,y), fe_b = folderrget(x,y-1);
                int fe_c = folderrget(x-1,y-1), fe_d = folderrget(x+1,y-1);
                int model    = contextModel(fe_a, fe_b, fe_c, fe_d);
                unsigned sig = compute_significant_Idx(q1,q2,q3, fe_a,fe_b,fe_c,fe_d, a,b,c,d);

                int k = 0;
                // for(k=0; (N[Q]<<k)<A[Q]; k++); 
                while ((ctx.N << k) < ctx.A) k++;

                eenc.encodeSample(static_cast<PGMImage::Sample>(e_fold), model, sig, fe_a, fe_b, fe_c, fe_d, q1, q2, q3, k);

                if (error_code == 0)           zero_residuals++;
                if (std::abs(error_code) <= 1) small_residuals++;

                if (debug == 1)
                    writeDebugSubtract(y*m_width+x, x, y, a, b, c, d,
                                       g1, g2, g3, q1, q2, q3, sign,
                                       pred_med, pred_corr, actual, error_code, ctxIdx, ctx);

                updateContextLocoI(ctxIdx, (sign < 0) ? -error_code : error_code);
            }
        }

        eenc.finish();
        std::cout << "\n=== SUBTRACTPREDICTION SUMMARY ===\n"
                  << "Total pixels: "            << total_pixels << "\n"
                  << "Bias corrections: "         << bias_used
                  << " (" << (100.0*bias_used/total_pixels) << "%)\n"
                  << "Zero residuals: "           << zero_residuals
                  << " (" << (100.0*zero_residuals/total_pixels) << "%)\n"
                  << "Small residuals (|e|<=1): " << small_residuals
                  << " (" << (100.0*small_residuals/total_pixels) << "%)\n"
                  << "Output dir: " << m_outputDir << "\n" << std::endl;
    }


    void addPrediction(EntropyDecoder& edec)
    {
        std::vector<PGMImage::Sample> tempData(m_data.size());
        std::vector<int16_t>          data(m_width * m_height);
        int debug = 0;

        for (int y = 0; y < m_height; ++y) {
            for (int x = 0; x < m_width; ++x) {
                auto get = [&](int xx, int yy) -> int {
                    if (xx < 0 || yy < 0 || xx >= m_width || yy >= m_height) return 0;
                    return int(tempData[yy * m_width + xx]);
                };
                auto folderrget = [&](int xx, int yy) -> int {
                    if (xx < 0 || yy < 0 || xx >= m_width || yy >= m_height) return 0;
                    return int(data[yy * m_width + xx]);
                };

                int a = get(x-1,y), b = get(x,y-1), c = get(x-1,y-1), d = get(x+1,y-1);
                int g1 = d-b, g2 = b-c, g3 = c-a;
                int q1 = quantize(g1), q2 = quantize(g2), q3 = quantize(g3);

                int sign = 1;
                if      (q1 < 0)                        sign = -1;
                else if (q1 == 0 && q2 < 0)             sign = -1;
                else if (q1 == 0 && q2 == 0 && q3 < 0)  sign = -1;

                int q1_n = q1, q2_n = q2, q3_n = q3;
                if (sign < 0) { q1_n = -q1; q2_n = -q2; q3_n = -q3; }

                int ctxIdx = q1_n*81 + q2_n*9 + q3_n;
                Context& ctx = m_ctxs[ctxIdx];

                int pred_med;
                if      (c >= std::max(a,b)) pred_med = std::min(a,b);
                else if (c <= std::min(a,b)) pred_med = std::max(a,b);
                else                         pred_med = a+b-c;

                int pred_corr = pred_med + (sign > 0 ? ctx.C : -ctx.C);
                pred_corr = std::clamp(pred_corr, 0, 255);

                int fe_a = folderrget(x-1,y), fe_b = folderrget(x,y-1);
                int fe_c = folderrget(x-1,y-1), fe_d = folderrget(x+1,y-1);
                int model    = contextModel(fe_a, fe_b, fe_c, fe_d);
                unsigned sig = compute_significant_Idx(q1,q2,q3, fe_a,fe_b,fe_c,fe_d, a,b,c,d);

                int k = 0;
                // for(k=0; (N[Q]<<k)<A[Q]; k++); 
                while ((ctx.N << k) < ctx.A) k++;

                int e_fold = edec.decodeSample(model, sig, fe_a, fe_b, fe_c, fe_d, q1, q2, q3, k);
                data[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(e_fold, -128, 127));

                int reconstructed = pred_corr + e_fold;
                if (reconstructed < 0)         reconstructed += 256;
                else if (reconstructed >= 256) reconstructed -= 256;
                tempData[y * m_width + x] = static_cast<PGMImage::Sample>(reconstructed);

                int error_code = reconstructed - pred_corr;
                if (debug == 1)
                    writeDebugAdd(y*m_width+x, x, y, a, b, c, d,
                                  g1, g2, g3, q1, q2, q3, sign,
                                  pred_med, pred_corr, reconstructed, error_code, ctxIdx, ctx);

                updateContextLocoI(ctxIdx, (sign < 0) ? -error_code : error_code);
            }
        }
        m_data = tempData;
        std::cout << "\n=== ADDPREDICTION SUMMARY ===\nReconstruction complete\n" << std::endl;
    }


private:
    struct Context { int A=0, B=0, N=1, C=0; };

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
        if (ctxIdx < 1 || ctxIdx >= 365) return;
        Context& ctx = m_ctxs[ctxIdx];
        ctx.B += error;
        ctx.A += std::abs(error);
        const int N_RESET = 64;
        if (ctx.N >= N_RESET) {
            ctx.A >>= 1;
            if (ctx.B >= 0) ctx.B >>= 1;
            else            ctx.B = -((1 - ctx.B) >> 1);
            ctx.N >>= 1;
        }
        ctx.N += 1;
        if (ctx.B <= -ctx.N) {
            ctx.C -= 1; ctx.B += ctx.N;
            if (ctx.B <= -ctx.N) ctx.B = -ctx.N + 1;
        } else if (ctx.B > 0) {
            ctx.C += 1; ctx.B -= ctx.N;
            if (ctx.B > 0) ctx.B = 0;
        }
        ctx.C = std::clamp(ctx.C, -128, 127);
    }

    inline int quantize(int grad) const
    {
        const int T1=3, T2=7, T3=21;
        if      (grad <= -T3) return -4;
        else if (grad <= -T2) return -3;
        else if (grad <= -T1) return -2;
        else if (grad <    0) return -1;
        else if (grad ==   0) return  0;
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
    bool          m_isDecoding;
    std::string   m_outputDir;
    int           bias_used;
    int           total_pixels;
    const int     m_alphabetSize = 256;
    int           m_width;
    int           m_height;

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