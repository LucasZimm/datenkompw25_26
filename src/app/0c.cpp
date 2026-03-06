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
struct cmdPars
{
    bool        decode;
    std::string inname;
    std::string outname;
    std::string outdir;   
    int         threshold1;  // NEU: Erste Schwelle für compute_significant_Idx
    int         threshold2;  // NEU: Zweite Schwelle für compute_significant_Idx
    int         threshold3;  // NEU: Dritte Schwelle für compute_significant_Idx
    int         threshold4;  // NEU: Vierte Schwelle für compute_significant_Idx
};


bool readCmdLine(int argc, char** argv, cmdPars& pars)
{
    std::stringstream err;
    int arg = 0;
    bool ok = true;

    ok = ok && (argc > ++arg);
    if (ok) {
        if (std::string(argv[arg]) == "-e")
            pars.decode = false;
        else if (std::string(argv[arg]) == "-d")
            pars.decode = true;
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
    if (ok) {
        pars.outname = argv[arg];
    }

    // ===== OPTIONAL: Output-Ordner =====
    if (argc > ++arg) {
        pars.outdir = argv[arg];
    } else {
        pars.outdir = "output";
    }

    // ===== OPTIONAL: Threshold1 =====
    if (argc > ++arg) {
        try {
            pars.threshold1 = std::stoi(argv[arg]);
            if (pars.threshold1 < 0 || pars.threshold1 > 1023) {
                err << "ERROR: Threshold1 must be between 0 and 1023." << std::endl;
                ok = false;
            }
        } catch (...) {
            err << "ERROR: Threshold1 must be an integer." << std::endl;
            ok = false;
        }
    } else {
        pars.threshold1 = 10;  // Default
    }

    // ===== OPTIONAL: Threshold2 =====
    if (argc > ++arg) {
        try {
            pars.threshold2 = std::stoi(argv[arg]);
            if (pars.threshold2 < 0 || pars.threshold2 > 1023) {
                err << "ERROR: Threshold2 must be between 0 and 1023." << std::endl;
                ok = false;
            }
        } catch (...) {
            err << "ERROR: Threshold2 must be an integer." << std::endl;
            ok = false;
        }
    } else {
        pars.threshold2 = 20;  // Default
    }

    // ===== OPTIONAL: Threshold3 =====
    if (argc > ++arg) {
        try {
            pars.threshold3 = std::stoi(argv[arg]);
            if (pars.threshold3 < 0 || pars.threshold3 > 1023) {
                err << "ERROR: Threshold3 must be between 0 and 1023." << std::endl;
                ok = false;
            }
        } catch (...) {
            err << "ERROR: Threshold3 must be an integer." << std::endl;
            ok = false;
        }
    } else {
        pars.threshold3 = 30;  // Default
    }

    // ===== OPTIONAL: Threshold4 =====
    if (argc > ++arg) {
        try {
            pars.threshold4 = std::stoi(argv[arg]);
            if (pars.threshold4 < 0 || pars.threshold4 > 1023) {
                err << "ERROR: Threshold4 must be between 0 and 1023." << std::endl;
                ok = false;
            }
        } catch (...) {
            err << "ERROR: Threshold4 must be an integer." << std::endl;
            ok = false;
        }
    } else {
        pars.threshold4 = 40;  // Default
    }

    // Validierung: threshold1 < threshold2 < threshold3 < threshold4
    if (ok && (pars.threshold1 >= pars.threshold2 || pars.threshold2 >= pars.threshold3 || pars.threshold3 >= pars.threshold4)) {
        err << "ERROR: Threshold1 < Threshold2 < Threshold3 < Threshold4 erforderlich." << std::endl;
        ok = false;
    }

    if (!ok) {
        std::string pname = argv[0];
        std::cerr << err.str() << std::endl;
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  decoding: " << pname << " -d inFile outFile [outputDir] [threshold1] [threshold2] [threshold3] [threshold4]" << std::endl;
        std::cerr << "  encoding: " << pname << " -e inFile outFile [outputDir] [threshold1] [threshold2] [threshold3] [threshold4]" << std::endl;
        std::cerr << "  threshold1-4: 0-1023 (optional, default: 10, 20, 30, 40) - threshold1 < threshold2 < threshold3 < threshold4" << std::endl;
        return false;
    }

    return true;
}






//======================================================
//
//   G L O B A L E   T E S T P A R A M E T E R
//
//======================================================
// Globale Schwellen für compute_significant_Idx (0-1023)
int g_threshold1 = 10;  // Erste Schwelle - wird beim Testen variiert
int g_threshold2 = 20;  // Zweite Schwelle - wird beim Testen variiert
int g_threshold3 = 30;  // Dritte Schwelle - wird beim Testen variiert
int g_threshold4 = 40;  // Vierte Schwelle - wird beim Testen variiert
std::ofstream g_bpp_results_file;  // Datei zum Speichern der BPP-Ergebnisse


//======================================================
//
//   E N T R O P Y   C O D I N G
//
//======================================================

class EntropyCoderBase
{
protected:
    static const unsigned N = 3;                     // max Bin Index für Magnitude
    static const unsigned NUM_CTX = 256 + 256 + 256 + 256; // Anzahl der Kontexte
    unsigned GROUP_SIZE;

    // PMF für Magnitude: [Kontext][Bin-Index]
    std::array<std::array<uint8_t, N>, NUM_CTX> m_pmfAbsCtx;

    // PMF für Sign: kontextabhängig
    std::array<uint8_t, NUM_CTX> m_pmfSignCtx;

    // PMF für Significance: kontextabhängig
    std::array<uint8_t, NUM_CTX> m_pmfSigCtx;

    inline unsigned mappedCtx(unsigned ctxIdx) const {
        // Für Testzwecke einfach 1 zurückgeben
        ctxIdx = 1;
        return ctxIdx;
    }

public:
    EntropyCoderBase(unsigned groupSize) 
        : GROUP_SIZE(groupSize)
    {
        // Initialisiere Magnitude-PMFs auf kleine Priors
        for (auto &ctx : m_pmfAbsCtx) {
            for (auto &p : ctx) p = 0;
        }

        // Initialisiere Sign-PMF auf kleinen Startwert
        for (auto &p : m_pmfSignCtx) p = 0;

        // Initialisiere Significance-PMF auf kleinen Startwert
        for (auto &p : m_pmfSigCtx) p = 0;
    }
};


class EntropyEncoder : protected EntropyCoderBase
{
public:
    EntropyEncoder(OBitstream &bs, unsigned groupSize)
        : EntropyCoderBase(groupSize), aenc(bs)
    {
        aenc.start();
    }
    void encodeSample(PGMImage::Sample s, unsigned ctxIdx, unsigned sigIdx);

    void finish() { aenc.finish(); }

private:
    ArithmeticEncoder aenc;
};


class EntropyDecoder : protected EntropyCoderBase
{
public:
    EntropyDecoder(IBitstream &bs, unsigned groupSize)
        : EntropyCoderBase(groupSize), adec(bs)
    {
        adec.start();
    }

    PGMImage::Sample decodeSample(unsigned ctxIdx, unsigned sigIdx);

private:
    ArithmeticDecoder adec;
};


// -----------------------------
// Encoder
// -----------------------------
void EntropyEncoder::encodeSample(PGMImage::Sample s,
                                  unsigned ctxIdx,
                                  unsigned sigIdx)
                                  
{
    unsigned absValue = unsigned(s < 0 ? -s : s);
    unsigned mIdx = mappedCtx(ctxIdx);

    // -----------------------------
    // 1) significance
    // -----------------------------
    bool nonZero = (absValue != 0);
    aenc.encBin(m_pmfSigCtx[sigIdx], nonZero);

    if (!nonZero)
        return;

    // -----------------------------
    // 2) sign (kontextabhängig!)
    // -----------------------------
    //unsigned signCtx = 1;
    aenc.encBit(s < 0);

    // -----------------------------
    // 3) magnitude (dein Code!)
    // -----------------------------
    unsigned rem = absValue - 1;
    unsigned binIdx = 1;

    while (rem--)
        aenc.encBin(
            m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)], 1);

    aenc.encBin(
        m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)], 0);
}



// -----------------------------
// Decoder
// -----------------------------
PGMImage::Sample EntropyDecoder::decodeSample(unsigned ctxIdx, unsigned sigIdx)
                                              
{
    unsigned mIdx = mappedCtx(ctxIdx);

    // -----------------------------
    // 1) significance
    // -----------------------------
    if (!adec.decBin(m_pmfSigCtx[sigIdx]))
        return 0;

    // -----------------------------
    // 2) sign
    // -----------------------------
    bool sign = adec.decBit();

    // -----------------------------
    // 3) magnitude
    // -----------------------------
    unsigned mag = 0;
    unsigned binIdx = 1;

    while (adec.decBin(
        m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)]))
    {
        mag++;
    }
    mag += 1;
    return sign ? -PGMImage::Sample(mag)
                :  PGMImage::Sample(mag);
}



inline unsigned compute_significant_Idx(
    int /*q1*/, int /*q2*/, int /*q3*/,
    int a_err, int b_err, int c_err, int d_err,
    int /*a*/, int /*b*/, int /*c*/, int /*d*/)
{
    // Testimplementierung: Vier Schwellen zur Einteilung in 5 Klassen
    int sum = std::abs(a_err) + std::abs(b_err) + std::abs(c_err) + std::abs(d_err);
    
    // Rückgabe: 0 wenn sum <= threshold1
    //          1 wenn threshold1 < sum <= threshold2
    //          2 wenn threshold2 < sum <= threshold3
    //          3 wenn threshold3 < sum <= threshold4
    //          4 wenn sum > threshold4
    // (5 Klassen)
    if (sum <= g_threshold1) return 0;
    else if (sum <= g_threshold2) return 1;
    else if (sum <= g_threshold3) return 2;
    else if (sum <= g_threshold4) return 3;
    else return 4;
}





//======================================================
//
//   LOCO-I PREDICTOR MIT DIAGNOSE
//
//======================================================
class Prediction
{
public:
    Prediction(int width, int height, std::vector<PGMImage::Sample>& img, bool isDecoding = false, const std::string& outputDir = "output")
        : m_width(width), m_height(height), m_data(img), m_org(m_data), m_isDecoding(isDecoding), m_outputDir(outputDir)
    {
        m_ctxs.resize(365);
        
        for (auto& ctx : m_ctxs) {
            ctx.N = 1;
            ctx.A = std::max(2, (256 + 32) / 64);
            ctx.B = 0;
            ctx.C = 0;
        }
        
        // Debug-Streams öffnen mit Pfad aus outputDir
        fs::create_directories(outputDir);
        
        if (isDecoding) {
            m_decDebugFile.open(outputDir + "/dec_debug.txt", std::ios::out);
        } else {
            m_encDebugFile.open(outputDir + "/enc_debug.txt", std::ios::out);
        }
        
        // Diagnose-Zähler
        bias_used = 0;
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
        PGMImage::Sample* data = m_data.data();
        std::vector<int16_t> residualRaw(m_width * m_height);
        int debug = 0;
        int zero_residuals = 0;
        int small_residuals = 0;
        gsize = gsize;  
        // --------------------
        // MiniBilder
        // --------------------

        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                total_pixels++;

                // STEP 0: Nachbarwerte laden (Rand = 0)
                auto get = [&](int xx, int yy) -> int {
                    if (xx < 0 || yy < 0) return 0;
                    if (xx >= m_width || yy >= m_height) return 0;
                    return int(orgData[yy * m_width + xx]);
                };
                auto folderrget = [&](int xx, int yy) -> int {
                    if (xx < 0 || yy < 0) return 0;
                    if (xx >= m_width || yy >= m_height) return 0;
                    return int(data[yy * m_width + xx]);
                };

                int a = get(x - 1, y);     // left
                int b = get(x,     y - 1); // up
                int c = get(x - 1, y - 1); // up-left
                int d = get(x + 1, y - 1); // up-right

                // STEP 1: Lokale Gradienten
                int g1 = d - b;
                int g2 = b - c;
                int g3 = c - a;

                // STEP 3 + 4: Gradienten quantisieren + Kontext
                int q1 = quantize(g1);
                int q2 = quantize(g2);
                int q3 = quantize(g3);

                int sign = 1;
                if (q1 < 0) sign = -1;
                else if (q1 == 0 && q2 < 0) sign = -1;
                else if (q1 == 0 && q2 == 0 && q3 < 0) sign = -1;

                int q1_norm = q1, q2_norm = q2, q3_norm = q3;
                if (sign < 0) {
                    q1_norm = -q1;
                    q2_norm = -q2;
                    q3_norm = -q3;
                }

                int ctxIdx = (q1_norm) * 81 + (q2_norm) * 9 + (q3_norm);
                //ctxIdx = std::clamp(ctxIdx, 0, 364);
                //if (ctxIdx == 0) ctxIdx = 1;



                Context& ctx = m_ctxs[ctxIdx];

                // STEP 5: MED-Prädiktion
                int pred_med;
                if (c >= std::max(a, b))
                    pred_med = std::min(a, b);
                else if (c <= std::min(a, b))
                    pred_med = std::max(a, b);
                else
                    pred_med = a + b - c;

                // STEP 6: Adaptive Korrektur
                int pred_corr = pred_med;
                if (sign > 0)
                    pred_corr += ctx.C;
                else
                    pred_corr -= ctx.C;
                pred_corr = std::clamp(pred_corr, 0, 255);

                if (std::abs(pred_corr - pred_med) > 0) {
                    bias_used++;
                }

                // STEP 7: Residuum berechnen
                int actual = int(orgData[y * m_width + x]);


                int error_code = actual - pred_corr;
                // Residuum speichern für Diagnose
                residualRaw[y * m_width + x] = static_cast<int16_t>(error_code);

                // Modulo-Faltung des Fehlers
                int e_fold = error_code;
                if (e_fold < 0) e_fold += 256;
                if (e_fold >= 128) e_fold -= 256;

                // Gefaltetes Residuum speichern
                data[y * m_width + x] = static_cast<PGMImage::Sample>(
                    std::clamp(e_fold, -128, 127)
                );

                // Kontextmodell für CABAC
                int a_folderrget = folderrget(x - 1, y);
                int b_folderrget = folderrget(x, y - 1);
                int c_folderrget = folderrget(x - 1, y - 1);
                int d_folderrget = folderrget(x + 1, y - 1);
                int model = contextModel(a_folderrget, b_folderrget, c_folderrget, d_folderrget);


                unsigned sigIdx = compute_significant_Idx(
                    q1, q2, q3,
                    folderrget(x - 1, y),
                    folderrget(x,     y - 1),
                    folderrget(x - 1, y - 1),
                    folderrget(x + 1, y - 1),
                    a, b, c, d
                );
                // Codieren mit CABAC
                eenc.encodeSample(static_cast<PGMImage::Sample>(e_fold), model, sigIdx);

                // Residuen-Statistiken
                if (error_code == 0) zero_residuals++;
                if (std::abs(error_code) <= 1) small_residuals++;


                // Debug-Ausgabe
                if (debug == 1) {
                    writeDebugSubtract(y * m_width + x, x, y, a, b, c, d, g1, g2, g3, 
                                       q1, q2, q3, sign, pred_med, pred_corr, actual, 
                                       error_code, ctxIdx, ctx);
                }

                // Context update
                int error_ctx = (sign < 0) ? -error_code : error_code;
                updateContextLocoI(ctxIdx, error_ctx);
            }
        }

        eenc.finish();

        std::cout << "\n=== SUBTRACTPREDICTION SUMMARY ===\n"
                  << "Total pixels: " << total_pixels << "\n"
                  << "Bias corrections: " << bias_used << " (" << (100.0*bias_used/total_pixels) << "%)\n"
                  << "Zero residuals: " << zero_residuals << " (" << (100.0*zero_residuals/total_pixels) << "%)\n"
                  << "Small residuals (|e|<=1): " << small_residuals << " (" << (100.0*small_residuals/total_pixels) << "%)\n"
                  << "Histogramme gespeichert im Ordner: " << m_outputDir << "\n"
                  << std::endl;
    }


    void addPrediction(EntropyDecoder& edec)
{
    // Temporären Puffer für rekonstruierte Daten erstellen
    std::vector<PGMImage::Sample> tempData(m_data.size());
    std::vector<int16_t> data(m_width * m_height);

    int debug = 0;
    for (int y = 0; y < m_height; ++y)
    {
        for (int x = 0; x < m_width; ++x)
        {
            // Für Prädiktion die bereits rekonstruierten Nachbarn verwenden
            auto get = [&](int xx, int yy) -> int {
                if (xx < 0 || yy < 0) return 0;
                if (xx >= m_width || yy >= m_height) return 0;
                return int(tempData[yy * m_width + xx]);
            };
            auto folderrget = [&](int xx, int yy) -> int {
                if (xx < 0 || yy < 0) return 0;
                if (xx >= m_width || yy >= m_height) return 0;
                return int(data[yy * m_width + xx]);
            };

            int a = get(x - 1, y);     // left
            int b = get(x,     y - 1); // up
            int c = get(x - 1, y - 1); // up-left
            int d = get(x + 1, y - 1); // up-right

            // STEP 1: Gradienten
            int g1 = d - b;
            int g2 = b - c;
            int g3 = c - a;

            // STEP 3 + 4: Quantisierung + Kontext
            int q1 = quantize(g1);
            int q2 = quantize(g2);
            int q3 = quantize(g3);

            int sign = 1;
            if (q1 < 0) sign = -1;
            else if (q1 == 0 && q2 < 0) sign = -1;
            else if (q1 == 0 && q2 == 0 && q3 < 0) sign = -1;

            int q1_norm = q1, q2_norm = q2, q3_norm = q3;
            if (sign < 0) {
                q1_norm = -q1;
                q2_norm = -q2;
                q3_norm = -q3;
            }

            int ctxIdx = (q1_norm) * 81 + (q2_norm) * 9 + (q3_norm);
            //ctxIdx = std::clamp(ctxIdx, 0, 364);
            //if (ctxIdx == 0) ctxIdx = 1;

            Context& ctx = m_ctxs[ctxIdx];

            // STEP 5: MED-Prädiktion
            int pred_med;
            if (c >= std::max(a, b))
                pred_med = std::min(a, b);
            else if (c <= std::min(a, b))
                pred_med = std::max(a, b);
            else
                pred_med = a + b - c;

            // STEP 6: Adaptive Korrektur
            int pred_corr = pred_med;
            if (sign > 0)
                pred_corr += ctx.C;
            else
                pred_corr -= ctx.C;
            pred_corr = std::clamp(pred_corr, 0, 255);

            // Kontextmodell für CABAC
            int a_folderrget = folderrget(x - 1, y);
            int b_folderrget = folderrget(x, y - 1);
            int c_folderrget = folderrget(x - 1, y - 1);
            int d_folderrget = folderrget(x + 1, y - 1);
            int model = contextModel(a_folderrget, b_folderrget, c_folderrget, d_folderrget);

            unsigned sigIdx = compute_significant_Idx(
                    q1, q2, q3,
                    folderrget(x - 1, y),
                    folderrget(x,     y - 1),
                    folderrget(x - 1, y - 1),
                    folderrget(x + 1, y - 1),
                    a, b, c, d
                );

            // STEP 7: Residuum direkt vom Decoder holen
            int e_fold = edec.decodeSample(model, sigIdx);
                      
            // Gefaltetes Residuum speichern
            data[y * m_width + x] = static_cast<PGMImage::Sample>(
                std::clamp(e_fold, -128, 127)
            );


            // Rekonstruktion mit modulo-Addition
            int reconstructed = (pred_corr + e_fold);
            if (reconstructed < 0) reconstructed += 256;
            else if (reconstructed >= 256) reconstructed -= 256;
            tempData[y * m_width + x] = static_cast<PGMImage::Sample>(reconstructed);   
            // echten Fehler für Kontextupdate neu berechnen
            int error_code = reconstructed - pred_corr;
            int error_ctx = (sign < 0) ? -error_code : error_code;


            // Debug
            if (debug ==1) {
             
            writeDebugAdd(y * m_width + x, x, y, a, b, c, d,
                          g1, g2, g3, q1, q2, q3, sign,
                          pred_med, pred_corr, reconstructed,
                          error_code, ctxIdx, ctx);
            }

            // STEP 11-12: Kontextupdate mit Residuum
            updateContextLocoI(ctxIdx, error_ctx);
        }
    }

    // Rekonstruierte Daten zurück in m_data kopieren
    m_data = tempData;

    std::cout << "\n=== ADDPREDICTION SUMMARY ===\n"
              << "Reconstruction complete\n" << std::endl;
}


private:
    // Struktur ZUERST definieren
    struct Context
    {
        int A = 0;  // Sum of absolute errors
        int B = 0;  // Sum of errors (for bias)
        int N = 1;  // Context count
        int C = 0;  // Bias correction value
    };

    // Debug-Methoden DANACH
    void writeDebugSubtract(int idx, int x, int y, int a, int b, int c, int d,
                           int g1, int g2, int g3, int q1, int q2, int q3,
                           int sign, int pred_med, int pred_corr, int actual,
                           int error, int ctxIdx, const Context& ctx)
    {
        if (!m_encDebugFile.is_open()) return;
        
        m_encDebugFile << "idx=" << idx
                      << " | (x,y)=(" << x << "," << y << ")"
                      << " | [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "]"
                      << " | grads={" << g1 << "," << g2 << "," << g3 << "}"
                      << " | q={" << q1 << "," << q2 << "," << q3 << "}"
                      << " | sign=" << sign
                      << " | pred_med=" << pred_med
                      << " | pred_corr=" << pred_corr
                      << " | actual=" << actual
                      << " | error=" << error
                      << " | ctx=" << ctxIdx
                      << " | ctx(A,B,N,C)=(" << ctx.A << "," << ctx.B << "," << ctx.N << "," << ctx.C << ")"
                      << std::endl;
    }

    void writeDebugAdd(int idx, int x, int y, int a, int b, int c, int d,
                      int g1, int g2, int g3, int q1, int q2, int q3,
                      int sign, int pred_med, int pred_corr, int reconstructed,
                      int error, int ctxIdx, const Context& ctx)
    {
        if (!m_decDebugFile.is_open()) return;
        
        m_decDebugFile << "idx=" << idx
                      << " | (x,y)=(" << x << "," << y << ")"
                      << " | [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "]"
                      << " | grads={" << g1 << "," << g2 << "," << g3 << "}"
                      << " | q={" << q1 << "," << q2 << "," << q3 << "}"
                      << " | sign=" << sign
                      << " | pred_med=" << pred_med
                      << " | pred_corr=" << pred_corr
                      << " | reconstructed=" << reconstructed
                      << " | error=" << error
                      << " | ctx=" << ctxIdx
                      << " | ctx(A,B,N,C)=(" << ctx.A << "," << ctx.B << "," << ctx.N << "," << ctx.C << ")"
                      << std::endl;
    }

    // LOCO-I konforme Context-Aktualisierung
    inline void updateContextLocoI(int ctxIdx, int error)
    {
        if (ctxIdx < 1 || ctxIdx >= 365) return;
        
        Context& ctx = m_ctxs[ctxIdx];
        
        // STEP 11: Zähler aktualisieren
        ctx.B += error;
        ctx.A += std::abs(error);

        // Reset-Schwellenwert
        const int N_RESET = 64;
        if (ctx.N >= N_RESET) {
            ctx.A >>= 1;
            if (ctx.B >= 0)
                ctx.B >>= 1;
            else
                ctx.B = -((1-ctx.B) >> 1);
            ctx.N >>= 1;
        }
        ctx.N += 1;

        // STEP 12: B und C aktualisieren nach LOCO-I Schema
        if (ctx.B <= -ctx.N) {
            ctx.C -= 1;
            ctx.B += ctx.N;
            if (ctx.B <= -ctx.N) ctx.B = -ctx.N + 1;
        } else if (ctx.B > 0) {
            ctx.C += 1;
            ctx.B -= ctx.N;
            if (ctx.B > 0) ctx.B = 0;
        }
        
        // C clampen
        ctx.C = std::clamp(ctx.C, -128, 127);
    }
    inline int quantize(int grad) const
    {
        int T1 = 3, T2 = 7, T3 = 21;
        if (grad <= -T3) return -4;
        else if (grad <= -T2) return -3;
        else if (grad <= -T1) return -2;
        else if (grad < 0)    return -1;
        else if (grad == 0)   return 0;
        else if (grad < T1)   return 1;
        else if (grad < T2)   return 2;
        else if (grad < T3)   return 3;
        else                  return 4;
    }

    int contextModel( int a, int b, int c, int d) const
    {
        return  (abs(a)+abs(b)+abs(c)+abs(d)); // Summe der Nachbarn
    }

private:
    // Debug-Streams
    std::ofstream m_encDebugFile;
    std::ofstream m_decDebugFile;
    bool m_isDecoding;
    std::string m_outputDir;  // NEU hinzufügen
    
    // Diagnose-Variablen
    int bias_used;
    int total_pixels;
    const int m_alphabetSize = 256;
    int m_width;
    int m_height;

    std::vector<PGMImage::Sample>& m_data;
    std::vector<PGMImage::Sample> m_org;

    std::vector<Context> m_ctxs;
};



//======================================================
//
//   M A I N ENCODING + DECODING
//
//======================================================
void encode(const std::string& inname, const std::string& outname, unsigned groupSize, const std::string& outputDir = "output") {
    PGMImage img;
    img.read(inname);
    // Ensure output directory exists and write the encoded file into it (so batch scripts find it)
    std::filesystem::create_directories(outputDir);
    const std::filesystem::path outpath = std::filesystem::path(outputDir) / outname;
    std::ofstream stream(outpath, std::ios::binary);
    OBitstream bs(stream);

    bs.addFixed<unsigned>(img.getWidth(), 16);
    bs.addFixed<unsigned>(img.getHeight(), 16);

    EntropyEncoder eenc(bs, groupSize);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData(), false, outputDir);
    pred.subtractPrediction(eenc, groupSize);
    bs.byteAlign();
    
    // BPP berechnen und speichern
    stream.close();
    unsigned long filesize = static_cast<unsigned long>(std::filesystem::file_size(outpath));
    
    unsigned long pixels = img.getWidth() * img.getHeight();
    unsigned long bits = filesize * 8;
    double bpp = (double)bits / pixels;
    
    if (g_bpp_results_file.is_open()) {
        g_bpp_results_file << "Threshold1=" << g_threshold1 
                          << ", Threshold2=" << g_threshold2
                          << ", Threshold3=" << g_threshold3
                          << ", Threshold4=" << g_threshold4
                          << " | File=" << inname 
                          << " | Size=" << filesize << " bytes"
                          << " | Pixels=" << pixels
                          << " | BPP=" << std::fixed << std::setprecision(4) << bpp << std::endl;
        g_bpp_results_file.flush();
    }
    
    std::cout << "  Thresholds [" << g_threshold1 << "," << g_threshold2 << "," << g_threshold3 << "," << g_threshold4 << "]: " << std::fixed << std::setprecision(4) << bpp << " BPP" << std::endl;
}

void decode(const std::string& inname, const std::string& outname, unsigned groupSize, const std::string& outputDir = "output") {
    std::ifstream stream(inname, std::ios::binary);
    IBitstream bs(stream);
    int width = bs.getFixed<unsigned>(16);
    int height = bs.getFixed<unsigned>(16);

    PGMImage img(width, height);
    EntropyDecoder edec(bs, groupSize);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData(), true, outputDir);
    pred.addPrediction(edec);
    // Ensure output dir exists and write decoded image into it
    std::filesystem::create_directories(outputDir);
    const std::filesystem::path outpath = std::filesystem::path(outputDir) / outname;
    img.write(outpath.string());
}


int main(int argc, char** argv)
{
    cmdPars pars;

    if (!readCmdLine(argc, argv, pars))
        return 1;

    fs::create_directories(pars.outdir);
    
    // Globale Schwellen setzen
    g_threshold1 = pars.threshold1;
    g_threshold2 = pars.threshold2;
    g_threshold3 = pars.threshold3;
    g_threshold4 = pars.threshold4;

    std::string fullOutputPath = pars.outdir + "/" + pars.outname;

    // BPP-Ergebnisse in separate Datei speichern
    std::string bpp_log_file = pars.outdir + "/threshold_bpp_results.txt";
    
    // Append mode - damit können mehrere Schwellentests hintereinander gemacht werden
    if (!g_bpp_results_file.is_open()) {
        g_bpp_results_file.open(bpp_log_file, std::ios::app);
        if (!g_bpp_results_file.is_open()) {
            std::cerr << "WARNING: Cannot open BPP results file: " << bpp_log_file << std::endl;
        }
    }

    unsigned groupSize = 1;

    try {
        if (pars.decode) {
            std::cout << "Dekodiere: " << pars.inname << " -> " << fullOutputPath << std::endl;
            decode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        } else {
            std::cout << "Kodiere: " << pars.inname << " -> " << fullOutputPath 
                     << " (Thresholds=[" << g_threshold1 << "," << g_threshold2 << "," << g_threshold3 << "," << g_threshold4 << "])" << std::endl;
            encode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        }
        std::cout << "✅ Erfolgreich!" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Fehler: " << e.what() << std::endl;
        if (g_bpp_results_file.is_open()) g_bpp_results_file.close();
        return 1;
    }

    if (g_bpp_results_file.is_open()) g_bpp_results_file.close();
    return 0;
}
