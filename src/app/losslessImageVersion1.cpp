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

    if (!ok) {
        std::string pname = argv[0];
        std::cerr << err.str() << std::endl;
        std::cerr << "Usage:" << std::endl;
        std::cerr << "  decoding: " << pname << " -d inFile outFile [outputDir]" << std::endl;
        std::cerr << "  encoding: " << pname << " -e inFile outFile [outputDir]" << std::endl;
        return false;
    }

    return true;
}






//======================================================
//
//   E N T R O P Y   C O D I N G
//
//======================================================
class EntropyCoderBase
{
protected:
    EntropyCoderBase(unsigned groupSize) : GROUP_SIZE(groupSize)
    {
        for (auto &ctx : m_pmfAbsCtx) {
            for (auto &p : ctx) p = 0;
        }
    }

protected:
    static const unsigned N = 3;
    static const unsigned NUM_CTX = 365;
    unsigned GROUP_SIZE; // jetzt variabel
    std::array<std::array<uint8_t, N>, NUM_CTX> m_pmfAbsCtx;

    inline unsigned mappedCtx(unsigned ctxIdx) const {
        if (ctxIdx >= NUM_CTX) ctxIdx = NUM_CTX - 1;
        return (ctxIdx / GROUP_SIZE) * GROUP_SIZE;
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
    void encodeSample(PGMImage::Sample s, unsigned ctxIdx);

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

    PGMImage::Sample decodeSample(unsigned ctxIdx);

private:
    ArithmeticDecoder adec;
};


// -----------------------------
// Encoder
// -----------------------------
void EntropyEncoder::encodeSample(PGMImage::Sample s, unsigned ctxIdx)
{
    unsigned absValue = unsigned(s < 0 ? -s : s);
    unsigned rem = absValue;
    unsigned binIdx = 0;

    unsigned mIdx = mappedCtx(ctxIdx); // gemappten Kontext verwenden

    while (rem--)
        aenc.encBin(m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)], 1);

    aenc.encBin(m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)], 0);

    if (absValue)
        aenc.encBit(s < 0);
}


// -----------------------------
// Decoder
// -----------------------------
PGMImage::Sample EntropyDecoder::decodeSample(unsigned ctxIdx)
{
    PGMImage::Sample s = 0;
    unsigned binIdx = 0;

    unsigned mIdx = mappedCtx(ctxIdx); // gemappten Kontext verwenden

    while (adec.decBin(m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)]))
        s++;

    if (s && adec.decBit())
        s = -s;

    return s;
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
            ctx.A = std::max(2, (width * height + 32) / 64);
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
        int debug = 1;
        int zero_residuals = 0;
        int small_residuals = 0;

        // --------------------
        // Histogramme initialisieren
        // --------------------
        // Kontext-Histogramme
        std::array<size_t, 256> histoA_ctx{};
        std::array<size_t, 256> histoB_ctx{};
        std::array<size_t, 256> histoC_ctx{};
        std::array<size_t, 256> histoD_ctx{};
        // Fehler-Histogramme
        std::array<size_t, 511> histoA_err{};
        std::array<size_t, 511> histoB_err{};
        std::array<size_t, 511> histoC_err{};
        std::array<size_t, 511> histoD_err{};
        std::array<size_t, 511> histoAll_err{};
        // --------------------
        // MiniBilder
        // --------------------
        struct ContextStats {
            uint64_t sumA = 0, sumB = 0, sumC = 0, sumD = 0, sumX = 0;
            uint64_t count = 0;
        };

        std::vector<ContextStats> ctxStats(366); // 0..364

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

                int ctxIdx = (q1_norm) * 25 + (q2_norm) * 5 + (q3_norm);
                ctxIdx = std::clamp(ctxIdx, 0, 364);
                if (ctxIdx == 0) ctxIdx = 1;



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

                ctxStats[ctxIdx].sumA += a;
                ctxStats[ctxIdx].sumB += b;
                ctxStats[ctxIdx].sumC += c;
                ctxStats[ctxIdx].sumD += d;
                ctxStats[ctxIdx].sumX += actual;
                ctxStats[ctxIdx].count++;


                int error_code = actual - pred_corr;
                residualRaw[y * m_width + x] = static_cast<int16_t>(error_code);

                // Modulo-Faltung des Fehlers
                int e_fold = ((error_code % 256) + 256) % 256;
                if (e_fold >= 128) e_fold -= 256;

                // Gefaltetes Residuum speichern
                data[y * m_width + x] = static_cast<PGMImage::Sample>(
                    std::clamp(e_fold, -128, 127)
                );

                // Codieren mit CABAC
                eenc.encodeSample(static_cast<PGMImage::Sample>(e_fold), ctxIdx);

                // Residuen-Statistiken
                if (error_code == 0) zero_residuals++;
                if (std::abs(error_code) <= 1) small_residuals++;

                // --------------------
                // Histogramme füllen
                // --------------------
                // Kontext-Histogramme
                histoA_ctx[a]++;
                histoB_ctx[b]++;
                histoC_ctx[c]++;
                histoD_ctx[d]++;
                int eA = (actual - a)+255; // Fehler nur relativ zum linken Nachbarn
                int eB = (actual - b)+255; // Fehler nur relativ zum oberen Nachbarn
                int eC = (actual - c)+255; // Fehler nur relativ zum oben-links Nachbarn
                int eD = (actual - d)+255; // Fehler nur relativ zum oben-rechts Nachbarn
                histoA_err[eA]++;
                histoB_err[eB]++;
                histoC_err[eC]++;
                histoD_err[eD]++;   
                histoAll_err[error_code + 255]++;               

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

        // --------------------
        // Histogramme abspeichern
        // --------------------
        auto saveHistogram = [&](const std::array<size_t,256>& hist, const std::string& name) {
            std::ofstream out(m_outputDir + "/" + name + ".txt");
            for (size_t i = 0; i < hist.size(); ++i)
                out << i << " " << hist[i] << "\n";
        };
        auto saveErrorHistogram = [&](const std::array<size_t,511>& hist, const std::string& name) {
            std::ofstream out(m_outputDir + "/" + name + ".txt");
            for (size_t i = 0; i < hist.size(); ++i)
                out << int(i-255) << " " << hist[i] << "\n";
        };
        
        // Kontext-Histogramme
        saveHistogram(histoA_ctx, "histoA_ctx");
        saveHistogram(histoB_ctx, "histoB_ctx");
        saveHistogram(histoC_ctx, "histoC_ctx");
        saveHistogram(histoD_ctx, "histoD_ctx");

        // Fehler-Histogramme
        saveErrorHistogram(histoA_err, "histoA_err");
        saveErrorHistogram(histoB_err, "histoB_err");
        saveErrorHistogram(histoC_err, "histoC_err");
        saveErrorHistogram(histoD_err, "histoD_err");
        saveErrorHistogram(histoAll_err, "histoAll_err");

        // --------------------
        // Residualbilder abspeichern
        // --------------------
        auto saveRawPGM = [&](const std::string& filename) {
            std::ofstream out(filename, std::ios::binary);
            if (!out.is_open()) return;

            out << "P5\n" << m_width << " " << m_height << "\n255\n";
            for (int i = 0; i < m_width * m_height; ++i) {
                int e = std::clamp<int>(residualRaw[i], -128, 127);
                uint8_t v = static_cast<uint8_t>(e + 128);
                out.write(reinterpret_cast<const char*>(&v), 1);
            }
        };
        auto saveRawPGM3color = [&](const std::string& filename) {
            std::ofstream out(filename, std::ios::binary);
            if (!out.is_open()) return;

            out << "P5\n" << m_width << " " << m_height << "\n255\n";
            for (int i = 0; i < m_width * m_height; ++i) {
                int e = std::clamp<int>(residualRaw[i], -128, 127);
                uint8_t v;
                if (abs(e) <= 1) v = 0;
                else if (abs(e) > 1 && abs(e) <= 127) v = 127;
                else v = 255;
                out.write(reinterpret_cast<const char*>(&v), 1);
            }
        };

        for (int i = 1; i < 365; ++i)
        {
            if (ctxStats[i].count == 0) continue;
        
            double A = double(ctxStats[i].sumA) / ctxStats[i].count;
            double B = double(ctxStats[i].sumB) / ctxStats[i].count;
            double C = double(ctxStats[i].sumC) / ctxStats[i].count;
            double D = double(ctxStats[i].sumD) / ctxStats[i].count;
            double X = double(ctxStats[i].sumX) / ctxStats[i].count;

            // Unterordner definieren
            std::string miniDir = m_outputDir + "/minibilder";
            std::filesystem::create_directories(miniDir);  
            std::string filename = miniDir + "/ctx3x2_" + std::to_string(i) + ".pgm";

            std::ofstream out(filename, std::ios::binary);
        
            out << "P5\n3 2\n255\n";

            uint8_t img[6];

            // Zeile 1
            img[0] = (uint8_t)std::clamp(C, 0.0, 255.0);  // C
            img[1] = (uint8_t)std::clamp(B, 0.0, 255.0);  // B
            img[2] = (uint8_t)std::clamp(D, 0.0, 255.0);  // D

            // Zeile 2
            img[3] = (uint8_t)std::clamp(A, 0.0, 255.0);  // A
            img[4] = (uint8_t)std::clamp(X, 0.0, 255.0);  // X (leer)
            img[5] = 0;                                   // leer

            out.write(reinterpret_cast<char*>(img), 6);
        }


        std::string residuals_raw = m_outputDir + "/residuals_raw_" + std::to_string(gsize) + ".pgm";
        std::string residuals_3color = m_outputDir + "/residuals_3color_" + std::to_string(gsize) + ".pgm";

        saveRawPGM(residuals_raw);
        saveRawPGM3color(residuals_3color);

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
    int debug = 1;
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

            int ctxIdx = (q1_norm) * 25 + (q2_norm) * 5 + (q3_norm);
            ctxIdx = std::clamp(ctxIdx, 0, 364);
            if (ctxIdx == 0) ctxIdx = 1;

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

            // STEP 7: Residuum direkt vom Decoder holen
            int e_fold = edec.decodeSample(ctxIdx);
                      
            // Rekonstruktion mit modulo-Addition
            int reconstructed = (pred_corr + e_fold) & 0xFF;
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
        ctx.N += 1;
        
        // Reset-Schwellenwert (Lmax-Äquivalent)
        const int N_RESET = 640;
        if (ctx.N >= N_RESET) {
            ctx.N >>= 1;
            ctx.A >>= 1;
            ctx.B >>= 1;
        }
        
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
    inline int quantize(int gradient) const
    {
        if (gradient == 0) return 0;
        int sign = (gradient < 0) ? -1 : 1;
        int a = std::abs(gradient);
        if (a <= 2)   return sign * 1;
        if (a <= 6)   return sign * 2;
        if (a <= 20)  return sign * 3;
        return sign * 4;
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
    std::ofstream stream(outname, std::ios::binary);
    OBitstream bs(stream);

    bs.addFixed<unsigned>(img.getWidth(), 16);
    bs.addFixed<unsigned>(img.getHeight(), 16);

    EntropyEncoder eenc(bs, groupSize);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData(), false, outputDir);
    pred.subtractPrediction(eenc, groupSize);
    bs.byteAlign();
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
    img.write(outname);
}


namespace fs = std::filesystem;


int main(int argc, char** argv)
{
    cmdPars pars;

    if (!readCmdLine(argc, argv, pars))
        return 1;

    fs::create_directories(pars.outdir);

    std::string fullOutputPath = pars.outdir + "/" + pars.outname;

    unsigned groupSize = 1;

    try {
        if (pars.decode) {
            std::cout << "Dekodiere: " << pars.inname << " -> " << fullOutputPath << std::endl;
            decode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        } else {
            std::cout << "Kodiere: " << pars.inname << " -> " << fullOutputPath << std::endl;
            encode(pars.inname, fullOutputPath, groupSize, pars.outdir);
        }
        std::cout << "✅ Erfolgreich!" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Fehler: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
