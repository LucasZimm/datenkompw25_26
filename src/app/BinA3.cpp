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
            for (auto &p : ctx) p = 1;
        }

        // Initialisiere Sign-PMF auf kleinen Startwert
        for (auto &p : m_pmfSignCtx) p = 1;

        // Initialisiere Significance-PMF auf kleinen Startwert
        for (auto &p : m_pmfSigCtx) p = 1;
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
    void encodeSample(PGMImage::Sample s, unsigned ctxIdx, unsigned signIdx);

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

    PGMImage::Sample decodeSample(unsigned ctxIdx, unsigned signIdx);

private:
    ArithmeticDecoder adec;
};


// -----------------------------
// Encoder
// -----------------------------
void EntropyEncoder::encodeSample(PGMImage::Sample s,
                                  unsigned ctxIdx, unsigned signIdx)
                                  
{
    unsigned absValue = unsigned(s < 0 ? -s : s);
    unsigned mIdx = mappedCtx(ctxIdx);

    // -----------------------------
    // 1) significance
    // -----------------------------
    bool nonZero = (absValue != 0);
    aenc.encBin(m_pmfSigCtx[mIdx], nonZero);

    if (!nonZero)
        return;

    // -----------------------------
    // 2) sign (kontextabhängig!)
    // -----------------------------
    unsigned signCtx = signIdx;
    aenc.encBin(m_pmfSignCtx[signCtx], s < 0);

    // -----------------------------
    // 3) magnitude (dein Code!)
    // -----------------------------
    unsigned rem = absValue;
    unsigned binIdx = 0;

    while (rem--)
        aenc.encBin(
            m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)], 1);

    aenc.encBin(
        m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)], 0);
}



// -----------------------------
// Decoder
// -----------------------------
PGMImage::Sample EntropyDecoder::decodeSample(unsigned ctxIdx, unsigned signIdx)
                                              
{
    unsigned mIdx = mappedCtx(ctxIdx);

    // -----------------------------
    // 1) significance
    // -----------------------------
    if (!adec.decBin(m_pmfSigCtx[mIdx]))
        return 0;

    // -----------------------------
    // 2) sign
    // -----------------------------
    unsigned signCtx = signIdx;
    bool sign = adec.decBin(m_pmfSignCtx[signCtx]);

    // -----------------------------
    // 3) magnitude
    // -----------------------------
    unsigned mag = 0;
    unsigned binIdx = 0;

    while (adec.decBin(
        m_pmfAbsCtx[mIdx][std::min<unsigned>(N - 1, binIdx++)]))
    {
        mag++;
    }

    return sign ? -PGMImage::Sample(mag)
                :  PGMImage::Sample(mag);
}


//======================================================
//
//   S I G N   I N D E X   B E R E C H N U N G E N
//
//======================================================

/**
 * Lineare signIdx-Berechnung (Binary: 0-1)
 * Input: Gradienten (q1, q2, q3) und Nachbarfehler (a_err, b_err)
 * Output: 0 oder 1 für Sign-Kontext
 */
inline unsigned computeSignIdx_Linear_Binary(int q1, int q2, int q3, 
                                              int a_err, int b_err)
{
    // Gradient-Vorzeichen mit Priorität
    int gradient_sign = (std::abs(q1) > 0) ? (q1 < 0 ? -1 : 1) :
                        (std::abs(q2) > 0) ? (q2 < 0 ? -1 : 1) : 
                                             (q3 < 0 ? -1 : 1);
    
    // Nachbarfehler-Vorzeichen
    int residual_sign = (a_err < 0 || b_err < 0) ? -1 : 1;
    
    // Lineare Kombination
    return ((gradient_sign + residual_sign) > 0) ? 1 : 0;
}


/**
 * Lineare signIdx-Berechnung (Ternär: 0-2)
 * Input: Gradienten (q1, q2, q3) und Nachbarfehler (a_err, b_err)
 * Output: 0, 1, oder 2 für 3 verschiedene Sign-Kontexte
 * 
 * Differenziert auch zwischen starken und schwachen Gradienten
 */
inline unsigned computeSignIdx_Linear_Ternary(int q1, int q2, int q3, 
                                               int a_err, int b_err)
{
    // Gradienten-Magnitude
    int max_grad = std::max({std::abs(q1), std::abs(q2), std::abs(q3)});
    
    // Gradienten-Vorzeichen
    int gradient_sign = (std::abs(q1) == max_grad) ? (q1 < 0 ? -1 : 1) :
                        (std::abs(q2) == max_grad) ? (q2 < 0 ? -1 : 1) : 
                                                     (q3 < 0 ? -1 : 1);
    
    // Nachbarfehler-Vorzeichen
    int residual_sign = (a_err < 0 || b_err < 0) ? -1 : 1;
    
    // Ternäre Klassifikation: schwach/neutral/stark
    if (max_grad > 10 && gradient_sign == residual_sign) {
        return 2;  // Starke Übereinstimmung
    } else if ((gradient_sign + residual_sign) > 0) {
        return 1;  // Schwache Übereinstimmung
    } else {
        return 0;  // Gegensätzlich
    }
}


/**
 * Nicht-lineare signIdx-Berechnung (4-bit: 0-15)
 * Input: Gradienten (q1, q2, q3) und Nachbarfehler (a_err, b_err)
 * Output: 0-15 für 16 verschiedene Sign-Kontexte
 * 
 * Nutzt 4 Bits: 2 für Gradient-Info, 2 für Fehler-Info
 */
inline unsigned computeSignIdx_Nonlinear_4bit(int q1, int q2, int q3,
                                               int a_err, int b_err)
{
    unsigned ctx = 0;
    
    // Bit 0-1: Gradient-Kategorie (2 Bits)
    int max_grad = std::max({std::abs(q1), std::abs(q2), std::abs(q3)});
    int gradient_sign = (std::abs(q1) == max_grad) ? (q1 < 0 ? -1 : 1) :
                        (std::abs(q2) == max_grad) ? (q2 < 0 ? -1 : 1) : 
                                                     (q3 < 0 ? -1 : 1);
    
    if (gradient_sign > 0) {
        ctx |= (max_grad > 10) ? 0b11 : 0b01;  // Stark positiv oder schwach positiv
    } else {
        ctx |= (max_grad > 10) ? 0b10 : 0b00;  // Stark negativ oder schwach negativ
    }
    
    // Bit 2-3: Fehler-Kategorie (2 Bits)
    int error_sum = a_err + b_err;
    //int error_diff = std::abs(a_err) - std::abs(b_err);
    
    unsigned error_ctx = 0;
    if (error_sum > 0) {
        error_ctx = (std::abs(error_sum) > 10) ? 0b11 : 0b01;
    } else {
        error_ctx = (std::abs(error_sum) > 10) ? 0b10 : 0b00;
    }
    
    ctx |= (error_ctx << 2);
    
    return ctx;
}


/**
 * Nicht-lineare signIdx-Berechnung (3-bit: 0-7)
 * Input: Gradienten (q1, q2, q3) und Nachbarfehler (a_err, b_err)
 * Output: 0-7 für 8 verschiedene Sign-Kontexte
 * 
 * Bit-Struktur: [error_sign:1 | error_mag:1 | gradient_sign:1]
 */
inline unsigned computeSignIdx_Nonlinear_3bit(int q1, int q2, int q3,
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


/**
 * Nicht-lineare signIdx-Berechnung (Maximum Entropy: 0-7)
 * Input: Gradienten mit Magnitude und Nachbarfehler
 * Output: 0-7 für optimale Kontextverteilung
 * 
 * Kodiert: [q1_sign:1 | q1_strong:1 | error_sign:1]
 */
inline unsigned computeSignIdx_Nonlinear_MaxEntropy(int q1, int q2, int q3,
                                                     int a_err, int b_err)
{
    unsigned ctx = 0;
    if (q2 < q3) {
        std::swap(q2, q3);
    }
    // Verwende q1 als primären Gradienten (meist aussagekräftig)
    ctx = (q1 < 0) ? 0b001 : 0b000;
    
    // Magnitude-Kategorie von q1
    if (std::abs(q1) > 15) {
        ctx |= 0b010;
    }
    
    // Kombiniertes Nachbar-Fehler-Vorzeichen
    int combined_error = a_err + b_err;
    if (combined_error < 0) {
        ctx |= 0b100;
    }
    
    return ctx;
}


/**
 * Verbesserte signIdx-Berechnung (alle Gradienten: 0-7)
 * Input: Alle 3 Gradienten (q1, q2, q3) und Nachbarfehler
 * Output: 0-7 für 8 verschiedene Sign-Kontexte
 * 
 * Nutzt ALLE Gradienten für bessere Vorhersage:
 * - Zählt wie viele Gradienten negativ sind (Voting)
 * - Berücksichtigt auch Magnitude-Unterschiede
 */
inline unsigned computeSignIdx_AllGradients_Voting(int q1, int q2, int q3,
                                                    int a_err, int b_err)
{
    unsigned ctx = 0;
    
    // Bit 0: Mehrheits-Vorzeichen der Gradienten (Voting)
    int negative_count = 0;
    if (q1 < 0) negative_count++;
    if (q2 < 0) negative_count++;
    if (q3 < 0) negative_count++;
    
    // Wenn mehr Gradienten negativ sind, setze Bit 0
    if (negative_count >= 2) {
        ctx |= 0b001;
    }
    
    // Bit 1: Magnitude-Konsistenz (sind die Gradienten ähnlich groß?)
    int max_grad = std::max({std::abs(q1), std::abs(q2), std::abs(q3)});
    int min_grad = std::min({std::abs(q1), std::abs(q2), std::abs(q3)});
    int mag_diff = max_grad - min_grad;
    
    if (mag_diff > 5) {
        ctx |= 0b010;  // Große Unterschiede
    }
    
    // Bit 2: Fehler-Vorzeichen
    int combined_error = a_err + b_err;
    if (combined_error < 0) {
        ctx |= 0b100;
    }
    
    return ctx;
}


/**
 * Gewichtete signIdx-Berechnung (alle Gradienten: 0-15)
 * Input: Alle 3 Gradienten (q1, q2, q3) und Nachbarfehler
 * Output: 0-15 für 16 verschiedene Sign-Kontexte (4-bit)
 * 
 * Nutzt gewichtete Summe aller 3 Gradienten:
 * - q1 bekommt Gewicht 4 (wichtigster Gradient)
 * - q2 bekommt Gewicht 2
 * - q3 bekommt Gewicht 1
 */
inline unsigned computeSignIdx_AllGradients_Weighted(int q1, int q2, int q3,
                                                      int a_err, int b_err)
{
    unsigned ctx = 0;
    
    // Gewichtete Summe aller Gradienten
    int gradient_sum = q1 * 4 + q2 * 2 + q3 * 1;
    
    // Bit 0-1: Gradient-Vorzeichen und Magnitude (2 Bits)
    if (gradient_sum < -20) {
        ctx |= 0b00;  // Stark negativ
    } else if (gradient_sum < 0) {
        ctx |= 0b01;  // Schwach negativ
    } else if (gradient_sum < 20) {
        ctx |= 0b10;  // Schwach positiv
    } else {
        ctx |= 0b11;  // Stark positiv
    }
    
    // Bit 2-3: Fehler-Kategorie (2 Bits)
    int error_sum = a_err + b_err;
    unsigned error_ctx = 0;
    
    if (error_sum < -10) {
        error_ctx = 0b00;
    } else if (error_sum < 0) {
        error_ctx = 0b01;
    } else if (error_sum < 10) {
        error_ctx = 0b10;
    } else {
        error_ctx = 0b11;
    }
    
    ctx |= (error_ctx << 2);
    
    return ctx;
}


/**
 * Kombinierte signIdx-Berechnung (alle Gradienten + Fehler: 0-7)
 * Input: Alle 3 Gradienten (q1, q2, q3) und Nachbarfehler (a_err, b_err)
 * Output: 0-7 für 8 verschiedene Sign-Kontexte
 * 
 * Kodiert: [max_grad_sign:1 | consistency:1 | error_alignment:1]
 * 
 * - max_grad_sign: Vorzeichen des dominantesten Gradienten
 * - consistency: Ob alle Gradienten ähnliche Vorzeichen haben
 * - error_alignment: Ob Fehler mit dominantem Gradienten aligned ist
 */
inline unsigned computeSignIdx_AllGradients_Combined(int q1, int q2, int q3,
                                                      int a_err, int b_err)
{
    unsigned ctx = 0;
    
    // Finde dominantesten Gradienten
    int max_grad = std::max({std::abs(q1), std::abs(q2), std::abs(q3)});
    int dominant_sign = 0;
    
    if (std::abs(q1) == max_grad) {
        dominant_sign = (q1 < 0) ? -1 : 1;
    } else if (std::abs(q2) == max_grad) {
        dominant_sign = (q2 < 0) ? -1 : 1;
    } else {
        dominant_sign = (q3 < 0) ? -1 : 1;
    }
    
    // Bit 0: Dominanter Gradienten-Vorzeichen
    if (dominant_sign < 0) ctx |= 0b001;
    
    // Bit 1: Konsistenz - stimmen mindestens 2 Gradienten überein?
    int matching = 0;
    if ((q1 < 0 && q2 < 0) || (q1 >= 0 && q2 >= 0)) matching++;
    if ((q2 < 0 && q3 < 0) || (q2 >= 0 && q3 >= 0)) matching++;
    if ((q1 < 0 && q3 < 0) || (q1 >= 0 && q3 >= 0)) matching++;
    
    if (matching >= 2) {
        ctx |= 0b010;  // Konsistent
    }
    
    // Bit 2: Error-Alignment mit dominantem Gradienten
    int combined_error = a_err + b_err;
    int error_sign = (combined_error < 0) ? -1 : 1;
    
    if (error_sign == dominant_sign) {
        ctx |= 0b100;  // Aligned mit Gradienten
    }
    
    return ctx;
}


/**
 * Weighted Non-linear signIdx-Berechnung (0-7)
 * Input: Gradienten (q1, q2, q3) und Nachbarfehler (a_err, b_err)
 * Output: 0-7 für 8 verschiedene Sign-Kontexte
 * 
 * Gewichtet q1 stärker als q2 und q3
 */
inline unsigned computeSignIdx_WeightedNonlinear(int q1, int q2, int q3, 
                                                  int a_err, int b_err)
{
    // Gewichtete Gradienten-Summe
    int gradient_weighted = q1 * 4 + q2 * 2 + q3 * 1;
    
    // Fehler-Summe
    int error_sum = a_err + b_err;
    
    // Klassifikation
    unsigned ctx = 0;
    
    // Bits 0-1: Gradient-Klasse
    if (gradient_weighted > 20) {
        ctx |= 0b11;
    } else if (gradient_weighted > 0) {
        ctx |= 0b01;
    } else if (gradient_weighted < -20) {
        ctx |= 0b10;
    }
    
    // Bits 2: Error-Vorzeichen
    if (error_sum < 0) {
        ctx |= 0b100;
    }
    
    return ctx;
}


/**
 * Magnitude-Dependent signIdx-Berechnung (0-7)
 * Input: Gradienten mit Magnitude und Nachbarfehler
 * Output: 0-7 für 8 verschiedene Sign-Kontexte
 * 
 * Codiert Magnitude-Information direkt in den Kontext
 */
inline unsigned computeSignIdx_MagnitudeDependent(int q1, int q2, int q3,
                                                   int a_err, int b_err)
{
    // Dominanter Gradient
    int max_grad = std::max({std::abs(q1), std::abs(q2), std::abs(q3)});
    
    int gradient_sign = (std::abs(q1) == max_grad) ? (q1 < 0 ? -1 : 1) :
                        (std::abs(q2) == max_grad) ? (q2 < 0 ? -1 : 1) : 
                                                     (q3 < 0 ? -1 : 1);
    
    // Magnitude-Kategorien: 0-2 (3 Werte)
    unsigned mag_cat = 0;
    if (max_grad > 20) mag_cat = 2;
    else if (max_grad > 5) mag_cat = 1;
    
    // Fehler-Vorzeichen
    int residual_sign = (a_err < 0 || b_err < 0) ? -1 : 1;
    
    // Kontext zusammensetzen: [magnitude:2 | error_sign:1 | gradient_sign:1]
    unsigned ctx = mag_cat << 2;
    if (residual_sign < 0) ctx |= 0b010;
    if (gradient_sign < 0) ctx |= 0b001;
    
    return ctx;
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
        // Kontext-Histogramme getrennt nach Error==0 / Error!=0
        std::array<size_t, 256> histoA_ctx_err0{};
        std::array<size_t, 256> histoB_ctx_err0{};
        std::array<size_t, 256> histoC_ctx_err0{};
        std::array<size_t, 256> histoD_ctx_err0{};
        std::array<size_t, 256> histoA_ctx_err1{};
        std::array<size_t, 256> histoB_ctx_err1{};
        std::array<size_t, 256> histoC_ctx_err1{};
        std::array<size_t, 256> histoD_ctx_err1{};
        std::array<size_t, 256> histoA_ctx_errN{};
        std::array<size_t, 256> histoB_ctx_errN{};
        std::array<size_t, 256> histoC_ctx_errN{};
        std::array<size_t, 256> histoD_ctx_errN{};
        // Fehler-Histogramme
        std::array<size_t, 511> histoA_err{};
        std::array<size_t, 511> histoB_err{};
        std::array<size_t, 511> histoC_err{};
        std::array<size_t, 511> histoD_err{};
        std::array<size_t, 511> histoAll_err{};
        // Fehler-Histogramme getrennt nach Error==0 / Error!=0




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

                ctxStats[ctxIdx].sumA += a;
                ctxStats[ctxIdx].sumB += b;
                ctxStats[ctxIdx].sumC += c;
                ctxStats[ctxIdx].sumD += d;
                ctxStats[ctxIdx].sumX += actual;
                ctxStats[ctxIdx].count++;


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

                // sign Index fuer CABAC - Kombination aus Gradienten + Residuen-Vorzeichen
                int sig_func = 3;  // 0=Linear_Binary, 1=Linear_Ternary, 2=Nonlinear_4bit, 3=Nonlinear_3bit, 4=MaxEntropy, 5=WeightedNonlinear, 6=MagnitudeDependent
                unsigned signIdx;
                
                if (sig_func == 0) {
                    signIdx = computeSignIdx_Linear_Binary(q1, q2, q3, a_folderrget, b_folderrget);
                } else if (sig_func == 1) {
                    signIdx = computeSignIdx_Linear_Ternary(q1, q2, q3, a_folderrget, b_folderrget);
                } else if (sig_func == 2) {
                    signIdx = computeSignIdx_Nonlinear_4bit(q1, q2, q3, a_folderrget, b_folderrget);
                } else if (sig_func == 3) {
                    signIdx = computeSignIdx_Nonlinear_3bit(q1, q2, q3, a_folderrget, b_folderrget);
                } else if (sig_func == 4) {
                    signIdx = computeSignIdx_Nonlinear_MaxEntropy(q1, q2, q3, a_folderrget, b_folderrget);
                } else if (sig_func == 5) {
                    signIdx = computeSignIdx_WeightedNonlinear(q1, q2, q3, a_folderrget, b_folderrget);
                } else if (sig_func == 6) {
                    signIdx = computeSignIdx_MagnitudeDependent(q1, q2, q3, a_folderrget, b_folderrget);
                } else {
                    signIdx = 0;  // Default fallback
                }

                // Codieren mit CABAC
                eenc.encodeSample(static_cast<PGMImage::Sample>(e_fold), model, signIdx);

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

                // Zusätzliche Histogramme: A,B,C,D getrennt nach (error == 0) und (error != 0)
                if (error_code == 0) {
                    histoA_ctx_err0[uint8_t(a)]++;
                    histoB_ctx_err0[uint8_t(b)]++;
                    histoC_ctx_err0[uint8_t(c)]++;
                    histoD_ctx_err0[uint8_t(d)]++;
                } else {
                    if (std::abs(error_code) <= 1) {
                        histoA_ctx_err1[uint8_t(a)]++;
                        histoB_ctx_err1[uint8_t(b)]++;
                        histoC_ctx_err1[uint8_t(c)]++;
                        histoD_ctx_err1[uint8_t(d)]++;
                    } else {
                        // für Fehler > 1
                    histoA_ctx_errN[uint8_t(a)]++;
                    histoB_ctx_errN[uint8_t(b)]++;
                    histoC_ctx_errN[uint8_t(c)]++;
                    histoD_ctx_errN[uint8_t(d)]++;
                }
                }

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

        // Kontext-Histogramme getrennt nach Error==0 / Error!=0 speichern
        saveHistogram(histoA_ctx_err0, "histoA_ctx_err0");
        saveHistogram(histoB_ctx_err0, "histoB_ctx_err0");
        saveHistogram(histoC_ctx_err0, "histoC_ctx_err0");
        saveHistogram(histoD_ctx_err0, "histoD_ctx_err0");
        saveHistogram(histoA_ctx_err1, "histoA_ctx_err1");
        saveHistogram(histoB_ctx_err1, "histoB_ctx_err1");
        saveHistogram(histoC_ctx_err1, "histoC_ctx_err1");
        saveHistogram(histoD_ctx_err1, "histoD_ctx_err1");
        saveHistogram(histoA_ctx_errN, "histoA_ctx_errN");
        saveHistogram(histoB_ctx_errN, "histoB_ctx_errN");
        saveHistogram(histoC_ctx_errN, "histoC_ctx_errN");
        saveHistogram(histoD_ctx_errN, "histoD_ctx_errN");
 
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
            // Write a simple P6 (RGB) image where each channel gets the same mapped value
            out << "P6\n" << m_width << " " << m_height << "\n255\n";
            for (int i = 0; i < m_width * m_height; ++i) {
                int e = std::clamp<int>(residualRaw[i], -128, 127);
                uint8_t v;
                if (std::abs(e) <= 1) v = 0;
                else if (std::abs(e) > 1 && std::abs(e) <= 127) v = 127;
                else v = 255;
                // write RGB triple
                out.put(static_cast<char>(v));
                out.put(static_cast<char>(v));
                out.put(static_cast<char>(v));
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
    std::vector<int16_t> data(m_width * m_height);

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

            // sign Index fuer CABAC - Kombination aus Gradienten + Residuen-Vorzeichen
            int sig_func = 3;  // 0=Linear_Binary, 1=Linear_Ternary, 2=Nonlinear_4bit, 3=Nonlinear_3bit, 4=MaxEntropy, 5=WeightedNonlinear, 6=MagnitudeDependent
            unsigned signIdx;
            
            if (sig_func == 0) {
                signIdx = computeSignIdx_Linear_Binary(q1, q2, q3, a_folderrget, b_folderrget);
            } else if (sig_func == 1) {
                signIdx = computeSignIdx_Linear_Ternary(q1, q2, q3, a_folderrget, b_folderrget);
            } else if (sig_func == 2) {
                signIdx = computeSignIdx_Nonlinear_4bit(q1, q2, q3, a_folderrget, b_folderrget);
            } else if (sig_func == 3) {
                signIdx = computeSignIdx_Nonlinear_3bit(q1, q2, q3, a_folderrget, b_folderrget);
            } else if (sig_func == 4) {
                signIdx = computeSignIdx_Nonlinear_MaxEntropy(q1, q2, q3, a_folderrget, b_folderrget);
            } else if (sig_func == 5) {
                signIdx = computeSignIdx_WeightedNonlinear(q1, q2, q3, a_folderrget, b_folderrget);
            } else if (sig_func == 6) {
                signIdx = computeSignIdx_MagnitudeDependent(q1, q2, q3, a_folderrget, b_folderrget);
            } else {
                signIdx = 0;  // Default fallback
            }

            // STEP 7: Residuum direkt vom Decoder holen
            int e_fold = edec.decodeSample(model, signIdx);
                      
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

    int signedContextModel( int a, int b, int c, int d) const
    {
        int sum = (a + b + c + d);
        return (sum >= 0) ? sum : -sum;
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
