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
};

bool readCmdLine( int argc, char** argv, cmdPars& pars )
{
  std::stringstream err;
  int               arg = 0;
  bool              ok  = true;
  ok = ok && ( argc > ++arg );
  if( ok ) {
    if( std::string( argv[ arg ] ) == std::string( "-e" ) ) {
      pars.decode = false;
    } else if( std::string( argv[ arg ] ) == std::string( "-d" ) ) {
      pars.decode = true;
    } else {
      err << "ERROR: First parameter must be '-e' or '-d'." << std::endl;
      ok = false;
    }
  }
  ok = ok && ( argc > ++arg );
  if( ok )  {
    pars.inname = std::string( argv[ arg ] );
    std::ifstream intest( pars.inname );
    if( !intest.good() ) { 
      err << "ERROR: Cannot open input file \"" << pars.inname << "\"." << std::endl;
      ok = false; 
    } 
  }
  ok = ok && ( argc > ++arg );
  if( ok )  {
    pars.outname = std::string( argv[ arg ] );
    std::ofstream outtest( pars.outname );
    if( !outtest.good() ) { 
      err << "ERROR: Cannot open output file \"" << pars.outname << "\"." << std::endl;
      ok = false; 
    } 
  }
  if( !ok) {
    std::string pname = std::string( argv[0] );
    std::cerr << err.str() << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << "  decoding: " << pname << " -d binFile recFile" << std::endl;
    std::cerr << "  encoding: " << pname << " -e orgFile binFile" << std::endl;
    std::cerr << std::endl;
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
  EntropyCoderBase() 
  {
    m_pmfAbs  = std::vector<uint8_t>(N, uint8_t{0});  
  }

protected:
  static const unsigned   N = 3;    // number of probability models for unary binarization of absolute values
  std::vector<uint8_t>    m_pmfAbs; // probability models for coding absolute values
};


class EntropyEncoder : protected EntropyCoderBase
{
public:
  EntropyEncoder( OBitstream& bs ) : EntropyCoderBase(), aenc(bs)
  {
    aenc.start();
  }
  void encodeSample  ( PGMImage::Sample s );
  void finish        ()    
  { 
    aenc.finish(); 
  }
private:
  ArithmeticEncoder aenc;
};

class EntropyDecoder : protected EntropyCoderBase
{
public:
  EntropyDecoder( IBitstream& bs ) : EntropyCoderBase(), adec(bs)
  {
    adec.start();
  }
  PGMImage::Sample  decodeSample  ();
private:
  ArithmeticDecoder adec;
};

void EntropyEncoder::encodeSample( PGMImage::Sample s )
{
  unsigned absValue = unsigned( s < 0 ? -s : s );
  unsigned rem      = absValue;
  unsigned binIdx   = 0;
  while( rem-- )
  {
    aenc.encBin ( m_pmfAbs[ std::min<unsigned>( N-1, binIdx++ ) ], 1 );
  }
  aenc.encBin   ( m_pmfAbs[ std::min<unsigned>( N-1, binIdx++ ) ], 0 );
  if( absValue )
    aenc.encBit ( s<0 );
}

PGMImage::Sample EntropyDecoder::decodeSample()
{
  PGMImage::Sample  s      = 0;
  unsigned          binIdx = 0;
  while( adec.decBin( m_pmfAbs[ std::min<unsigned>( N-1, binIdx++ ) ] ) )
  {
    s++;
  }
  if( s && adec.decBit() )
    s = -s;
  return s;
}






//======================================================
//
//   LOCO-I PREDICTOR (JPEG-LS) - KORRIGIERTE IMPLEMENTIERUNG
//
//======================================================
//======================================================
//
//   LOCO-I PREDICTOR MIT DIAGNOSE
//
//======================================================
class Prediction
{
public:
    Prediction(int width, int height, std::vector<PGMImage::Sample>& img)
        : m_width(width), m_height(height), m_data(img), m_org(m_data)
    {
        m_ctxs.resize(365);
        
        for (auto& ctx : m_ctxs) {
            ctx.N = 1;
            ctx.A = std::max(2, (width * height + 32) / 64);
            ctx.B = 0;
            ctx.C = 0;
        }
        
        // Diagnose-Zähler
        bias_used = 0;
        total_pixels = 0;
        
    }

    void subtractPrediction()
    {
        const PGMImage::Sample* orgData = m_org.data();
        PGMImage::Sample* data = m_data.data();
        
        int zero_residuals = 0;
        int small_residuals = 0;
        int dbg_count = 0;
        constexpr int DBG_PRINT_PIXELS = 16; // Weniger für Übersichtlichkeit
        
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
                
                // STEP 1: Lokale Gradienten berechnen
                int g1 = d - b;  // db
                int g2 = b - c;  // bc
                int g3 = c - a;  // ca
                
                // STEP 2: Check für Run-Mode (alle Gradienten = 0?)
                // Für diese Implementierung: Wir behandeln alle Pixel im "regular mode"
                // (Run-Mode Implementation wäre optional)
                
                // STEP 3 + 4: Gradienten quantisieren und Kontext bestimmen
                int q1 = quantize(g1);
                int q2 = quantize(g2);
                int q3 = quantize(g3);
                
                // Vorzeichen bestimmen (erste nicht-Null Komponente)
                int sign = 1;
                if (q1 < 0) sign = -1;
                else if (q1 == 0 && q2 < 0) sign = -1;
                else if (q1 == 0 && q2 == 0 && q3 < 0) sign = -1;
                
                // Gradienten normalisieren für Kontext-Index
                int q1_norm = q1, q2_norm = q2, q3_norm = q3;
                if (sign < 0) {
                    q1_norm = -q1;
                    q2_norm = -q2;
                    q3_norm = -q3;
                }
                
                // Kontext-Index berechnen (1..364)
                int ctxIdx = (q1_norm) * 25 + (q2_norm) * 5 + (q3_norm);
                ctxIdx = std::clamp(ctxIdx, 0, 364);
                if (ctxIdx == 0) ctxIdx = 1;

                size_t safeIdx = static_cast<size_t>(ctxIdx);
                Context& ctx = m_ctxs[safeIdx];
                
                // STEP 5: Feste MED-Prädiktion berechnen
                int pred_med;
                if (c >= std::max(a, b))
                    pred_med = std::min(a, b);
                else if (c <= std::min(a, b))
                    pred_med = std::max(a, b);
                else
                    pred_med = a + b - c;
                
                // STEP 6: Adaptive Korrektur anwenden
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

                int error = actual - pred_corr;
                
                // Negatives Vorzeichen anwenden falls nötig
                if (sign < 0)
                    error = -error;
                
                // Modulo-Reduktion auf [-128, 127]
                int error_mod = error;
                while (error_mod > 127) error_mod -= 256;
                while (error_mod < -128) error_mod += 256;
                
                // STEP 8-10: Golomb Parameter k (nicht implementiert, verwende nur error_mod)
                // Hier würde Golomb-Codierung folgen
                
                // Speichern des Residuums
                data[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(error_mod, -128, 127));
                
                // Residuen-Statistiken
                if (error_mod == 0) zero_residuals++;
                if (std::abs(error_mod) <= 1) small_residuals++;
                
                // Debug-Ausgabe
                if (dbg_count < DBG_PRINT_PIXELS) {
                    std::cout << "[DBG SUB] idx=" << (y*m_width + x)
                              << " (" << x << "," << y << ")"
                              << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "]"
                              << " grads={" << g1 << "," << g2 << "," << g3 << "}"
                              << " q={" << q1 << "," << q2 << "," << q3 << "}"
                              << " sign=" << sign
                              << " pred_med=" << pred_med
                              << " pred_corr=" << pred_corr
                              << " actual=" << actual
                              << " error=" << error
                              << " error_mod=" << error_mod
                              << " ctx=" << ctxIdx
                              << " ctx(A,B,N,C)=(" << ctx.A << "," << ctx.B << "," << ctx.N << "," << ctx.C << ")"
                              << std::endl;
                    ++dbg_count;
                }
                if (dbg_count < DBG_PRINT_PIXELS) std::cout << "ctxIdx: " << ctxIdx << "error" << error << std::endl;
                // STEP 11-12: Context-Zähler aktualisieren
                updateContextLocoI(ctxIdx, error);
            }
        }

        std::cout << "\n=== SUBTRACTPREDICTION SUMMARY ===\n"
                  << "Total pixels: " << total_pixels << "\n"
                  << "Bias corrections: " << bias_used << " (" << (100.0*bias_used/total_pixels) << "%)\n"
                  << "Zero residuals: " << zero_residuals << " (" << (100.0*zero_residuals/total_pixels) << "%)\n"
                  << "Small residuals (|e|<=1): " << small_residuals << " (" << (100.0*small_residuals/total_pixels) << "%)\n"
                  << std::endl;
    }

void addPrediction()
{
    //PGMImage::Sample* recData = m_data.data(); // Verwende die Residuumsdaten als Basis
    const PGMImage::Sample* residualData = m_data.data(); // Residuen-Daten
    
    int dbg_count = 0;
    constexpr int DBG_PRINT_PIXELS = 16;

    // Temporären Puffer für rekonstruierte Daten erstellen
    std::vector<PGMImage::Sample> tempData(m_data.size());
    
    for (int y = 0; y < m_height; ++y)
    {
        for (int x = 0; x < m_width; ++x)
        {
            // WICHTIG: Für Prädiktion die bereits rekonstruierten Nachbarn verwenden
            auto get = [&](int xx, int yy) -> int {
                if (xx < 0 || yy < 0) return 0;
                if (xx >= m_width || yy >= m_height) return 0;
                // Verwende tempData für bereits rekonstruierte Pixel
                return int(tempData[yy * m_width + xx]);
            };
            
            int a = get(x - 1, y);     // left
            int b = get(x,     y - 1); // up
            int c = get(x - 1, y - 1); // up-left
            int d = get(x + 1, y - 1); // up-right
            
            // STEP 1: Lokale Gradienten berechnen
            int g1 = d - b;
            int g2 = b - c;
            int g3 = c - a;
            
            // STEP 3 + 4: Gradienten quantisieren und Kontext bestimmen
            int q1 = quantize(g1);
            int q2 = quantize(g2);
            int q3 = quantize(g3);
            
            // Vorzeichen bestimmen
            int sign = 1;
            if (q1 < 0) sign = -1;
            else if (q1 == 0 && q2 < 0) sign = -1;
            else if (q1 == 0 && q2 == 0 && q3 < 0) sign = -1;
            
            // Gradienten normalisieren für Kontext-Index
            int q1_norm = q1, q2_norm = q2, q3_norm = q3;
            if (sign < 0) {
                q1_norm = -q1;
                q2_norm = -q2;
                q3_norm = -q3;
            }
            
            // Kontext-Index berechnen
            int ctxIdx = (q1_norm) * 25 + (q2_norm) * 5 + (q3_norm);
            ctxIdx = std::clamp(ctxIdx, 0, 364);
            if (ctxIdx == 0) ctxIdx = 1;

            size_t safeIdx = static_cast<size_t>(ctxIdx);
            Context& ctx = m_ctxs[safeIdx];
            
            // STEP 5: Feste MED-Prädiktion berechnen
            int pred_med;
            if (c >= std::max(a, b))
                pred_med = std::min(a, b);
            else if (c <= std::min(a, b))
                pred_med = std::max(a, b);
            else
                pred_med = a + b - c;
            
            // STEP 6: Adaptive Korrektur anwenden
            int pred_corr = pred_med;
            if (sign > 0)
                pred_corr += ctx.C;
            else
                pred_corr -= ctx.C;
            pred_corr = std::clamp(pred_corr, 0, 255);
            
            // STEP 7: Residuum aus den decodierten Daten laden
            int error_stored = int(residualData[y * m_width + x]);
            
            // Negatives Vorzeichen anwenden falls nötig
            int error = error_stored;
            if (sign < 0)
                error = -error;
            
            // Rekonstruktion: Prädiktion + Residuum
            int reconstructed = pred_corr + error;
            
            // Auf gültigen Bereich beschränken
            reconstructed = std::clamp(reconstructed, 0, 255);
            
            // In temporären Puffer speichern
            tempData[y * m_width + x] = static_cast<PGMImage::Sample>(reconstructed);
            
            // Debug-Ausgabe
            if (dbg_count < DBG_PRINT_PIXELS) {
                std::cout << "[DBG ADD] idx=" << (y*m_width + x)
                          << " (" << x << "," << y << ")"
                          << " [a,b,c,d]=[" << a << "," << b << "," << c << "," << d << "]"
                          << " grads={" << g1 << "," << g2 << "," << g3 << "}"
                          << " q={" << q1 << "," << q2 << "," << q3 << "}"
                          << " sign=" << sign
                          << " pred_med=" << pred_med
                          << " pred_corr=" << pred_corr
                          << " error_stored=" << error_stored
                          << " error=" << error
                          << " reconstructed=" << reconstructed
                          << " ctx=" << ctxIdx
                          << " ctx(A,B,N,C)=(" << ctx.A << "," << ctx.B << "," << ctx.N << "," << ctx.C << ")"
                          << std::endl;
                ++dbg_count;
            }
            if (dbg_count < DBG_PRINT_PIXELS) std::cout << "ctxIdx: " << ctxIdx << "error" << error << std::endl;
            // STEP 11-12: Context-Zähler aktualisieren
            updateContextLocoI(ctxIdx, error_stored);
        }
    }
    
    // Rekonstruierte Daten zurück in m_data kopieren
    m_data = tempData;

    std::cout << "\n=== ADDPREDICTION SUMMARY ===\n"
              << "Reconstruction complete\n" << std::endl;
}

private:
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
    struct Context
    {
        int A = 0;  // Sum of absolute errors
        int B = 0;  // Sum of errors (for bias)
        int N = 1;  // Context count
        int C = 0;  // Bias correction value
    };
private:
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
void encode(const std::string& inname, const std::string& outname)
{
    PGMImage img;
    img.read(inname);

    std::ofstream stream(outname, std::ios::out | std::ios::binary);
    OBitstream bs(stream);

    bs.addFixed<unsigned>(img.getWidth(), 16);
    bs.addFixed<unsigned>(img.getHeight(), 16);
    
    EntropyEncoder eenc(bs);
    Prediction pred(img.getWidth(), img.getHeight(), img.getData());

    // Prädiktion anwenden
    pred.subtractPrediction();
    
    // Residuen codieren
    PGMImage::Sample* data = img.getData().data();
    for (int k = 0; k < img.getSize(); k++) {
        eenc.encodeSample(data[k]);
    }
    eenc.finish();

    bs.byteAlign();
}


void decode( const std::string& inname, const std::string& outname )
{
  // open bitstream and read header
  std::ifstream stream( inname, std::ios::in|std::ios::binary );
  IBitstream    bs( stream );
  assert( stream.good() );
  int width  = bs.getFixed<unsigned>( 16 );
  int height = bs.getFixed<unsigned>( 16 );

  // create image and decode color planes block by block
  PGMImage        img( width, height );
  EntropyDecoder  edec( bs );
  Prediction      pred( img.getWidth(), img.getHeight(), img.getData() );

  // decode prediction error signal
  PGMImage::Sample* data = img.getData().data();
  for( int k = 0; k < img.getSize(); k++ )
    data[k] = edec.decodeSample();

  // apply prediction
  pred.addPrediction();

  // output reconstructed image
  img.write( outname );
}

namespace fs = std::filesystem;

int main() {
    std::string inputFolder  = "kodak-pgm";
    std::string outputFolder = "output";

    // Aktuelles Arbeitsverzeichnis ausgeben
    std::cout << "Current path: " << fs::current_path() << std::endl;

    // Prüfen, ob Input-Ordner existiert
    if (!fs::exists(inputFolder) || !fs::is_directory(inputFolder)) {
        std::cerr << "Fehler: Input-Ordner \"" << inputFolder << "\" existiert nicht!" << std::endl;
        return 1;
    }

    // Output-Ordner erstellen, falls nicht vorhanden
    fs::create_directories(outputFolder);

    // Bitraten-Datei öffnen
    std::ofstream report(outputFolder + "/bitrates.txt", std::ios::out);
    if (!report.is_open()) {
        std::cerr << "Fehler: Konnte bitrates.txt nicht erstellen!" << std::endl;
        return 1;
    }

    int fileCount = 0;

    // Alle PGM-Dateien im Input-Ordner verarbeiten
    for (const auto& entry : fs::directory_iterator(inputFolder)) {
        if (entry.path().extension() == ".pgm") {
            fileCount++;
            std::string inputFile = entry.path().string();
            std::string baseName  = entry.path().stem().string();

            std::string encodedFile = outputFolder + "/" + baseName + ".bin";
            std::string decodedFile = outputFolder + "/" + baseName + "_decoded.pgm";

            std::cout << "=== " << baseName << " ===" << std::endl;

            // Encoding
            std::cout << "Encoding " << inputFile << " -> " << encodedFile << std::endl;
            encode(inputFile, encodedFile);
            std::cout << "Tests ";
            // Bitrate berechnen
            uintmax_t fileSize = fs::file_size(encodedFile);
            PGMImage img;
            img.read(inputFile);
            double bitrate = (fileSize * 8.0) / img.getSize();

            std::cout << "Bitrate: " << std::fixed << std::setprecision(3) 
                      << bitrate << " bpp" << std::endl;
            report << baseName << ": " << std::fixed << std::setprecision(3) 
                   << bitrate << " bpp" << std::endl;

            // Decoding
            std::cout << "Decoding " << encodedFile << " -> " << decodedFile << std::endl;
            decode(encodedFile, decodedFile);

            std::cout << std::endl;
        }
    }

    if (fileCount == 0) {
        std::cout << "Keine .pgm-Dateien im Ordner \"" << inputFolder << "\" gefunden!" << std::endl;
    } else {
        std::cout << "Fertig! " << fileCount << " Dateien verarbeitet." << std::endl;
    }

    return 0;
}