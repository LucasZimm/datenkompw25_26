
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
//   LOCO-I PREDICTOR (JPEG-LS) 
//
//======================================================
class Prediction
{
public:
    Prediction(int width, int height, std::vector<PGMImage::Sample>& img)
        : m_width(width), m_height(height), m_data(img), m_org(m_data)
    {
        // Kontexte initialisieren
        m_ctxs.resize(365);
        for (auto& ctx : m_ctxs) {
            ctx.N = 1;
            ctx.A = std::max(2, (width * height + 32) / 64);
            ctx.B = 0;
            ctx.C = 0;
        }


        bias_used = 0;
        total_pixels = 0;
    }

    // --------------------------------------------------
    // Vorhersage subtrahieren (Steps 1-12)
    // --------------------------------------------------
    void subtractPrediction()
    {
        const PGMImage::Sample* orgData = m_org.data();
        PGMImage::Sample* data = m_data.data();

        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                // Step 1: lokale Gradienten berechnen
                int g1, g2, g3;
                std::tie(g1, g2, g3) = computeGradients(orgData, x, y);

                // Step 3: Gradienten quantisieren
                int q1 = quantize(g1);
                int q2 = quantize(g2);
                int q3 = quantize(g3);

                // Step 4: Kontextindex und Sign bestimmen
                int sign;
                int ctxIdx = indexBestimmung(q1, q2, q3, sign);

                // Step 5: Feste MED-Vorhersage
                int fixedPred = getFixedPrediction(orgData, x, y);

                // Step 6: Adaptive Bias-Korrektur
                int correctedPred = MedianCorrection(fixedPred, ctxIdx, sign);

                // Step 7: Prediction Residuum berechnen
                int actual = orgData[y * m_width + x];
                int predictionResiduum = CalcPredResiduum(actual, correctedPred, sign);

                // Residuum speichern
                data[y * m_width + x] = static_cast<PGMImage::Sample>(
                    std::clamp(predictionResiduum, -128, 127)
                );

                // Step 11: Kontext-Statistiken aktualisieren
                updateContext(ctxIdx, predictionResiduum);

                // Step 12: Feinjustierung von C
                updateCumoContext(ctxIdx, predictionResiduum);
            }
        }
    }

    void addPrediction()
    {
        const PGMImage::Sample* residualData = m_data.data(); // Hier sind die Residuen
        PGMImage::Sample* decodedData = m_data.data();         // Wir schreiben die rekonstruierten Pixel hierhin
    
        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                // Step 1: lokale Gradienten berechnen
                int g1, g2, g3;
                std::tie(g1, g2, g3) = computeGradients(decodedData, x, y);
            
                // Step 3: Gradienten quantisieren
                int q1 = quantize(g1);
                int q2 = quantize(g2);
                int q3 = quantize(g3);
            
                // Step 4: Kontextindex und Sign bestimmen
                int sign;
                int ctxIdx = indexBestimmung(q1, q2, q3, sign);
            
                // Step 5: Feste MED-Vorhersage
                int fixedPred = getFixedPrediction(decodedData, x, y);
            
                // Step 6: Adaptive Bias-Korrektur
                int correctedPred = MedianCorrection(fixedPred, ctxIdx, sign);
            
                // Step 7: Residuum holen
                int residual = residualData[y * m_width + x];
            
                // Sign berücksichtigen
                int e = (sign > 0) ? residual : -residual;
            
                // Pixelwert rekonstruieren
                int decodedPixel = correctedPred + e;
                decodedPixel = std::clamp(decodedPixel, 0, 255);
            
                decodedData[y * m_width + x] = static_cast<PGMImage::Sample>(decodedPixel);
            
                // Kontext-Update
                updateContext(ctxIdx, e);
                updateCumoContext(ctxIdx, e);
            }
        }
    }


private:
    // ------------------------------
    // 1. Gradientenberechnung
    // ------------------------------
    std::tuple<int,int,int> computeGradients(const PGMImage::Sample* data, int x, int y)
    {
        auto get = [&](int xx, int yy) -> PGMImage::Sample {
            if (xx < 0) xx = 0;
            if (yy < 0) yy = 0;
            if (xx >= m_width) xx = m_width - 1;
            if (yy >= m_height) yy = m_height - 1;
            return data[yy * m_width + xx];
        };

        PGMImage::Sample a = get(x - 1, y);     // links
        PGMImage::Sample b = get(x, y - 1);     // oben
        PGMImage::Sample c = get(x - 1, y - 1); // oben-links
        PGMImage::Sample d = get(x + 1, y - 1); // oben-rechts

        int g1 = int(d) - int(b);
        int g2 = int(b) - int(c);
        int g3 = int(c) - int(a);

        return {g1, g2, g3};
    }

    // ------------------------------
    // 3. Gradienten-Quantisierung
    // ------------------------------
    inline int quantize(int g)
    {
        const int T1 = 3, T2 = 7, T3 = 21;

        if (g <= -T3) return -4;
        else if (g <= -T2) return -3;
        else if (g <= -T1) return -2;
        else if (g < 0) return -1;
        else if (g == 0) return 0;
        else if (g < T1) return 1;
        else if (g < T2) return 2;
        else if (g < T3) return 3;
        else return 4;
    }

    // ------------------------------
    // 4. Kontextindex bestimmen
    // ------------------------------
    inline int indexBestimmung(int q1, int q2, int q3, int& sign)
    {
        sign = 1;
        if (q1 < 0 || (q1 == 0 && q2 < 0) || (q1 == 0 && q2 == 0 && q3 < 0))
        {
            sign = -1;
            q1 = -q1; q2 = -q2; q3 = -q3;
        }

        int ctxIndex = ((q1 + 4) * 9 + (q2 + 4)) * 9 + (q3 + 4);
        return std::clamp(ctxIndex, 1, 364);
    }

    // ------------------------------
    // 5. MED Predictor
    // ------------------------------
    int getFixedPrediction(const PGMImage::Sample* data, int x, int y) const
    {
        auto get = [&](int xx, int yy) -> PGMImage::Sample {
            if (xx < 0) xx = 0;
            if (yy < 0) yy = 0;
            if (xx >= m_width) xx = m_width - 1;
            if (yy >= m_height) yy = m_height - 1;
            return data[yy * m_width + xx];
        };

        PGMImage::Sample a = get(x - 1, y);
        PGMImage::Sample b = get(x, y - 1);
        PGMImage::Sample c = get(x - 1, y - 1);

        if (c >= std::max(a,b)) return std::min(a,b);
        else if (c <= std::min(a,b)) return std::max(a,b);
        else return a + b - c;
    }

    // ------------------------------
    // 6. Adaptive Bias-Korrektur
    // ------------------------------
    inline int MedianCorrection(int fixedPred, int ctxIdx, int sign)
    {
        int C = m_ctxs[ctxIdx].C; // direkt auf Member zugreifen
        int correctedPred = (sign > 0) ? fixedPred + C : fixedPred - C;
        return std::clamp(correctedPred, 0, 255);
    }

    // ------------------------------
    // 7. Prediction Residuum
    // ------------------------------
    inline int CalcPredResiduum(int actual, int correctedPred, int sign)
    {
        int error = actual - correctedPred;
        if (sign < 0) error = -error;
        return (error + 128) % 256 - 128;
    }

    // ------------------------------
    // 11. Kontext-Update
    // ------------------------------
    inline void updateContext(int ctxIdx, int e)
    {
        auto& ctx = m_ctxs[ctxIdx];

        ctx.B += e;
        ctx.A += std::abs(e);

        if (ctx.N == 64)
        {
            ctx.A /= 2;
            ctx.B /= 2;
            ctx.N /= 2;
        }

        ctx.N += 1;
        ctx.C = (ctx.B >= 0) ? (ctx.B + (ctx.N/2)) / ctx.N : (ctx.B - (ctx.N/2)) / ctx.N;
    }

    // ------------------------------
    // 12. Feinjustierung von C
    // ------------------------------
    inline void updateCumoContext(int ctxIdx, int)
    {
        auto& ctx = m_ctxs[ctxIdx];

        if (ctx.B <= -ctx.N) { 
          ctx.C -= 1;
          ctx.B += ctx.N;
          if (ctx.B <= -ctx.N) ctx.B = -ctx.N + 1; 
        }
        else if (ctx.B > 0){ 
          ctx.C += 1;
          ctx.B -= ctx.N;
          if (ctx.B > 0) ctx.B = 0;
        }

        ctx.C = std::clamp(ctx.C, -128, 127);
    }

    // ------------------------------
    // Kontextstruktur
    // ------------------------------
    struct Context
    {
        int A = 0;
        int B = 0;
        int N = 1;
        int C = 0;
    };

private:
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





//======================================================
//
//   M A I N
//
//======================================================
//int main( int argc, char** argv )
//{
//  cmdPars pars;
//  if( !readCmdLine( argc, argv, pars ) )
//    return 1;
//
//  if( pars.decode )
//    decode( pars.inname, pars.outname );
//  else
//    encode( pars.inname, pars.outname );
//
//  return 0;
//}
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