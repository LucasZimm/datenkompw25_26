
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
        m_ctxs.resize(365); // 9^3 = 729 / 2 (für Vorzeichen) → 365 Kontexte
    }

    // ------------------------------
    // Subtract prediction (encoding)
    // ------------------------------
    void subtractPrediction()
    {
        const PGMImage::Sample* orgData = m_org.data();
        PGMImage::Sample* data = m_data.data();

        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                int ctx;
                int pred = getPrediction(orgData, x, y, ctx);

                int e = int(orgData[y * m_width + x]) - pred;
                if (ctx < 0) e = -e;

                // modulo correction for residuals
                int halfA = a / 2;
                e = ((e + halfA) % a) - halfA;

                data[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(e, -128, 127));

                updateContext(ctx, e);
            }
        }
    }

    // ------------------------------
    // Add prediction (decoding)
    // ------------------------------
    void addPrediction()
    {
        PGMImage::Sample* recData = m_org.data();
        const PGMImage::Sample* data = m_data.data();

        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                int ctx;
                int pred = getPrediction(recData, x, y, ctx); // auf rekonstruierten Nachbarn

                int e = int(data[y * m_width + x]);
                if (ctx < 0) e = -e;

                int val = (pred + e + a) % a;
                recData[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(val, 0, 255));

                updateContext(ctx, e); // wichtig: hier den Fehlerwert verwenden
            }
        }

        m_data = m_org; // optional, falls du rekonstruiertes Bild zurück in m_data willst
    }

private:
    // ------------------------------
    // Median Predictor + Context
    // ------------------------------
    PGMImage::Sample getPrediction(const PGMImage::Sample* data, int x, int y, int& ctx)
    {
        // Hilfsfunktion: clamp auf gültige Koordinaten (Randbehandlung)
        auto get = [&](int xx, int yy) -> PGMImage::Sample {
            if (xx < 0) xx = 0;
            if (yy < 0) yy = 0;
            if (xx >= m_width) xx = m_width - 1;
            if (yy >= m_height) yy = m_height - 1;
            return data[yy * m_width + xx];
        };

        PGMImage::Sample A = get(x - 1, y);     // links
        PGMImage::Sample B = get(x, y - 1);     // oben
        PGMImage::Sample C = get(x - 1, y - 1); // oben-links
        PGMImage::Sample D = get(x + 1, y - 1); // oben-rechts

        ctx = computeContext(A, B, C, D);
        size_t ctxIdx = static_cast<size_t>(std::clamp(std::abs(ctx), 0, static_cast<int>(m_ctxs.size() - 1)));
        Context& c = m_ctxs[ctxIdx];

        // LOCO-I Median Predictor
        int pred;
        if (C >= std::max(A, B))
            pred = std::min(A, B);
        else if (C <= std::min(A, B))
            pred = std::max(A, B);
        else
            pred = A + B - C;

        int corr = c.C;
        if (ctx < 0) pred -= corr;
        else pred += corr;

        return static_cast<PGMImage::Sample>(std::clamp(pred, 0, 255));
    }

    // ------------------------------
    // Gradient Quantisierung
    // ------------------------------
    inline int quantize(int g) const
    {
        const int T1 = 3, T2 = 7, T3 = 21;
        if (g <= -T3) return -4;
        if (g <= -T2) return -3;
        if (g <= -T1) return -2;
        if (g < 0) return -1;
        if (g == 0) return 0;
        if (g < T1) return 1;
        if (g < T2) return 2;
        if (g < T3) return 3;
        return 4;
    }

    // ------------------------------
    // Contextnummer berechnen
    // ------------------------------
    inline int computeContext(PGMImage::Sample A, PGMImage::Sample B, PGMImage::Sample C, PGMImage::Sample D)
    {
        int g1 = D - B;
        int g2 = B - C;
        int g3 = C - A; // korrigiert

        int q1 = quantize(g1);
        int q2 = quantize(g2);
        int q3 = quantize(g3);

        int sign = 1;
        if (q1 < 0) sign = -1;
        else if (q1 == 0 && q2 < 0) sign = -1;
        else if (q1 == 0 && q2 == 0 && q3 < 0) sign = -1;

        if (sign < 0) { q1 = -q1; q2 = -q2; q3 = -q3; }

        int idx = (q1 + 4) * 81 + (q2 + 4) * 9 + (q3 + 4);
        return (sign < 0) ? -idx : idx;
    }

    // ------------------------------
    // Kontext aktualisieren
    // ------------------------------
    inline void updateContext(int ctxId, int e)
    {
        if (ctxId == 0) return;

        size_t idx = static_cast<size_t>(std::clamp(std::abs(ctxId), 0, static_cast<int>(m_ctxs.size() - 1)));
        Context& ctx = m_ctxs[idx];

        ctx.B += e;
        ctx.N += 1;

        if (ctx.B <= -ctx.N)
        {
            ctx.C -= 1;
            ctx.B += ctx.N;
            if (ctx.B <= -ctx.N) ctx.B = -ctx.N + 1;
        }
        else if (ctx.B > 0)
        {
            ctx.C += 1;
            ctx.B -= ctx.N;
            if (ctx.B > 0) ctx.B = 0;
        }

        if (ctx.N >= 64)
        {
            ctx.N >>= 1;
            ctx.B >>= 1;
        }
    }

    // ------------------------------
    // Kontextstruktur
    // ------------------------------
    struct Context
    {
        int B = 0;
        int N = 1;
        int C = 0;
    };

private:
    const int a = 256;
    int m_width;
    int m_height;

    std::vector<PGMImage::Sample>& m_data;
    std::vector<PGMImage::Sample> m_org;

    std::vector<Context> m_ctxs;
};




//======================================================
//
//   M A I N    E N C O D I N G   +   D E C O D I N G
//
//======================================================
void encode( const std::string& inname, const std::string& outname )
{

  // read original image
  PGMImage img;

  img.read( inname );

  // create bitstream + write header
  std::ofstream stream( outname, std::ios::out|std::ios::binary );
  OBitstream    bs( stream );
  assert( stream.good() );
  bs.addFixed<unsigned>( img.getWidth(),  16 );
  bs.addFixed<unsigned>( img.getHeight(), 16 );
  
  // code image block by block
  EntropyEncoder  eenc( bs );
  Prediction      pred( img.getWidth(), img.getHeight(), img.getData() );

  // apply prediction
  pred.subtractPrediction();
  {
    std::ofstream residFile("residues.txt");
    PGMImage::Sample* data = img.getData().data();
    for (int k = 0; k < img.getSize(); k++) {
        residFile << int(data[k]) << "\n";  // int cast für negative Werte
    }
}
  // encode prediction error signal
  PGMImage::Sample* data = img.getData().data();
  for( int k = 0; k < img.getSize(); k++ )
    eenc.encodeSample( data[k] );
  eenc.finish();

  bs.byteAlign();
  assert( stream.good() );
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