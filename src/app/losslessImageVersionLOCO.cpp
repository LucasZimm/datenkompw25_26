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
//   LOCO-PREDICTOR (JPEG-LS) - 1-STUFIG (NUR MED)
//
//======================================================
class Prediction
{
public:
    Prediction(int width, int height, std::vector<PGMImage::Sample>& img, const std::string& outputDir = "output")
        : m_width(width), m_height(height), m_data(img), m_org(m_data), m_outputDir(outputDir)
    {
        // Keine Kontexte needed für 1-stufig
    }

    // ------------------------------
    // Prädiktion subtrahieren (Encoding)
    // ------------------------------
    void subtractPrediction()
    {
        const PGMImage::Sample* orgData = m_org.data();
        PGMImage::Sample* data = m_data.data();
        // --------------------
        // Histogramme initialisieren
        // --------------------
        // Kontext-Histogramme
        std::array<size_t, 256> histoA_ctx{};
        std::array<size_t, 256> histoB_ctx{};
        std::array<size_t, 256> histoC_ctx{};
        // Fehler-Histogramme
        std::array<size_t, 511> histoA_err{};
        std::array<size_t, 511> histoB_err{};
        std::array<size_t, 511> histoC_err{};
        std::array<size_t, 511> histoAll_err{};


        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                // Nur feste MED-Prädiktion
                int a = 0, b = 0, c = 0;
                int pred = getFixedPrediction(orgData, x, y, a,b,c);

                // Prädiktionsfehler berechnen
                int actual = orgData[y * m_width + x];

                              // --------------------
                // Histogramme füllen
                // --------------------
                // Kontext-Histogramme
                histoA_ctx[a]++;
                histoB_ctx[b]++;
                histoC_ctx[c]++;
                int eA = (actual - a)+255; // Fehler nur relativ zum linken Nachbarn
                int eB = (actual - b)+255; // Fehler nur relativ zum oberen Nachbarn
                int eC = (actual - c)+255; // Fehler nur relativ zum oben-links Nachbarn
                histoA_err[eA]++;
                histoB_err[eB]++;
                histoC_err[eC]++;
 



                int error = actual - pred;
                histoAll_err[error + 255]++;       

                // Modulo-Korrektur für Residuen
                //error = moduloCorrection(error);
                
                // In Daten speichern
                data[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(error, -255, 255));
            }
        }
      
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

        // Fehler-Histogramme
        saveErrorHistogram(histoA_err, "histoA_err");
        saveErrorHistogram(histoB_err, "histoB_err");
        saveErrorHistogram(histoC_err, "histoC_err");
        saveErrorHistogram(histoAll_err, "histoAll_err");



    }

    // ------------------------------
    // Prädiktion addieren (Decoding)
    // ------------------------------
    void addPrediction()
    {
        PGMImage::Sample* recData = m_org.data();
        const PGMImage::Sample* data = m_data.data();
    
        for (int y = 0; y < m_height; ++y)
        {
            for (int x = 0; x < m_width; ++x)
            {
                // Nur feste MED-Prädiktion
                int a = 0, b = 0, c = 0;
                int pred = getFixedPrediction(recData, x, y, a, b ,c);
                
                // Fehlerwert aus Daten
                int error = data[y * m_width + x];
                // Rekonstruierten Wert berechnen
                int val = pred + error;
                // Modulo-Korrektur rückgängig machen
                //if (val < 0) val += m_alphabetSize;
                //if (val >= m_alphabetSize) val -= m_alphabetSize;
                
                recData[y * m_width + x] = static_cast<PGMImage::Sample>(std::clamp(val, 0, 255));
            }
        }
      
        m_data = m_org; // Rekonstruierte Daten zurückkopieren
    }

private:
    // ------------------------------
    // Feste MED-Prädiktion (Median Edge Detector)
    // ------------------------------
    int getFixedPrediction(const PGMImage::Sample* data, int x, int y, int& a, int& b, int& c) const
    {
        auto get = [&](int xx, int yy) -> PGMImage::Sample {
                    if (xx < 0 || yy < 0) return 0;
                    if (xx >= m_width || yy >= m_height) return 0;
                    return data[yy * m_width + xx];
                };
      
        PGMImage::Sample aI = get(x - 1, y);     // links
        PGMImage::Sample bI = get(x, y - 1);     // oben  
        PGMImage::Sample cI = get(x - 1, y - 1); // oben-links
        a = aI;
        b = bI;
        c = cI;
        // LOCO-I Median Predictor (MED)
        if (c >= std::max(aI, bI))
            return std::min(aI, bI);
        else if (c <= std::min(aI, bI))
            return std::max(aI, bI);
        else
            return aI + bI - cI;
    }

    // ------------------------------
    // Modulo-Korrektur für Residuen
    // ------------------------------
    inline int moduloCorrection(int error) const
    {
        const int halfA = m_alphabetSize / 2;
        error = ((error + halfA) % m_alphabetSize) - halfA;
        if (error < -halfA) error += m_alphabetSize;
        return error;
    }

private:
    const int m_alphabetSize = 256;
    int m_width;
    int m_height;
    std::string m_outputDir;

    std::vector<PGMImage::Sample>& m_data;
    std::vector<PGMImage::Sample> m_org;
};



//======================================================
//
//   M A I N    E N C O D I N G   +   D E C O D I N G
//
//======================================================
void encode(const std::string& inname, const std::string& outname, const std::string& outdir)

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
  Prediction pred(img.getWidth(), img.getHeight(), img.getData(), outdir);


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

int main(int argc, char** argv)
{
    cmdPars pars;

    if (!readCmdLine(argc, argv, pars))
        return 1;

    fs::create_directories(pars.outdir);

    std::string fullOutputPath = pars.outdir + "/" + pars.outname;

    try {
        if (pars.decode) {
            std::cout << "Dekodiere: " << pars.inname << " -> " << fullOutputPath << std::endl;
            decode(pars.inname, fullOutputPath);
        } else {
            std::cout << "Kodiere: " << pars.inname << " -> " << fullOutputPath << std::endl;
            encode(pars.inname, fullOutputPath, pars.outdir);

        }
        std::cout << "✅ Erfolgreich!" << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Fehler: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}