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
//   P R E D I C T I O N
//
//======================================================
class Prediction
{
public:
  Prediction( int width, int height, std::vector<PGMImage::Sample>& img, const std::string& outputDir = "output" )
    : m_width(width), m_height(height), m_data(img), m_org(m_data)
    , m_outputDir(outputDir)
  {}
  void subtractPrediction()
  {
    const PGMImage::Sample* orgLine = m_org .data();
    PGMImage::Sample*       line    = m_data.data();
    for( int row = 0; row < m_height; row++, line+=m_width, orgLine+=m_width )
      for( int col = 0; col < m_width; col++ )
        line[col] -= getPrediction( orgLine, col, row );
  }
  void addPrediction()
  {
    PGMImage::Sample*       recLine = m_org .data();
    const PGMImage::Sample* line    = m_data.data();
    for( int row = 0; row < m_height; row++, line+=m_width, recLine+=m_width )
      for( int col = 0; col < m_width; col++ )
        recLine[col] = line[col] + getPrediction( recLine, col, row );
    m_data = m_org;
  }

private:
  PGMImage::Sample getPrediction( const PGMImage::Sample* line, int x, int y ) const
  {
    PGMImage::Sample X = ( x > 0 ? line[x-1] : y > 0 ? line[x-m_width]   : 128 ); // fallback
    PGMImage::Sample L = ( x > 0                     ? line[x-1]         : X );
    PGMImage::Sample A = ( y > 0                     ? line[x-m_width]   : X );
    return (L+A+1)>>1;
  }

private:
  int                             m_width;
  int                             m_height;
  std::vector<PGMImage::Sample>&  m_data;
  std::vector<PGMImage::Sample>   m_org;
  std::string                     m_outputDir;
};




//======================================================
//
//   M A I N    E N C O D I N G   +   D E C O D I N G
//
//======================================================
void encode( const std::string& inname, const std::string& outname, const std::string& outdir = "output" )
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
  Prediction      pred( img.getWidth(), img.getHeight(), img.getData(), outdir );

  // apply prediction
  pred.subtractPrediction();

  // encode prediction error signal
  PGMImage::Sample* data = img.getData().data();
  for( int k = 0; k < img.getSize(); k++ )
    eenc.encodeSample( data[k] );
  eenc.finish();

  bs.byteAlign();
  assert( stream.good() );
}


void decode( const std::string& inname, const std::string& outname, const std::string& outdir = "output" )
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
  Prediction      pred( img.getWidth(), img.getHeight(), img.getData(), outdir );

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
            decode(pars.inname, fullOutputPath, pars.outdir);
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