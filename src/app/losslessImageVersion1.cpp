
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
class BinaryPmfFixed : public IPmf
{
public:
  BinaryPmfFixed( unsigned Vbits ) : V( Vbits ) {}
  uint64_t operator[] ( unsigned ) const  { return uint64_t(1)<<(V-1); }
  void     update     ( unsigned )        {}
private:
  const unsigned   V;
};


class BinaryPmfAdaptive : public IPmf
{
public:
  BinaryPmfAdaptive( unsigned Vbits, unsigned winSize ) : V( Vbits ), W( winSize ), p0( 1<<(V-1) )
  { 
    if( W >= V-1 ) throw 1; 
  }
  uint64_t operator[] ( unsigned index ) const  { return index ? (1<<V)-p0 : p0; }
  void     update     ( unsigned index )
  {
    if( index )   p0 -= ( p0 >> W );
    else          p0 += ( ( (1<<V) - p0 ) >> W );
  }
private:
  const unsigned   V;
  const unsigned   W; // window size
  unsigned         p0;
};


class EntropyCoderBase
{
protected:
  EntropyCoderBase() : m_pmfFixed(V)
  {
    // initialization of binary pmfs for absolute values
    for( unsigned k = 0; k < N; k++ )
      m_pmfAbs.push_back( BinaryPmfAdaptive( V, W ) );
  }
protected:
  static const unsigned     V = 11;     // number of bits for probability masses
  static const unsigned     U = 20;     // number of bits for interval width
  static const unsigned     W = 5;      // window size for pmf update
  static const unsigned     N = 3;      // number of probability models for unary binarization of absolute values
  std::vector<BinaryPmfAdaptive>  m_pmfAbs;   // probability models for coding absolute values
  BinaryPmfFixed                  m_pmfFixed; // fixed probability model (for sign, mode bits)
};


class EntropyEncoder : protected EntropyCoderBase
{
public:
  EntropyEncoder( OBitstream& bs ) : EntropyCoderBase(), aenc( U, V, bs )
  {}
  void encodeSample  ( PGMImage::Sample s );
  void finish        ()    { aenc.terminate(); }
private:
  ArithEnc aenc;
};

class EntropyDecoder : protected EntropyCoderBase
{
public:
  EntropyDecoder( IBitstream& bs ) : EntropyCoderBase(), adec( U, V, bs )
  {}
  PGMImage::Sample  decodeSample  ();
private:
  ArithDec adec;
};

void EntropyEncoder::encodeSample( PGMImage::Sample s )
{
  unsigned absValue = unsigned( s < 0 ? -s : s );
  unsigned rem      = absValue;
  unsigned binIdx   = 0;
  while( rem-- )
  {
    aenc.encode ( 1, m_pmfAbs[ std::min<unsigned>( N-1, binIdx++ ) ] );
  }
  aenc.encode   ( 0, m_pmfAbs[ std::min<unsigned>( N-1, binIdx++ ) ] );
  if( absValue )
    aenc.encode ( s<0, m_pmfFixed );
}

PGMImage::Sample EntropyDecoder::decodeSample()
{
  PGMImage::Sample  s      = 0;
  unsigned          binIdx = 0;
  while( adec.decode( m_pmfAbs[ std::min<unsigned>( N-1, binIdx++ ) ] ) )
  {
    s++;
  }
  if( s && adec.decode( m_pmfFixed ) )
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
  Prediction( int width, int height, std::vector<PGMImage::Sample>& img )
    : m_width(width), m_height(height), m_data(img), m_org(m_data)
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

  // encode prediction error signal
  PGMImage::Sample* data = img.getData().data();
  for( int k = 0; k < img.getSize(); k++ )
    eenc.encodeSample( data[k] );
  eenc.finish();

  bs.flush();
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
    std::string inputFolder = "kodak-pgm";
    std::string outputFolder = "output";

    std::cout << fs::current_path() << std::endl;
    std::ofstream report("bitrates.txt", std::ios::out);
    if (!report.is_open()) {
        std::cerr << "Fehler: Konnte bitrates.txt nicht öffnen!" << std::endl;
        return 1;
    }

    // sicherstellen, dass der Output-Ordner existiert
    fs::create_directories(outputFolder);

    // durch alle Dateien im Eingabeordner gehen
    for (const auto& entry : fs::directory_iterator(inputFolder)) {
        if (entry.path().extension() == ".pgm") {
            std::string inputFile = entry.path().string();

            // Name ohne Endung:
            std::string baseName = entry.path().stem().string();

            // Ausgabe-Dateinamen
            std::string encodedFile = outputFolder + "/" + baseName + ".bin";
            std::string decodedFile = outputFolder + "/" + baseName + "_decoded.pgm";

            std::cout << "Encoding " << inputFile << " -> " << encodedFile << std::endl;
            encode(inputFile, encodedFile);

            uintmax_t fileSize = fs::file_size(encodedFile);
            PGMImage img;
            img.read(inputFile);
            size_t numPixels = img.getSize();
            double bitrate = (fileSize * 8.0) / numPixels;
                    
            // Ausgabe in Konsole
            std::cout << baseName << ": Bitrate = " 
                      << std::fixed << std::setprecision(3) 
                      << bitrate << " bpp" << std::endl;
                    
            // Ausgabe in Datei
            report << baseName << ": " << std::fixed << std::setprecision(3) 
                   << bitrate << " bpp" << std::endl;



            std::cout << "Decoding " << encodedFile << " -> " << decodedFile << std::endl;
            decode(encodedFile, decodedFile);
        }
    }

    std::cout << "Fertig!" << std::endl;
    return 0;
}