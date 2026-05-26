
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
//   P R E D I C T I O N
//
//======================================================
class Prediction
{
public:
  Prediction( int width, int height, std::vector<PGMImage::Sample>& img )
    : m_width(width), m_height(height), m_data(img), m_org(m_data)
  {

  }
  void subtractPrediction()
  {
    const PGMImage::Sample* orgLine = m_org .data();
    PGMImage::Sample*       line    = m_data.data();

    std::vector<Context> contexts(729);

    
    for( int row = 0; row < m_height; row++, line+=m_width, orgLine+=m_width )
      for( int col = 0; col < m_width; col++ ) {
        //line[col] -= getPrediction( orgLine, col, row );
        PGMImage::Sample Ra = getNeighbor(-1, 0, orgLine, col, row);
        PGMImage::Sample Rb = getNeighbor(0, -1, orgLine, col, row);
        PGMImage::Sample Rc = getNeighbor(-1, -1, orgLine, col, row);
        PGMImage::Sample Rd = getNeighbor(1, -1, orgLine, col, row);

        int D1 = Rd - Rb;
        int D2 = Rb - Rc;
        int D3 = Rc - Ra;

        int Q1 = codeSegmentA4(D1);
        int Q2 = codeSegmentA4(D2);
        int Q3 = codeSegmentA4(D3);

        int SIGN = findSIGN(Q1, Q2, Q3);

        int c_Nr = findContextNumber(Q1, Q2, Q3);

        int Px;
        Px = getPrediction(orgLine, col, row); 

        int C = contexts[c_Nr].C;
        codeSegmentA6(SIGN, Px, C);

        int Errval;
        codeSegmentA7(Errval, orgLine[col], Px, SIGN);

        codeSegmentA9(Errval);
 
        line[col] = static_cast<PGMImage::Sample>(Errval);

        codeSegmentA12(contexts[c_Nr].N, contexts[c_Nr].A, contexts[c_Nr].B, Errval);

        codeSegmentA13(contexts[c_Nr].N, contexts[c_Nr].B, contexts[c_Nr].C);
      }
  }
  void addPrediction()
  {
    PGMImage::Sample*       recLine = m_org .data();
    const PGMImage::Sample* line    = m_data.data();

    std::vector<Context> contexts(729);

    for( int row = 0; row < m_height; row++, line+=m_width, recLine+=m_width )
      for( int col = 0; col < m_width; col++ ) {
        //recLine[col] = line[col] + getPrediction( recLine, col, row );
        PGMImage::Sample Ra = getNeighbor(-1, 0, recLine, col, row);
        PGMImage::Sample Rb = getNeighbor(0, -1, recLine, col, row);
        PGMImage::Sample Rc = getNeighbor(-1, -1, recLine, col, row);
        PGMImage::Sample Rd = getNeighbor(1, -1, recLine, col, row);

        int D1 = Rd - Rb;
        int D2 = Rb - Rc;
        int D3 = Rc - Ra;

        int Q1 = codeSegmentA4(D1);
        int Q2 = codeSegmentA4(D2);
        int Q3 = codeSegmentA4(D3);

        int SIGN = findSIGN(Q1, Q2, Q3);

        int c_Nr = findContextNumber(Q1, Q2, Q3);

        int Px;
        Px = getPrediction(recLine, col, row); // codeSegmentA5 

        int C = contexts[c_Nr].C;
        codeSegmentA6(SIGN, Px, C);

        int Errval = line[col];

        codeSegmentA12(contexts[c_Nr].N, contexts[c_Nr].A, contexts[c_Nr].B, Errval);

        Errval = (SIGN == 1) ? Errval : -Errval;

        int Rx;
        Rx = Px + Errval; 

        Rx %= 256; // Modulo-Operation, um sicherzustellen, dass Rx im Bereich [0..255] liegt

        if (Rx < 0) Rx += 256; // Falls Rx negativ ist, korrigieren
        else if (Rx > 255) Rx -= 256; // Falls Rx über 255 liegt, korrigiere

        if (Rx < 0) Rx = 0;
        else if (Rx > 255) Rx = 255;

        recLine[col] = static_cast<PGMImage::Sample>(Rx);

        codeSegmentA13(contexts[c_Nr].N, contexts[c_Nr].B, contexts[c_Nr].C);

    
    }
    m_data = m_org;
  }

private:
  PGMImage::Sample getPrediction(const PGMImage::Sample* line,
                                 int x,
                                 int y) const
  {
      constexpr PGMImage::Sample fallback = 128;
  
      PGMImage::Sample L =
          (x > 0)
          ? line[x - 1]
          : fallback;
  
      PGMImage::Sample A =
          (y > 0)
          ? line[x - m_width]
          : fallback;
  
      PGMImage::Sample C =
          (x > 0 && y > 0)
          ? line[x - m_width - 1]
          : fallback;
  
      if (C >= std::max(L, A))
          return std::min(L, A);
  
      if (C <= std::min(L, A))
          return std::max(L, A);
  
      return L + A - C;
  }

PGMImage::Sample getNeighbor(int dx, int dy, const PGMImage::Sample* line, int col, int row) const {
    int nx = col + dx;
    int ny = row + dy;
    
    if (nx < 0 || nx >= m_width || ny < 0 || ny >= m_height)
        return 128;
    
    return line[dy * m_width + nx];  // dy*m_width verschiebt die Zeile, nx die Spalte
}


  int codeSegmentA4(int Di) {
    int T1 = 3, T2 = 7, T3 = 21;
    if (Di <= -T3)      return -4;
    else if (Di <= -T2) return -3;
    else if (Di <= -T1) return -2;
    else if (Di < 0)    return -1;
    else if (Di <= 0)   return 0;
    else if (Di < T1)   return 1;
    else if (Di < T2)   return 2;
    else if (Di < T3)   return 3;
    else 
        return 4;
  }


  void codeSegmentA6(int SIGN, int& Px, int C) {
    if (SIGN == 1) Px = Px + C;
    else           Px = Px - C;
    if (Px > MAXVAL) Px = MAXVAL;
    else if (Px < 0) Px = 0;
  }

  void codeSegmentA7(int& Errval, PGMImage::Sample Ix, int& Px, int SIGN) {
    Errval = Ix - Px;
    if (SIGN == -1) Errval = -Errval;
  }

  void codeSegmentA9(int& Errval) {
    if (Errval<0) Errval = Errval + Range;
    if (Errval >= (Range+1)/2) Errval = Errval - Range;
  }

  void codeSegmentA9inverse(int& Errval) {
    if (Errval >= (Range+1)/2) Errval = Errval - Range;
    if (Errval < 0) Errval = Errval + Range;
  }

  int codeSegmentA10(int N, int A) {
    int k = 0;
    for (k=0; (N << k) < A; k++);
    return k;
  }

  void codeSegmentA12(int& N, int& A, int& B, int Errval) {
    B = B + Errval;
    A = A + std::abs(Errval);
    if (N == 64) {
      A = A >> 1;
      if (B >= 0) B = B >> 1;
      else B= -((1-B) >> 1);
      N = N >> 1;
    }
    N=N+1;
  }

  void codeSegmentA13(int& N, int& B, int& C) {
    if (B <= -N) {
      B = B + N;
      if (C > MIN_C) C = C - 1;
      if (B <= -N) B = -N + 1;
      
    }
    else if (B > 0) {
      B = B - N;
      if (C < MAX_C) C = C + 1;
      if (B > 0) B = 0;
    }
  }

  int findSIGN(int& Q1, int& Q2, int& Q3) {
      int sign = 1;

      if (Q1 != 0) {
          if (Q1 < 0) sign = -1;
      } 
      else if (Q1 == 0 && Q2 != 0) {
          if (Q2 < 0) sign = -1;
      } 
      else if ((Q1 == 0) && Q2 == 0 && Q3 != 0) {
          if (Q3 < 0) sign = -1;
      }

      if (sign == -1) {
          Q1 = -Q1;
          Q2 = -Q2;
          Q3 = -Q3;
      }
      return sign;
  }

  int findContextNumber(int Q1, int Q2, int Q3) {
    //3d array for context numbers
    Q1 +=4; Q2 +=4; Q3 +=4; // shift to positive range
    int contextNumbers[9][9][9];
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        for (int k = 0; k < 9; k++) {
          contextNumbers[i][j][k] = i * 81 + j * 9 + k;
        }
      }
    }
    return contextNumbers[Q1][Q2][Q3];
  }
  
private:
  int MAXVAL = 255;
  int Range = MAXVAL + 1;
  int MIN_C = -128;
  int MAX_C = 127;

  struct Context {
    int N = 1;
    int A = (256+32)/64; // initial value for cumulative absolute error
    int B = 0;
    int C = 0;
  };

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
    std::string inputFolder = "kodak-pgm";
    std::string outputFolder = "output\\losslessImageJPEG-LS";

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