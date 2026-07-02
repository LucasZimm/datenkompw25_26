
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
#include <cstdlib>

#include <math.h>

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
    m_pmfAbs    = std::vector<uint8_t>(N, uint8_t{0}); 
    m_pmfSigCtx = std::vector<uint8_t>(7, uint8_t{0});
    m_pmfSignCtx= std::vector<uint8_t>(128, uint8_t{0});
    m_pmfGt1Ctx = std::vector<uint8_t>(7, uint8_t{0});
    m_pmfGt2Ctx = std::vector<uint8_t>(7, uint8_t{0});
    m_pmfGt3Ctx = std::vector<uint8_t>(7, uint8_t{0});
    m_pmfGt4Ctx = std::vector<uint8_t>(7, uint8_t{0});
    m_pmfExpCtx = std::vector<uint8_t>(1000, uint8_t{0});

  }
protected:
  static const unsigned   N = 8;
  std::vector<uint8_t>    m_pmfAbs;
  std::vector<uint8_t>    m_pmfSigCtx;
  std::vector<uint8_t>    m_pmfGt1Ctx;
  std::vector<uint8_t>    m_pmfGt2Ctx;
  std::vector<uint8_t>    m_pmfGt3Ctx;
  std::vector<uint8_t>    m_pmfGt4Ctx;
  std::vector<uint8_t>    m_pmfExpCtx;
  std::vector<uint8_t>    m_pmfSignCtx;
};

class EntropyEncoder : protected EntropyCoderBase
{
public:
  EntropyEncoder( OBitstream& bs ) : EntropyCoderBase(), aenc(bs)
  {
    aenc.start();
  }
  void encodeSampleLG( PGMImage::Sample s, int k );
  void encodeSample(PGMImage::Sample s, unsigned sigIdx, int a_err, int b_err, int c_err, int d_err, int q1, int q2, int q3);
  void encodeSampleRice(
    PGMImage::Sample s,
    unsigned sigIdx,
    int a_err, int b_err,
    int q1, int q2, int q3,
    int k);
  void encodeSampleSigRiceMapped(
    int Errval,
    unsigned sigIdx,
    int k);
  void finish()    
  { 
    aenc.finish(); 
  }
private:
  void encFixed ( unsigned val, int r );
  void encUnary ( unsigned pre );
  void encodeUnary(int x);
  ArithmeticEncoder aenc;
};

class EntropyDecoder : protected EntropyCoderBase
{
public:
  EntropyDecoder( IBitstream& bs ) : EntropyCoderBase(), adec(bs)
  {
    adec.start();
  }
  int decodeSampleLG( int k );
  PGMImage::Sample decodeSample(unsigned sigIdx, int a_err, int b_err, int c_err, int d_err, int q1, int q2, int q3);
  PGMImage::Sample decodeSampleRice(
    unsigned sigIdx,
    int a_err, int b_err,
    int q1, int q2, int q3,
    int k);
  int decodeSampleSigRiceMapped(
    unsigned sigIdx,
    int k);
private:
  unsigned decFixed ( int r );
  unsigned decUnary ();
    int decodeUnary();
  ArithmeticDecoder adec;
};

//------------------------------------------------------
// Encoder helpers
//------------------------------------------------------
static inline unsigned get_gt1_index(int a_err, int b_err, int c_err, int d_err)
  {
    int sum = std::abs(a_err) +
              std::abs(b_err) +
              std::abs(c_err) +
              std::abs(d_err);

    if      (sum <= 1) return 0;
    else if (sum <= 5) return 1;
    else if (sum <= 10) return 2;
    else if (sum <= 17) return 3;
    else if (sum <= 29) return 4;
    else if (sum <= 61) return 5;
    else                return 6;
  }
    static inline unsigned get_gtN_index(int a_err, int b_err, int c_err, int d_err, unsigned binIdx)
    {
        int N = 3;
        int sum = std::abs(a_err) + std::abs(b_err) + std::abs(c_err) + std::abs(d_err);
        unsigned bin = std::min<unsigned>(N-1, binIdx);
        if      (sum <= 1)  return 0  + bin;
        else if (sum <= 2)  return 1+N  + bin;  // Abstand muss >= N sein!
        else if (sum <= 4)  return 1+2*N  + bin;
        else if (sum <= 8)  return 1+3*N  + bin;
        else if (sum <= 16) return 1+4*N  + bin;  // Abstand muss >= N sein!
        else if (sum <= 32) return 1+5*N  + bin;
        else                return 1+6*N  + bin;
    }


    inline unsigned get_sign_index(int q1, int q2, int q3,
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

        static inline unsigned get_gt2_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);
        if      (sum <= 2)  return 0;
        else if (sum <= 7)  return 1;
        else if (sum <= 14) return 2;
        else if (sum <= 23) return 3;
        else if (sum <= 38) return 4;
        else if (sum <= 72) return 5;
        else                return 6;
    }

    static inline unsigned get_gt3_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);
        if      (sum <= 3)  return 0;
        else if (sum <= 9)  return 1;
        else if (sum <= 18) return 2;
        else if (sum <= 30) return 3;
        else if (sum <= 50) return 4;
        else if (sum <= 90) return 5;
        else                return 6;
    }

    static inline unsigned get_gt4_index(int a_err, int b_err, int c_err, int d_err)
    {
        int sum = std::abs(a_err) +
                  std::abs(b_err) +
                  std::abs(c_err) +
                  std::abs(d_err);
        if      (sum <= 4)  return 0;
        else if (sum <= 11) return 1;
        else if (sum <= 22) return 2;
        else if (sum <= 35) return 3;
        else if (sum <= 60) return 4;
        else if (sum <= 110) return 5;
        else                return 6;
    }

// Encode 'r' fixed-length bits of 'val' (LSB-first order after shift).
void EntropyEncoder::encFixed( unsigned val, int r )
{
  while( r-- )
  {
    Bit b = Bit( val >> r ) & 1;
    aenc.encBit( b );
  }
}

// Encode a unary code: 'pre' zeros followed by a terminating one.
void EntropyEncoder::encUnary( unsigned pre )
{
  while( pre-- )
    aenc.encBit( 0 );
  aenc.encBit( 1 );
}

//------------------------------------------------------
// Decoder helpers
//------------------------------------------------------

// Decode 'r' fixed-length bits and return the value.
unsigned EntropyDecoder::decFixed( int r )
{
  unsigned val = 0;
  while( r-- )
  {
    val <<= 1;
    val += adec.decBit();
  }
  return val;
}

// Decode a unary code: count zeros until the terminating one.
unsigned EntropyDecoder::decUnary()
{
  unsigned pre = 0;
  while( !adec.decBit() )
    pre++;
  return pre;
}



void EntropyEncoder::encodeUnary(int x)
{
    for (int n = 0; n < x; n++)
        aenc.encBit(0);
    aenc.encBit(1);
}

int EntropyDecoder::decodeUnary()
{
    int x = 0;
    while (!adec.decBit())
        x++;
    return x;
}
//------------------------------------------------------
// encodeSample / decodeSample
//------------------------------------------------------

void EntropyEncoder::encodeSampleLG(PGMImage::Sample s, int k)
{
    //int bpp = 8;
    int LIMIT = 32;
    int QBPP = 8;
    int val    = (int)s;
    int prefix = val >> k;
    if (prefix < LIMIT - QBPP - 1) {
        encodeUnary(prefix);
        aenc.encBits(val - (prefix << k), k);
    }
    else {
        encodeUnary(LIMIT - QBPP - 1);
        aenc.encBits(val - 1, QBPP);
    }
}


int EntropyDecoder::decodeSampleLG(int k)
{
    //int bpp = 8;
    int LIMIT = 32;
    int QBPP = 8;

    int val    = 0;
    int prefix = decodeUnary();
    if (prefix < LIMIT - QBPP - 1)
        val = (prefix << k) + adec.decBits(k);
    else
        val = 1 + adec.decBits(QBPP);
    return static_cast<PGMImage::Sample>(val);
}

void EntropyEncoder::encodeSample(PGMImage::Sample s, unsigned sigIdx, int a_err, int b_err, int c_err, int d_err, int q1, int q2, int q3)
{
    unsigned binIdx  = 0;
    unsigned absValue = unsigned(s < 0 ? -s : s);
    bool nonZero = (absValue != 0);
    aenc.encBin(m_pmfSigCtx[sigIdx], nonZero);   // ← eigener Kontext
    if (!nonZero) return;


    unsigned signCtx = get_sign_index(q1, q2, q3, a_err, b_err);
    aenc.encBin(m_pmfSignCtx[signCtx], s < 0);



    bool gt1 = (absValue > 1);
    unsigned gt1Idx = get_gt1_index(a_err, b_err, c_err, d_err);
    aenc.encBin(m_pmfGt1Ctx[gt1Idx], gt1);

    if (gt1){
      bool gt2 = (absValue > 2);
      aenc.encBin(m_pmfGt2Ctx[get_gt2_index(a_err, b_err, c_err, d_err)], gt2);
      if (gt2) {
        bool gt3 = (absValue > 3);
        aenc.encBin(m_pmfGt3Ctx[get_gt3_index(a_err, b_err, c_err, d_err)], gt3);
        if (gt3) {
            bool gt4 = (absValue > 4);
            aenc.encBin(m_pmfGt4Ctx[get_gt4_index(a_err, b_err, c_err, d_err)], gt4);
            if (gt4) {
              
              unsigned rem = absValue - 5;
              while (rem--) {
                  aenc.encBin(m_pmfExpCtx[get_gtN_index(a_err, b_err, c_err, d_err, binIdx)], 1);
                  binIdx++;
              }
              aenc.encBin(m_pmfExpCtx[get_gtN_index(a_err, b_err, c_err, d_err, binIdx)], 0);
            }
        }
    }
  }
}

PGMImage::Sample EntropyDecoder::decodeSample(unsigned sigIdx, int a_err, int b_err, int c_err, int d_err, int q1, int q2, int q3)
{
    unsigned binIdx = 0;
    if (!adec.decBin(m_pmfSigCtx[sigIdx]))       // ← eigener Kontext
        return 0;

    unsigned signCtx = get_sign_index(q1, q2, q3, a_err, b_err);
    bool sign = adec.decBin(m_pmfSignCtx[signCtx]);

    unsigned mag = 1;

    unsigned gt1Idx = get_gt1_index(a_err, b_err, c_err, d_err);
    bool gt1 = adec.decBin(m_pmfGt1Ctx[gt1Idx]);
    
    if (gt1) {
      mag++;
      bool gt2 = adec.decBin(m_pmfGt2Ctx[get_gt2_index(a_err, b_err, c_err, d_err)]);
      if (gt2) {
          mag++;
          bool gt3 = adec.decBin(m_pmfGt3Ctx[get_gt3_index(a_err, b_err, c_err, d_err)]);
          if (gt3) {
            mag++;
            bool gt4 = adec.decBin(m_pmfGt4Ctx[get_gt4_index(a_err, b_err, c_err, d_err)]);
            if (gt4) {
              mag++;

              while (adec.decBin(m_pmfExpCtx[get_gtN_index(a_err, b_err, c_err, d_err, binIdx)])) {
                binIdx++;
                mag++;
              }
            }
          }
      }
    }
    return sign ? -PGMImage::Sample(mag) : PGMImage::Sample(mag);
}

void EntropyEncoder::encodeSampleRice(
    PGMImage::Sample s,
    unsigned sigIdx,
    int a_err, int b_err,
    int q1, int q2, int q3,
    int k)
{
    unsigned absValue = unsigned(s < 0 ? -s : s);

    bool nonZero = (absValue != 0);
    aenc.encBin(m_pmfSigCtx[sigIdx], nonZero);

    if (!nonZero)
        return;

    unsigned signCtx = get_sign_index(q1, q2, q3, a_err, b_err);
    aenc.encBin(m_pmfSignCtx[signCtx], s < 0);

    // Rice auf Betrag-1
    encodeSampleLG(PGMImage::Sample(absValue - 1), k);
}

PGMImage::Sample EntropyDecoder::decodeSampleRice(
    unsigned sigIdx,
    int a_err, int b_err,
    int q1, int q2, int q3,
    int k)
{
    if (!adec.decBin(m_pmfSigCtx[sigIdx]))
        return 0;

    unsigned signCtx = get_sign_index(q1, q2, q3, a_err, b_err);
    bool sign = adec.decBin(m_pmfSignCtx[signCtx]);

    unsigned riceVal = decodeSampleLG(k);

    unsigned absValue = riceVal + 1;

    return sign
        ? -static_cast<PGMImage::Sample>(absValue)
        :  static_cast<PGMImage::Sample>(absValue);
}

// EntropyEncoder
void EntropyEncoder::encodeSampleSigRiceMapped(int MErrval, unsigned sigIdx, int k)
{
    bool nonZero = (MErrval != 0);
    aenc.encBin(m_pmfSigCtx[sigIdx], nonZero);
    if (!nonZero)
        return;

    // MErrval ist bereits >=0 durch das Mapping; kein Sign-Bit, kein -1 nötig
    // ABER: wenn nonZero garantiert MErrval>=1, können wir wie bei Rice 1 abziehen,
    // um den Code dichter zu halten:
    encodeSampleLG(static_cast<PGMImage::Sample>(MErrval - 1), k);
}

// EntropyDecoder
int EntropyDecoder::decodeSampleSigRiceMapped(unsigned sigIdx, int k)
{
    if (!adec.decBin(m_pmfSigCtx[sigIdx]))
        return 0;

    unsigned val = decodeSampleLG(k);
    return static_cast<int>(val + 1);
}
//==========================
//
//Hilfsfunktion
//
//=================
inline unsigned compute_significant_Idx(
    int /*q1*/, int /*q2*/, int /*q3*/,
    int a_err, int b_err, int c_err, int d_err,
    int /*a*/, int /*b*/, int /*c*/, int /*d*/)
{
    const int T1 = 1;
    const int T2 = 3;
    const int T3 = 8;
    const int T4 = 14;
    const int T5 = 29;
    const int T6 = 61;

    int sum = std::abs(a_err) +
              std::abs(b_err) +
              std::abs(c_err) +
              std::abs(d_err);

    if (sum <= T1) return 0;
    else if (sum <= T2) return 1;
    else if (sum <= T3) return 2;
    else if (sum <= T4) return 3;
    else if (sum <= T5) return 4;
    else if (sum <= T6) return 5;
    else return 6;
}

//constexpr int k_val = 0;
constexpr int k_swap = 1;
//======================================================
//
//   P R E D I C T I O N
//
//======================================================
class Prediction
{
public:
  Prediction( int width, int height, std::vector<PGMImage::Sample>& img, int blockSize = 16, const std::vector<PGMImage::Sample>* referenceData = nullptr )
    : m_width(width), m_height(height), m_data(img), m_org(m_data), m_blockSize(std::max(1, blockSize)), m_reference(referenceData)
  {

  }

  void setBlockSize(int blockSize)
  {
    m_blockSize = std::max(1, blockSize);
  }
  void subtractPrediction(EntropyEncoder& eenc, const std::string& histPath = "", const std::string& csvPath = "")
  {
    std::vector<int>        errorStorage(m_width * m_height, 0);
    std::vector<int>        merrorStorage(m_width * m_height, 0);
    std::ofstream csvOut;
    if (!csvPath.empty()) {
        csvOut.open(csvPath);
        csvOut << "row,col,k_jpls,k_neighbor,k_optimal,a_err,b_err,c_err,d_err\n";
    }
    std::vector<Context> contexts(729);

    int counter = 0;
    
    for( int blockRow = 0; blockRow < m_height; blockRow += m_blockSize )
      for( int blockCol = 0; blockCol < m_width; blockCol += m_blockSize ) {
        const int rowEnd = std::min(blockRow + m_blockSize, m_height);
        const int colEnd = std::min(blockCol + m_blockSize, m_width);
        for( int row = blockRow; row < rowEnd; ++row ) {
          const PGMImage::Sample* orgLine = m_org.data() + row * m_width;
          int*                    errorLine = errorStorage.data() + row * m_width;
          int*                    merrorLine = merrorStorage.data() + row * m_width;

          for( int col = blockCol; col < colEnd; ++col ) {
            //line[col] -= getPrediction( orgLine, col, row );
            PGMImage::Sample Ra = getNeighbor(-1, 0, orgLine, col, row);
            PGMImage::Sample Rb = getNeighbor(0, -1, orgLine, col, row);
            PGMImage::Sample Rc = getNeighbor(-1, -1, orgLine, col, row);
            PGMImage::Sample Rd = getNeighbor(1, -1, orgLine, col, row);
            int a_err = getError(-1,  0, errorLine,     col, row);
            int b_err = getError( 0, -1, errorLine,     col, row);
            int c_err = getError(-1, -1, errorLine,     col, row);
            int d_err = getError( 1, -1, errorLine,     col, row);
            int a_merr = getError(-1,  0, merrorLine, col, row);
            int b_merr = getError( 0, -1, merrorLine, col, row);
            int c_merr = getError(-1, -1, merrorLine, col, row);
            int d_merr = getError( 1, -1, merrorLine, col, row);

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
 
            //line[col] = static_cast<PGMImage::Sample>(Errval);
            //int k = getNeighborBasedK(Ra,Rb,Rc,Rd);
            int k = 0 ;// codeSegmentA10(contexts[c_Nr].N, contexts[c_Nr].A);

            int MErrval;
            if constexpr (k_swap==0) {
              k = getNeighborBasedK(norm_err(a_err),norm_err(b_err),norm_err(c_err),norm_err(d_err));
            }
            else if constexpr (k_swap==1) {
              k = getNeighborBasedK(std::abs(a_merr),std::abs(b_merr),std::abs(c_merr),std::abs(d_merr));
            }
            else {
              // unten
            }
            codeSegmentA11(k, contexts[c_Nr].B, contexts[c_Nr].N, Errval, MErrval);

            errorLine[col] = Errval;
            merrorLine[col] = MErrval;


            unsigned sigIdx = compute_significant_Idx(
              Q1, Q2, Q3,
              a_err, b_err, c_err, d_err,
              Ra, Rb, Rc, Rd);
            

            if constexpr (k_swap==0) {
              eenc.encodeSampleRice(static_cast<PGMImage::Sample>(Errval), sigIdx, a_err, b_err, Q1, Q2, Q3,k);
            }
            else if constexpr (k_swap==1){
              eenc.encodeSampleLG(static_cast<PGMImage::Sample>(MErrval),k);
              int k_jpls = codeSegmentA10(contexts[c_Nr].N, contexts[c_Nr].A);
              int k_neighbor = k; //getNeighborBasedK(norm_err(a_err), norm_err(b_err), norm_err(c_err), norm_err(d_err));
              int k_optimal  = getOptimalK(MErrval);
              m_histJpls[k_jpls]++;
              m_histNeighbor[k_neighbor]++;
              m_histOptimal[k_optimal]++;
              if (csvOut.is_open())
                 csvOut << row << "," << col << "," << k_jpls << "," << k_neighbor << "," << k_optimal << "," << norm_err(a_err) << "," << norm_err(b_err) << "," << norm_err(c_err) << "," << norm_err(d_err) << "\n";
            }
            else {
            int k_jpls = codeSegmentA10(contexts[c_Nr].N, contexts[c_Nr].A);
            int MErrval_jpls;
            codeSegmentA11(k_jpls, contexts[c_Nr].B, contexts[c_Nr].N, Errval, MErrval_jpls);
                    
            sigIdx = compute_significant_Idx(
                Q1, Q2, Q3, a_err, b_err, c_err, d_err, Ra, Rb, Rc, Rd);
            
            eenc.encodeSampleSigRiceMapped(MErrval_jpls, sigIdx, k_jpls);
            }
            

            if (m_debugCounter < 120) {
                std::cerr << "[ENC-DBG] row=" << row << " col=" << col
                          << " orig=" << static_cast<int>(orgLine[col])
                          << " Px=" << Px
                          << " Q1=" << Q1 << " Q2=" << Q2 << " Q3=" << Q3
                          << " SIGN=" << SIGN << " cNr=" << c_Nr
                          << " k=" << k
                          << " Errval=" << Errval
                          << " MErrval=" << MErrval
                          << " absErr=" << std::abs(Errval)
                          << std::endl;
            }
            ++m_debugCounter;
            ++counter;

            codeSegmentA12(contexts[c_Nr].N, contexts[c_Nr].A, contexts[c_Nr].B, Errval);

            codeSegmentA13(contexts[c_Nr].N, contexts[c_Nr].B, contexts[c_Nr].C);
          }
        }
      }
    eenc.finish();
    if (!histPath.empty())
        writeKHistograms(histPath);
    if (csvOut.is_open())
        csvOut.close();
  }
  void addPrediction(EntropyDecoder& edec)
  {
    std::vector<int>        errorStorage(m_width * m_height, 0);
    std::vector<int>        merrorStorage(m_width * m_height, 0);

    std::vector<Context> contexts(729);

    int counter = 0;

    for( int blockRow = 0; blockRow < m_height; blockRow += m_blockSize )
      for( int blockCol = 0; blockCol < m_width; blockCol += m_blockSize ) {
        const int rowEnd = std::min(blockRow + m_blockSize, m_height);
        const int colEnd = std::min(blockCol + m_blockSize, m_width);
        for( int row = blockRow; row < rowEnd; ++row ) {
          PGMImage::Sample*       recLine = m_org.data() + row * m_width;
          int*                    errorLine    = errorStorage.data() + row * m_width;
          int*                    merrorLine   = merrorStorage.data() + row * m_width;
          for( int col = blockCol; col < colEnd; ++col ) {
            //recLine[col] = line[col] + getPrediction( recLine, col, row );
            PGMImage::Sample Ra = getNeighbor(-1, 0, recLine, col, row);
            PGMImage::Sample Rb = getNeighbor(0, -1, recLine, col, row);
            PGMImage::Sample Rc = getNeighbor(-1, -1, recLine, col, row);
            PGMImage::Sample Rd = getNeighbor(1, -1, recLine, col, row);
            int a_err = getError(-1,  0, errorLine,     col, row);
            int b_err = getError( 0, -1, errorLine,     col, row);
            int c_err = getError(-1, -1, errorLine,     col, row);
            int d_err = getError( 1, -1, errorLine,     col, row);
            int a_merr = getError(-1,  0, merrorLine, col, row);
            int b_merr = getError( 0, -1, merrorLine, col, row);
            int c_merr = getError(-1, -1, merrorLine, col, row);
            int d_merr = getError( 1, -1, merrorLine, col, row);

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

            //int k = getNeighborBasedK(Ra,Rb,Rc,Rd);
            int k = 0; //codeSegmentA10(contexts[c_Nr].N, contexts[c_Nr].A);
            if constexpr (k_swap==0) {
              k = getNeighborBasedK(norm_err(a_err),norm_err(b_err),norm_err(c_err),norm_err(d_err));
            }
            else if constexpr (k_swap==1) {
              k = getNeighborBasedK(std::abs(a_merr),std::abs(b_merr),std::abs(c_merr),std::abs(d_merr));
            }
            else {
              // unten
            }
            int Errval;
            int MErrval;
            //Errval = line[col];

            unsigned sigIdx = compute_significant_Idx(
                Q1, Q2, Q3,
                a_err, b_err, c_err, d_err,
                Ra, Rb, Rc, Rd);
              //MErrval = edec.decodeSample(sigIdx, a_err, b_err, c_err, d_err, Q1, Q2, Q3);

           

            if constexpr (k_swap==0) {
              MErrval = edec.decodeSampleRice(sigIdx, a_err, b_err, Q1, Q2, Q3,k);
              Errval = MErrval;
            }
            else if constexpr (k_swap==1){
              MErrval = edec.decodeSampleLG(k);
              inversecodeSegmentA11(k, contexts[c_Nr].B, contexts[c_Nr].N, Errval, MErrval);
            } else {
              int k_jpls = codeSegmentA10(contexts[c_Nr].N, contexts[c_Nr].A);

              sigIdx = compute_significant_Idx(
                  Q1, Q2, Q3, a_err, b_err, c_err, d_err, Ra, Rb, Rc, Rd);
              
              MErrval = edec.decodeSampleSigRiceMapped(sigIdx, k_jpls);
              inversecodeSegmentA11(k_jpls, contexts[c_Nr].B, contexts[c_Nr].N, Errval, MErrval);
            }

            errorLine[col] = Errval;  
            merrorLine[col] = MErrval;

            codeSegmentA12(contexts[c_Nr].N, contexts[c_Nr].A, contexts[c_Nr].B, Errval);

            Errval = (SIGN == 1) ? Errval : -Errval;

            int Rx;
            Rx = Px + Errval; 

            Rx %= 256; // Modulo-Operation, um sicherzustellen, dass Rx im Bereich [0..255] liegt

            if (Rx < 0) Rx += 256; // Falls Rx negativ ist, korrigieren
            else if (Rx > 255) Rx -= 256; // Falls Rx über 255 liegt, korrigiere

            if (Rx < 0) Rx = 0;
            else if (Rx > 255) Rx = 255;

            if (m_debugCounter < 12000) {
                std::cerr << "[DEC-DBG] row=" << row << " col=" << col
                          << " Px=" << Px
                          << " Q1=" << Q1 << " Q2=" << Q2 << " Q3=" << Q3
                          << " SIGN=" << SIGN << " cNr=" << c_Nr
                          << " k=" << k
                          << " MErrval=" << MErrval
                          << " Errval=" << Errval
                          << " Rx=" << Rx;
                if (m_reference != nullptr && !m_reference->empty()) {
                    const int refVal = static_cast<int>((*m_reference)[row * m_width + col]);
                    std::cerr << " ref=" << refVal << " diff=" << (Rx - refVal);
                }
                std::cerr << std::endl;
            }
            if (m_reference != nullptr && !m_reference->empty()) {
                const int refVal = static_cast<int>((*m_reference)[row * m_width + col]);
                if (Rx != refVal) {
                    std::cerr << "[DEC-MISMATCH] first mismatch at row=" << row << " col=" << col
                              << " decoded=" << Rx << " expected=" << refVal
                              << " after " << m_debugCounter << " samples"
                              << " positionClass=" << ((row == 0 || row == m_height - 1 || col == 0 || col == m_width - 1) ? "edge" : "inner")
                              << " prevRow=" << (row > 0 ? static_cast<int>(recLine[col - 1]) : -1)
                              << " prevCol=" << (col > 0 ? static_cast<int>(recLine[col - 1]) : -1)
                              << std::endl;
                    std::abort();
                }
            }
            ++m_debugCounter;
            counter++;

            if (m_debugCounter < 12000) {
                std::cerr << "[DEC-WRITE] row=" << row << " col=" << col
                          << " beforeWrite=" << static_cast<int>(recLine[col])
                          << " computed=" << Rx
                          << " finalStored=" << static_cast<int>(static_cast<PGMImage::Sample>(Rx))
                          << " positionClass=" << ((row == 0 || row == m_height - 1 || col == 0 || col == m_width - 1) ? "edge" : "inner")
                          << std::endl;
            }
            recLine[col] = static_cast<PGMImage::Sample>(Rx);

            codeSegmentA13(contexts[c_Nr].N, contexts[c_Nr].B, contexts[c_Nr].C);

        
          }
        }
      }
    m_data = m_org;
  }

private:
  PGMImage::Sample getPrediction(const PGMImage::Sample* line,
                                 int x,
                                 int y) const
  {
      constexpr PGMImage::Sample fallback = 0;
  
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

int norm_err(int err) {
  if (err <= 1) return 0;
  return std::abs(err)-1;
}

int getError(int dx, int dy, const int* errorLine, int col, int row) const {
    int nx = col + dx;
    int ny = row + dy;
    if (nx < 0 || nx >= m_width || ny < 0 || ny >= m_height)
        return 0;
    return errorLine[dy * m_width + nx];
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

  void codeSegmentA11( int k, int B, int N, int Errval, int& MErrval)
  {
      if ((k == 0) && (2 * B <= -N)) {
          // prediction correction case
          if (Errval >= 0) MErrval = 2 * Errval + 1;
          else             MErrval = -2 * (Errval + 1);
      } else {
          // normal case
          if (Errval >= 0) MErrval = 2 * Errval;
          else             MErrval = -2 * Errval - 1;
      }
  }

  void inversecodeSegmentA11(int k, int B, int N, int& Errval, int MErrval)
  {
      if ((k == 0) && (2 * B <= -N)) {
          // prediction correction case
          if (MErrval % 2 == 1) Errval =  (MErrval - 1) / 2;
          else                  Errval = -(MErrval / 2 + 1);
      } else {
          // normal case
          if (MErrval % 2 == 0) Errval =  MErrval / 2;
          else                  Errval = -(MErrval + 1) / 2;
      }
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

int bitsForValue(int val, int k)
{
    const int LIMIT = 32;
    const int QBPP  = 8;

    int prefix = val >> k;

    if (prefix < LIMIT - QBPP - 1) {
        // Unary(prefix) + k Suffixbits
        return (prefix + 1) + k;
    } else {
        // Escape-Code
        return (LIMIT - QBPP - 1 + 1) + QBPP;
    }
}
int getOptimalK(int val) 
{
    int best_k = 0;
    int best_bits = INT_MAX;
    for (int k = 0; k < 8; k++) {
        int bits = bitsForValue(val, k);
        if (bits < best_bits) {
            best_bits = bits;
            best_k = k;
        }
    }
    return best_k;
}
int getNeighborBasedK(int a, int b, int c, int d)
{
    int neighbors[4] = { a, b, c, d };

    int best_k = 0;
    int best_bits = INT_MAX;

    for (int k = 0; k < 8; k++)
    {
        int bits = 0;

        for (int i = 0; i < 4; i++)
            bits += bitsForValue(neighbors[i], k);

        if (bits < best_bits)
        {
            best_bits = bits;
            best_k = k;
        }
    }

    return best_k;
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
  int                             m_blockSize;
  std::vector<PGMImage::Sample>&  m_data;
  std::vector<PGMImage::Sample>   m_org;
  const std::vector<PGMImage::Sample>* m_reference;
  int                             m_debugCounter = 0;

public:
  std::array<long long, 8> m_histJpls{};
  std::array<long long, 8> m_histNeighbor{};
  std::array<long long, 8> m_histOptimal{};

  void writeKHistograms(const std::string& path) const
  {
      std::ofstream out(path);
      out << "k\tjpls\tneighbor\toptimal\n";
      for (int k = 0; k < 8; k++) {
          out << k << "\t" << m_histJpls[k] << "\t" << m_histNeighbor[k] << "\t" << m_histOptimal[k] << "\n";
      }
  }
};





//======================================================
//
//   M A I N    E N C O D I N G   +   D E C O D I N G
//
//======================================================
void encode( const std::string& inname, const std::string& outname, const std::string& histPath = "", const std::string& csvPath = "" )
{
  PGMImage img;
  img.read( inname );

  std::ofstream stream( outname, std::ios::out|std::ios::binary );
  OBitstream    bs( stream );
  assert( stream.good() );
  bs.addFixed<unsigned>( img.getWidth(),  16 );
  bs.addFixed<unsigned>( img.getHeight(), 16 );

  EntropyEncoder  eenc( bs );
  Prediction      pred( img.getWidth(), img.getHeight(), img.getData() );

  pred.subtractPrediction(eenc, histPath, csvPath);

  bs.byteAlign();
  assert( stream.good() );
}


void decode( const std::string& inname, const std::string& outname, const std::string& referencePath = "" )
{
  // open bitstream and read header
  std::ifstream stream( inname, std::ios::in|std::ios::binary );
  IBitstream    bs( stream );
  assert( stream.good() );
  int width  = bs.getFixed<unsigned>( 16 );
  int height = bs.getFixed<unsigned>( 16 );

  const std::vector<PGMImage::Sample>* referenceData = nullptr;
  PGMImage referenceImg;
  if( !referencePath.empty() ) {
    referenceImg.read( referencePath );
    referenceData = &referenceImg.getData();
    std::cerr << "[DEBUG] Using reference image for comparison: " << referencePath << std::endl;
  }

  // create image and decode color planes block by block
  PGMImage        img( width, height );
  EntropyDecoder  edec( bs );
  Prediction      pred( img.getWidth(), img.getHeight(), img.getData(), 16, referenceData );

  // decode prediction error signal
  //PGMImage::Sample* data = img.getData().data();
  //for( int k = 0; k < img.getSize(); k++ )
  //  data[k] = edec.decodeSample();

  // apply prediction
  pred.addPrediction(edec);
  std::cout << "Decoding finished." << std::endl;
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
    std::string outputFolder = "output\\losslessImageJPEG-LSTEST";

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
            std::string histFile = outputFolder + "/" + baseName + "_khist.txt";
            std::string csvFile  = outputFolder + "/" + baseName + "_kpixels.csv";

            std::cout << "Encoding " << inputFile << " -> " << encodedFile << std::endl;
            encode(inputFile, encodedFile, histFile, csvFile);

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
            decode(encodedFile, decodedFile, inputFile);
        }
    }

    std::cout << "Fertig!" << std::endl;
    return 0;
}


