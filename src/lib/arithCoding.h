
#pragma once


#include <bit>
#include "bitstream.h"



//----------------------------------------------------------------------
// 
//   Base class for arithmetic encoder and arithmetic decoder
// 
//----------------------------------------------------------------------
class ArithmeticCoderBase
{
protected:
  static const uint8_t next_state_mps [ 128 ];
  static const uint8_t next_state_lps [ 128 ];
  static const int     entropy_bits   [ 128 ];
  static const uint8_t lps_table      [  64 ][ 4 ];
  static const uint8_t renorm_table   [  32 ];
};





//----------------------------------------------------------------------
// 
//   Arithmetic encoder
// 
//----------------------------------------------------------------------
class ArithmeticEncoder : public ArithmeticCoderBase
{
public:
  ArithmeticEncoder (OBitstream& bitstream) : m_bitstream(bitstream) {}
  ~ArithmeticEncoder() {}

  // start arithmetic codeword
  void start ()
  {
    m_range           = 510u;
    m_low             = 0u;
    m_numOutstanding  = 0u;
  }

  // finalize arithmetic codeword
  void finish()
  {
    m_low        <<= 9u;
    _outputBits     (9u);
    _putOutstanding ();
  }

  // encode a bin in regular mode (using adaptive probability)
  void encBin(uint8_t& prob, bool bin)
  {
    if (mode != Mode::CODING)
    {
      sumScaledBits += entropy_bits[prob ^ uint8_t(bin)];
      if (mode == Mode::EST_WITH_UPDATE)
      {
        if (bin == bool(prob & 1))
          prob = next_state_mps[prob];
        else
          prob = next_state_lps[prob];
      }
      return;
    }

    unsigned lps_range = lps_table[ prob >> 1 ][ ( m_range >> 6 ) & 3 ];
    unsigned num_bits  = 0;
    m_range           -= lps_range;

    if( bin != bool(prob & 1) ) 
    { // lps
      m_low      += m_range;
      m_range     = lps_range;
      prob        = next_state_lps[prob];
      num_bits    = renorm_table[ m_range >> 3 ];
    } 
    else
    { // mps
      prob        = next_state_mps[prob];
      num_bits    = int(m_range < 256);
    }
    m_range <<= num_bits;
    m_low   <<= num_bits;
    _outputBits(num_bits);
  }

  // encode a bin in bypass mode (with probability 1/2)
  void encBit(bool bit)
  {
    if (mode != Mode::CODING)
    {
      sumScaledBits += (1u << 15);
      return;
    }

    m_low  <<= 1;
    if (bit)
    {
      m_low += m_range;
    }
    _outputBits(1u);
  }

  // encode multiple bins in regular mode (with same probability model)
  void encBins(uint8_t& prob, unsigned pattern, unsigned num)
  {
    while (num--)
    {
      encBin(prob, (pattern >> num) & 1u);
    }
  }

  // encode multiple bins in bypass mode
  void encBits(unsigned pattern, unsigned num)
  {
    while (num--)
    {
      encBit((pattern >> num) & 1u);
    }
  }

  // set estimation mode (with or without probability update)
  void set_est_mode(bool withUpdate = true)
  {
    mode = (withUpdate ? Mode::EST_WITH_UPDATE : Mode::EST_WITHOUT_UPDATE);
    sumScaledBits = 0;
  }

  // unset estimation mode and get estimate for number of coded bits
  float unset_est_mode()
  {
    static constexpr float rateScale = 1.0f / (float)(1 << 15);
    mode = Mode::CODING;
    return rateScale * float(sumScaledBits);
  }

private:
  void _outputBits(unsigned num)
  {
    const unsigned carry    = 1u << num;
    unsigned       pattern  = m_low >> 9;
    m_low                  &= 511u;
    if (pattern >= carry)
    {
      _putCarry();
      pattern -= carry;
    }
    const unsigned trailingOnes = (unsigned)std::countr_one(pattern);
    if (trailingOnes < num)
    {
      _putOutstanding();
      m_numOutstanding = trailingOnes + 1u;
      m_bitstream.addFixed(pattern >> m_numOutstanding, num - m_numOutstanding);
    }
    else if (m_numOutstanding)
    {
      m_numOutstanding += num;
    }
    else
    {
      m_bitstream.addFixed(pattern, num);
    }
  }

  void _putCarry()
  {
    m_bitstream.addBit(1);
    while (--m_numOutstanding > 1u)
    {
      m_bitstream.addBit(0);
    }
  }

  void _putOutstanding()
  {
    if (m_numOutstanding)
    {
      m_bitstream.addBit(0);
      while (--m_numOutstanding)
      {
        m_bitstream.addBit(1);
      }
    }
  }

private:
  OBitstream&   m_bitstream;
  unsigned      m_range           = 510u;
  unsigned      m_low             = 0u;
  unsigned      m_numOutstanding  = 0u;

  enum Mode{ CODING, EST_WITHOUT_UPDATE, EST_WITH_UPDATE};
  Mode          mode              = CODING;
  int64_t       sumScaledBits     = 0;
};





//----------------------------------------------------------------------
// 
//   Arithmetic decoder
// 
//----------------------------------------------------------------------
class ArithmeticDecoder : public ArithmeticCoderBase
{
public:
  ArithmeticDecoder (IBitstream& bitstream) : m_bitstream(bitstream) {}
  ~ArithmeticDecoder() {}

  // start decoding of arithmetic codeword
  void start ()
  {
    m_range = 510u;
    m_value = m_bitstream.getFixed<unsigned>(9u);
  }

  // finalize decoding of arithmetic codeword
  void finish(bool byteAlign = false)
  {
    if (byteAlign)
    {
      m_bitstream.byteAlign();
    }
  }

  // decode a bin in regular mode (using adaptive probability)
  bool decBin(uint8_t& prob)
  {
    bool     bin       = bool(prob & 1u);
    unsigned lps_range = lps_table[ prob >> 1 ][ ( m_range >> 6 ) & 3 ];
    unsigned num_bits  = 0;
    m_range           -= lps_range;

    if (m_value >= m_range)
    { // lps
      m_value    -= m_range;
      m_range     = lps_range;
      prob        = next_state_lps[prob];
      num_bits    = renorm_table[ m_range >> 3 ];
      bin         = !bin;
    }
    else
    { // mps
      prob        = next_state_mps[prob];
      num_bits    = int(m_range < 256);
    }
    m_range <<= num_bits;
    m_value <<= num_bits;
    m_value  += m_bitstream.getFixed<unsigned>(num_bits);
    return bin;
  }

  // decode a bin in bypass mode (with probability 1/2)
  bool decBit()
  {
    m_value    <<= 1;
    m_value     += m_bitstream.getFixed<unsigned>(1);
    if (m_value >= m_range)
    {
      m_value   -= m_range;
      return true;
    }
    return false;
  }

  // decode multiple bins in regular mode (with same probability model)
  unsigned decBins(uint8_t& prob, unsigned num)
  {
    unsigned pattern = 0;
    while (num--)
    {
      pattern = (pattern << 1) + (unsigned)decBin(prob);
    }
    return pattern;
  }

  // decode multiple bins in bypass mode
  unsigned decBits(unsigned num)
  {
    unsigned pattern = 0;
    while (num--)
    {
      pattern = (pattern << 1) + (unsigned)decBit();
    }
    return pattern;
  }

private:
  IBitstream&   m_bitstream;
  unsigned      m_range     = 510u;
  unsigned      m_value     = 0u;
};
























//===== interface for quantized pmf =====
class IPmf
{
public:
  virtual ~IPmf()  {}
  virtual uint64_t operator[] ( unsigned index ) const = 0;
  virtual void     update     ( unsigned index )       = 0;
};



//===== ARITHMETIC ENCODER =====
class ArithEnc
{
public:
  ArithEnc ( unsigned Ubits, unsigned Vbits, OBitstream& outBS );
  ~ArithEnc();

  void encode( unsigned index, IPmf& pmf );
  void terminate();

private:
  // helper functions
  uint64_t getCv            ( unsigned index, const IPmf& pmf ); // get value of modified cmf
  unsigned leadingZerosUV   ( uint64_t value );                  // get number of leading zeros
                                                                 //   in (U+V)-bit integer
  unsigned trailingOnesUV   ( uint64_t value, unsigned lz );     // get number of trailing ones
                                                                 //   in first lz bits of (U+V)-bit integer
  void     putOutstanding   ();                                  // output outstanding bits
                                                                 //   and update counter for these
  void     putInvOutstanding();                                  // output inverted outstanding bits
                                                                 //   and update counter for these
  void     putFirstOfUVInt  ( uint64_t value, unsigned num );    // output first num bits of (U+V)-bit integer

private:
  // configuration
  const unsigned U;     // number of bits for interval width
  const unsigned V;     // number of bits for probability masses
  const uint64_t carry; // mask for carry bit at position U+V+1
  const uint64_t mask;  // mask for all bits of an (U+V)-bit integer

  // for actual coding
  OBitstream&    bs;    // output bitstream
  uint64_t       A;     // U-bit integer representing interval width
  uint64_t       B;     // (U+V)-bit integer representing lower interval boundary
  int64_t        c;     // number of outstanding bits
};



//===== ARITHMETIC DECODER =====
class ArithDec
{
public:
  ArithDec ( unsigned Ubits, unsigned Vbits, IBitstream& outBS );
  ~ArithDec();

  unsigned decode( IPmf& pmf );

private:
  // helper functions
  unsigned leadingZerosUV( uint64_t value );  // get number of leading zeros in (U+V)-bit integer

private:
  // configuration
  const unsigned U;     // number of bits for interval width
  const unsigned V;     // number of bits for probability masses

  // for actual coding
  IBitstream&    bs;    // output bitstream
  uint64_t       A;     // U-bit integer representing interval width
  uint64_t       v;     // (U+V)-bit integer representing value minus lower interval boundary
};

