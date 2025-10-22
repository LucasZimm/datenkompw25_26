
#pragma once


#include "bitstream.h"



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

