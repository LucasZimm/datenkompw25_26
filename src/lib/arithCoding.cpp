
#include <cassert>
#include <algorithm>

#include "arithCoding.h"



//===== ARITHMETIC ENCODER =====
ArithEnc::ArithEnc ( unsigned Ubits, unsigned Vbits, OBitstream& outBS )
  : U    ( Ubits )
  , V    ( Vbits )
  , carry( uint64_t(1) << ( U + V ) )
  , mask ( carry - 1 )
  , bs   ( outBS )
  , A    ( ( uint64_t(1) << U ) - 1 )
  , B    ( 0 )
  , c    ( 0 )
{
  assert( U+V < 64 );
}


ArithEnc::~ArithEnc()
{}


// get modified cmf[ index ]
uint64_t ArithEnc::getCv( unsigned index, const IPmf& pmf )
{
  uint64_t cv = 0;
  for( unsigned k = 0; k < index; cv += pmf[ k++ ] );
  return   cv;
}


// get number of leading zeros in (U+V)-bit integer
unsigned ArithEnc::leadingZerosUV( uint64_t value )
{
  unsigned num = 0;
  for( unsigned bitPos = U+V-1; ; bitPos--, num++ )
  {
    if( value >> bitPos )
      return num;
  }
}


// get number of trailing ones in first lz bits of (U+V)-bit integer
unsigned ArithEnc::trailingOnesUV( uint64_t value, unsigned lz )
{
  unsigned num = 0;
  for( unsigned bitPos = U+V-lz; ; bitPos++, num++ )
  {
    if( ( ( value >> bitPos ) & 1 ) == 0 )
      return num;
  }
}


// output outstanding bits and update counter for these
void ArithEnc::putOutstanding()
{
  if( c )
  {
    bs.addBit( 0 );
    c--;
  }
  while( c )
  {
    bs.addBit( 1 );
    c--;
  }
}


// output inverted outstanding bits and update counter for these
void ArithEnc::putInvOutstanding()
{
  if( c )
  {
    bs.addBit( 1 );
    c--;
  }
  while( c > 1 )
  {
    bs.addBit( 0 );
    c--;
  }
}


// output first num bits of (U+V)-bit integer
void ArithEnc::putFirstOfUVInt( uint64_t value, unsigned num )
{
  for( unsigned k = 1; k <= num; k++ )
  {
    bs.addBit( ( value >> ( U+V-k ) ) & 1 );
  }
}



//=====  MAIN ENCODING FUNCTION  =====
void ArithEnc::encode( unsigned index, IPmf& pmf )
{
  uint64_t Ax =     A * pmf   [ index ];       //   (U+V)-bit integer [intermediate new interval width]
  uint64_t Bx = B + A * getCv ( index, pmf );  // (U+V+1)-bit integer [intermediate new lower boundary]
  unsigned lz = leadingZerosUV( Ax );          // number of leading zeros in (U+V)-bit representation of Ax

  // check for carry
  if( Bx >= carry )
  {
    // carry converts "0111..1" pattern in outstanding bits into "1000..0"
    // we have to output the first (c-1) bits and keep the last
    putInvOutstanding(); // output inverted outstanding bits (except last 0) and update counter c
    Bx -= carry;         // remove carry bit from register --> Bx is now a (U+V) integer
  }

  // output settled bits and update outstanding bits
  if( lz ) // do nothing if there aren't any leading zeros in Ax (in this case, we don't shift around)
  {
    unsigned t1 = trailingOnesUV( Bx, lz ); // get number of trailing ones in first lz bits of Bx
    if( t1 < lz )
    {
      // CASE (1): There is a "0" at the most significant position in (U+V)-bit integer Bx.
      //           --> output all outstanding bits "011..1" (need to check whether there is at least one)
      //           --> output (lz-t1-1) first bits in (U+V)-integer Bx
      //           --> the next 2 bits with pattern "01" in Bx become the new outstanding bits
      putOutstanding ();              // output all outstanding bits and update counter c
      putFirstOfUVInt( Bx, lz-t1-1 ); // output first (lz-t1-1) bits of (U+V)-bit integer Bx
      c = t1 + 1;                     // set new outstanding bit counter c
                                      //   (one leading "0" + lz bits equal to "1")
    }
    else if( c )
    {
      // CASE (2): There only "1"s in the first lz bits of Bx,
      //           and we already have at least one outstanding bit.
      //           --> increase counter for outstanding bits
      c += t1;
    }
    else
    {
      // CASE (3): There only "1"s in the first lz bits of Bx, but we don't have any outstanding bits.
      //           (this case is rare, but can happen at the beginning of a sequence)
      //           --> output all lz ones at the start of Bx (still no outstanding bits)
      putFirstOfUVInt( Bx, t1 );
    }	
  }

  // update integers for interval width and lower boundary
  A = ( Ax << lz ) >> V;     // keep first U-bits with pattern "1..." of Ax
  B = ( Bx << lz ) &  mask;  // masking out remaining bits of Bx and align it with A

  // update pmf (for non-adaptive pmf implementations, this has no effect)
  pmf.update( index );
}



//=====  CODEWORD TERMINATION  =====
void ArithEnc::terminate()
{
  const unsigned numBitsOfB  = 1;                                     // (use 2 for prefix-free version)
  const uint64_t maskLastBit = uint64_t(1) << ( U + V - numBitsOfB ); // mask for last bit

  // round of lower interval boundary
  if( B & ( maskLastBit - 1 ) ) // check whether there are any "1"s after the numBitsOfB bits
  {
    B += maskLastBit;           // round up
  }

  // check for carry (same as in encode(..) )
  if( B >= carry )
  {
    putInvOutstanding(); // output inverted outstanding bits (except last 0) and update counter c
    B  -= carry;         // remove carry bit from register --> Bx is now a (U+V) integer
  }

  // output all outstanding bits and first numBitsOfB of B
  putOutstanding ();                // output all outstanding bits
  putFirstOfUVInt( B, numBitsOfB ); // output first numBitsOfB bits of (U+V)-bit integer B

                                    // reset arithmetic coder
  A = ( uint64_t(1) << U ) - 1;
  B = 0;
  c = 0;
}





//===== ARITHMETIC DECODER =====
ArithDec::ArithDec ( unsigned Ubits, unsigned Vbits, IBitstream& inBS )
  : U    ( Ubits )
  , V    ( Vbits )
  , bs   ( inBS )
  , A    ( ( uint64_t(1) << U ) - 1 )
  , v    ( 0 )
{
  assert( U+V < 64 );
  for( unsigned k = 0; k < U+V; k++ ) // fill with first U+V bits from bitstream
  {
    v = ( v << 1 ) + bs.getBit();
  }
}


ArithDec::~ArithDec()
{}


// get number of leading zeros in (U+V)-bit integer
unsigned ArithDec::leadingZerosUV( uint64_t value )
{
  unsigned num = 0;
  for( unsigned bitPos = U+V-1; ; bitPos--, num++ )
  {
    if( value >> bitPos )
      return num;
  }
}


//=====  MAIN DECODING FUNCTION  =====
unsigned ArithDec::decode( IPmf& pmf )
{
  uint64_t Uk = 0;
  for( unsigned index = 0; ; index++ ) // loop over symbol alphabet
  {
    uint64_t Ak = A * pmf[ index ];    // candidate for new interval width
    uint64_t Ck = Uk;                  // lower boundary = last upper boundary
    Uk         += Ak;                  // new upper boundary
    if( v < Uk )
    {
      // symbol found -> proceed with update
      unsigned lz = leadingZerosUV( Ak ); // number of leading zeros in (U+V)-bit representation of Ak
      A           = ( Ak << lz ) >> V;    // update interval width
      v          -= Ck;                   // subtract lower interval boundary
      for( unsigned k = 0; k < lz; k++ )  // fill v with next lz bits from bitstream
      {
        v = ( v << 1 ) + bs.getBit();
      }

      // update pmf (for non-adaptive pmf implementations, this has no effect)
      pmf.update( index );
      
      // return index
      return index;
    }
  }

}




