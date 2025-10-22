
#pragma once

#include <iostream>


typedef uint8_t Bit;


//===== class for output bitstreams =====
class OBitstream
{
public:
  // constructor and destructor
  OBitstream ( std::ostream& _os );
  ~OBitstream();

  // check for errors
  bool good   () const { return m_outStream.good(); }

  // add single bit to bitstream
  void        addBit    ( const Bit _bit );
  OBitstream& operator<<( const Bit _bit );  // for convenience

  // fill last byte (required at end)
  void flush  ();

  // add bit pattern for representing integer values using fixed-length code
  template< typename T >
  void addFixed( T _value, unsigned _numBits = 8*sizeof(T) )
  {
    while( _numBits-- )
      addBit( Bit( ( _value >> _numBits ) & T(1) ) );
  }

private:
  std::ostream&   m_outStream;
  uint8_t         m_nextByte;
  uint8_t         m_freeBits;
};



//===== class for output bitstreams =====
class IBitstream
{
public:
  // constructor and destructor
  IBitstream ( std::istream& _is ); 
  ~IBitstream();

  // check for eof or errors
  bool eof    () const { return m_eof; }
  bool good   () const { return m_inStream.good(); }

  // get next single bit from bitstream
  Bit         getBit    ();
  IBitstream& operator>>( Bit& _bit );  // for convenience

  // get bit pattern for reading integer values using fixed-length code
  template< typename T > 
  T getFixed( unsigned _numBits = 8*sizeof(T) )
  {
    T value = T(0);
    while( _numBits-- )
      value = ( value << 1 ) + (T)getBit();
    return value;
  }

private:
  std::istream&   m_inStream;
  uint8_t         m_currByte;
  uint8_t         m_availBits;
  bool            m_eof;
};

