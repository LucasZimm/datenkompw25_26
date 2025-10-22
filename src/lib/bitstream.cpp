
#include "bitstream.h"


OBitstream::OBitstream( std::ostream& _os ) 
  : m_outStream ( _os )
  , m_nextByte  ( 0 )
  , m_freeBits  ( 8 )
{}

OBitstream::~OBitstream() 
{ 
  // fill last byte (if not done already)
  flush(); 
}

void OBitstream::addBit( const Bit _bit )
{
  // fill byte (order: from msb to lsb)
  m_nextByte = ( m_nextByte << 1 ) + uint8_t( !!_bit );
  
  // write full bytes to ostream
  if( --m_freeBits == 0 ) {
    m_outStream.put( m_nextByte );
    m_freeBits = 8;
  }
}

OBitstream& OBitstream::operator<<( const Bit _bit )
{
  addBit( _bit );
  return *this;
}

void OBitstream::flush()
{
  // add bits equal to 0 until last byte is completed
  while( m_freeBits < 8 )
    addBit( Bit(0) );
  m_outStream.flush();
}



IBitstream::IBitstream( std::istream& _is ) 
  : m_inStream  ( _is )
  , m_currByte  ( 0 )
  , m_availBits ( 0 )
  , m_eof       ( false )
{}

IBitstream::~IBitstream() 
{}

Bit IBitstream::getBit()
{
  // read next byte from istream (if no bits available in current byte)
  if( m_availBits == 0 )
  {
    m_currByte  = 0;
    m_availBits = 8;
    // if eof is reached, set eof flag and set byte equal to 0
    if( !m_eof )
    {
      int byte = m_inStream.get();
      if( byte == std::char_traits<char>::eof() )
        m_eof = true;
      else
        m_currByte = uint8_t( byte );
    }
  }
  // return bits (order: from msb to lsb)
  return Bit( ( m_currByte >> --m_availBits ) & 1 );
}

IBitstream& IBitstream::operator>>( Bit& _bit )
{
  _bit = getBit();
  return *this;
}

