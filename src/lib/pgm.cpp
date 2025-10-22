
#include "pgm.h"


PGMImage::PGMImage()
  : m_width ( 0 )
  , m_height( 0 )
{}


PGMImage::PGMImage( int _width, int _height )
  : m_width ( _width )
  , m_height( _height )
  , m_data  ( m_width*m_height, Sample(0) )
{}


PGMImage::PGMImage( int _width, int _height, const std::vector<Sample>& _data )
  : m_width ( _width )
  , m_height( _height )
  , m_data  ( _data )
{
  if( m_width * m_height != (int)m_data.size() ) {
    std::cerr << "ERROR: PGMImage::PGMImage(..): Wrong array dimensions!" << std::endl;
    throw 1;
  }
}


bool PGMImage::read( const std::string& _filename )
{
  std::ifstream stream( _filename, std::ios::in|std::ios::binary );
  if( !stream.is_open() || !stream.good() ) {
    std::cerr << "ERROR: Cannot open file \"" << _filename << "\"" << std::endl;
    return false;
  }

  std::string imheader;
  int         imwidth;
  int         imheight;
  int         immaxval;
  stream >> imheader;
  stream >> imwidth;
  stream >> imheight;
  stream >> immaxval;
  stream.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );    
  if( !stream.good() || imheader != m_header || imwidth <= 0 || imheight <= 0 || immaxval != (int)m_maxValue ) {
    std::cerr << "ERROR: File \"" << _filename << "\" has no correct PGM header" << std::endl;
    return false;
  }

  int byte = 0;
  std::vector<Sample> data;  // temporary array
  for( int y = 0; y < imheight; y++ )
  {
    for( int x = 0; x < imwidth; x++ )
    {
      byte = stream.get();
      data.push_back( Sample( byte ) );
    }
  }
  if( !stream.good() ) {
    std::cerr << "ERROR: Could not correctly read file \"" << _filename << "\"" << std::endl;
    return false;
  }

  m_width   = imwidth;
  m_height  = imheight;
  m_data    = data;
  return true;
}


bool PGMImage::write( const std::string& _filename ) const
{
  std::ofstream stream( _filename, std::ios::out|std::ios::binary );
  if( !stream.is_open() || !stream.good() ) {
    std::cerr << "ERROR: Cannot open file \"" << _filename << "\" for writing" << std::endl;
    return false;
  }

  stream << m_header << std::endl;
  stream << m_width << " " << m_height << std::endl;
  stream << m_maxValue << std::endl;

  int imsize = getSize();
  for( int k = 0; k < imsize; k++ )
  {
    uint8_t byte = (uint8_t)std::max<Sample>( Sample(0), std::min<Sample>( m_maxValue, m_data[k] ) );
    stream.put( byte );
  }
  if( !stream.good() ) {
    std::cerr << "ERROR: Could not correctly write file \"" << _filename << "\"" << std::endl;
    return false;
  }
  return true;
}

