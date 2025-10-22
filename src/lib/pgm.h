
#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>


class PGMImage
{
public:
  typedef int16_t Sample;

public:
  PGMImage();
  PGMImage( int _width, int _height );
  PGMImage( int _width, int _height, const std::vector<Sample>& _data );

  int                         getWidth  ()  const { return m_width; }
  int                         getHeight ()  const { return m_height; }
  int                         getSize   ()  const { return m_width*m_height; }
  const std::vector<Sample>&  getData   ()  const { return m_data; }
  std::vector<Sample>&        getData   ()        { return m_data; }

  bool  read  ( const std::string& _filename );
  bool  write ( const std::string& _filename )  const;

private:
  const std::string   m_header    = "P5";
  const Sample        m_maxValue  = 255;
  int                 m_width;
  int                 m_height;
  std::vector<Sample> m_data;
};

