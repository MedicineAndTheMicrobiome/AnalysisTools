#ifndef taxmatrix_cfareader__hpp
#define taxmatrix_cfareader__hpp

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <cctype>
#include <assert.h>

#include "of-scopes.hpp"


BEGIN_OLIGOFAR_SCOPES

class CFaReader
{
public:
    CFaReader( const std::string& name );
    CFaReader( std::istream& in );
    const std::string& GetIds() const { return m_ids; }
    const std::string& GetDefline() const { return m_defline; }
    size_t GetOrdinal() const { return m_ordinal; }
    const std::string& GetContents() const { return m_content; }
    bool next();
    bool good() const { return m_good; }
protected:
    std::auto_ptr<std::ifstream> m_fin;
    std::istream * m_in;
    std::string m_ids;
    std::string m_defline;
    std::string m_content;
    std::string m_buffer;
    size_t m_ordinal;
    bool m_good;
};

inline CFaReader::CFaReader( const std::string& fname ) : m_fin( new std::ifstream( fname.c_str() ) ), m_in( m_fin.get() ), m_ordinal(0) 
{ 
    m_good = std::getline( *m_fin, m_buffer );
    next(); 
}

inline CFaReader::CFaReader( std::istream& in ) : m_in( &in ), m_ordinal(0) 
{ 
    m_good = std::getline( *m_fin, m_buffer );
    next(); 
}

inline bool CFaReader::next() 
{
    if( !m_good ) return false;
    if( m_fin->eof() ) return m_good = false;
    ++m_ordinal;
    assert( m_buffer[0] == '>' );
    const char * x = m_buffer.c_str() + 1;
    while( std::isspace( *x ) ) ++x;
    const char * y = x;
    while( *y && !std::isspace( *y ) ) ++y;
    m_ids = std::string( x, y-x );
    while( std::isspace( *y ) ) ++y;
    m_defline = y;
    while( m_defline.length() && std::isspace( m_defline[m_defline.length()-1] ) ) m_defline.resize( m_defline.length()-1 );
    m_content.clear();
    while( std::getline( *m_fin, m_buffer ) ) {
        if( m_buffer[0] == '>' ) return true;
        size_t x = m_buffer.size();
        while( x && std::isspace( m_buffer[x-1] ) ) --x;
        m_content += m_buffer.substr( 0, x );
    }
    return true;
}
    

END_OLIGOFAR_SCOPES

#endif
