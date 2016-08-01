#ifndef of_util__hpp
#define of_util__hpp

#include "of-scopes.hpp"
#include "of-types.hpp"
#include "of-debug.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <stdexcept>

BEGIN_OLIGOFAR_SCOPES

#define THROW( exception, message ) \
    do { \
        ostringstream __msg__; \
        __msg__ << __FILE__ << ":" << __LINE__ << "\t" << #exception << "\t" << message; \
        throw exception( __msg__.str() ); \
    } while(0)

#define ITERATE( Type, var, container ) \
    for( Type::const_iterator var = container.begin(), _end = container.end(); var != _end; ++var )

enum EOnError { eOnError_throw, eOnError_ignore, eOnError_default = eOnError_throw };

template<class T, class Fun>
T RealFromStr( const char * str, Fun fun, EOnError e ) 
{
    char * x = 0;
    T t = fun( str, &x );
    if( e == eOnError_throw && ( x == 0 || *x != 0 ) ) 
        THROW( invalid_argument, "real" << (sizeof(T)*8) << "_t from \"" << str << "\"" );
    return t;
}
template<class T, class Fun>
T IntFromStr( const char * str, Fun fun, EOnError e ) 
{
    char * x = 0;
    T t = fun( str, &x, 0 );
    if( e == eOnError_throw && ( x == 0 || *x != 0 ) ) 
        THROW( invalid_argument, "int" << (sizeof(T)*8) << "_t from \"" << str << "\"" );
    return t;
}

template <class Type> Type FromStr( const char * str, EOnError e );
template <class Type> Type FromStr( const char * str ) { return FromStr<Type>( str, eOnError_throw ); }
template <class Type> Type FromStrX( const char * str ) { return FromStr<Type>( str, eOnError_throw ); }

template<> inline string FromStr<string>( const char * str, EOnError ) { return str; }
template<> inline const char * FromStr<const char *>( const char * str, EOnError ) { return str; }
template<> inline double FromStr<double>( const char * str, EOnError e ) { return RealFromStr<double>( str, strtod, e ); }
template<> inline float FromStr<float>( const char * str, EOnError e ) { return RealFromStr<float>( str, strtof, e ); }
template<> inline Int2 FromStr<Int2>( const char * str, EOnError e ) { return IntFromStr<Int2>( str, strtol, e ); }
template<> inline Int4 FromStr<Int4>( const char * str, EOnError e ) { return IntFromStr<Int4>( str, strtol, e ); }
template<> inline Int8 FromStr<Int8>( const char * str, EOnError e ) { return IntFromStr<Int8>( str, strtoll, e ); }
template<> inline Uint2 FromStr<Uint2>( const char * str, EOnError e ) { return IntFromStr<Uint2>( str, strtoul, e ); }
template<> inline Uint4 FromStr<Uint4>( const char * str, EOnError e ) { return IntFromStr<Uint4>( str, strtoul, e ); }
template<> inline Uint8 FromStr<Uint8>( const char * str, EOnError e ) { return IntFromStr<Uint8>( str, strtoull, e ); }
template<class Type> Type FromStr( const string& str ) { return FromStr<Type>( str.c_str() ); }

enum ERealFmt { eReal_any, eReal_fixed, eReal_exponental };
enum EIntBase { eBase_bin = 2, eBase_oct = 8, eBase_dec = 10, eBase_hex = 16 };

inline std::string ToStr( double value, ERealFmt f = eReal_any, int prec = 0 ) 
{
    char buffer[1024]; 
    if( prec ) {
        static const char * fmt[] = { "%.*lg", "%.*lf", "%.*le" };
        std::snprintf( buffer, sizeof( buffer ), fmt[f], prec, value );
    } else {
        static const char * fmt[] = { "%lg", "%lf", "%le" };
        std::snprintf( buffer, sizeof( buffer ), fmt[f], value );
    }
    return buffer;
}

inline std::string ToStr( float value, ERealFmt f = eReal_any, int prec = 0 ) { return ToStr( (double)value, f, prec ); }

inline std::string ToStr( Uint8 value, int base = eBase_dec, int fill = 0 ) 
{
    char buffer[1024];
    char * x = buffer + sizeof(buffer);
    *--x = 0;
    char * y = x;
    if( fill > sizeof( buffer ) - 5 ) fill = sizeof( buffer ) - 5;
    static const char * digits = "0123456789abcdefghijklmnopqrstuvwxyz";
    if( value == 0 ) *--x = '0'; 
    else {
        if( base < 2 ) base = 2;
        if( base > 36 ) base = 36;
        while( value ) {
            int d = value%base;
            value /= base;
            *--x = digits[d];
        }
    }
    while( y - x < fill ) *--x = '0';
    switch( base ) {
        case eBase_bin: *--x = 'b'; *--x = '0'; break;
        case eBase_hex: *--x = 'x'; 
        case eBase_oct: *--x = '0'; 
        case eBase_dec: break;
        default: *--x = ')'; *--x = digits[base-1]; *--x = '('; break;
    }
    return x;
}

inline std::string ToStr( Int8 value, int base = eBase_dec, int fill = 0 ) 
{
    if( value < 0 ) return "-" + ToStr( Uint8( -value ), base, fill );
    else return ToStr( Uint8( value ), base, fill );
}
inline std::string ToStr( Uint2 value, int base = eBase_dec, int fill = 0 ) { return ToStr( Uint8( value ), base, fill ); }
inline std::string ToStr( Uint4 value, int base = eBase_dec, int fill = 0 ) { return ToStr( Uint8( value ), base, fill ); }
inline std::string ToStr( Int2 value, int base = eBase_dec, int fill = 0 ) { return ToStr( Int8( value ), base, fill ); }
inline std::string ToStr( Int4 value, int base = eBase_dec, int fill = 0 ) { return ToStr( Int8( value ), base, fill ); }

inline pair<int,int> ParseRange( const char * str, const char * delim ) 
{
    const char * x = 0;
    pair<int,int> ret(0,0);
    ret.first = strtol( str, const_cast<char**>( &x ), 10 );
    if( x == 0 || x == str )
        THROW( runtime_error, "Integer or integer pair expected, got [" << str << "]" );
    if( *x == 0 ) { ret.second = ret.first; return ret; }
    if( strchr( delim, *x ) ) {
        const char * s = x+1;
        ret.second = strtol( s, const_cast<char**>( &x ), 10 );
        if( x != 0 && x != s && *x == 0 ) 
            return ret;
    }
    THROW( runtime_error, "Integer or integer pair expected, got trailing characters (" << x << ") in [" << str << "]" );
}

class CFromStr
{
public:
    typedef const char * cptr;
    typedef const string& cstrref;
    CFromStr( const char * str ) : m_str( str ) {}
    template<class To>
    operator To () const { return FromStr<To>( m_str ); }
    operator cptr () const { return m_str; }
protected:
    const char * m_str;
};

template<class iterator, class converter> 
inline iterator Split(const string& str, const string& delim,
                      iterator i, const converter& conv, bool tokenize = true) 
{
	const char * c = str.c_str();
	int cnt = 0;
    char Delim[256]; 
    memset( Delim, 0, 256 ); 
    for( const char * d = delim.c_str(); *d; ++d ) Delim[(unsigned char)*d] =1;
	while( const char * cc = strpbrk( c, delim.c_str() ) ) {
		if( cc > c || !tokenize ) { string s( c, cc ); *i++ = conv( s.c_str() );  }
		c = cc + 1;
		if( tokenize ) while( *c && Delim[(unsigned char)*c] ) ++c;
		++cnt;
	}
    *i++ = conv( c );
	return i;
}

template<class iterator> 
inline iterator Split(const string& str, const string& delim,
                      iterator i, bool tokenize = true) 
{ return Split( str, delim, i, FromStr<const char *>, tokenize ); }

template<class container>
inline string Join(const string& delim, const container& c) 
{
	ostringstream ret;
	typedef typename container::const_iterator iterator;
	for( iterator q = c.begin(); q != c.end(); ++q ) {
		if( q != c.begin() ) ret << delim;
		ret << *q;
	}
	return ret.str();
}

template<class iterator>
inline string Join(const string& delim, iterator b, iterator e) 
{
	ostringstream ret;
	for( iterator i = b; i != e; ++i ) {
		if( i != b ) ret << delim;
		ret << *i;
	}
	return ret.str();
}

END_OLIGOFAR_SCOPES

#endif
