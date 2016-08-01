#ifndef of_debug__hpp
#define of_debug__hpp

#include <stdexcept>
#include <iostream>
#include <sstream>

#include "of-util.hpp"

#define VT100(c) "\x1b["c

#define TRACE( what ) \
    cerr << VT100("31m") << __FILE__ << VT100("32m") << ":" << VT100( "33m" ) << __LINE__ \
    << VT100( "34m" ) << "\t" << __FUNCTION__ \
    << VT100("0m") << "\t" << what \
    << VT100("0m") << "\n"

#define DISPLAY( what ) \
    "\t" << VT100("35m") << #what << VT100("32m") "=" << VT100("36m") << "[" << VT100("37m") << (what) << VT100("36m") << "]" << VT100("0m")
#define DISPLAYh( what ) \
    "\t" << VT100("35m") << #what << VT100("32m") "=" << VT100("36m") << "0x" << hex << VT100("37m") << (what) << dec << VT100("0m")

#define ASSERT( what ) \
    if( ! (what) ) THROW( std::logic_error, "Assertion " << #what << " failed" )

#define TESTVAL(type,expr,val)                                          \
    do {                                                                \
        type a_ = expr;                                                 \
        type b_ = val;                                                  \
        if( a_ != b_ ) {                                                \
            std::ostringstream o;                                       \
            o << "Assertion " << #expr << " == " << #val << " failed:\n"; \
            o << " : " << (a_) << " (" #expr ")\n";                      \
            o << " : " << (b_) << " (" #val ")\n";                      \
            string x_ = o.str();                                        \
            THROW( logic_error, x_ );                                   \
        }                                                               \
    }                                                                   \
    while(0)

#endif
