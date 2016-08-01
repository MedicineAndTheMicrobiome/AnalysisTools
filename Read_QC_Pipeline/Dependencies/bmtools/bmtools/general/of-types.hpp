#ifndef of_types__hpp
#define of_types__hpp

#include "of-scopes.hpp"
#include <stdint.h>

BEGIN_OLIGOFAR_SCOPES

#ifdef STANDALONE
typedef int8_t Int1;
typedef int16_t Int2;
typedef int32_t Int4;
typedef int64_t Int8;
typedef uint8_t Uint1;
typedef uint16_t Uint2;
typedef uint32_t Uint4;
typedef uint64_t Uint8;
#endif

END_OLIGOFAR_SCOPES

#endif
