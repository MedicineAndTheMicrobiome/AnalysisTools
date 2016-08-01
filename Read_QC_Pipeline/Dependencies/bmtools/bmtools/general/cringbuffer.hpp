#ifndef OLIGOFAR__CRINGBUFFER__HPP
#define OLIGOFAR__CRINGBUFFER__HPP

#include "of-scopes.hpp"

#include <stdexcept>
#include <algorithm>
#include <cstring>
#include <assert.h>

BEGIN_OLIGOFAR_SCOPES

#define OLIGOFAR__NOINCLUDES
#include "cringbuffer-impl.hpp"
#undef OLIGOFAR__NOINCLUDES

END_OLIGOFAR_SCOPES

#endif
