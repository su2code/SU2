#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <cerrno>
#include <sstream>
#include <iostream>
#include <cctype>  // For tolower() and toupper()

#include "./shunting-yard.h"
#include "./shunting-yard-exceptions.h"

/* * * * * Built-in Features: * * * * */

/**
 * Force the Startup classes to run
 * at static initialization time:
 */
#ifdef STATIC_CPARSE_STARTUP
#define __CPARSE_STARTUP __Startup__
#else
#define __CPARSE_STARTUP
#endif

#include "./builtin-features/functions.inc"

#include "./builtin-features/operations.inc"

#include "./builtin-features/reservedWords.inc"

#include "./builtin-features/typeSpecificFunctions.inc"

/**
 * If STATIC_CPARSE_STARTUP is undefined
 * declare a startup function to initialize
 * all built-in features:
 */
#ifndef STATIC_CPARSE_STARTUP

void cparse_startup() {
  builtin_functions::Startup();
  builtin_operations::Startup();
  builtin_reservedWords::Startup();
  builtin_typeSpecificFunctions::Startup();
}

#endif
