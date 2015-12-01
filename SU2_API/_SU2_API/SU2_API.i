// SWIG input file of the 'SU2Solver' API module

%feature("autodoc","1");

%module(docstring=
"'SU2Solver' module",
directors="1",
threads="1"
) SU2Solver
%{

#include "SU2_API.h"

%}

// ----------- USED MODULES ------------
%include "std_string.i"

// ----------- API CLASSES ----------------
%include "SU2_API.h"
