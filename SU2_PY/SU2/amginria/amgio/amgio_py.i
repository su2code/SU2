/*  Example of wrapping cos function from math.h using SWIG. */

%module amgio
%{
    /* the resulting C file should be built as a python extension */
    #define SWIG_FILE_WITH_INIT
    /*  Includes the header in the wrapper code */
	#include "amgio_py.h" 
%}

/*  Parse the header file to generate wrappers */
%include "amgio_py.h"

