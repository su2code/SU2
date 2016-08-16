// SWIG input file of the 'SU2Solver' API module

%feature("autodoc","1");

%module(docstring=
"'SU2Solver' module",
directors="1",
threads="1"
) SU2Solver
%{

#include "../../SU2_CFD/include/driver_structure.hpp"

%}

// ----------- USED MODULES ------------
%import "../../Common/include/datatypes/primitive_structure.hpp"
%include "std_string.i"
%include "typemaps.i"
//%include "numpy.i"
#ifdef HAVE_MPI			//Need mpi4py only for a parallel build of the wrapper.
  %include "mpi4py/mpi4py.i"
  %mpi4py_typemap(Comm, MPI_Comm)
#endif
// ----------- API CLASSES ----------------
%include "../../SU2_CFD/include/driver_structure.hpp"
