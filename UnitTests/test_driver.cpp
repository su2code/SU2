/*--- This tells Catch to provide a main().  This should only be done in
 * one cpp file. ---*/
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "../../../Common/include/mpi_structure.hpp"
#include "../../../Common/include/option_structure.hpp"

int main(int argc, char *argv[]) {

  /*--- Startup MPI, if supported ---*/
#ifdef HAVE_MPI
  int  buffsize;
  char *buffptr;
#ifdef HAVE_OMP
  int provided;
  SU2_MPI::Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
#else
  SU2_MPI::Init(&argc, &argv);
#endif
  SU2_MPI::Buffer_attach( malloc(BUFSIZE), BUFSIZE );
  SU2_Comm MPICommunicator(MPI_COMM_WORLD);
#else
  SU2_Comm MPICommunicator(0);
#endif

  int result = Catch::Session().run(argc, argv);

  /*--- Finalize MPI parallelization ---*/
#ifdef HAVE_MPI
  SU2_MPI::Buffer_detach(&buffptr, &buffsize);
  free(buffptr);
  SU2_MPI::Finalize();
#endif

  return result;
}
