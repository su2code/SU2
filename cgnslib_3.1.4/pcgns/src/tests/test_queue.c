//! @file test_queue.c
//! @author Kyle Horne <horne.kyle@gmail.com>
//! @version 0.2
//!
//! @section LICENSE
//! BSD style license
//!
//! @section DESCRIPTION
//! Test program for pcgns library

#include "pcgnslib.h"

#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"

int main(int argc, char* argv[]) {
	int err;
	int comm_size;
	int comm_rank;
	MPI_Info info;
	int fn;
	int B;
	int Z;
	char basename[100];
	char zonename[100];
	int cell_dim = 3;
	int phys_dim = 3;
	int nijk[3][3];

	nijk[0][0] = 10;
	nijk[0][1] = 10;
	nijk[0][2] = 1;
	nijk[1][0] = nijk[0][0]-1;
	nijk[1][1] = nijk[0][1]-1;
	nijk[1][2] = nijk[0][2]-1;
	nijk[2][0] = 0;
	nijk[2][1] = 0;
	nijk[2][2] = 0;

	err = MPI_Init(&argc,&argv);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Info_create(&(info));
	if(err!=MPI_SUCCESS) cgp_doError;

	err = cgp_open("test_queue.cgns", 0, MPI_COMM_WORLD, &info, &fn);
	if(err!=0) cgp_doError;
	err = cgp_base_write(fn, "Base 1", cell_dim, phys_dim, &B);
	if(err!=0) cgp_doError;
	err = cgp_zone_write(fn, B, "Zone 1", &(nijk[0][0]), Structured, &Z);
	if(err!=0) cgp_doError;
	err = cgp_zone_read(fn, B, Z, zonename, &(nijk[0][0]));
	if(err!=0) cgp_doError;

	int min = 0;
	int max = 0;
	double data = 0.0;
	err = queue_slice_write(Empty, fn, B, Z, NULL,1, &min, &max, &data);
	err = queue_flush();

	err = cgp_close(fn);
	if(err!=0) cgp_doError;

	err = MPI_Finalize();
	if(err!=MPI_SUCCESS) cgp_doError;
	return 0;
	}
