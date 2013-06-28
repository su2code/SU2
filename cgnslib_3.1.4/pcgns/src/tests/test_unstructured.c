//! @file test_unstructured.c
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
	int S;
	char basename[100];
	char zonename[100];
	int cell_dim = 3;
	int phys_dim = 3;
	int nijk[3][1];
	int k;

	nijk[0][0] = 10;
	nijk[1][0] = nijk[0][0]-1;
	nijk[2][0] = 0;

	err = MPI_Init(&argc,&argv);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	if(err!=MPI_SUCCESS) cgp_doError;
	err = MPI_Info_create(&(info));
	if(err!=MPI_SUCCESS) cgp_doError;

	err = cgp_open("test_unstructured.cgns", 0, MPI_COMM_WORLD, &info, &fn);
	if(err!=0) cgp_doError;
	err = cgp_base_write(fn, "Base 1", cell_dim, phys_dim, &B);
	if(err!=0) cgp_doError;
	err = cgp_zone_write(fn, B, "Zone 1", &(nijk[0][0]), Unstructured, &Z);
	if(err!=0) cgp_doError;

	double x[10/comm_size];
	double y[10/comm_size];
	double z[10/comm_size];

	int min[1] = {10/comm_size*comm_rank};
	int max[1] = {10/comm_size*(comm_rank+1)-1};

	for(k=0;k<10/comm_size;k++) {
		x[k] = (double) (min[0]+k);
		y[k] = 0.0;
		z[k] = 0.0;
		}

	int Cx,Cy,Cz;

	err = cgp_coord_write(fn,B,Z,Unstructured,"CoordinateX",&Cx);
	//~ err = cgp_coord_write_data(fn,B,Z,Cz,min,max,x);

	err = cgp_coord_write(fn,B,Z,Unstructured,"CoordinateY",&Cy);
	//~ err = cgp_coord_write_data(fn,B,Z,Cy,min,max,y);

	err = cgp_coord_write(fn,B,Z,Unstructured,"CoordinateZ",&Cz);
	//~ err = cgp_coord_write_data(fn,B,Z,Cz,min,max,z);

	err = queue_slice_write(Coords, fn, B, Z, &Cx, 1, min, max, x);
	err = queue_slice_write(Coords, fn, B, Z, &Cy, 1, min, max, y);
	err = queue_slice_write(Coords, fn, B, Z, &Cz, 1, min, max, z);
	err = queue_flush();

	int start = 1;
	int end = 9;
	err = cgp_section_write(fn,B,Z,"Elements",BAR_2,start,end,0,&S);

	int nelems = (comm_rank!=comm_size-1)?9/comm_size:9-(9/comm_size)*(comm_size-1);
	printf("%d:%d\n",comm_rank,nelems);
	int emin = (9/comm_size)*comm_rank+1;
	int emax = (comm_rank!=comm_size-1)?(9/comm_size)*(comm_rank+1):9;
	int elements[nelems*2];
	for(k=0;k<nelems;k++) {
		elements[2*k] = k+emin;
		elements[2*k+1] = k+emin+1;
		}
	printf("%d:%d %d %d\n",comm_rank,nelems,emin,emax);
   //~ err =  cgp_section_write_data(fn,B,Z,S,emin,emax,&(elements[0]));

	err = queue_slice_write(Elements, fn, B, Z, &S, 1, &emin, &emax, elements);
	err = queue_flush();

	err = cgp_close(fn);
	if(err!=0) cgp_doError;

	err = MPI_Finalize();
	if(err!=MPI_SUCCESS) cgp_doError;
	return 0;
	}
