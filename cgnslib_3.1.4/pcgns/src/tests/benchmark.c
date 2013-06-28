//! @file benchmark.c
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
#include "string.h"
#include "math.h"

#define MEGA_BYTES 256
#define BUF_LENGTH (MEGA_BYTES*1024*1024/sizeof(double))
//#define BUF_LENGTH (9)

#define N ((int) sqrt((double) BUF_LENGTH))

int comm_size;
int comm_rank;
MPI_Info info;
int nijk[3][3];

double* x;
double* y;
double* z;

int min[3];
int max[3];

int fn;
int B;
int Z;
int C;

double t0;
double t1;
double ta;

int initialize(int* argc, char** argv[]) {
	MPI_Init(argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Info_create(&(info));

	nijk[0][0] = N;
	nijk[0][1] = N;
	nijk[0][2] = comm_size;
	nijk[1][0] = nijk[0][0]-1;
	nijk[1][1] = nijk[0][1]-1;
	nijk[1][2] = nijk[0][2]-1;
	nijk[2][0] = 0;
	nijk[2][1] = 0;
	nijk[2][2] = 0;

	x = (double*) malloc(BUF_LENGTH*sizeof(double));
	y = (double*) malloc(BUF_LENGTH*sizeof(double));
	z = (double*) malloc(BUF_LENGTH*sizeof(double));

	int i,j;
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			x[i*N+j] = (double) (i);
			y[i*N+j] = (double) (j);
			z[i*N+j] = (double) (comm_rank);
			}
		}

	min[0] = comm_rank;
	min[1] = 0;
	min[2] = 0;
	max[0] = comm_rank;
	max[1] = N-1;
	max[2] = N-1;


	return 0;
	}

int finalize() {
	free(x);
	free(y);
	free(z);

	MPI_Finalize();

	return 0;
	}

int doTimer(const char* msg, double time) {
	double min;
	double max;
	double avg;
	MPI_Reduce(&time, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&time, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	avg = avg/((double) comm_size);
	if(comm_rank==0) printf("%20s Time = { min: %-20f max: %-20f avg: %-20f} s\n",msg,min,max,avg);
	return 0;
	}

int doBandwidth(const char* msg, double time) {
	double min;
	double max;
	double avg;

	double MB = ((double) BUF_LENGTH*sizeof(double))/(1024.0*1024.0);
	MPI_Reduce(&time, &max, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	max = MB/max;
	MPI_Reduce(&time, &min, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	min = MB/min;
	MPI_Reduce(&time, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	avg = avg/((double) comm_size);
	avg = MB/avg;
	if(comm_rank==0) printf("%20s Band = { min: %-20f max: %-20f avg: %-20f} MB/s (local)\n",msg,min,max,avg);
	return 0;
	}

int doBandwidthAgg(const char* msg, double time) {
	double min;
	double max;
	double avg;

	double MB = ((double) BUF_LENGTH*sizeof(double))/(1024.0*1024.0)*((double) comm_size);
	MPI_Reduce(&time, &max, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	max = MB/max;
	MPI_Reduce(&time, &min, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	min = MB/min;
	MPI_Reduce(&time, &avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	avg = avg/((double) comm_size);
	avg = MB/avg;
	if(comm_rank==0) printf("%20s Band = { min: %-20f max: %-20f avg: %-20f} MB/s (aggregate)\n",msg,min,max,avg);
	return 0;
	}

int main(int argc, char* argv[]) {
	// Initialize varaibles
	initialize(&argc,&argv);

	// Time the creation of a file
	t0 = MPI_Wtime();
	cgp_open("benchmark.cgns", 0, MPI_COMM_WORLD, &info, &fn);
	t1 = MPI_Wtime();
	doTimer("File Open", t1-t0);

	// Time the creation of a base
	t0 = MPI_Wtime();
	cgp_base_write(fn, "Base 1", 3, 3, &B);
	t1 = MPI_Wtime();
	doTimer("Base Write", t1-t0);

	// Time the creation of a zone
	t0 = MPI_Wtime();
	cgp_zone_write(fn, B, "Zone 1", &(nijk[0][0]), 0, &Z);
	t1 = MPI_Wtime();
	doTimer("Zone Write", t1-t0);

	// Time the creation of coordinates X
	t0 = MPI_Wtime();
	cgp_coord_write(fn,B,Z,0,"CoordinateX",&C);
	t1 = MPI_Wtime();
	doTimer("Coord X Write", t1-t0);

	// Time the write speed of coordinates X
	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	cgp_coord_write_data(fn,B,Z,C,min,max,x);
	t1 = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	ta = MPI_Wtime();

	doTimer("Coord X Write Data", t1-t0);
	doBandwidth("Coord X Write Data", t1-t0);
	doBandwidthAgg("Coord X Write Data", ta-t0);

	// Time the creation of coordinates Y
	t0 = MPI_Wtime();
	cgp_coord_write(fn,B,Z,0,"CoordinateY",&C);
	t1 = MPI_Wtime();
	doTimer("Coord Y Write", t1-t0);

	// Time the write speed of coordinates Y
	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	cgp_coord_write_data(fn,B,Z,C,min,max,y);
	t1 = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	ta = MPI_Wtime();

	doTimer("Coord Y Write Data", t1-t0);
	doBandwidth("Coord Y Write Data", t1-t0);
	doBandwidthAgg("Coord Y Write Data", ta-t0);

	// Time the creation of coordinates Z
	t0 = MPI_Wtime();
	cgp_coord_write(fn,B,Z,0,"CoordinateZ",&C);
	t1 = MPI_Wtime();
	doTimer("Coord Z Write", t1-t0);

	// Time the write speed of coordinates Z
	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	cgp_coord_write_data(fn,B,Z,C,min,max,z);
	t1 = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	ta = MPI_Wtime();

	doTimer("Coord Z Write Data", t1-t0);
	doBandwidth("Coord Z Write Data", t1-t0);
	doBandwidthAgg("Coord Z Write Data", ta-t0);

	// Time closing of the file
	t0 = MPI_Wtime();
	cgp_close(fn);
	t1 = MPI_Wtime();
	doTimer("File Close", t1-t0);

	finalize();

	return 0;
	}
