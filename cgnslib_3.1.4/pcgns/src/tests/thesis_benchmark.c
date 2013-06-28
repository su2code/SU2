#include "pcgnslib.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mpi.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

int comm_size;
int comm_rank;
MPI_Info info;

double data_size;
int N;
int Nl;
int pc;
int zpp;
int ppz;
int zc;

int* zones;
int* subzones;

double* x;
double* y;
double* z;

int* e;

double* u;
double* v;
double* w;

double* h;

int read_inputs(int* argc, char*** argv) {
	int k;
	if(comm_rank==0) {
		if(*argc<7) exit(1);
		for(k=1;k<*argc;k++) {
			if(strcmp((*argv)[k],"-ds")==0) {
				k++;
				sscanf((*argv)[k],"%lf",&data_size);
				printf("data_size=%lf\n",data_size);
				}
			if(strcmp((*argv)[k],"-zpp")==0) {
				k++;
				sscanf((*argv)[k],"%d",&zpp);
				printf("zpp=%d\n",zpp);
				}
			if(strcmp((*argv)[k],"-ppz")==0) {
				k++;
				sscanf((*argv)[k],"%d",&ppz);
				printf("ppz=%d\n",ppz);
				}
			}
		}
	MPI_Bcast(&data_size,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&zpp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&ppz,1,MPI_INT,0,MPI_COMM_WORLD);
	return 0;
	}

int initialize(int* argc, char*** argv) {
	int j,k;

	MPI_Init(argc,argv);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
	MPI_Info_create(&info);

	read_inputs(argc,argv);

	N = (data_size*1024*1024)/((double) sizeof(double));
	pc = comm_size;
	zc = (pc*zpp)/ppz;
	Nl = N/comm_size/zpp;

	zones = malloc(zc*sizeof(int));
	subzones = malloc(zpp*sizeof(int));
	for(k=0;k<zpp;k++) {
		zones[k] = comm_rank/ppz*zpp+k;
		subzones[k] = comm_rank%ppz;
		}

	// Initialize Arrays
	x = malloc(Nl*sizeof(double));
	y = malloc(Nl*sizeof(double));
	z = malloc(Nl*sizeof(double));

	e = malloc(Nl*sizeof(int));

	u = malloc(Nl*sizeof(double));
	v = malloc(Nl*sizeof(double));
	w = malloc(Nl*sizeof(double));

	h = malloc(Nl*sizeof(double));

	double theta;
	double r;
	for(k=0;k<Nl;k++) {
		j = Nl*subzones[0]+k;
		theta = ((double) j)/((double) Nl*zpp);
		r = theta;
		x[k] = r*cos(theta);
		y[k] = r*sin(theta);
		z[k] = r;
		e[k] = j+1;
		u[k] = x[k];
		v[k] = y[k];
		w[k] = z[k];
		h[k] = r;
		}

	//~ printf("%d: Nl %d\n",comm_rank,Nl);
	//~ for(k=0;k<zpp;k++) printf("%d: Z%d.%d\n",comm_rank,zones[k],subzones[k]);

	return 0;
	}

int finalize(void) {
	free(zones);
	free(subzones);

	free(x);
	free(y);
	free(z);
	free(e);
	free(u);
	free(v);
	free(w);

	MPI_Finalize();
	return 0;
	}

int main(int argc, char* argv[]) {
	int k;
	int F;
	int B;

	int nijk[3][1];

	initialize(&argc,&argv);

	cgp_open("thesis_benchmark.cgns",0,MPI_COMM_WORLD, &info, &F);
	cgp_base_write(F,"Base",3,3,&B);

	nijk[0][0] = Nl*ppz;
	nijk[1][0] = Nl*ppz;
	nijk[2][0] = 0;

	default_pio_mode = CGP_INDEPENDENT;

	int Z[zc];
	int Cx[zc];
	int Cy[zc];
	int Cz[zc];
	int E[zc];
	int Su[zc];
	int Sv[zc];
	int Sw[zc];

	for(k=0;k<zc;k++) {
		char zonename[100+1];
		sprintf(zonename,"%s %d","Zone",k);
		cgp_zone_write(F,B,zonename,&(nijk[0][0]),Unstructured,&(Z[k]));
		cgp_coord_write(F,B,Z[k],0,"CoordinateX",&(Cx[k]));
		cgp_coord_write(F,B,Z[k],0,"CoordinateY",&(Cy[k]));
		cgp_coord_write(F,B,Z[k],0,"CoordinateZ",&(Cz[k]));
		cgp_section_write(F,B,Z[k],"Elements",NODE,1,Nl*ppz,0,&(E[k]));
		cgp_sol_write(F,B,Z[k],"MomentumX",Vertex,&(Su[k]));
		cgp_sol_write(F,B,Z[k],"MomentumY",Vertex,&(Sv[k]));
		cgp_sol_write(F,B,Z[k],"MomentumZ",Vertex,&(Sw[k]));
		cgp_array_write(F,B,Z[k],"phi",Vertex);
		}

	double T0,T1;
	double Tw0,Tw1;
	double Tr0,Tr1;
	double t0,t1;

	MPI_Barrier(MPI_COMM_WORLD);
	T0 = MPI_Wtime();
	Tw0 = MPI_Wtime();

	// Writes

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		cgp_coord_write_data(F,B,Z[zones[k]],Cx[zones[k]],&min,&max,x);
		cgp_coord_write_data(F,B,Z[zones[k]],Cy[zones[k]],&min,&max,y);
		cgp_coord_write_data(F,B,Z[zones[k]],Cz[zones[k]],&min,&max,z);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Coords Write\n");
		printf("\tTime=%lf\n",comm_rank,t1-t0);
		printf("\tBandwidth=%lf\n",comm_rank,3.0*data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		cgp_sol_write_data(F,B,Z[zones[k]],Su[zones[k]],&min,&max,u);
		cgp_sol_write_data(F,B,Z[zones[k]],Sv[zones[k]],&min,&max,v);
		cgp_sol_write_data(F,B,Z[zones[k]],Sw[zones[k]],&min,&max,w);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Solutions Write\n");
		printf("\tTime=%lf\n",t1-t0);
		printf("\tBandwidth=%lf\n",3.0*data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		cgp_array_write_data(F,B,Z[zones[k]],"phi",&min,&max,h);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Arrays Write\n");
		printf("\tTime=%lf\n",t1-t0);
		printf("\tBandwidth=%lf\n",data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		min++;
		max++;
		cgp_section_write_data(F,B,Z[zones[k]],E[zones[k]],min,max,e);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Elements Write\n");
		printf("\tTime=%lf\n",t1-t0);
		printf("\tBandwidth=%lf\n",((double) sizeof(int))/((double) sizeof(double))*data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	Tw1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Total Write Time=%lf\n",Tw1-Tw0);
		printf("Total Write Bandwidth=%lf\n",(6.0+((double) sizeof(int))/((double) sizeof(double)))*data_size/(Tw1-Tw0));
		}

	//=======//
	//=Reads=//
	//=======//

	MPI_Barrier(MPI_COMM_WORLD);
	Tr0 = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		cgp_coord_read_data(F,B,Z[zones[k]],Cx[zones[k]],&min,&max,x);
		cgp_coord_read_data(F,B,Z[zones[k]],Cy[zones[k]],&min,&max,y);
		cgp_coord_read_data(F,B,Z[zones[k]],Cz[zones[k]],&min,&max,z);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Coords Read\n");
		printf("\tTime=%lf\n",comm_rank,t1-t0);
		printf("\tBandwidth=%lf\n",comm_rank,3.0*data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		cgp_sol_read_data(F,B,Z[zones[k]],Su[zones[k]],&min,&max,u);
		cgp_sol_read_data(F,B,Z[zones[k]],Sv[zones[k]],&min,&max,v);
		cgp_sol_read_data(F,B,Z[zones[k]],Sw[zones[k]],&min,&max,w);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Solutions Read\n");
		printf("\tTime=%lf\n",t1-t0);
		printf("\tBandwidth=%lf\n",3.0*data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		cgp_array_read_data(F,B,Z[zones[k]],"phi",&min,&max,h);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Arrays Read\n");
		printf("\tTime=%lf\n",t1-t0);
		printf("\tBandwidth=%lf\n",data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for(k=0;k<zpp;k++) {
		int min = subzones[k]*Nl;
		int max = (subzones[k]+1)*Nl-1;

		min++;
		max++;
		cgp_section_read_data(F,B,Z[zones[k]],E[zones[k]],min,max,e);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Elements Read\n");
		printf("\tTime=%lf\n",t1-t0);
		printf("\tBandwidth=%lf\n",((double) sizeof(int))/((double) sizeof(double))*data_size/(t1-t0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	Tr1 = MPI_Wtime();
	if(comm_rank==0) {
		printf("Total Read Time=%lf\n",Tr1-Tr0);
		printf("Total Read Bandwidth=%lf\n",(6.0+((double) sizeof(int))/((double) sizeof(double)))*data_size/(Tr1-Tr0));
		}

	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	cgp_close(F);
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(comm_rank==0) printf("Close_Time=%lf\n",t1-t0);

	MPI_Barrier(MPI_COMM_WORLD);
	T1 = MPI_Wtime();

	if(comm_rank==0) {
		printf("Total Time=%lf\n",T1-T0);
		}

	finalize();

	return 0;
	}
