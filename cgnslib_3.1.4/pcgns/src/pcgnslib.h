//! @file pcgnslib.h
//! @author Kyle Horne <horne.kyle@gmail.com>
//! @version 0.2
//!
//! @section LICENSE
//! BSD style license
//!
//! @section DESCRIPTION
//! Header file for all public functions of the pcgns library

#ifndef PCGNSLIB_H_
#define PCGNSLIB_H_

#ifdef _DEBUG
#define cgp_doError {printf("Error at %s:%u\n",__FILE__, __LINE__); return 1;}
#else
#define cgp_doError ;
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

//#define printTime printf("Time at %s:%u = %f\n",__FILE__,__LINE__,MPI_Wtime())
#define printTime

#include "mpi.h"
#include "hdf5.h"

//=====================//
//== Begin Datatypes ==//
//=====================//

typedef enum {Structured, Unstructured} ZoneType_t;
typedef enum {
	ElementTypeNull, ElementTypeUserDefined, // 0 1
	NODE,                                    // 2
	BAR_2, BAR_3,                            // 3 4
	TRI_3, TRI_6,                            // 5 6
	QUAD_4, QUAD_8, QUAD_9,                  // 7 8 9
	TETRA_4, TETRA_10,                       // 10 11
	PYRA_5, PYRA_13, PYRA_14,                // 12 13 14
	PENTA_6, PENTA_15, PENTA_18,             // 15 16 17
	HEXA_8, HEXA_20, HEXA_27,                // 18 19 20
	MIXED,                                   // 21
	NGON_n, NFACE_n                          // 22 23
	} ElementType_t;
typedef int DataType_t;
typedef enum {Vertex, CellCenter} GridLocation_t;
typedef enum {Empty, Coords, Elements, Solution, Array} SliceType_t;

extern int preallocate;

enum {CGP_COLLECTIVE=TRUE,CGP_INDEPENDENT=FALSE};

extern int default_pio_mode;

//===============================//
//== Begin Function Prototypes ==//
//===============================//

//= File IO Prototypes =//
//! Open a file for reading and writing
//! @param filename [in]: Name of the file to open
//! @param mode [in]: IO mode (read/write)
//! @param comm [in]: MPI communicator on which to open the file
//! @param info [in]: MPI info object to allow hints passed to MPI-IO
//! @param fn [out]: Handle of the opened file
//! @return Error code
int cgp_open(const char* filename, int mode, MPI_Comm comm, MPI_Info* info, int* fn);

//! Close a previously opened file
//! @param fn [in]: Handle of the file to close
//! @return Error code
int cgp_close(int fn);

//= Base IO Prototypes =//
//! Read info about a base
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param basename [out]: Name of the base
//! @param cell_dim [out]: Cell dimensions of the base
//! @param phys_dim [out]: Physical dimensions of the base
//! @return Error code
int cgp_base_read(int fn, int B, char* basename, int* cell_dim, int* phys_dim);

//! Write a base to a file
//! @param fn int[in]: Handle f the file
//! @param basename [in]: Name of the base to write
//! @param cell_dim [in]: Cell dimensions of the base
//! @param phys_dim [in]: Physical dimensions of the base
//! @param B [out]: Index of the base
//! @return Error code
int cgp_base_write(int fn, char const* basename, int cell_dim, int phys_dim, int* B);

//! Read the number of bases in a file
//! @param fn [in]: Handle of the file
//! @param nbases [out]: Number of bases in the specified file
//! @return Error code
int cgp_nbases(int fn, int *nbases);

//= Zone IO Prototypes =//
//! Read info about a zone
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone to read
//! @param zonename [out]: Name of the zone
//! @param nijk [out]: Dimensions of the zone
//! @return Error code
int cgp_zone_read(int fn, int B, int Z, char* zonename, int* nijk);

//! Read the type of a zone
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param zonetype [out]: Type of zone
//! @return Error code
int cgp_zone_type(int fn, int B, int Z, ZoneType_t *zonetype);

//! Write a zone to a base
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param zonename [in]: Name of the zone to write
//! @param nijk [in]: Dimensions of the zone
//! @param type [in]: Type of zone
//! @param Z [out]: Index of the zone
//! @return Error code
int cgp_zone_write(int fn, int B, const char* zonename, const int* nijk, ZoneType_t type, int* Z);

//! Read the number of zones in a base
//! @param fn [in]: Handle of file
//! @param B [in]: Index of base
//! @param nzones [out]: Number of zones in the specified base
//! @return Error code
int cgp_nzones(int fn, int B, int *nzones);

//= Grid IO Prototypes =//

//! Write coords group, but not data, to a grid
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param datatype [in]: Type of floats stored
//! @param coordname [in]: Name of the coords
//! @param C [out]: Index of the coords
int cgp_coord_write(int fn, int B, int Z, DataType_t type, const char* coordname, int* C);

//! Write coords to a grid in parallel
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param C [in]: Index of the coords
//! @param range_min [in]: Array of lower bound index
//! @param range_max [in]: Array of upper bound index
//! @param coord_array [in]: Pointer to the data
//! @return Error code
int cgp_coord_write_data(int fn, int B, int Z, int C, int* min, int* max, void* coord_array);

//! Write coords to a grid in parallel
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param C [in]: Index of the coords
//! @param range_min [in]: Array of lower bound index
//! @param range_max [in]: Array of upper bound index
//! @param coord_array [out]: Pointer to the data
//! @return Error code
int cgp_coord_read_data(int fn, int B, int Z, int C, int* min, int* max, void* coord_array);

//= Solution IO Prototypes =//
//! Write a solution to a zone
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param solname [in]: Name of solution
//! @param location [in]: Location of solution within each cell
//! @param S [out]: Index of solution
//! @return Error code
int cgp_sol_write(int fn, int B, int Z, char *solname, GridLocation_t location, int *S);

//! Write a solution's data in parallel
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param S [in]: Index of soltution
//! @param min [in]: Lower bound array for data
//! @param max [in]: Upper bound array for data
//! @param data [in]: Data to be written
//! @return Error code
int cgp_sol_write_data(int fn, int B, int Z, int S, int* min, int* max, void* data);

//! Write a solution's data in parallel
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param S [in]: Index of soltution
//! @param min [in]: Lower bound array for data
//! @param max [in]: Upper bound array for data
//! @param data [out]: Data to be written
//! @return Error code
int cgp_sol_read_data(int fn, int B, int Z, int S, int* min, int* max, void* data);

//! Read the number of solutions in a zone
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param nsols [out]: Number of solutions in the specified zone
//! @return Error code
int cgp_nsols(int fn, int B, int Z, int* nsols);

//= Unstructured Grid Prototypes =//
//! Write the element connectivity groups for a section
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param C [in]: Index of the coords
//! @param sectionname [in]: Name of element section
//! @param type [in]: Type of element data
//! @param start [in]: Element lower bound index
//! @param end [in]: Element upper bound index
//! @param nbndry [in]: Number of boundary elements (unused)
//! @param S [out]: Section index
//! @return Error code
int cgp_section_write(int fn, int B, int Z, char* sectionname, ElementType_t type,
	int start, int end, int nbndry, int* S);

//! Write the element connectivity data for a section
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param C [in]: Index of the coords
//! @param S [in]: Section index
//! @param min [in]: Output array lower bound index
//! @param max [in]: Output array  upper bound index
//! @param elements [in]: Pointer to the data
//! @return Error code
int cgp_section_write_data(int fn, int B, int Z, int S, int min, int max, int *elements);

//! Write the element connectivity data for a section
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param C [in]: Index of the coords
//! @param S [in]: Section index
//! @param min [in]: Output array lower bound index
//! @param max [in]: Output array  upper bound index
//! @param elements [out]: Pointer to the data
//! @return Error code
int cgp_section_read_data(int fn, int B, int Z, int S, int min, int max, int *elements);

//= Array IO Prototypes =//
//! Write an array to a zone
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param arrayname [in]: Name of array
//! @param location [in]: Location of solution within each cell
//! @return Error code
int cgp_array_write(int fn, int B, int Z, char *arrayname, GridLocation_t location);

//! Write an array's data in parallel
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param arrayname [in]: Name of array
//! @param min [in]: Lower bound array for data
//! @param max [in]: Upper bound array for data
//! @param data [in]: Data to be written
//! @return Error code
int cgp_array_write_data(int fn, int B, int Z, char* arrayname, int* min, int* max, void* data);

//! Read an array's data in parallel
//! @param fn [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param arrayname [in]: Name of array
//! @param min [in]: Lower bound array for data
//! @param max [in]: Upper bound array for data
//! @param data [out]: Data to be written
//! @return Error code
int cgp_array_read_data(int fn, int B, int Z, char* arrayname, int* min, int* max, void* data);

//= Queue IO Prototypes =//
//! Queue an IO write operation for flushing later
//! @param type [in]: Type of operation to queue
//! @param F [in]: Handle of the file
//! @param B [in]: Index of the base
//! @param Z [in]: Index of the zone
//! @param SN [in]: Pointer to array locator, which is an int for coordinates, solutions and sections, but a string for arrays
//! @param rank [in]: Rank of data to be written
//! @param min [in]: Pointer to the minumum location array
//! @param max [in]: Pointer to the maximum location array
//! @param data [in]: Pointer to the data to be written
//! @return Error code
int queue_slice_write(SliceType_t type, int F, int B, int Z, void* SN, int rank,
	int* min, int* max, void* data);

//! Flush all the IO operations waiting in the queue
//! @return Error code
int queue_flush(void);
#endif
