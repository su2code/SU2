//! @file pcgns_util.h
//! @author Kyle Horne {horne.kyle@gmail.com}
//! @version 0.2
//!
//! @section LICENSE
//! BSD style license
//!
//! @section DESCRIPTION
//! Header file for utility functions

#ifndef PCGNS_UTIL_H_
#define PCGNS_UTIL_H_

#include "pcgnslib.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


//=====================//
//== Begin Datatypes ==//
//=====================//

//! Struct to describe a coords node in the file
typedef struct coords_s {
	//! Index of the node
	int idx;
	//! Name of the node
	char coordname[100+1];
	//! HDF5 handle to the node
	hid_t group_id;
	} coords_t;

//! Struct to describe a section node in the file
typedef struct section_s {
	//! Index of the node
	int idx;
	//! Name of the node
	char sectionname[100+1];
	//! HDF5 handle to the node
	hid_t group_id;
	//! Type of elements in the section
	ElementType_t type;
	//! Range of element ids in the section
	int range[2];
	//! HDF5 handle to ElementConnectivity
	hid_t connectivity_id;
	//! HDF5 handle to ElementRange
	hid_t range_id;
	} section_t;

//! Struct to describe a solution node in the file
typedef struct sol_s {
	//! Index of the node
	int idx;
	//! Name of the node
	char solname[100+1];
	//! HDF5 handle to the node
	hid_t group_id;
	} sol_t;

//! Struct to describe a zone node in the file
typedef struct zone_s {
	//! Index of the node
	int idx;
	//! Name of the node
	char zonename[100+1];
	//! Array of coord_t objects which describe the coords in the zone
	coords_t* coords;
	//! Length of the array {coords}
	int ncoords;
	//! Array of section_s objects which describe the element sections in the zone
	section_t* sections;
	//! Length of the array {sections}
	int nsections;
	//! Array of sol_t objects which describe the solutions in the zone
	sol_t* sols;
	//! Length of the array {sols}
	int nsols;
	//! HDF5 handle to the node
	hid_t group_id;
	//! Zone dimensions
	int* nijk;
	//! Zone type
	ZoneType_t type;
	//! HDF5 handle to the grid node
	hid_t grid_id;
	//! HDF5 handle to the flow node
	hid_t flow_id;
	} zone_t;

//! Struct to describe a base node in the file
typedef struct base_s {
	//! Index of the node
	int idx;
	//! Name of the node
	char basename[100+1];
	//! Array of zone_t objects which describe the zones in the base
	zone_t *zones;
	//! Length of the array {zones}
	int nzones;
	//! HDF5 handle to the node
	hid_t group_id;
	//! Cell dimensions
	int cell_dim;
	//! Physical dimensions
	int phys_dim;
	} base_t;

//! Struct to describe a CGNS file
typedef struct file_s {
	//! The index of this file in {files}
	int idx;
	//! The name of this file
	char filename[100+1];
	//! Flag to tell if this file is open
	int isOpen;
	//! Array of base_t objects which describe the bases in the file
	base_t* bases;
	//! Length of the array {bases}
	int nbases;
	//! HDF5 handle to the file
	hid_t file_id;
	//! HDF5 property list of the file
	hid_t plist_id;
	//! Flag to set collective/independent IO
	int collective;
	//! MPI comm on which the file was opened
	MPI_Comm comm;
	//! MPI info
	MPI_Info info;
	//! MPI rank of this process in this comm
	int rank;
	//! MPI size of this comm
	int size;
	} file_t;

typedef struct slice_s {
	SliceType_t type;
	int rank;
	int* min;
	int* max;
	void* data;
	int F;
	int B;
	int Z;
	int Selector;
	char name[100+1];
	} slice_t;

//==================================//
//== Begin Global Data Prototypes ==//
//=================================//

//! Array of file_t's which decribe the open files
extern file_t* files;

//! Internal count of used slots in  {files}
extern int files_count;

//! Internal count of the slots in {files}
extern int files_size;

//! Queue of IO write operations
extern slice_t* write_queue;

//! Length of write queue
extern int write_queue_len;

//===============================//
//== Begin Function Prototypes ==//
//===============================//

//! Return the number of the next availible file_t in {files}
//! @param fn [out]: Index of next file
//! @return Error code
int next_file(int* fn);

//! Free the memory of a file_t object
//! @param file [in]: Pointer the file
//! @return Error code
int free_file(file_t* file);

//! Free the memory of a base_t object
//! @param base [in]: Pointer the base
//! @return Error code
int free_base(base_t* base);

//! Free the memory of a zone_t object
//! @param zone [in]: Pointer the zone
//! @return Error code
int free_zone(zone_t* zone);

//! Free the memory of a coords_t object
//! @param coords [in]: Pointer the coords
//! @return Error code
int free_coords(coords_t* coords);

//! Free the memory of a section_t object
//! @param section [in]: Pointer the section
//! @return Error code
int free_section(section_t* section);

//! Free the memory of a sol_t object
//! @param sol [in]: Pointer the sol
//! @return Error code
int free_sol(sol_t* sol);

//==============================================================================

//! Create a new string with parent {pid}, name {name}, and value {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of string
//! @param value [in]: Contents of string
//! @param len [in]: Length of string in file
//! @return Error code
int new_str(hid_t pid, const char* name, const char* value, int len);

//! Read a string with parent {pid}, name {name}, and value {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of string
//! @param len [in]: Length of string in file
//! @param value [out]: Contents of string
//! @return Error code
int get_str(hid_t pid, const char* name, int len, char* value);

//! Create a new string attribute at location {pid} with {name} with {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of string
//! @param value [in]: Contents of string
//! @param len [in]: Length of string in file
//! @return Error code
int new_str_attb(hid_t pid, const char* name, const char* value, int len);

//! Read a string attribute at location {pid} with {name} with {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of string
//! @param len [in]: Length of string in file
//! @param value [out]: Contents of string
//! @return Error code
int get_str_attb(hid_t pid, const char* name, int len, char* value);

//==============================================================================

//! Create a new integer with parent {pid}, name {name}, and value {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of integer
//! @param value [in]: Integer to write
//! @return Error code
int new_int(hid_t pid, const char* name, const int* value);

//! Read an integer with parent {pid}, name {name}, and value {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of integer
//! @param value [out]: Integer to read
//! @return Error code
int get_int(hid_t pid, const char* name, int* value);

//! Create a new integer attribute at location {pid} with {name} with {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of integer
//! @param value [in]: Integer to write
//! @return Error code
int new_int_attb(hid_t pid, const char* name, const int* value);

//! Read an integer attribute at location {pid} with {name} with {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of integer
//! @param value [out]: Integer to read
//! @return Error code
int get_int_attb(hid_t pid, const char* name, int* value);

//==============================================================================

//! Create a new float with parent {pid}, name {name}, and value {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of float
//! @param value [in]: Float to write
//! @return Error code
int new_float(hid_t pid, const char* name, const float* value);

//! Read an float with parent {pid}, name {name}, and value {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of float
//! @param value [out]: Float to read
//! @return Error code
int get_float(hid_t pid, const char* name, float* value);

//! Create a new float attribute at location {pid} with {name} with {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of float
//! @param value [in]: Float to write
//! @return Error code
int new_float_attb(hid_t pid, const char* name, const float* value);

//! Read an float attribute at location {pid} with {name} with {value}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of float
//! @param value [out]: Float to read
//! @return Error code
int get_float_attb(hid_t pid, const char* name, float* value);

//==============================================================================

//! Create a new node at with parent {pid}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of node
//! @param label [in]: Label of node
//! @param type [in]: Type of node
//! @return Error code
int new_node(hid_t pid, const char* name, const char* label, const char* type);

//! Delete a node with parent {pid}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of node
//! @return Error code
int del_node(hid_t pid, const char* name);

//! Check if a node exists with name {name} and parent {pid}
//! @param pid [in]: HDF5 locator for parent
//! @param name [in]: Name of node
//! @return {1=Exists, 0=Does not exist}
int node_exists(hid_t pid, const char* name);

//==============================================================================

//! Count the number of nodes of type {label} with parent {pid}
//! @param pid [in]: HDF5 locator for parent
//! @param label [in]: Label of nodes
//! @param num [out]: Number of nodes
//! @return Error code
int num_nodes(hid_t pid, const char* label, int* num);

//! Convert a node name to a node id
//! @param pid [in]: HDF5 locator for parent
//! @param label [in]: Label of node
//! @param name [in]: Name of node
//! @param idx [out]: Index of node
//! @return Error code
int node_name2idx(hid_t pid, const char* label, const char* name, int* idx);

//! Convert a node id to a node name
//! @param pid [in]: HDF5 locator for parent
//! @param label [in]: Label of node
//! @param idx [in]: Index of node
//! @param name [out]: Name of node
int node_idx2name(hid_t pid, const char* label, int idx, char* name);

//==============================================================================

//! Fill a buffer with the HDF5 version number
//! @param len [in]: Length of buf
//! @param buf [out]: String with HDF5 version number
//! @return Error code
int hdf5_version_str(int len, char* buf);

//! Fill a buffer with the HDF5 format string
//! @param len [in]: Length of buf
//! @param buf [out]: String with the HDF5 format
//! @return Error code
int hdf5_format_str(int len, char* buf);

#endif
