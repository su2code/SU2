//! @file pcgns_util.c
//! @author Kyle Horne <horne.kyle@gmail.com>
//! @version 0.2
//!
//! @section LICENSE
//! BSD style license
//!
//! @section DESCRIPTION
//! Implementation of utility functions

#include "pcgns_util.h"

#include "stdlib.h"
#include "string.h"

//=======================//
//== Begin Global Data ==//
//=======================//


file_t* files = NULL;
int files_count = 0;
int files_size = 0;

slice_t* write_queue = NULL;
int write_queue_len = 0;

//================================//
//== Begin Function Definitions ==//
//================================//

//! Function to free the array {files} at exit
void cleanup_files(void) {
	if(files!=NULL) free(files);
	}

int next_file(int* fn) {
	int err = 0;
	// Loop index
	int k;
	// Return the index of the next slot
	*fn = files_count;
	// Increment the number of used slots
	files_count++;
	// If this is the first execution, register to cleanup files at exit
	if(files_size==0) atexit(cleanup_files);
	// If the returned slot does not exist, extend the list so that it does
	if(files_size<=files_count) {
		// Keep track of how many slots were used
		int old_files_size = files_size;
		// Set the new size of the array
		files_size = 2*files_count;
		// Allocate a new array to replace the old one
		file_t* new_files = malloc(2*files_count*sizeof(file_t));
		if(new_files==NULL) cgp_doError;
		// Copy the old array's data to the new array
		for(k=0;k<old_files_size;k++) {
			new_files[k] = files[k];
			}
		// Free the old array
		free(files);
		// Point to the new array
		files = new_files;
		}
	// Assume nothing bad happened
	return 0;
	}

int free_file(file_t* file) {
	printTime;
	int err = 0;
	herr_t herr;
	// Loop index
	int k;
	// Free all the bases in this file
	printTime;
	for(k=0;k<file->nbases;k++) err = free_base(&(file->bases[k]));
	printTime;
	if(err!=0) cgp_doError;
	// If bases was allocated, free it
	if(file->bases!=NULL) free(file->bases);
	// Free the MPI info object
	err = MPI_Info_free(&(file->info));
	if(err!=MPI_SUCCESS) cgp_doError;
	// Free the property list
	printTime;
	herr = H5Pclose(file->plist_id);
	printTime;
	if(herr<0) cgp_doError;
	// Free the HDF5 file
	printTime;
	herr = H5Fclose(file->file_id);
	printTime;
	if(herr<0) cgp_doError;
	// Set the file's status to closed
	file->isOpen = FALSE;
	// Zero out HDF5 pointers
	file->plist_id = 0;
	file->file_id = 0;
	return 0;
	}

int free_base(base_t* base) {
	int err = 0;
	herr_t herr;
	// Loop index
	int k;
	// Free the zones in this base
	for(k=0;k<base->nzones;k++) err = free_zone(&(base->zones[k]));
	if(err!=0) cgp_doError;
	// If zones was allocated, free it
	if(base->zones!=NULL) free(base->zones);
	// Free the HDF5 group
	herr = H5Gclose(base->group_id);
	if(herr<0) cgp_doError;
	// Zero out the HDF5 pointer
	base->group_id = 0;
	return 0;
	}

int free_zone(zone_t* zone) {
	int err = 0;
	herr_t herr;
	// Loop index
	int k;
	// Free the zone dimensions
	if(zone->nijk!=NULL) free(zone->nijk);
	// Free the coords in this zone
	for(k=0;k<zone->ncoords;k++) err = free_coord(&(zone->coords[k]));
	if(err!=0) cgp_doError;
	// If coords was allocated, free it
	if(zone->coords!=NULL) free(zone->coords);
	// Free the sections in this zone
	for(k=0;k<zone->nsections;k++) err = free_section(&(zone->sections[k]));
	if(zone->sections!=NULL) free(zone->sections);
	// Free the sols in this zone
	for(k=0;k<zone->nsols;k++) err = free_sol(&(zone->sols[k]));
	if(err!=0) cgp_doError;
	// If sols was allocated, free it
	if(zone->sols!=NULL) free(zone->sols);
	// Free the HDF5 groups
	herr = H5Gclose(zone->group_id);
	if(herr<0) cgp_doError;
	herr = H5Gclose(zone->grid_id);
	if(herr<0) cgp_doError;
	herr = H5Gclose(zone->flow_id);
	if(herr<0) cgp_doError;
	// Zero out the HDF5 pointer
	zone->group_id = 0;
	return 0;
	}

int free_coord(coords_t* coords) {
	herr_t herr;
	// Free the HDF5 group
	herr = H5Gclose(coords->group_id);
	if(herr<0) cgp_doError;
	// Zero out the HDF5 pointer
	coords->group_id = 0;
	return 0;
	}

int free_section(section_t* section) {
	herr_t herr;
	// Free the HDF5 groups
	herr = H5Gclose(section->group_id);
	if(herr<0) cgp_doError;
	herr = H5Gclose(section->connectivity_id);
	if(herr<0) cgp_doError;
	herr = H5Gclose(section->range_id);
	if(herr<0) cgp_doError;
	// Zero out the HDf5 pointers
	section->group_id = 0;
	section->connectivity_id = 0;
	section->range_id = 0;
	return 0;
	}

int free_sol(sol_t* sol) {
	herr_t herr;
	// Free the HDF5 group
	herr = H5Gclose(sol->group_id);
	if(herr<0) cgp_doError;
	// Zero out the HDF5 pointer
	sol->group_id = 0;
	return 0;
	}


//==============================================================================

int new_str(hid_t pid, const char* name, const char* value, int len) {
	herr_t herr;
	// Set the dimension of the data to the length of the string, including the terminator
	hsize_t dim = len+1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Create the data in the file
	hid_t data_id = H5Dcreate2(pid, name, H5T_NATIVE_CHAR, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Write the data to the file
	herr = H5Dwrite(data_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int get_str(hid_t pid, const char* name, int len, char* value) {
	herr_t herr;
	// Open the data in the file
	hid_t data_id = H5Dopen2(pid, name, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Get the shape description
	hid_t shape_id = H5Dget_space(data_id);
	// Read the dimensions
	hsize_t dim;
	H5Sget_simple_extent_dims(shape_id, &dim, NULL);
	if(dim>len) cgp_doError;
	// Zero out the string
	memset(value,'\0',len+1);
	// Read the data from the file
	herr = H5Dread(data_id, H5T_NATIVE_CHAR, shape_id, shape_id, H5P_DEFAULT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int new_str_attb(hid_t pid, const char* name, const char* value, int len) {
	herr_t herr;
	// Create a shape object for the data
	hid_t shape_id = H5Screate(H5S_SCALAR);
	if(shape_id<0) cgp_doError;
	// Create a type for the data
	hid_t type_id = H5Tcopy(H5T_C_S1);
	if(type_id<0) cgp_doError;
	// Set the dimension of the data to the length of the string, including the terminator
	herr = H5Tset_size(type_id, len+1);
	if(herr<0) cgp_doError;
	// Create the data in the file
	hid_t attb_id = H5Acreate(pid, name, type_id, shape_id, H5P_DEFAULT, H5P_DEFAULT);
	if(attb_id<0) cgp_doError;
	// Write the data to the file
	herr = H5Awrite(attb_id, type_id, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Aclose(attb_id);
	if(herr<0) cgp_doError;
	// Close the type
	herr = H5Tclose(type_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

// This probably needs to be fixed to be more like get_str()
int get_str_attb(hid_t pid, const char* name, int len, char* value) {
	herr_t herr;
	// Create a shape object for the data
	hid_t shape_id = H5Screate(H5S_SCALAR);
	if(shape_id<0) cgp_doError;
	// Create a type for the data
	hid_t type_id = H5Tcopy(H5T_C_S1);
	if(type_id<0) cgp_doError;
	// Set the dimension of the data to the length of the string, including the terminator
	herr = H5Tset_size(type_id, len+1);
	if(herr<0) cgp_doError;
	// Open the data in the file
	hid_t attb_id = H5Aopen(pid, name, H5P_DEFAULT);
	if(attb_id<0) cgp_doError;
	// Read the data from the file
	herr = H5Aread(attb_id, type_id, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Aclose(attb_id);
	if(herr<0) cgp_doError;
	// Close the type
	herr = H5Tclose(type_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

//==============================================================================

int new_int(hid_t pid, const char* name, const int* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Create the data in the file
	hid_t data_id = H5Dcreate2(pid, name, H5T_NATIVE_INT, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Write the data to the file
	herr = H5Dwrite(data_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int get_int(hid_t pid, const char* name, int* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Open the data in the file
	hid_t data_id = H5Dopen2(pid, name, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Read the data from the file
	herr = H5Dread(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int new_int_attb(hid_t pid, const char* name, const int* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Create the data in the file
	hid_t attb_id = H5Acreate(pid, name, H5T_NATIVE_INT, shape_id, H5P_DEFAULT, H5P_DEFAULT);
	if(attb_id<0) cgp_doError;
	// Write the data to the file
	herr = H5Awrite(attb_id, H5T_NATIVE_INT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Aclose(attb_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int get_int_attb(hid_t pid, const char* name, int* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Open the data in the file
	hid_t attb_id = H5Aopen(pid, name, H5P_DEFAULT);
	if(attb_id<0) cgp_doError;
	// Read the data from the file
	herr = H5Aread(attb_id, H5T_NATIVE_INT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Aclose(attb_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

//==============================================================================

int new_float(hid_t pid, const char* name, const float* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Create the data in the file
	hid_t data_id = H5Dcreate2(pid, name, H5T_NATIVE_FLOAT, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Write the data to the file
	herr = H5Dwrite(data_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int get_float(hid_t pid, const char* name, float* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Open the data in the file
	hid_t data_id = H5Dopen2(pid, name, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Read the data from the file
	herr = H5Dread(data_id, H5T_NATIVE_FLOAT, shape_id, shape_id, H5P_DEFAULT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int new_float_attb(hid_t pid, const char* name, const float* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Create the data in the file
	hid_t attb_id = H5Acreate(pid, name, H5T_NATIVE_FLOAT, shape_id, H5P_DEFAULT, H5P_DEFAULT);
	if(attb_id<0) cgp_doError;
	// Write the data to the file
	herr = H5Awrite(attb_id, H5T_NATIVE_FLOAT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Aclose(attb_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int get_float_attb(hid_t pid, const char* name, float* value) {
	herr_t herr;
	// Set the dimension to be one element
	hsize_t dim = 1;
	// Create a shape object for the data
	hid_t shape_id = H5Screate_simple(1,&dim, NULL);
	if(shape_id<0) cgp_doError;
	// Open the data in the file
	hid_t attb_id = H5Aopen(pid, name, H5P_DEFAULT);
	if(attb_id<0) cgp_doError;
	// Read the data from the file
	herr = H5Aread(attb_id, H5T_NATIVE_FLOAT, value);
	if(herr<0) cgp_doError;
	// Close the data
	herr = H5Aclose(attb_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

//==============================================================================

int new_node(hid_t pid, const char* name, const char* label, const char* type) {
	int err = 0;
	herr_t herr;
	// Create the node in the file
	hid_t group_id = H5Gcreate2(pid, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(group_id<0) cgp_doError;
	// Set the label attribute; length is 32 for compatability
	err = new_str_attb(group_id, "label", label, 32);
	if(err!=0) cgp_doError;
	// Set the name attribute; length is 32 for compatability
	// The name should always correspond to the HDF5 group name
	err = new_str_attb(group_id, "name", name, 32);
	if(err!=0) cgp_doError;
	// Set the type attribute; length is 2 for compatability
	err = new_str_attb(group_id, "type", type, 2);
	if(err!=0) cgp_doError;
	// Set the flags attribure
	int val = 1;
	err = new_int_attb(group_id, "flags", &val);
	if(err!=0) cgp_doError;
	// Close the group
	herr = H5Gclose(group_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int del_node(hid_t pid, const char* name) {
	herr_t herr;
	// Delete the node
	herr = H5Ldelete(pid, name, H5P_DEFAULT);
	if(herr<0) cgp_doError;
	return 0;
	}

int node_exists(hid_t pid, const char* name) {
	// Check for the existance of a node
	return H5Lexists(pid, name, H5P_DEFAULT)?TRUE:FALSE;
	}

//==============================================================================

//! Structure to describe nodes during iteration
struct iter_s {int pos; int counter; char* name; char* label;};

//! Call-back function used to count the nodes with a particular label
herr_t node_counter(hid_t gid, const char* name, const H5L_info_t* linfo, void* vdata) {
	int err = 0;
	herr_t herr;
	// Cast vdata as an iter_s
	struct iter_s* data = (struct iter_s*) vdata;
	// info object to hold data
	H5O_info_t info;
	// Get the info for the current object
	herr = H5Oget_info_by_name(gid, name, &info, H5P_DEFAULT);
	if(herr<0) cgp_doError;
	// If the object is a group
	if(info.type==H5O_TYPE_GROUP) {
		// Buffer to store the label
		char label[100+1];
		// Open the group
		hid_t group_id = H5Gopen2(gid, name, H5P_DEFAULT);
		if(group_id<0) cgp_doError;
		// Read the label
		err = get_str_attb(group_id, "label", 100, label);
		if(err!=0) cgp_doError;
		// Close the group
		herr = H5Gclose(group_id);
		if(herr<0) cgp_doError;
		// If the label is the same as the given label, increment the counter
		if(strcmp(label,data->label)==0) data->counter++;
		}
	return 0;
	}

int num_nodes(hid_t pid, const char* label, int* num) {
	herr_t herr;
	// An iter_s to hold the iteration results
	struct iter_s data;
	// Initialize position to be -1
	data.pos = -1;
	// Initailize the counter to zero
	data.counter = 0;
	// Point the label to the given label
	data.label = (char*) label;
	// Iteration though the children of pid
	herr = H5Literate(pid, H5_INDEX_NAME, H5_ITER_INC, NULL, node_counter, &data);
	if(herr<0) cgp_doError;
	// Return the number of children whos label's match label
	*num = data.counter;
	return 0;
	}

//! Call-back function used to find the index of a nodes with a particular label
herr_t node_idx_finder(hid_t gid, const char* name, const H5L_info_t* linfo, void* vdata) {
	int err = 0;
	herr_t herr;
	// Cast vdata as an iter_s
	struct iter_s* data = (struct iter_s*) vdata;
	// Info for each object
	H5O_info_t info;
	// Get the info for the current object
	herr = H5Oget_info_by_name(gid, name, &info, H5P_DEFAULT);
	if(herr<0) cgp_doError;
	// If the object is a group
	if(info.type==H5O_TYPE_GROUP) {
		// Buffer for label
		char label[100+1];
		// Buffer for name
		char nname[100+1];
		// Open the group
		hid_t group_id = H5Gopen2(gid, name, H5P_DEFAULT);
		if(group_id<0) cgp_doError;
		// Read the label
		err = get_str_attb(group_id, "label", 100, label);
		if(err!=0) cgp_doError;
		// Read the name
		err = get_str_attb(group_id, "name", 100, nname);
		if(err!=0) cgp_doError;
		// Close the group
		herr = H5Gclose(group_id);
		if(herr<0) cgp_doError;
		// If the labels are the same
		if(strcmp(label,data->label)==0) {
			// If the names are the same, set the position to the counter
			if(strcmp(nname,data->name)==0) data->pos=data->counter;
			// Increment the counter
			data->counter++;
			}
		}
	return 0;
	}

int node_name2idx(hid_t pid, const char* label, const char* name, int* idx) {
	herr_t herr;
	// Data for iteration
	struct iter_s data;
	// Initialize position
	data.pos = -1;
	// Initialize counter
	data.counter = 0;
	// Point the label to the given label
	data.label = (char*) label;
	// Point the name to the given name
	data.name = (char*) name;
	// Iterate though all objects in this group
	herr = H5Literate(pid, H5_INDEX_NAME, H5_ITER_INC, NULL, node_idx_finder, &data);
	if(herr<0) cgp_doError;
	// Set the index to pos
	*idx = data.pos;
	return 0;
	}

//! Call-back function used to find the name of a node at index {idx}
herr_t node_name_finder(hid_t gid, const char* name, const H5L_info_t* linfo, void* vdata) {
	int err = 0;
	herr_t herr;
	// Cast vdata to an iter_s object
	struct iter_s* data = (struct iter_s*) vdata;
	// Info about the object
	H5O_info_t info;
	// Read the info
	herr = H5Oget_info_by_name(gid, name, &info, H5P_DEFAULT);
	if(herr<0) cgp_doError;
	// If the object is a group
	if(info.type==H5O_TYPE_GROUP) {
		// Buffer for the label
		char label[100+1];
		// Buffer for the name
		char nname[100+1];
		// Open the group
		hid_t group_id = H5Gopen2(gid, name, H5P_DEFAULT);
		if(group_id<0) cgp_doError;
		// Read the label
		err = get_str_attb(group_id, "label", 100, label);
		if(err!=0) cgp_doError;
		// Read the name
		err = get_str_attb(group_id, "name", 100, nname);
		if(err!=0) cgp_doError;
		// Close the group
		herr = H5Gclose(group_id);
		if(herr<0) cgp_doError;
		// If the label is the same as the given label
		if(strcmp(label,data->label)==0) {
			// If the counter is the same as the given coutner, copy the name
			if(data->counter==data->pos) strcpy(data->name, nname);
			// Increment the counter
			data->counter++;
			}
		}
	return 0;
	}

int node_idx2name(hid_t pid, const char* label, int idx, char* name) {
	herr_t herr;
	// Data for iteration
	struct iter_s data;
	// Initalize the position
	data.pos = idx;
	// Initialize the counter
	data.counter = 0;
	// Point the label
	data.label = (char*) label;
	// Point the name
	data.name = name;
	// Iterate through the objects
	herr = H5Literate(pid, H5_INDEX_NAME, H5_ITER_INC, NULL, node_name_finder, &data);
	if(herr<0) cgp_doError;
	return 0;
	}

//==============================================================================

int hdf5_version_str(int len, char* buf) {
	herr_t herr;
	// Three components of the HDF5 version
	unsigned int maj, min, rel;
	// Read the version
	herr = H5get_libversion(&maj, &min, &rel);
	if(herr<0) cgp_doError;
	// Preset the buffer to zeros
	memset(buf,'\0',len+1);
	// Write the verions string to the buffer
	sprintf(buf, "HDF5 Version %d.%d.%d", maj, min, rel);
	return 0;
	}

int hdf5_format_str(int len, char* buf) {
	herr_t herr;
	// Get a copy of the native float type
	hid_t type_id = H5Tcopy(H5T_NATIVE_FLOAT);
	if(type_id<0) cgp_doError;
	// Preset the buffer to zeros
	memset(buf,'\0',len+1);
	// Write the appropriate formate string to the buffer
	if      (H5Tequal(type_id, H5T_IEEE_F32BE)) strcpy(buf, "IEEE_BIG_32");
	else if (H5Tequal(type_id, H5T_IEEE_F32LE)) strcpy(buf, "IEEE_LITTLE_32");
	else if (H5Tequal(type_id, H5T_IEEE_F64BE)) strcpy(buf, "IEEE_BIG_64");
	else if (H5Tequal(type_id, H5T_IEEE_F64LE)) strcpy(buf, "IEEE_LITTLE_64");
	else sprintf(buf, "NATIVE_%d", H5Tget_precision(type_id));
	// Close the type
	herr = H5Tclose(type_id);
	if(herr<0) cgp_doError;
	return 0;
	}
