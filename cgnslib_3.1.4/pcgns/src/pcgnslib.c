//! @file pcgnslib.c
//! @author Kyle Horne <horne.kyle@gmail.com>
//! @version 0.2
//!
//! @section LICENSE
//! BSD style license
//!
//! @section DESCRIPTION
//! Implimentation of functions provided by the pcgns library

#include "pcgnslib.h"
#include "pcgns_util.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

// Number of nodes per element
// Index corresponds to the ElementType_t enum
// Negative numbers indicate non-supported types
int node_counts[24] = {
	-1, -1,
	1,
	2,3,
	3,6,
	4,8,9,
	4,10,
	5,13,14,
	6,15,18,
	8,20,27,
	-1,
	-1,-1
	};

int preallocate = 0;

int default_pio_mode = FALSE;

//================================//
//== Begin Function Definitions ==//
//================================//

//= File IO Prototypes =//
int cgp_open(const char* filename, int mode, MPI_Comm comm, MPI_Info* info, int* fn) {
	int err;
	herr_t herr;
	// Test file to check for file's existance
	FILE* test_file = fopen(filename, "rb");
	// Flag for file's existance
	int file_exists = (test_file==NULL)?FALSE:TRUE;
	// Flag for file's HDF5 status
	int file_isHDF5 = FALSE;
	// If the file exists
	if(file_exists) {
		// Close the file
		fclose(test_file);
		// Test to see if it is an HDF5 file
		file_isHDF5 = H5Fis_hdf5(filename)?TRUE:FALSE;
		}
	// Create a new file reference in the global files array
	next_file(fn);
	// Create a pointer to this file in that array
	file_t* file = &(files[*fn]);
	// Set the file's index
	file->idx = *fn;
	// Set the file's state to open
	file->isOpen = TRUE;
	// Set the parallel IO mode
	file->collective = default_pio_mode;
	// Set the file's name
	strcpy(file->filename,filename);
	// Set the file's communicator
	file->comm = comm;
	// Set the MPI_Info for the file
	err = MPI_Info_dup(*info,&(file->info));
	if(err!=MPI_SUCCESS) cgp_doError;
	// Set the rank on this communicator
	err = MPI_Comm_rank(comm,&(file->rank));
	if(err!=MPI_SUCCESS) cgp_doError;
	// Set the size of this communicator
	err = MPI_Comm_size(comm,&(file->size));
	if(err!=MPI_SUCCESS) cgp_doError;
	// Set the access property list
	file->plist_id = H5Pcreate(H5P_FILE_ACCESS);
	if(file->plist_id<0) cgp_doError;
	// Set the access property list to use MPI
	herr = H5Pset_fapl_mpio(file->plist_id, file->comm, file->info);
	if(herr<0) cgp_doError;
	// If the file exists and it is an HDF5 file, try to read is as CGNS
	if(file_exists&&file_isHDF5) {
		// Open the file with HDF5 and set the file_id
		file->file_id = H5Fopen(filename, H5F_ACC_RDWR, file->plist_id);
		if(file->file_id<0) cgp_doError;
		// Read the number of bases in the file
		num_nodes(file->file_id, "CGNSBase_t",&(file->nbases));
		// Allocate space to store the bases
		file->bases = malloc(file->nbases*sizeof(base_t));
		if(file->bases==NULL) cgp_doError;
		// Loop counter
		int k;
		// For each base in bases...
		for(k=0;k<file->nbases;k++) {
			// Buffer to store the name of the base
			char basename[100+1];
			// Find the name of the k'th base in bases
			node_idx2name(file->file_id, "CGNSBase_t", k, basename);
			// Open that base
			file->bases[k].group_id = H5Gopen2(file->file_id, basename, H5P_DEFAULT);
			if(file->bases[k].group_id<0) cgp_doError;
			// Set the base's index
			file->bases[k].idx = k;
			// Copy the base's name
			strcpy(file->bases[k].basename,basename);
			// Iniialize it to have no zones
			file->bases[k].zones = NULL;
			file->bases[k].nzones = 0;
			}
		}
	// If the file does not exist, create an empty CGNS file
	else if((!file_exists)||(!file_isHDF5)) {
		// Create a new HDF5 file
		file->file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, file->plist_id);
		if(file->file_id<0) cgp_doError;
		// Create the needed attributes to describe the root node
		new_str_attb(file->file_id, "label", "Root Node of ADF File", 32);
		new_str_attb(file->file_id, "name", "HDF5 MotherNode", 32);
		new_str_attb(file->file_id, "type", "MT", 2);
		// Buffer for version information
		char version[100+1];
		// Buffer for format information
		char format[100+1];
		// Get the HDF5 format string
		hdf5_format_str(100,format);
		// Get the HDF5 version string
		hdf5_version_str(100,version);
		// Write the HDF5 format string to the file
		new_str(file->file_id, " format", format, strlen(format));
		// Write the HDF5 version string to the file
		new_str(file->file_id, " hdf5version", version, 32);
		// Create a CGNS library version number
		new_node(file->file_id, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4");
		// Open that version number group
		hid_t group_id = H5Gopen2(file->file_id, "CGNSLibraryVersion", H5P_DEFAULT);
		if(group_id<0) cgp_doError;
		// Set the files to be compatible with CGNS version 3.0
		float libversion = 3.0;
		// Write the version number to the file
		new_float(group_id, " data", &libversion);
		// Close the version number group
		herr = H5Gclose(group_id);
		if(herr<0) cgp_doError;

		// set the file to have no bases
		file->bases = NULL;
		file->nbases = 0;
		}
	int dummy = 0;
	err = new_int(file->file_id, " dummy",&dummy);
	if(err!=0) cgp_doError;
	return 0;
	}

int cgp_close(int fn) {
	// Free the file referenced by fn
	printTime;
	del_node(files[fn].file_id, " dummy");
	free_file(&(files[fn]));
	printTime;
	return 0;
	}

//= Base IO Prototypes =//
int cgp_base_read(int fn, int B, char* basename, int* cell_dim, int* phys_dim) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	herr_t herr;
	// Pointer to the base being read
	base_t* base = &(files[fn].bases[B-1]);
	// Get the name of the base
	strcpy(basename, base->basename);

	// Size of array to read
	hsize_t dim = 2;
	// Create shape object for shape of array
	hid_t shape_id = H5Screate_simple(1,&dim,NULL);
	if(shape_id<0) cgp_doError;
	// Open the data array in the file
	hid_t data_id = H5Dopen2(base->group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Buffer to store the data
	int data[2];
	// Read the data array condaining this base's dimensions
	herr = H5Dread(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, data);
	if(herr<0) cgp_doError;
	// Copy the data out of the array [The order here needs verification]
	*cell_dim = data[0];
	base->cell_dim = data[0];
	*phys_dim = data[1];
	base->phys_dim = data[1];
	// Close the data array in the file
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Destroy the shape object
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;

	// Read Zone info
	// Read the number of zones in the base
	num_nodes(base->group_id, "Zone_t", &(base->nzones));
	// Allocate space for zone descriptors
	base->zones = (zone_t*) malloc(base->nzones*sizeof(zone_t));
	if(base->zones==NULL) cgp_doError;
	// Loop index
	int k;
	// For each zone in zones
	for(k=0;k<base->nzones;k++) {
		// Buffer to store the zone's anme
		char zonename[100+1];
		// Get the name of the k'th zone
		node_idx2name(base->group_id, "Zone_t", k, zonename);
		// Open the group for this zone
		base->zones[k].group_id = H5Gopen2(base->group_id,zonename,H5P_DEFAULT);
		if(base->zones[k].group_id<0) cgp_doError;
		// Set the zone's index
		base->zones[k].idx = k;
		// Copy the zone's name
		strcpy(base->zones[k].zonename,zonename);
		// Default to no coords or solutions
		base->zones[k].coords = NULL;
		base->zones[k].ncoords = 0;
		base->zones[k].sols = NULL;
		base->zones[k].nsols = 0;
		}
	return 0;
	}

int cgp_base_write(int fn, char const* basename, int cell_dim, int phys_dim, int* B) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	herr_t herr;
	// Pointer to this file
	file_t* file = &(files[fn]);
	// Index
	int idx;
	// If this base already exists, replace the old one
	if(node_exists(file->file_id, basename)) {
		// Loop index
		int k;
		// Find this base in the list of bases for this file
		for(k=0;k<file->nbases;k++) if(strcmp(basename,file->bases[k].basename)==0) idx = k;
		// Free the old base in memory
		free_base(&(file->bases[idx]));
		// Delete the old base in the file
		del_node(file->file_id, basename);
		}
	// If this base does not yet exist, create one
	else {
		// Loop index
		int k;
		// Allocate space for a longer list of bases
		base_t* bases = malloc((file->nbases+1)*sizeof(base_t));
		if(bases==NULL) cgp_doError;
		// Copy all the old bases to the new list
		for(k=0;k<file->nbases;k++) bases[k] = file->bases[k];
		// Free the old list
		if(file->bases!=NULL) free(file->bases);
		// Point the file to use the new list
		file->bases = bases;
		// Increment the number of bases in the file
		file->nbases++;
		// Set the index to the last base in the list
		idx = file->nbases-1;
		}
	// Pointer to the current base
	base_t* base = &(file->bases[idx]);
	// Copy the name of the base
	strcpy(file->bases[idx].basename, basename);
	// Set the  index of the base
	base->idx = idx;
	// Create no zones in the base
	base->zones = NULL;
	base->nzones = 0;
	// Set the base's dimensions
	base->cell_dim = cell_dim;
	base->phys_dim = phys_dim;
	// Return the index of the base
	*B = idx+1;
	// Create the base in the file
	new_node(file->file_id, basename, "CGNSBase_t", "I4");
	// Open the base
	hid_t group_id = H5Gopen2(file->file_id, basename, H5P_DEFAULT);
	if(group_id<0) cgp_doError;
	// Copy the group_id for later use
	file->bases[idx].group_id = group_id;
	// Set the size of the array to write
	hsize_t dim = 2;
	// Create a shape object for the array to write
	hid_t shape_id = H5Screate_simple(1,&dim,NULL);
	if(shape_id<0) cgp_doError;
	// Create the data in the file
	hid_t data_id = H5Dcreate2(group_id, " data", H5T_NATIVE_INT, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Buffer to hold that data array
	int data[2];
	// Fill the buffer [Order needs verification]
	data[0] = cell_dim;
	data[1] = phys_dim;
	// Write the data array to the file
	herr = H5Dwrite(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, data);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Destroy the shape object
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;

	return 0;
	}

int cgp_nbases(int fn, int *nbases) {
	// Read the number of bases in this file
	*nbases = files[fn].nbases;
	return 0;
	}

//= Zone IO Prototypes =//
int cgp_zone_read(int fn, int B, int Z, char* zonename, int* nijk) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	int err;
	// Pointer to the current base
	base_t* base = &(files[fn].bases[B-1]);
	// Pointer to the zone to be read
	zone_t* zone = &(files[fn].bases[B-1].zones[Z-1]);
	// Copy the zone's name
	strcpy(zonename, zone->zonename);
	// Loop index
	int k;
	// If the zone already has a coordinates sub-group, read it
	if(node_exists(zone->group_id, "GridCoordinates")) {
		// Test if the existing HDF5 handle for the grid_id is valid
		// If it is, use it, otherwise open a new handle
		if(!H5Iis_valid(zone->grid_id)) zone->grid_id = H5Gopen2(zone->group_id, "GridCoordinates", H5P_DEFAULT);
		if(zone->grid_id<0) cgp_doError;
		// Count the number of data arrays in the group
		num_nodes(zone->grid_id, "DataArray_t", &(zone->ncoords));
		// Allocate memory to describe the data arrays
		zone->coords = (coords_t*) malloc(zone->ncoords*sizeof(coords_t));
		if(zone->coords==NULL) cgp_doError;
		// For each ncoord in ncoords
		for(k=0;k<zone->ncoords;k++) {
			// Buffer for the name of the coords
			char coordname[100+1];
			// Read the name of the coords
			node_idx2name(zone->grid_id, "DataArray_t", k, coordname);
			// Open the data array
			zone->coords[k].group_id = H5Gopen2(zone->grid_id, coordname, H5P_DEFAULT);
			if(zone->coords[k].group_id<0) cgp_doError;
			// Set the index of the data
			zone->coords[k].idx = k;
			// Copy the name of the coords
			strcpy(zone->coords[k].coordname,coordname);
			}
		}
	// If the zone does not have a coordinates sub-group, create one for it
	else {
		// Create a new group to hold the coords
		new_node(zone->group_id, "GridCoordinates", "GridCoordinates_t", "MT");
		// Open the new group
		zone->grid_id = H5Gopen2(zone->group_id, "GridCoordinates", H5P_DEFAULT);
		if(zone->grid_id<0) cgp_doError;
		// Default to hold no coords data
		zone->ncoords = 0;
		zone->coords = NULL;
		}

	// If the zone already has a flow sub-group, read it
	if(node_exists(zone->group_id, "GridCoordinates")) {
		// Test if the existing HDF5 handle for the flow_id is valid
		// If it is, use it, otherwise open a new handle
		if(!H5Iis_valid(zone->flow_id)) zone->flow_id = H5Gopen2(zone->flow_id, "FlowSolution", H5P_DEFAULT);
		if(zone->flow_id<0) cgp_doError;
		// Count the number of data arrays in the group
		num_nodes(zone->flow_id, "DataArray_t", &(zone->ncoords));
		// Allocate memory to describe the data arrays
		zone->sols = (sol_t*) malloc(zone->nsols*sizeof(sol_t));
		if(zone->sols==NULL) cgp_doError;
		// For each sol in nsols
		for(k=0;k<zone->nsols;k++) {
			// Buffer for the name of the coords
			char solname[100+1];
			// Read the name of the coords
			node_idx2name(zone->flow_id, "DataArray_t", k, solname);
			// Open the data array
			zone->sols[k].group_id = H5Gopen2(zone->flow_id, solname, H5P_DEFAULT);
			if(zone->sols[k].group_id<0) cgp_doError;
			// Set the index of the data
			zone->sols[k].idx = k;
			// Copy the name of the coords
			strcpy(zone->sols[k].solname,solname);
			}
		}
	// If the zone does not have a flow sub-group, create one for it
	else {
		// Create a new group to hold the coords
		new_node(zone->group_id, "FlowSolution", "FlowSolution_t", "MT");
		// Open the new group
		zone->flow_id = H5Gopen2(zone->group_id, "FlowSolution", H5P_DEFAULT);
		if(zone->flow_id<0) cgp_doError;
		// Default to hold no coords data
		zone->nsols = 0;
		zone->sols = NULL;
		}

	// Count the number of element sections
	num_nodes(zone->group_id, "Elements_t", &(zone->nsections));
	// Allocate memory to describe the element sections
	zone->sections = (section_t*) malloc(zone->nsections*sizeof(section_t));
	// For each section in nsections
	for(k=0;k<zone->nsections;k++) {
		// Buffer for the section's name
		char sectionname[100+1];
		// Read the name of the k'th section
		node_idx2name(zone->group_id, "Elements_t", k, sectionname);
		// Open the section
		zone->sections[k].group_id = H5Gopen2(zone->group_id, sectionname, H5P_DEFAULT);
		if(zone->sections[k].group_id<0) cgp_doError;
		// Open the section's connectivity
		zone->sections[k].connectivity_id = H5Gopen2(zone->sections[k].group_id, "ElementConnectivity", H5P_DEFAULT);
		if(zone->sections[k].connectivity_id<0) cgp_doError;
		// Open the section's ranges
		zone->sections[k].range_id = H5Gopen2(zone->sections[k].group_id, "ElementRange", H5P_DEFAULT);
		if(zone->sections[k].range_id<0) cgp_doError;
		// Set the index of the section
		zone->sections[k].idx = k;
		// Copy the name of the section
		strcpy(zone->sections[k].sectionname,sectionname);
		// Read the element type 
		{
			// Dimensions of data
			hsize_t dims[1] = {2};
			// Shape of data in memory and the file
			hid_t shape_id = H5Screate_simple(1,dims,NULL);
			// Open the array
			hid_t data_id = H5Dopen2(zone->sections[k].group_id, " data", H5P_DEFAULT);
			// Buffer to hold results
			int data[2];
			// Read array
			H5Dread(data_id,H5T_NATIVE_INT, shape_id,shape_id,H5P_DEFAULT,data);
			// Set the type
			zone->sections[k].type = data[0];
			// Clsoe the array
			H5Dclose(data_id);
			// Close the shape
			H5Sclose(shape_id);
			}
		// Read the element range
		{
			// Dimensions of data
			hsize_t dims[1] = {2};
			// Shape of data in memory and the file
			hid_t shape_id = H5Screate_simple(1,dims,NULL);
			// Open the array
			hid_t data_id = H5Dopen2(zone->sections[k].range_id, " data", H5P_DEFAULT);
			// Buffer to hold results
			int data[2];
			// Read array
			H5Dread(data_id,H5T_NATIVE_INT, shape_id,shape_id,H5P_DEFAULT,data);
			// Set the range
			zone->sections[k].range[0] = data[0];
			zone->sections[k].range[1] = data[1];
			// Clsoe the array
			H5Dclose(data_id);
			// Close the shape
			H5Sclose(shape_id);
			}
		}

	// Count the number of flow solutions
	num_nodes(zone->group_id, "FlowSolution_t", &(zone->nsols));
	// Allocate memory to describe the flow solutions
	zone->sols = (sol_t*) malloc(zone->nsols*sizeof(sol_t));
	if(zone->sols==NULL) cgp_doError;
	// For each solution in nsols
	for(k=0;k<zone->nsols;k++) {
		// Buffer for the solution's name
		char solname[100+1];
		// Read the name of the k'th solution
		node_idx2name(zone->group_id, "FlowSolution_t", k, solname);
		// Open the solution
		zone->sols[k].group_id = H5Gopen2(zone->group_id, solname, H5P_DEFAULT);
		if(zone->sols[k].group_id<0) cgp_doError;
		// Set the index of the solution
		zone->sols[k].idx = k;
		// Copy the name of the solution
		strcpy(zone->sols[k].solname,solname);
		}

	// Read the zone type
	char ztype[100];
	err = get_str(zone->group_id, "ZoneType/ data", 99, ztype);
	if(strcmp(ztype,"Structured")) zone->type = Structured;
	else if(strcmp(ztype,"Untructured")) zone->type = Unstructured;
	else cgp_doError;

	// Read nijk
	// Set the size of the array to read
	int cols = (zone->type==Structured)?base->cell_dim:1;
	hsize_t dim[2] = {3, cols};
	// Create the shape object for the array
	hid_t shape_id = H5Screate_simple(2,dim, NULL);
	if(shape_id<0) cgp_doError;
	// Open the array in the file
	hid_t data_id = H5Dopen2(zone->group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Read the array from the file to the local nijk
	herr = H5Dread(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, nijk);
	if(herr<0) cgp_doError;
	// Allocate space to store the data
	zone->nijk = malloc(3*cols*sizeof(int));
	if(nijk==NULL) cgp_doError;
	// Read the array from the file to the global nijk
	herr = H5Dread(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, zone->nijk);
	if(herr<0) cgp_doError;
	// Close the array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Destroy the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int cgp_zone_type(int fn, int B, int Z, ZoneType_t *zonetype) {
	// Read Zonetype from file and return in enum
	file_t* file = &(files[fn]);
	base_t* base = &(file->bases[B-1]);
	zone_t* zone = &(base->zones[Z-1]);
	*zonetype = zone->type;
	return 0;
	}

int cgp_zone_write(int fn, int B, const char* zonename, const int* nijk, ZoneType_t type, int* Z) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	herr_t herr;
	// Pointer to the base
	base_t* base = &(files[fn].bases[B-1]);
	// Index
	int idx;
	// Loop index
	int k;
	// If the zone exists, replace it
	if(node_exists(base->group_id, zonename)) {
		// Find the zone in this base
		for(k=0;k<base->nzones;k++) if(strcmp(zonename,base->zones[k].zonename)==0) idx = k;
		// Free the memory for the zone
		free_zone(&(base->zones[idx]));
		// Delete the zone from the file
		del_node(base->group_id,zonename);
		}
	// If the zone does not exist, initialize memory for it
	else {
		// Allocate a bigger list of zones
		zone_t* zones = (zone_t*) malloc((base->nzones+1)*sizeof(zone_t));
		if(zones==NULL) cgp_doError;
		// Copy the old zones to the new list
		for(k=0;k<base->nzones;k++) zones[k] = base->zones[k];
		// Deallocate the old list
		if(base->zones!=NULL) free(base->zones);
		// Point the base's list to the new list
		base->zones = zones;
		// Increment the number of zones
		base->nzones++;
		// Set the index of the zone
		idx = base->nzones-1;
		}
	// Pointer to the current zone
	zone_t* zone= &(base->zones[idx]);
	// Copy the name of the zone
	strcpy(zone->zonename, zonename);
	// Set the index of the zone
	zone->idx = idx;
	// Default to no coords or solutions
	zone->coords = NULL;
	zone->ncoords = 0;
	zone->sections = NULL;
	zone->nsections = 0;
	zone->sols = NULL;
	zone->nsols = 0;
	// Return the index of the zone
	*Z = idx+1;

	// Create the zone in the file
	new_node(base->group_id, zonename, "Zone_t", "I4");
	// Open the new zone
	zone->group_id = H5Gopen2(base->group_id, zonename, H5P_DEFAULT);
	if(zone->group_id<0) cgp_doError;

	// Write ZoneType
	zone->type = type;
	// Create a new node for the zone type
	new_node(zone->group_id, "ZoneType", "ZoneType_t", "C1");
	// Open the new node
	hid_t group_id = H5Gopen2(zone->group_id, "ZoneType", H5P_DEFAULT);
	if(group_id<0) cgp_doError;
	// Write a string as the data for the new node
	if(type==Structured) new_str(group_id, " data", "Structured", strlen("Structured")-1);
	else if(type==Unstructured) new_str(group_id, " data", "Unstructured", strlen("Unstructured")-1);
	else cgp_doError;
	// Close the node
	herr = H5Gclose(group_id);
	if(herr<0) cgp_doError;

	// Write nijk
	// Size of array to write
	int cols = (type==Structured)?base->cell_dim:1;
	hsize_t dim[2] = {3, cols};
	// Create shape object for array to write
	hid_t shape_id = H5Screate_simple(2,dim, NULL);
	if(shape_id<0) cgp_doError;
	// Create the data array in the file
	hid_t data_id = H5Dcreate2(zone->group_id, " data", H5T_NATIVE_INT, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Write the data to the array
	herr = H5Dwrite(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, nijk);
	if(herr<0) cgp_doError;
	// Allocate space to store the data
	zone->nijk = malloc(3*cols*sizeof(int));
	if(zone->nijk==NULL) cgp_doError;
	// Read the data to the zone in memory
	herr = H5Dread(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, zone->nijk);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Destroy the shape object
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;

	// Write GridCoordinates
	// Create a node for the grid coordinates
	new_node(zone->group_id, "GridCoordinates", "GridCoordinates_t", "MT");
	// Open the node
	zone->grid_id = H5Gopen2(zone->group_id, "GridCoordinates", H5P_DEFAULT);
	if(zone->grid_id<0) cgp_doError;

	// Write FlowSolution
	// Create a node for the grid coordinates
	new_node(zone->group_id, "FlowSolution", "FlowSolution_t", "MT");
	// Open the node
	zone->flow_id = H5Gopen2(zone->group_id, "FlowSolution", H5P_DEFAULT);
	if(zone->grid_id<0) cgp_doError;

	return 0;
	}

int cgp_nzones(int fn, int B, int *nzones) {
	// Read the number of zones in this base
	*nzones = files[fn].bases[B-1].nzones;
	return 0;
	}

//= Grid IO Prototypes =//
int cgp_coord_write(int fn, int B, int Z, DataType_t type, const char* coordname, int* C) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Index
	int idx;
	// Loop index
	int k;
	// If the coordinates already exist, replace them
	if(node_exists(zone->grid_id, coordname)) {
		// Find the coordinates to be replaced
		for(k=0;k<zone->ncoords;k++) if(strcmp(coordname,zone->coords[k].coordname)==0) idx = k;
		// Delete them
		del_node(zone->grid_id,coordname);
		}
	// If the coordinates do not exist, prepare the memory structures
	else {
		// New list of coords
		coords_t* coords = (coords_t*) malloc((zone->ncoords+1)*sizeof(coords_t));
		if(coords==NULL) cgp_doError;
		// Copy old list to new list
		for(k=0;k<zone->ncoords;k++) coords[k] = zone->coords[k];
		// Free the old list
		if(zone->coords!=NULL) free(zone->coords);
		// Point the zone to the new list
		zone->coords = coords;
		// Incement the number of coords
		zone->ncoords++;
		// Set the coords index
		idx = zone->ncoords-1;
		}
	// Pointer to coords
	coords_t* coords = &(zone->coords[idx]);
	// Copy the name
	strcpy(coords->coordname, coordname);
	// Set the index
	coords->idx = idx;
	// Return the index
	*C = idx+1;

	// Create a new node in the file
	new_node(zone->grid_id, coordname, "DataArray_t", "R8");
	// Open the node
	coords->group_id = H5Gopen2(zone->grid_id, coordname, H5P_DEFAULT);
	if(coords->group_id<0) cgp_doError;

	// Set the rank of the dimensions
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Dimensions of the array
	hsize_t dims[rank];
	// Set these to correspond with the size of the zone
	for(k=0;k<rank;k++) dims[k] = zone->nijk[(rank-1)-k];

	// Create a shape for the array
	hid_t shape_id = H5Screate_simple(rank,dims,NULL);
	if(shape_id<0) cgp_doError;
	// Property list for dataset
	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
	if(preallocate) herr = H5Pset_alloc_time(plist_id,H5D_ALLOC_TIME_EARLY);
	if(preallocate) herr = H5Pset_fill_time(plist_id, H5D_FILL_TIME_ALLOC);
	// Create the array in the file
	hid_t data_id = H5Dcreate2(coords->group_id, " data", H5T_NATIVE_DOUBLE, shape_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Close the array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;	
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int cgp_coord_write_data(int fn, int B, int Z, int C, int* min, int* max, void* coord_array) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	if(C>files[fn].bases[B-1].zones[Z-1].ncoords||C<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Pointer to the current coords
	coords_t* coords = &(zone->coords[C-1]);

	// Open the data
	hid_t data_id = H5Dopen2(coords->group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Set the start position for the data write
	hsize_t start[rank];
	for(k=0;k<rank;k++) start[k] = min[k];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	for(k=0;k<rank;k++) dims[k] = max[k]-min[k]+1;
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Write the data in collective parallel I/O
	herr = H5Dwrite(data_id, H5T_NATIVE_DOUBLE, mem_shape_id, data_shape_id, plist_id, coord_array);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;

	return 0;
	}

int cgp_coord_read_data(int fn, int B, int Z, int C, int* min, int* max, void* coord_array) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	if(C>files[fn].bases[B-1].zones[Z-1].ncoords||C<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Pointer to the current coords
	coords_t* coords = &(zone->coords[C-1]);

	// Open the data
	hid_t data_id = H5Dopen2(coords->group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Set the start position for the data write
	hsize_t start[rank];
	for(k=0;k<rank;k++) start[k] = min[k];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	for(k=0;k<rank;k++) dims[k] = max[k]-min[k]+1;
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Read the data in collective parallel I/O
	herr = H5Dread(data_id, H5T_NATIVE_DOUBLE, mem_shape_id, data_shape_id, plist_id, coord_array);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;

	return 0;
	}

int cgp_sol_write_data(int fn, int B, int Z, int S, int* min, int* max, void* data) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	if(S>files[fn].bases[B-1].zones[Z-1].nsols||S<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Pointer to the current coords
	sol_t* sol = &(zone->sols[S-1]);

	// Open the data
	hid_t data_id = H5Dopen2(sol->group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Set the start position for the data write
	hsize_t start[rank];
	for(k=0;k<rank;k++) start[k] = min[k];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	for(k=0;k<rank;k++) dims[k] = max[k]-min[k]+1;
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Write the data in collective parallel I/O
	herr = H5Dwrite(data_id, H5T_NATIVE_DOUBLE, mem_shape_id, data_shape_id, plist_id, data);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;

	return 0;
	}

int cgp_sol_read_data(int fn, int B, int Z, int S, int* min, int* max, void* data) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	if(S>files[fn].bases[B-1].zones[Z-1].nsols||S<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Pointer to the current coords
	sol_t* sol = &(zone->sols[S-1]);

	// Open the data
	hid_t data_id = H5Dopen2(sol->group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Set the start position for the data write
	hsize_t start[rank];
	for(k=0;k<rank;k++) start[k] = min[k];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	for(k=0;k<rank;k++) dims[k] = max[k]-min[k]+1;
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Read the data in collective parallel I/O
	herr = H5Dread(data_id, H5T_NATIVE_DOUBLE, mem_shape_id, data_shape_id, plist_id, data);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;

	return 0;
	}

int cgp_sol_write(int fn, int B, int Z, char *solname, GridLocation_t location, int *S) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Index
	int idx;
	// Loop index
	int k;
	// If the solution already exist, replace it
	if(node_exists(zone->flow_id, solname)) {
		// Find the coordinates to be replaced
		for(k=0;k<zone->nsols;k++) if(strcmp(solname,zone->sols[k].solname)==0) idx = k;
		// Delete them
		del_node(zone->flow_id,solname);
		}
	// If the coordinates do not exist, prepare the memory structures
	else {
		// New list of coords
		sol_t* sols = (sol_t*) malloc((zone->nsols+1)*sizeof(sol_t));
		if(sols==NULL) cgp_doError;
		// Copy old list to new list
		for(k=0;k<zone->nsols;k++) sols[k] = zone->sols[k];
		// Free the old list
		if(zone->sols!=NULL) free(zone->sols);
		// Point the zone to the new list
		zone->sols = sols;
		// Incement the number of coords
		zone->nsols++;
		// Set the coords index
		idx = zone->nsols-1;
		}
	// Pointer to coords
	sol_t* sol = &(zone->sols[idx]);
	// Copy the name
	strcpy(sol->solname, solname);
	// Set the index
	sol->idx = idx;
	// Return the index
	*S = idx+1;

	// Create a new node in the file
	new_node(zone->flow_id, solname, "DataArray_t", "R8");
	// Open the node
	sol->group_id = H5Gopen2(zone->flow_id, solname, H5P_DEFAULT);
	if(sol->group_id<0) cgp_doError;

	// Set the rank of the dimensions
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Dimensions of the array
	hsize_t dims[rank];
	// Set these to correspond with the size of the zone
	if(location==Vertex) for(k=0;k<rank;k++) dims[k] = zone->nijk[(rank-1)-k];
	else if(location==CellCenter) for(k=0;k<rank;k++) dims[k] = zone->nijk[rank+(rank-1)-k];

	// Create a shape for the array
	hid_t shape_id = H5Screate_simple(rank,dims,NULL);
	if(shape_id<0) cgp_doError;
	// Property list for dataset
	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
	if(preallocate) herr = H5Pset_alloc_time(plist_id,H5D_ALLOC_TIME_EARLY);
	if(preallocate) herr = H5Pset_fill_time(plist_id, H5D_FILL_TIME_ALLOC);
	// Create the array in the file
	hid_t data_id = H5Dcreate2(sol->group_id, " data", H5T_NATIVE_DOUBLE, shape_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Close the array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int cgp_nsols(int fn, int B, int Z, int* nsols) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	*nsols = zone->nsols;
	return 0;
	}

int cgp_section_write(int fn, int B, int Z, char* sectionname, ElementType_t type,
	int start, int end, int nbndry, int* S) {

	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);

	// Index
	int idx;
	// Loop index
	int k;

	// If the section already exist, replace it
	if(node_exists(zone->grid_id, sectionname)) {
		// Find the section to be replaced
		for(k=0;k<zone->nsections;k++) if(strcmp(sectionname,zone->sections[k].sectionname)==0) idx = k;
		// Delete it
		del_node(zone->grid_id,sectionname);
		}
	// If the section does not exist, prepare the memory structures
	else {
		// New list of sections
		section_t* sections = (section_t*) malloc((zone->nsections+1)*sizeof(section_t));
		if(sections==NULL) cgp_doError;
		// Copy old list to new list
		for(k=0;k<zone->nsections;k++) sections[k] = zone->sections[k];
		// Free the old list
		if(zone->sections!=NULL) free(zone->sections);
		// Point the zone to the new list
		zone->sections = sections;
		// Incement the number of coords
		zone->nsections++;
		// Set the coords index
		idx = zone->nsections-1;
		}
	// Pointer to coords
	section_t* section = &(zone->sections[idx]);
	// Copy the name
	strcpy(section->sectionname, sectionname);
	// Set the index
	section->idx = idx;
	// Set the type
	section->type = type;
	// Set the range
	section->range[0] = start;
	section->range[1] = end;
	// Return the index
	*S = idx+1;

	// Create the section in the file
	new_node(zone->group_id, sectionname, "Elements_t", "I4");
	// Open the new zone
	section->group_id = H5Gopen2(zone->group_id, sectionname, H5P_DEFAULT);
	if(section->group_id<0) cgp_doError;

	// Write data
	{
		// Size of array to write
		hsize_t dim[1] = {2};
		// Create shape object for array to write
		hid_t shape_id = H5Screate_simple(1,dim, NULL);
		if(shape_id<0) cgp_doError;
		// Create the data array in the file
		hid_t data_id = H5Dcreate2(section->group_id, " data", H5T_NATIVE_INT, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(data_id<0) cgp_doError;
		// Write the data to the array
		int data[2] = {type,0};
		herr = H5Dwrite(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, data);
		if(herr<0) cgp_doError;
		// Close the data array
		herr = H5Dclose(data_id);
		if(herr<0) cgp_doError;
		// Destroy the shape object
		herr = H5Sclose(shape_id);
		if(herr<0) cgp_doError;
		}

	// Create the ElementConnectivity in the section
	new_node(section->group_id, "ElementConnectivity", "DataArray_t", "I4");
	// Open the new zone
	section->connectivity_id = H5Gopen2(section->group_id, "ElementConnectivity", H5P_DEFAULT);
	if(section->connectivity_id<0) cgp_doError;

	// Create array for data
	// We choose not to support the MIXED type
	// The array size is (number of elements)*(nodes per element)
	// Write data
	{
		// Size of array to write
		hsize_t dim[1] = {(end-start+1)*node_counts[type]};
		// Create shape object for array to write
		hid_t shape_id = H5Screate_simple(1,dim, NULL);
		if(shape_id<0) cgp_doError;
		// Property list for dataset
		hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
		if(preallocate) herr = H5Pset_alloc_time(plist_id,H5D_ALLOC_TIME_EARLY);
		if(preallocate) herr = H5Pset_fill_time(plist_id, H5D_FILL_TIME_ALLOC);
		// Create the data array in the file
		hid_t data_id = H5Dcreate2(section->connectivity_id, " data", H5T_NATIVE_INT, shape_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		if(data_id<0) cgp_doError;
		// Close the data array
		herr = H5Dclose(data_id);
		if(herr<0) cgp_doError;
		// Close the property list
		herr = H5Pclose(plist_id);
		if(herr<0) cgp_doError;
		// Destroy the shape object
		herr = H5Sclose(shape_id);
		if(herr<0) cgp_doError;
		}

	// Create the ElementRange in the section
	new_node(section->group_id, "ElementRange", "IndexRange_t", "I4");
	// Open the new zone
	section->range_id = H5Gopen2(section->group_id, "ElementRange", H5P_DEFAULT);
	if(section->range_id<0) cgp_doError;

	// Write data
	{
		// Size of array to write
		hsize_t dim[1] = {2};
		// Create shape object for array to write
		hid_t shape_id = H5Screate_simple(1,dim, NULL);
		if(shape_id<0) cgp_doError;
		// Create the data array in the file
		hid_t data_id = H5Dcreate2(section->range_id, " data", H5T_NATIVE_INT, shape_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if(data_id<0) cgp_doError;
		// Write the data to the array
		int data[2] = {start,end};
		herr = H5Dwrite(data_id, H5T_NATIVE_INT, shape_id, shape_id, H5P_DEFAULT, data);
		if(herr<0) cgp_doError;
		// Close the data array
		herr = H5Dclose(data_id);
		if(herr<0) cgp_doError;
		// Destroy the shape object
		herr = H5Sclose(shape_id);
		if(herr<0) cgp_doError;
		}

	return 0;
	}

int cgp_section_write_data(int fn, int B, int Z, int S, int min, int max, int *elements) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	if(S>files[fn].bases[B-1].zones[Z-1].nsections||S<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Pointer to the current coords
	section_t* section = &(zone->sections[S-1]);

	// Open the data
	hid_t data_id = H5Dopen2(section->connectivity_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = 1;
	// Set the start position for the data write
	hsize_t start[rank];
	start[0] = (min-section->range[0])*node_counts[section->type];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	dims[0] = (max-min+1)*node_counts[section->type];
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Write the data in collective parallel I/O
	herr = H5Dwrite(data_id, H5T_NATIVE_INT, mem_shape_id, data_shape_id, plist_id, elements);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;

	return 0;
	}

int cgp_section_read_data(int fn, int B, int Z, int S, int min, int max, int *elements) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	if(S>files[fn].bases[B-1].zones[Z-1].nsections||S<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Pointer to the current coords
	section_t* section = &(zone->sections[S-1]);

	// Open the data
	hid_t data_id = H5Dopen2(section->connectivity_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = 1;
	// Set the start position for the data write
	hsize_t start[rank];
	start[0] = (min-section->range[0])*node_counts[section->type];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	dims[0] = (max-min+1)*node_counts[section->type];
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Read the data in collective parallel I/O
	herr = H5Dread(data_id, H5T_NATIVE_INT, mem_shape_id, data_shape_id, plist_id, elements);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;

	return 0;
	}

//= Array IO Prototypes =//

int cgp_array_write(int fn, int B, int Z, char *arrayname, GridLocation_t location) {
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);
	// Index
	int idx;
	// Loop index
	int k;
	// If the solution already exist, replace it
	if(node_exists(zone->group_id, arrayname)) {
		// Delete them
		del_node(zone->group_id,arrayname);
		}

	// Create a new node in the file
	new_node(zone->group_id, arrayname, "DataArray_t", "R8");
	// Open the node
	hid_t group_id = H5Gopen2(zone->group_id, arrayname, H5P_DEFAULT);
	if(group_id<0) cgp_doError;

	// Set the rank of the dimensions
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Dimensions of the array
	hsize_t dims[rank];
	// Set these to correspond with the size of the zone
	if(location==Vertex) for(k=0;k<rank;k++) dims[k] = zone->nijk[(rank-1)-k];
	else if(location==CellCenter) for(k=0;k<rank;k++) dims[k] = zone->nijk[rank+(rank-1)-k];

	// Create a shape for the array
	hid_t shape_id = H5Screate_simple(rank,dims,NULL);
	if(shape_id<0) cgp_doError;
	// Property list for dataset
	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
	if(preallocate) herr = H5Pset_alloc_time(plist_id,H5D_ALLOC_TIME_EARLY);
	if(preallocate) herr = H5Pset_fill_time(plist_id, H5D_FILL_TIME_ALLOC);
	// Create the array in the file
	hid_t data_id = H5Dcreate2(group_id, " data", H5T_NATIVE_DOUBLE, shape_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	if(data_id<0) cgp_doError;
	// Close the array
	herr = H5Dclose(data_id);
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape
	herr = H5Sclose(shape_id);
	if(herr<0) cgp_doError;
	// Close the group
	herr = H5Gclose(group_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int cgp_array_write_data(int fn, int B, int Z, char* arrayname, int* min, int* max, void* data) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);

	// Open the group
	hid_t group_id = H5Gopen2(zone->group_id, arrayname, H5P_DEFAULT);
	// Open the data
	hid_t data_id = H5Dopen2(group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Set the start position for the data write
	hsize_t start[rank];
	for(k=0;k<rank;k++) start[k] = min[k];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	for(k=0;k<rank;k++) dims[k] = max[k]-min[k]+1;
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Write the data in collective parallel I/O
	herr = H5Dwrite(data_id, H5T_NATIVE_DOUBLE, mem_shape_id, data_shape_id, plist_id, data);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	herr = H5Gclose(group_id);
	if(herr<0) cgp_doError;
	return 0;
	}

int cgp_array_read_data(int fn, int B, int Z, char* arrayname, int* min, int* max, void* data) {
	int k;
	if(fn>=files_count||fn<0) cgp_doError;
	if(!files[fn].isOpen) cgp_doError;
	if(B>files[fn].nbases||B<=0) cgp_doError;
	if(Z>files[fn].bases[B-1].nzones||Z<=0) cgp_doError;
	herr_t herr;
	// Pointer to the current file
	file_t* file = &(files[fn]);
	// Pointer to the current base
	base_t* base = &(file->bases[B-1]);
	// Pointer to the current zone
	zone_t* zone = &(base->zones[Z-1]);

	// Open the group
	hid_t group_id = H5Gopen2(zone->group_id, arrayname, H5P_DEFAULT);
	// Open the data
	hid_t data_id = H5Dopen2(group_id, " data", H5P_DEFAULT);
	if(data_id<0) cgp_doError;

	// Set the rank of the data
	hsize_t rank = (zone->type==Structured)?base->cell_dim:1;
	// Set the start position for the data write
	hsize_t start[rank];
	for(k=0;k<rank;k++) start[k] = min[k];
	// Compute the counts in each dimension
	hsize_t dims[rank];
	for(k=0;k<rank;k++) dims[k] = max[k]-min[k]+1;
	// Create a shape for the data in memory
	hid_t mem_shape_id = H5Screate_simple(rank,dims,NULL);
	if(mem_shape_id<0) cgp_doError;
	// Create a shape for the data in the file
	hid_t data_shape_id = H5Dget_space(data_id);
	if(data_shape_id<0) cgp_doError;
	// Select a section of the array in the file
	herr = H5Sselect_hyperslab(data_shape_id, H5S_SELECT_SET, start, NULL, dims, NULL);
	if(herr<0) cgp_doError;
	// Set the access property list for data transfer
	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	if(plist_id<0) cgp_doError;
	// Set MPI-IO collective communication
	herr = H5Pset_dxpl_mpio(plist_id, file->collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
	if(herr<0) cgp_doError;
	// Read the data in collective parallel I/O
	herr = H5Dread(data_id, H5T_NATIVE_DOUBLE, mem_shape_id, data_shape_id, plist_id, data);
	herr = H5Fflush(file->file_id, H5F_SCOPE_GLOBAL);
	if(herr<0) cgp_doError;
	// Close the property list
	herr = H5Pclose(plist_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in the file
	herr = H5Sclose(data_shape_id);
	if(herr<0) cgp_doError;
	// Close the shape of the data in memory
	herr = H5Sclose(mem_shape_id);
	if(herr<0) cgp_doError;
	// Close the data array
	herr = H5Dclose(data_id);
	if(herr<0) cgp_doError;
	herr = H5Gclose(group_id);
	if(herr<0) cgp_doError;
	return 0;
	}

// The mallocs and memory copies in here need to be verified
int queue_slice_write(SliceType_t type, int F, int B, int Z, void* SN, int rank,
	int* min, int* max, void* data) {
	int k;
	slice_t* slices = malloc((write_queue_len+1)*sizeof(slice_t));
	for(k=0;k<write_queue_len;k++) slices[k] = write_queue[k];
	slice_t* slice = &(slices[write_queue_len]);

	slice->type = type;
	slice->F = F;
	slice->B = B;
	slice->Z = Z;
	slice->rank = rank;

	slice->min = (int*) malloc(rank*sizeof(int));
	for(k=0;k<rank;k++) slice->min[k] = min[k];
	slice->max = (int*) malloc(rank*sizeof(int));
	for(k=0;k<rank;k++) slice->max[k] = max[k];
	slice->data = data;
	if(type!=Array&&type!=Empty) slice->Selector = *((int*) SN);
	else if(type==Array) strcpy(slice->name,(char*) SN);
	else {}

	if(write_queue!=NULL) free(write_queue);
	write_queue = slices;
	write_queue_len++;
	return 0;
	}

// The mallocs and memory copies in here need to be verified
int queue_flush(void) {
	int err;
	herr_t herr;
	int i,j,k;
	int max_queue_len = 0;
	int world_rank = 0;
	int world_size = 0;

	err = MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	err = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	err = MPI_Allreduce(&write_queue_len, &max_queue_len, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	if(max_queue_len==0) max_queue_len = 1;

	slice_t queue[max_queue_len];
	for(k=0;k<write_queue_len;k++) {
		queue[k] = write_queue[k];
		if(queue[k].type==Array) strcpy(queue[k].name,write_queue[k].name);
		}
	for(k=write_queue_len;k<max_queue_len;k++) {
		queue[k].type = Empty;
		queue[k].F = -1;
		}

	int fn = queue[0].F;
	int Fs[world_size];
	err = MPI_Allgather(&fn, 1,MPI_INT, Fs, world_size, MPI_INT, MPI_COMM_WORLD);

	for(k=0;k<max_queue_len;k++) {
		hid_t data_id;
		switch(queue[k].type) {
			case(Empty):
				data_id = H5Dopen2(files[fn].file_id, " dummy", H5P_DEFAULT);
				hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
				herr = H5Pset_dxpl_mpio(plist_id, files[fn].collective==TRUE?H5FD_MPIO_COLLECTIVE:H5FD_MPIO_INDEPENDENT);
				hsize_t dim = 1;
				hid_t shape_id = H5Screate_simple(1,&dim, NULL);
				int buf = 0;
				H5Dwrite(data_id, H5T_NATIVE_INT, shape_id, shape_id, plist_id, &buf);
				H5Sclose(shape_id);
				H5Pclose(plist_id);
				H5Dclose(data_id);
				break;
			case(Coords):
				err = cgp_coord_write_data(queue[k].F,queue[k].B,queue[k].Z,queue[k].Selector,queue[k].min,queue[k].max,queue[k].data);
				break;
			case(Elements):
				err = cgp_section_write_data(queue[k].F,queue[k].B,queue[k].Z,queue[k].Selector,queue[k].min[0],queue[k].max[0],queue[k].data);
				break;
			case(Solution):
				err = cgp_sol_write_data(queue[k].F,queue[k].B,queue[k].Z,queue[k].Selector,queue[k].min,queue[k].max,queue[k].data);
				break;
			case(Array):
				err = cgp_array_write_data(queue[k].F,queue[k].B,queue[k].Z,queue[k].name,queue[k].min,queue[k].max,queue[k].data);
				break;
			}
		}
	for(k=0;k<write_queue_len;k++) {
		if(write_queue[k].min!=NULL) free(write_queue[k].min);
		if(write_queue[k].max!=NULL) free(write_queue[k].max);
		}
	if(write_queue!=NULL) free(write_queue);
	write_queue = NULL;
	write_queue_len = 0;
	return 0;
	}
