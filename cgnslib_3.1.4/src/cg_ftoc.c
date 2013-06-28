/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "fortran_macros.h"
#include "cgnslib.h"
#include "cgns_header.h"
#include <string.h>
#include "cgns_io.h"
#ifdef MEM_DEBUG
#include "cg_malloc.h"
#endif

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Convert between Fortran and C strings                            *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

static void string_2_C_string(char *string, int string_length,
	char *c_string, int max_len, cgsize_t *ierr)
{
    int i, iend;

    if (string == NULL || c_string == NULL) {
        cgi_error ("NULL string pointer");
        *ierr = CG_ERROR;
        return;
    }

    /** Skip and trailing blanks **/
    for (iend = string_length-1; iend >= 0; iend--) {
        if (string[iend] != ' ') break;
    }
    if (iend >= max_len) iend = max_len - 1;

    /** Copy the non-trailing blank portion of the string **/
    for (i = 0; i <= iend; i++)
        c_string[i] = string[i];

    /** NULL terminate the C string **/
    c_string[i] = '\0';
    *ierr = CG_OK;
}

/*-----------------------------------------------------------------------*/

static void string_2_F_string(char *c_string, char *string,
	int string_length, cgsize_t *ierr)
{
    int i, len;

    if (c_string == NULL || string == NULL) {
        cgi_error ("NULL string pointer");
        *ierr = CG_ERROR;
        return;
    }
    len = (int)strlen(c_string);
    if (len > string_length) len = string_length;

    for (i = 0; i < len; i++)
        string[i] = c_string[i];
    while (i < string_length)
        string[i++] = ' ';
    *ierr = CG_OK;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      LIBRARY FUNCTIONS                                                *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_is_cgns_f, CG_IS_CGNS_F) (STR_PSTR(filename),
	cgsize_t *file_type, cgsize_t *ier STR_PLEN(filename))
{
    int length, i_file_type;
    char *c_name;

    length = (int) STR_LEN(filename);
    c_name = CGNS_NEW(char, length+1);

    string_2_C_string(STR_PTR(filename), STR_LEN(filename), c_name, length, ier);
    if (*ier == 0) {
        *ier = cg_is_cgns(c_name, &i_file_type);
        *file_type = i_file_type;
    }
    CGNS_FREE(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_open_f, CG_OPEN_F) (STR_PSTR(filename), cgsize_t *mode,
	cgsize_t *fn, cgsize_t *ier STR_PLEN(filename))
{
    int length, i_fn;
    char *c_name;

    length = (int) STR_LEN(filename);
    c_name = CGNS_NEW(char, length+1);

    string_2_C_string(STR_PTR(filename), STR_LEN(filename), c_name, length, ier);
    if (*ier == 0) {
#if DEBUG_FTOC
        printf("filename='%s'\n",c_name);
#endif
        *ier = cg_open(c_name, (int)*mode, &i_fn);
        *fn  = i_fn;
    }
    free(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_version_f, CG_VERSION_F) (cgsize_t *fn,
	float *FileVersion, cgsize_t *ier)
{
    *ier = cg_version((int)*fn, FileVersion);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_precision_f, CG_PRECISION_F) (cgsize_t *fn,
	cgsize_t *precision, cgsize_t *ier)
{
    int i_precision;

    *ier = cg_precision((int)*fn, &i_precision);
    *precision = i_precision;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_close_f, CG_CLOSE_F) (cgsize_t *fn, cgsize_t *ier)
{
    *ier = cg_close((int)*fn);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_save_as_f, CG_SAVE_AS_F) (cgsize_t *fn,
	STR_PSTR(filename), cgsize_t *file_type, cgsize_t *follow_links,
	cgsize_t *ier STR_PLEN(filename))
{
    int length;
    char *c_name;

    length = (int) STR_LEN(filename);
    c_name = CGNS_NEW(char, length+1);

    string_2_C_string(STR_PTR(filename), STR_LEN(filename), c_name, length, ier);
    if (*ier == 0)
        *ier = cg_save_as((int)*fn, c_name, (int)*file_type, (int)*follow_links);
    free(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_set_file_type_f, CG_SET_FILE_TYPE_F) (
	cgsize_t *ft, cgsize_t *ier)
{
    *ier = cg_set_file_type((int)*ft);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_get_file_type_f, CG_GET_FILE_TYPE_F) (
	cgsize_t *fn, cgsize_t *ft, cgsize_t *ier)
{
    int i_ft;

    *ier = cg_get_file_type((int)*fn, &i_ft);
    *ft = i_ft;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_set_compress_f, CG_SET_COMPRESS_F) (
	cgsize_t *cmpr, cgsize_t *ier)
{
    *ier = cg_set_compress((int)*cmpr);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_get_compress_f, CG_GET_COMPRESS_F) (
	cgsize_t *cmpr, cgsize_t *ier)
{
    int i_cmpr;

    *ier = cg_get_compress(&i_cmpr);
    *cmpr = i_cmpr;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_set_path_f, CG_SET_PATH_F) (STR_PSTR(pathname),
	cgsize_t *ier STR_PLEN(pathname))
{
    int length;
    char *c_name;

    length = (int) STR_LEN(pathname);
    c_name = CGNS_NEW(char, length+1);

    string_2_C_string(STR_PTR(pathname), STR_LEN(pathname), c_name, length, ier);
    if (*ier == 0)
        *ier = cg_set_path(c_name);
    free(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_add_path_f, CG_ADD_PATH_F) (STR_PSTR(pathname),
	cgsize_t *ier STR_PLEN(pathname))
{
    int length;
    char *c_name;

    length = (int) STR_LEN(pathname);
    c_name = CGNS_NEW(char, length+1);

    string_2_C_string(STR_PTR(pathname), STR_LEN(pathname), c_name, length, ier);
    if (*ier == 0)
        *ier = cg_add_path(c_name);
    free(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_get_cgio_f, CG_GET_CGIO_F) (cgsize_t *fn,
	cgsize_t *cgio_num, cgsize_t *ier)
{
    int i_cgio_num;

    *ier = cg_get_cgio((int)*fn, &i_cgio_num);
    *cgio_num = i_cgio_num;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_root_id_f, CG_ROOT_ID_F) (cgsize_t *fn,
	double *rootid, cgsize_t *ier)
{
    *ier = cg_root_id((int)*fn, rootid);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write CGNSBase_t Nodes                                  *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nbases_f, CG_NBASES_F) (cgsize_t *fn,
	cgsize_t *nbases, cgsize_t *ier)
{
    int i_nbases;

    *ier = cg_nbases((int)*fn, &i_nbases);
    *nbases = i_nbases;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_base_read_f, CG_BASE_READ_F) (cgsize_t *fn, cgsize_t *B,
	STR_PSTR(basename), cgsize_t *cell_dim, cgsize_t *phys_dim,
	cgsize_t *ier STR_PLEN(basename))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_cell_dim, i_phys_dim;

    *ier = cg_base_read((int)*fn, (int)*B, c_name, &i_cell_dim, &i_phys_dim);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(basename), STR_LEN(basename), ier);
    *cell_dim = i_cell_dim;
    *phys_dim = i_phys_dim;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_base_id_f, CG_BASE_ID_F) (cgsize_t *fn, cgsize_t *B,
	double *base_id, cgsize_t *ier)
{
    *ier = cg_base_id((int)*fn, (int)*B, base_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_base_write_f, CG_BASE_WRITE_F) (cgsize_t *fn,
	STR_PSTR(basename), cgsize_t *cell_dim, cgsize_t *phys_dim,
	cgsize_t *B, cgsize_t *ier STR_PLEN(basename))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_B;

    string_2_C_string(STR_PTR(basename), STR_LEN(basename),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("\nbasename='%s'\n", c_name);
    printf("cell_dim=%d\n",*cell_dim);
    printf("phys_dim=%d\n",*phys_dim);
#endif
    *ier = cg_base_write((int)*fn, c_name, (int)*cell_dim, (int)*phys_dim, &i_B);
    *B = i_B;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_cell_dim_f, CG_CELL_DIM_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *dim, cgsize_t *ier)
{
    int i_dim;

    *ier = cg_cell_dim((int)*fn, (int)*B, &i_dim);
    *dim = i_dim;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write Zone_t Nodes                                      *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nzones_f, CG_NZONES_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *nzones, cgsize_t *ier)
{
    int i_nzones;

    *ier = cg_nzones((int)*fn, (int)*B, &i_nzones);
    *nzones = i_nzones;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zone_type_f, CG_ZONE_TYPE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *type, cgsize_t *ier)
{
    CGNS_ENUMT(ZoneType_t) i_type;

    *ier = cg_zone_type((int)*fn, (int)*B, (int)*Z, &i_type);
    *type = i_type;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zone_read_f, CG_ZONE_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(zonename), cgsize_t *size,
	cgsize_t *ier STR_PLEN(zonename))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_zone_read((int)*fn, (int)*B, (int)*Z, c_name, size);
    if (*ier == 0)
        string_2_F_string(c_name, STR_PTR(zonename), STR_LEN(zonename), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zone_id_f, CG_ZONE_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, double *zone_id, cgsize_t *ier)
{
    *ier = cg_zone_id((int)*fn, (int)*B, (int)*Z, zone_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zone_write_f, CG_ZONE_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	STR_PSTR(zonename), cgsize_t *size, CGNS_ENUMT(ZoneType_t) *type,
	cgsize_t *Z, cgsize_t *ier STR_PLEN(zonename))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_Z;

    string_2_C_string(STR_PTR(zonename), STR_LEN(zonename),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("\n  zonename='%s'\n", c_name);
#endif
    *ier = cg_zone_write((int)*fn, (int)*B, c_name, size,
               (CGNS_ENUMT(ZoneType_t))*type, &i_Z);
    *Z = i_Z;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_index_dim_f, CG_INDEX_DIM_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *dim, cgsize_t *ier)
{
    int i_dim;

    *ier = cg_index_dim((int)*fn, (int)*B, (int)*Z, &i_dim);
    *dim = i_dim;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write Family_t Nodes                                    *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nfamilies_f, CG_NFAMILIES_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *nfamilies, cgsize_t *ier)
{
    int i_nfamilies;

    *ier = cg_nfamilies((int)*fn, (int)*B, &i_nfamilies);
    *nfamilies = i_nfamilies;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_family_read_f, CG_FAMILY_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *F, STR_PSTR(family_name), cgsize_t *nboco, cgsize_t *ngeos,
	cgsize_t *ier STR_PLEN(family_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_nboco, i_ngeos;

    *ier = cg_family_read((int)*fn, (int)*B, (int)*F, c_name, &i_nboco, &i_ngeos);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(family_name), STR_LEN(family_name), ier);
    *nboco = i_nboco;
    *ngeos = i_ngeos;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_family_write_f, CG_FAMILY_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	STR_PSTR(family_name), cgsize_t *F, cgsize_t *ier STR_PLEN(family_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_F;

    string_2_C_string(STR_PTR(family_name), STR_LEN(family_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_family_write((int)*fn, (int)*B, c_name, &i_F);
    *F = i_F;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write FamBC_t Nodes                                     *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_fambc_read_f, CG_FAMBC_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *F, cgsize_t *BC, STR_PSTR(fambc_name), cgsize_t *bocotype,
	cgsize_t *ier STR_PLEN(fambc_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(BCType_t) i_bocotype;

    *ier = cg_fambc_read((int)*fn, (int)*B, (int)*F, (int)*BC,
               c_name, &i_bocotype);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(fambc_name), STR_LEN(fambc_name), ier);
    *bocotype = i_bocotype;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_fambc_write_f, CG_FAMBC_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *F, STR_PSTR(fambc_name), cgsize_t *bocotype,
	cgsize_t *BC, cgsize_t *ier STR_PLEN(fambc_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_BC;

    string_2_C_string(STR_PTR(fambc_name), STR_LEN(fambc_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_fambc_write((int)*fn, (int)*B, (int)*F, c_name,
               (CGNS_ENUMT(BCType_t))*bocotype, &i_BC);
    *BC = i_BC;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GeometryReference_t Nodes                         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_geo_read_f, CG_GEO_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *F, cgsize_t *G, STR_PSTR(geo_name), STR_PSTR(geo_file),
	STR_PSTR(CAD_name), cgsize_t *npart, cgsize_t *ier STR_PLEN(geo_name)
	STR_PLEN(geo_file) STR_PLEN(CAD_name))
{
    char c_geo_name[CGIO_MAX_NAME_LENGTH+1];
    char c_CAD_name[CGIO_MAX_NAME_LENGTH+1];
    char *c_geo_file;
    int i_npart;

    *ier = cg_geo_read((int)*fn, (int)*B, (int)*F, (int)*G, c_geo_name,
               &c_geo_file, c_CAD_name, &i_npart);
    if (*ier) return;
    *npart = i_npart;
    string_2_F_string(c_geo_file, STR_PTR(geo_file), STR_LEN(geo_file), ier);
    free(c_geo_file);
    if (*ier) return;
    string_2_F_string(c_geo_name, STR_PTR(geo_name), STR_LEN(geo_name), ier);
    if (*ier) return;
    string_2_F_string(c_CAD_name, STR_PTR(CAD_name), STR_LEN(CAD_name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_geo_write_f, CG_GEO_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *F, STR_PSTR(geo_name), STR_PSTR(geo_file), STR_PSTR(CAD_name),
	cgsize_t *G, cgsize_t *ier STR_PLEN(geo_name) STR_PLEN(geo_file)
	STR_PLEN(CAD_name))
{
    char c_geo_name[CGIO_MAX_NAME_LENGTH+1];
    char c_CAD_name[CGIO_MAX_NAME_LENGTH+1];
    char *c_geo_file;
    int length, i_G;

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(geo_name), STR_LEN(geo_name),
        c_geo_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    string_2_C_string(STR_PTR(CAD_name), STR_LEN(CAD_name),
        c_CAD_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;

    length = STR_LEN(geo_file);
    c_geo_file = CGNS_NEW(char, length+1);
    string_2_C_string(STR_PTR(geo_file), STR_LEN(geo_file),
        c_geo_file, length, ier);
    if (*ier == 0) {
        *ier = cg_geo_write((int)*fn, (int)*B, (int)*F, c_geo_name,
                   c_geo_file, c_CAD_name, &i_G);
        *G = i_G;
    }
    free(c_geo_file);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GeometryEntity_t Nodes                            *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_part_read_f, CG_PART_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *F, cgsize_t *G, cgsize_t *P, STR_PSTR(part_name),
	cgsize_t *ier STR_PLEN(part_name))
{
    char c_part_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_part_read((int)*fn, (int)*B, (int)*F, (int)*G, (int)*P, c_part_name);
    if (*ier == 0)
        string_2_F_string(c_part_name, STR_PTR(part_name), STR_LEN(part_name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_part_write_f, CG_PART_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *F, cgsize_t *G, STR_PSTR(part_name), cgsize_t *P,
	cgsize_t *ier STR_PLEN(part_name))
{
    char c_part_name[CGIO_MAX_NAME_LENGTH+1];
    int i_P;

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(part_name), STR_LEN(part_name),
        c_part_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_part_write((int)*fn, (int)*B, (int)*F, (int)*G, c_part_name, &i_P);
    *P = i_P;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write DiscreteData_t Nodes                              *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_ndiscrete_f, CG_NDISCRETE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *ndiscrete, cgsize_t *ier)
{
    int i_ndiscrete;

    *ier = cg_ndiscrete((int)*fn, (int)*B, (int)*Z, &i_ndiscrete);
    *ndiscrete = i_ndiscrete;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_discrete_read_f, CG_DISCRETE_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *D, STR_PSTR(discrete_name),
	cgsize_t *ier STR_PLEN(discrete_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_discrete_read((int)*fn, (int)*B, (int)*Z, (int)*D, c_name);
    if (*ier == 0)
        string_2_F_string(c_name, STR_PTR(discrete_name), STR_LEN(discrete_name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_discrete_write_f, CG_DISCRETE_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(discrete_name), cgsize_t *D,
	cgsize_t *ier STR_PLEN(discrete_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_D;

    string_2_C_string(STR_PTR(discrete_name), STR_LEN(discrete_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("    discrete_name='%s'\n", c_name);
#endif
    *ier = cg_discrete_write((int)*fn, (int)*B, (int)*Z, c_name, &i_D);
    *D = i_D;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_discrete_size_f, CG_DISCRETE_SIZE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *D, cgsize_t *ndim,
	cgsize_t *dims, cgsize_t *ier)
{
    int i_ndim;

    *ier = cg_discrete_size((int)*fn, (int)*B, (int)*Z, (int)*D,
                &i_ndim, dims);
    *ndim = i_ndim;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_discrete_ptset_info_f, CG_DISCRETE_PTSET_INFO_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	cgsize_t *ptype, cgsize_t *npnts, cgsize_t *ier)
{
    CGNS_ENUMT(PointSetType_t) i_ptype;

    *ier = cg_discrete_ptset_info((int)*fn, (int)*B, (int)*Z,
               (int)*S, &i_ptype, npnts);
    *ptype = i_ptype;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_discrete_ptset_read_f, CG_DISCRETE_PTSET_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_discrete_ptset_read((int)*fn, (int)*B, (int)*Z,
               (int)*S, pnts);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_discrete_ptset_write_f, CG_DISCRETE_PTSET_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(name),
	cgsize_t *location, cgsize_t *ptype, cgsize_t *npnts,
	cgsize_t *pnts, cgsize_t *S, cgsize_t *ier STR_PLEN(name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_S;

    string_2_C_string(STR_PTR(name), STR_LEN(name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;

    *ier = cg_discrete_ptset_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(GridLocation_t))*location,
               (CGNS_ENUMT(PointSetType_t))*ptype, *npnts, pnts, &i_S);
    *S = i_S;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GridCoordinates_t/DataArray_t Nodes               *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_ncoords_f, CG_NCOORDS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *ncoords, cgsize_t *ier)
{
    int i_ncoords;

    *ier = cg_ncoords((int)*fn, (int)*B, (int)*Z, &i_ncoords);
    *ncoords = i_ncoords;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_coord_info_f, CG_COORD_INFO_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *C, cgsize_t *type, STR_PSTR(coordname),
	cgsize_t *ier STR_PLEN(coordname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(DataType_t) i_type;

    *ier = cg_coord_info((int)*fn, (int)*B, (int)*Z, (int)*C, &i_type, c_name);
    if (*ier) return;
    *type = i_type;
    string_2_F_string(c_name, STR_PTR(coordname), STR_LEN(coordname), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_coord_read_f, CG_COORD_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(coordname), cgsize_t *type, cgsize_t *rmin,
	cgsize_t *rmax, void *coord, cgsize_t *ier STR_PLEN(coordname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(coordname), STR_LEN(coordname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("coordname='%s'\n", c_name);
#endif
    *ier = cg_coord_read((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(DataType_t))*type, rmin, rmax, coord);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_coord_id_f, CG_COORD_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *C, double *coord_id, cgsize_t *ier)
{
    *ier = cg_coord_id((int)*fn, (int)*B, (int)*Z, (int)*C, coord_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_coord_write_f, CG_COORD_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *type, STR_PSTR(coordname), void *coord, cgsize_t *C,
	cgsize_t *ier STR_PLEN(coordname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_C;

    string_2_C_string(STR_PTR(coordname), STR_LEN(coordname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("    coordname='%s'\n", c_name);
#endif
    *ier = cg_coord_write((int)*fn, (int)*B, (int)*Z,
               (CGNS_ENUMT(DataType_t))*type, c_name, coord, &i_C);
    *C = i_C;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_coord_partial_write_f, CG_COORD_PARTIAL_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *type, STR_PSTR(coordname),
	cgsize_t *rmin, cgsize_t *rmax, void *coord, cgsize_t *C,
	cgsize_t *ier STR_PLEN(coordname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_C;

    string_2_C_string(STR_PTR(coordname), STR_LEN(coordname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("    coordname='%s'\n", c_name);
#endif
    *ier = cg_coord_partial_write((int)*fn, (int)*B, (int)*Z,
               (CGNS_ENUMT(DataType_t))*type, c_name, rmin, rmax,
               coord, &i_C);
    *C = i_C;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write Elements_t Nodes                                  *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nsections_f, CG_NSECTIONS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *nsections, cgsize_t *ier)
{
    int i_nsections;

    *ier = cg_nsections((int)*fn, (int)*B, (int)*Z, &i_nsections);
    *nsections = i_nsections;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_section_read_f, CG_SECTION_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *E, STR_PSTR(section_name),
	cgsize_t *type, cgsize_t *start, cgsize_t *end, cgsize_t *nbndry,
	cgsize_t *parent_flag, cgsize_t *ier STR_PLEN(section_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(ElementType_t) i_type;
    int i_nbndry, i_parent_flag;

    *ier = cg_section_read((int)*fn, (int)*B, (int)*Z, (int)*E, c_name,
               &i_type, start, end, &i_nbndry, &i_parent_flag);
    if (*ier) return;
    *type = i_type;
    *nbndry = i_nbndry;
    *parent_flag = i_parent_flag;
    string_2_F_string(c_name, STR_PTR(section_name), STR_LEN(section_name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_elements_read_f, CG_ELEMENTS_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *E, cgsize_t *elements,
	cgsize_t *parent_data, cgsize_t *ier)
{
    *ier = cg_elements_read((int)*fn, (int)*B, (int)*Z, (int)*E,
               elements, parent_data);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_elementdatasize_f, CG_ELEMENTDATASIZE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *E, cgsize_t *ElementDataSize,
	cgsize_t *ier)
{
    *ier = cg_ElementDataSize((int)*fn, (int)*B, (int)*Z, (int)*E,
               ElementDataSize);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_elementpartialsize_f, CG_ELEMENTPARTIALSIZE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *E, cgsize_t *start, cgsize_t *end,
	cgsize_t *ElementDataSize, cgsize_t *ier)
{
    *ier = cg_ElementPartialSize((int)*fn, (int)*B, (int)*Z, (int)*E,
               *start, *end, ElementDataSize);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_section_write_f, CG_SECTION_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(section_name), cgsize_t *type,
	cgsize_t *start, cgsize_t *end, cgsize_t *nbndry, cgsize_t *elements,
	cgsize_t *S, cgsize_t *ier STR_PLEN(section_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_S;

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(section_name), STR_LEN(section_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_section_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(ElementType_t))*type, *start, *end,
               (int)*nbndry, elements, &i_S);
    *S = i_S;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_parent_data_write_f, CG_PARENT_DATA_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *S, cgsize_t *parent_data, cgsize_t *ier)
{
    *ier = cg_parent_data_write((int)*fn, (int)*B, (int)*Z, (int)*S, parent_data);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_section_partial_write_f, CG_SECTION_PARTIAL_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(section_name),
	cgsize_t *type, cgsize_t *start, cgsize_t *end, cgsize_t *nbndry,
	cgsize_t *S, cgsize_t *ier STR_PLEN(section_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_S;

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(section_name), STR_LEN(section_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_section_partial_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(ElementType_t))*type, *start,
               *end, (int)*nbndry, &i_S);
    *S = i_S;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_elements_partial_write_f, CG_ELEMENTS_PARTIAL_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S, cgsize_t *rmin,
	cgsize_t *rmax, cgsize_t *elements, cgsize_t *ier)
{
    *ier = cg_elements_partial_write((int)*fn, (int)*B, (int)*Z, (int)*S,
               *rmin, *rmax, elements);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_parent_data_partial_write_f, CG_PARENT_DATA_PARTIAL_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S, cgsize_t *rmin,
	cgsize_t *rmax, cgsize_t *parent_data, cgsize_t *ier)
{
    *ier = cg_parent_data_partial_write((int)*fn, (int)*B, (int)*Z, (int)*S,
               *rmin, *rmax, parent_data);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_elements_partial_read_f, CG_ELEMENTS_PARTIAL_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S, cgsize_t *rmin,
	cgsize_t *rmax, cgsize_t *elements, cgsize_t *parent, cgsize_t *ier)
{
    *ier = cg_elements_partial_read((int)*fn, (int)*B, (int)*Z, (int)*S,
               *rmin, *rmax, elements, parent);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write FlowSolution_t Nodes                              *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nsols_f, CG_NSOLS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *nsols, cgsize_t *ier)
{
    int i_nsols;

    *ier = cg_nsols((int)*fn, (int)*B, (int)*Z, &i_nsols);
    *nsols = i_nsols;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_info_f, CG_SOL_INFO_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, STR_PSTR(solname), cgsize_t *location,
	cgsize_t *ier STR_PLEN(solname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(GridLocation_t) i_location;

    *ier = cg_sol_info((int)*fn, (int)*B, (int)*Z, (int)*S, c_name, &i_location);
    if (*ier) return;
    *location = i_location;
    string_2_F_string(c_name, STR_PTR(solname), STR_LEN(solname), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_id_f, CG_SOL_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, double *sol_id, cgsize_t *ier)
{
    *ier = cg_sol_id((int)*fn, (int)*B, (int)*Z, (int)*S, sol_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_write_f, CG_SOL_WRITE_F)(cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(solname), cgsize_t *location, cgsize_t *S,
	cgsize_t *ier STR_PLEN(solname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_S;

    string_2_C_string(STR_PTR(solname), STR_LEN(solname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("\n    solname='%s'\n", c_name);
#endif
    *ier = cg_sol_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(GridLocation_t))*location, &i_S);
    *S = i_S;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_size_f, CG_SOL_SIZE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *S, cgsize_t *ndim,
	cgsize_t *dims, cgsize_t *ier)
{
    int i_ndim;

    *ier = cg_sol_size((int)*fn, (int)*B, (int)*Z, (int)*S,
                &i_ndim, dims);
    *ndim = i_ndim;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_ptset_info_f, CG_SOL_PTSET_INFO_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	cgsize_t *ptype, cgsize_t *npnts, cgsize_t *ier)
{
    CGNS_ENUMT(PointSetType_t) i_ptype;

    *ier = cg_sol_ptset_info((int)*fn, (int)*B, (int)*Z,
               (int)*S, &i_ptype, npnts);
    *ptype = i_ptype;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_ptset_read_f, CG_SOL_PTSET_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_sol_ptset_read((int)*fn, (int)*B, (int)*Z,
               (int)*S, pnts);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_sol_ptset_write_f, CG_SOL_PTSET_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(name),
	cgsize_t *location, cgsize_t *ptype, cgsize_t *npnts,
	cgsize_t *pnts, cgsize_t *S, cgsize_t *ier STR_PLEN(name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_S;

    string_2_C_string(STR_PTR(name), STR_LEN(name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;

    *ier = cg_sol_ptset_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(GridLocation_t))*location,
               (CGNS_ENUMT(PointSetType_t))*ptype, *npnts, pnts, &i_S);
    *S = i_S;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write solution DataArray_t Nodes                        *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nfields_f, CG_NFIELDS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, cgsize_t *nfields, cgsize_t *ier)
{
    int i_nfields;

    *ier = cg_nfields((int)*fn, (int)*B, (int)*Z, (int)*S, &i_nfields);
    *nfields = i_nfields;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_field_info_f, CG_FIELD_INFO_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, cgsize_t *F, cgsize_t *type, STR_PSTR(fieldname),
	cgsize_t *ier STR_PLEN(fieldname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(DataType_t) i_type;

    *ier = cg_field_info((int)*fn, (int)*B, (int)*Z, (int)*S, (int)*F,
               &i_type, c_name);
    if (*ier) return;
    *type = i_type;
    string_2_F_string(c_name, STR_PTR(fieldname), STR_LEN(fieldname), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_field_read_f, CG_FIELD_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, STR_PSTR(fieldname), cgsize_t *type, cgsize_t *rmin,
	cgsize_t *rmax, void *field_ptr, cgsize_t *ier STR_PLEN(fieldname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(fieldname), STR_LEN(fieldname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("fieldname='%s'\n", c_name);
#endif
    *ier = cg_field_read((int)*fn, (int)*B, (int)*Z, (int)*S, c_name,
               (CGNS_ENUMT(DataType_t))*type, rmin, rmax, field_ptr);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_field_id_f, CG_FIELD_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, cgsize_t *F, double *field_id, cgsize_t *ier)
{
    *ier = cg_field_id((int)*fn, (int)*B, (int)*Z, (int)*S, (int)*F, field_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_field_write_f, CG_FIELD_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *S, cgsize_t *type, STR_PSTR(fieldname), void *field_ptr,
	cgsize_t *F, cgsize_t *ier STR_PLEN(fieldname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_F;

    string_2_C_string(STR_PTR(fieldname), STR_LEN(fieldname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("      fieldname='%s'\n", c_name);
#endif
    *ier = cg_field_write((int)*fn, (int)*B, (int)*Z, (int)*S,
               (CGNS_ENUMT(DataType_t))*type, c_name, field_ptr, &i_F);
    *F = i_F;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_field_partial_write_f, CG_FIELD_PARTIAL_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *S, cgsize_t *type, STR_PSTR(fieldname),
	cgsize_t *rmin, cgsize_t *rmax, void *field_ptr, cgsize_t *F,
	cgsize_t *ier STR_PLEN(fieldname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_F;

    string_2_C_string(STR_PTR(fieldname), STR_LEN(fieldname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("      fieldname='%s'\n", c_name);
#endif
    *ier = cg_field_partial_write((int)*fn, (int)*B, (int)*Z, (int)*S,
               (CGNS_ENUMT(DataType_t))*type, c_name,
               rmin, rmax, field_ptr, &i_F);
    *F = i_F;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write ZoneSubRegion_t Nodes  			         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nsubregs_f, CG_NSUBREGS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *nsubreg, cgsize_t *ier)
{
    int i_nsub;

    *ier = cg_nsubregs((int)*fn, (int)*B, (int)*Z, &i_nsub);
    *nsubreg = i_nsub;
}

CGNSDLL void FMNAME(cg_subreg_info_f, CG_SUBREG_INFO_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *S, STR_PSTR(regname),
	cgsize_t *dimension, cgsize_t *location, cgsize_t *ptset_type,
	cgsize_t *npnts, cgsize_t *bcname_len, cgsize_t *gcname_len,
	cgsize_t *ier STR_PLEN(regname))
{
    char c_regname[CGIO_MAX_NAME_LENGTH+1];
    int i_dimension, i_bcname_len, i_gcname_len;
    CGNS_ENUMT(GridLocation_t) i_location;
    CGNS_ENUMT(PointSetType_t) i_ptset_type;

    *ier = cg_subreg_info((int)*fn, (int)*B, (int)*Z, (int)*S, c_regname,
               &i_dimension, &i_location, &i_ptset_type, npnts,
               &i_bcname_len, &i_gcname_len);
    if (*ier) return;
    string_2_F_string(c_regname, STR_PTR(regname), STR_LEN(regname), ier);
    *dimension = i_dimension;
    *location = i_location;
    *ptset_type = i_ptset_type;
    *bcname_len = i_bcname_len;
    *gcname_len = i_gcname_len;
}

CGNSDLL void FMNAME(cg_subreg_ptset_read_f, CG_SUBREG_PTSET_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_subreg_ptset_read((int)*fn, (int)*B, (int)*Z, (int)*S, pnts);
}

CGNSDLL void FMNAME(cg_subreg_bcname_read_f, CG_SUBREG_BCNAME_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	STR_PSTR(bcname), cgsize_t *ier STR_PLEN(bcname))
{
    char *name = 0;
    char regname[CGIO_MAX_NAME_LENGTH+1];
    int dimension, bclen, gclen;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptset_type;
    cgsize_t npnts;

    *ier = cg_subreg_info((int)*fn, (int)*B, (int)*Z, (int)*S, regname,
               &dimension, &location, &ptset_type, &npnts, &bclen, &gclen);
    if (*ier) return;
    if (bclen) name = CGNS_NEW(char, bclen+1);
    *ier = cg_subreg_bcname_read((int)*fn, (int)*B, (int)*Z, (int)*S, name);
    if (!*ier && name)
        string_2_F_string(name, STR_PTR(bcname), STR_LEN(bcname), ier);
    if (name) free(name);
}

CGNSDLL void FMNAME(cg_subreg_gcname_read_f, CG_SUBREG_GCNAME_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *S,
	STR_PSTR(gcname), cgsize_t *ier STR_PLEN(gcname))
{
    char *name = 0;
    char regname[CGIO_MAX_NAME_LENGTH+1];
    int dimension, bclen, gclen;
    CGNS_ENUMT(GridLocation_t) location;
    CGNS_ENUMT(PointSetType_t) ptset_type;
    cgsize_t npnts;

    *ier = cg_subreg_info((int)*fn, (int)*B, (int)*Z, (int)*S, regname,
               &dimension, &location, &ptset_type, &npnts, &bclen, &gclen);
    if (*ier) return;
    if (gclen) name = CGNS_NEW(char, gclen+1);
    *ier = cg_subreg_gcname_read((int)*fn, (int)*B, (int)*Z, (int)*S, name);
    if (!*ier && name)
        string_2_F_string(name, STR_PTR(gcname), STR_LEN(gcname), ier);
    if (name) free(name);
}

CGNSDLL void FMNAME(cg_subreg_ptset_write_f, CG_SUBREG_PTSET_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(regname),
	cgsize_t *dimension, cgsize_t *location, cgsize_t *ptset_type,
	cgsize_t *npnts, cgsize_t *pnts, cgsize_t *S,
	cgsize_t *ier STR_PLEN(regname))
{
    char c_regname[CGIO_MAX_NAME_LENGTH+1];
    int i_S;

    string_2_C_string(STR_PTR(regname), STR_LEN(regname),
        c_regname, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;

    *ier = cg_subreg_ptset_write((int)*fn, (int)*B, (int)*Z, c_regname,
	       (int)*dimension, (CGNS_ENUMT(GridLocation_t))*location,
	       (CGNS_ENUMT(PointSetType_t))*ptset_type, *npnts,
	       pnts, &i_S);
    *S = i_S;
}

CGNSDLL void FMNAME(cg_subreg_bcname_write_f, CG_SUBREG_BCNAME_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(regname),
	cgsize_t *dimension, STR_PSTR(bcname), cgsize_t *S,
	cgsize_t *ier STR_PLEN(regname) STR_PLEN(bcname))
{
    char c_regname[CGIO_MAX_NAME_LENGTH+1];
    char *name;
    int length, i_S;

    string_2_C_string(STR_PTR(regname), STR_LEN(regname),
        c_regname, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    length = (int)STR_LEN(bcname);
    name = CGNS_NEW(char, length+1);
    string_2_C_string(STR_PTR(bcname), STR_LEN(bcname), name, length, ier);
    if (!*ier) {
        *ier = cg_subreg_bcname_write((int)*fn, (int)*B, (int)*Z, c_regname,
	           (int)*dimension, name, &i_S);
	*S = i_S;
    }
    free(name);
}

CGNSDLL void FMNAME(cg_subreg_gcname_write_f, CG_SUBREG_GCNAME_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(regname),
	cgsize_t *dimension, STR_PSTR(gcname), cgsize_t *S,
	cgsize_t *ier STR_PLEN(regname) STR_PLEN(gcname))
{
    char c_regname[CGIO_MAX_NAME_LENGTH+1];
    char *name;
    int length, i_S;

    string_2_C_string(STR_PTR(regname), STR_LEN(regname),
        c_regname, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    length = (int)STR_LEN(gcname);
    name = CGNS_NEW(char, length+1);
    string_2_C_string(STR_PTR(gcname), STR_LEN(gcname), name, length, ier);
    if (!*ier) {
        *ier = cg_subreg_gcname_write((int)*fn, (int)*B, (int)*Z, c_regname,
	           (int)*dimension, name, &i_S);
	*S = i_S;
    }
    free(name);
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write ZoneGridConnectivity_t Nodes  			 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nzconns_f, CG_NZCONNS_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *nzconns, cgsize_t *ier)
{
    int i_nzconns;

    *ier = cg_nzconns((int)*fn, (int)*B, (int)*Z, &i_nzconns);
    *nzconns = i_nzconns;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zconn_read_f, CG_ZCONN_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *C, STR_PSTR(name),
	cgsize_t *ier STR_PLEN(name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_zconn_read((int)*fn, (int)*B, (int)*Z, (int)*C, c_name);
    if (!*ier)
        string_2_F_string(c_name, STR_PTR(name), STR_LEN(name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zconn_write_f, CG_ZCONN_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(name), cgsize_t *C,
	cgsize_t *ier STR_PLEN(name))
{
    int i_C;
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(name), STR_LEN(name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_zconn_write((int)*fn, (int)*B, (int)*Z, c_name, &i_C);
    *C = i_C;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zconn_get_f, CG_ZCONN_GET_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *C, cgsize_t *ier)
{
    int i_C;

    *ier = cg_zconn_get((int)*fn, (int)*B, (int)*Z, &i_C);
    *C = i_C;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_zconn_set_f, CG_ZCONN_SET_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *C, cgsize_t *ier)
{
    *ier = cg_zconn_set((int)*fn, (int)*B, (int)*Z, (int)*C);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write OversetHoles_t Nodes                              *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nholes_f, CG_NHOLES_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *nholes, cgsize_t *ier)
{
    int i_nholes;

    *ier = cg_nholes((int)*fn, (int)*B, (int)*Z, &i_nholes);
    *nholes = i_nholes;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_hole_info_f, CG_HOLE_INFO_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, STR_PSTR(holename), cgsize_t *location,
	cgsize_t *ptset_type, cgsize_t *nptsets, cgsize_t *npnts,
	cgsize_t *ier STR_PLEN(holename))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(GridLocation_t) i_location;
    CGNS_ENUMT(PointSetType_t) i_ptset_type;
    int i_nptsets;

    *ier = cg_hole_info((int)*fn, (int)*B, (int)*Z, (int)*I, c_name,
               &i_location, &i_ptset_type, &i_nptsets, npnts);
    if (*ier) return;
    *location = i_location;
    *ptset_type = i_ptset_type;
    *nptsets = i_nptsets;
    string_2_F_string(c_name, STR_PTR(holename), STR_LEN(holename), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_hole_read_f, CG_HOLE_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_hole_read((int)*fn, (int)*B, (int)*Z, (int)*I, pnts);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_hole_id_f, CG_HOLE_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, double *hole_id, cgsize_t *ier)
{
    *ier = cg_hole_id((int)*fn, (int)*B, (int)*Z, (int)*I, hole_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_hole_write_f, CG_HOLE_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(holename), cgsize_t *location,
	cgsize_t *ptset_type, cgsize_t *nptsets, cgsize_t *npnts,
	cgsize_t *pnts, cgsize_t *I, cgsize_t *ier STR_PLEN(holename))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_I;

    string_2_C_string(STR_PTR(holename), STR_LEN(holename),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("holename='%s'\n", c_name);
#endif
    *ier = cg_hole_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(GridLocation_t))*location,
               (CGNS_ENUMT(PointSetType_t))*ptset_type,
               (int)*nptsets, *npnts, pnts, &i_I);
    *I = i_I;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GridConnectivity_t Nodes                          *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nconns_f, CG_NCONNS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *nconns, cgsize_t *ier)
{
    int i_nconns;

    *ier = cg_nconns((int)*fn, (int)*B, (int)*Z, &i_nconns);
    *nconns = i_nconns;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_info_f, CG_CONN_INFO_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, STR_PSTR(connectname), cgsize_t *location,
	cgsize_t *type, cgsize_t *ptset_type, cgsize_t *npnts, STR_PSTR(donorname),
	cgsize_t *donor_zonetype, cgsize_t *donor_ptset_type, cgsize_t *donor_datatype,
	cgsize_t *ndata_donor, cgsize_t *ier STR_PLEN(connectname) STR_PLEN(donorname)) 
{
    char cc_name[CGIO_MAX_NAME_LENGTH+1], dc_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(GridLocation_t) i_location;
    CGNS_ENUMT(GridConnectivityType_t) i_type;
    CGNS_ENUMT(PointSetType_t) i_ptset_type;
    CGNS_ENUMT(ZoneType_t) i_donor_zonetype;
    CGNS_ENUMT(PointSetType_t) i_donor_ptset_type;
    CGNS_ENUMT(DataType_t) i_donor_datatype;

    *ier = cg_conn_info((int)*fn, (int)*B, (int)*Z, (int)*I, cc_name, &i_location,
               &i_type, &i_ptset_type, npnts, dc_name,
               &i_donor_zonetype, &i_donor_ptset_type,
               &i_donor_datatype, ndata_donor);
    if (*ier) return;
    string_2_F_string(cc_name, STR_PTR(connectname), STR_LEN(connectname), ier);
    if (*ier) return;
    string_2_F_string(dc_name, STR_PTR(donorname), STR_LEN(donorname), ier);
    if (*ier) return;
    *location = i_location;
    *type = i_type;
    *ptset_type = i_ptset_type;
    *donor_zonetype = i_donor_zonetype;
    *donor_ptset_type = i_donor_ptset_type;
    *donor_datatype = i_donor_datatype;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_read_f, CG_CONN_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, cgsize_t *pnts, cgsize_t *donor_datatype,
	cgsize_t *donor_data, cgsize_t *ier)
{
    *ier = cg_conn_read((int)*fn, (int)*B, (int)*Z, (int)*I, pnts,
               (CGNS_ENUMT(DataType_t))*donor_datatype, donor_data);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_read_short_f, CG_CONN_READ_SHORT_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *I, cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_conn_read_short((int)*fn, (int)*B, (int)*Z, (int)*I, pnts);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_id_f, CG_CONN_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, double *conn_id, cgsize_t *ier)
{
    *ier = cg_conn_id((int)*fn, (int)*B, (int)*Z, (int)*I, conn_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_write_f, CG_CONN_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(connectname), cgsize_t *location, cgsize_t *type,
	cgsize_t *ptset_type, cgsize_t *npnts, cgsize_t *pnts,
	STR_PSTR(donorname), cgsize_t *donor_zonetype, cgsize_t *donor_ptset_type,
	cgsize_t *donor_datatype, cgsize_t *ndata_donor, cgsize_t *donor_data,
	cgsize_t *I, cgsize_t *ier STR_PLEN(connectname) STR_PLEN(donorname))
{
    char cc_name[CGIO_MAX_NAME_LENGTH+1], dc_name[CGIO_MAX_NAME_LENGTH+1];
    int i_I;

    string_2_C_string(STR_PTR(connectname), STR_LEN(connectname),
        cc_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    string_2_C_string(STR_PTR(donorname), STR_LEN(donorname),
        dc_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("connectname='%s'\n", cc_name);
    printf("donorname='%s'\n", dc_name);
#endif
    *ier = cg_conn_write((int)*fn, (int)*B, (int)*Z, cc_name,
               (CGNS_ENUMT(GridLocation_t))*location,
               (CGNS_ENUMT(GridConnectivityType_t))*type,
               (CGNS_ENUMT(PointSetType_t))*ptset_type,
               *npnts, pnts, dc_name,
               (CGNS_ENUMT(ZoneType_t))*donor_zonetype,
               (CGNS_ENUMT(PointSetType_t))*donor_ptset_type,
               (CGNS_ENUMT(DataType_t))*donor_datatype,
               *ndata_donor, donor_data, &i_I);
    *I = i_I;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_write_short_f, CG_CONN_WRITE_SHORT_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(connectname), cgsize_t *location,
	cgsize_t *type, cgsize_t *ptset_type, cgsize_t *npnts,
	cgsize_t *pnts, STR_PSTR(donorname), cgsize_t *I,
	cgsize_t *ier STR_PLEN(connectname) STR_PLEN(donorname))
{
    char cc_name[CGIO_MAX_NAME_LENGTH+1], dc_name[CGIO_MAX_NAME_LENGTH+1];
    int i_I;

    string_2_C_string(STR_PTR(connectname), STR_LEN(connectname),
        cc_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    string_2_C_string(STR_PTR(donorname), STR_LEN(donorname),
        dc_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("connectname='%s'\n", cc_name);
    printf("donorname='%s'\n", dc_name);
#endif
    *ier = cg_conn_write_short((int)*fn, (int)*B, (int)*Z, cc_name,
               (CGNS_ENUMT(GridLocation_t))*location,
               (CGNS_ENUMT(GridConnectivityType_t))*type,
               (CGNS_ENUMT(PointSetType_t))*ptset_type,
               *npnts, pnts, dc_name, &i_I);
    *I = i_I;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GridConnectivity1to1_t Nodes in a zone            *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_n1to1_f, CG_N1TO1_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *n1to1, cgsize_t *ier)
{
    int i_n1to1;

    *ier = cg_n1to1((int)*fn, (int)*B, (int)*Z, &i_n1to1);
    *n1to1 = i_n1to1;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_read_f, CG_1TO1_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, STR_PSTR(connectname), STR_PSTR(donorname),
	cgsize_t *range, cgsize_t *donor_range, cgsize_t *transform,
	cgsize_t *ier STR_PLEN(connectname) STR_PLEN(donorname))
{
    char cc_name[CGIO_MAX_NAME_LENGTH+1], dc_name[CGIO_MAX_NAME_LENGTH+1];
    int n, index_dim, i_transform[3];

    *ier = cg_index_dim((int)*fn, (int)*B, (int)*Z, &index_dim);
    if (*ier) return;
    *ier = cg_1to1_read((int)*fn, (int)*B, (int)*Z, (int)*I, cc_name, dc_name,
               range, donor_range, i_transform);
    if (*ier) return;
    string_2_F_string(cc_name, STR_PTR(connectname), STR_LEN(connectname), ier);
    if (*ier) return;
    string_2_F_string(dc_name, STR_PTR(donorname), STR_LEN(donorname), ier);
    if (*ier) return;
    for (n = 0; n < index_dim; n++)
        transform[n] = i_transform[n];
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_id_f, CG_1TO1_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *I, double *one21_id, cgsize_t *ier)
{
    *ier = cg_1to1_id((int)*fn, (int)*B, (int)*Z, (int)*I, one21_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_write_f, CG_1TO1_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(connectname), STR_PSTR(donorname), cgsize_t *range,
	cgsize_t *donor_range, cgsize_t *transform, cgsize_t *I,
	cgsize_t *ier STR_PLEN(connectname) STR_PLEN(donorname))
{
    char cc_name[CGIO_MAX_NAME_LENGTH+1], dc_name[CGIO_MAX_NAME_LENGTH+1];
    int n, index_dim, i_I, i_transform[3];

    *ier = cg_index_dim((int)*fn, (int)*B, (int)*Z, &index_dim);
    if (*ier) return;
    string_2_C_string(STR_PTR(connectname), STR_LEN(connectname),
        cc_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    string_2_C_string(STR_PTR(donorname), STR_LEN(donorname),
        dc_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("connectname='%s'\n", cc_name);
    printf("donorname='%s'\n", dc_name);
#endif
    for (n = 0; n < index_dim; n++)
        i_transform[n] = (int)transform[n];
    *ier = cg_1to1_write((int)*fn, (int)*B, (int)*Z, cc_name, dc_name, range,
               donor_range, i_transform, &i_I);
    *I = i_I;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read all GridConnectivity1to1_t Nodes of a base                  *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_n1to1_global_f, CG_N1TO1_GLOBAL_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *n1to1_global, cgsize_t *ier)
{
    int i_n1to1_global;

    *ier = cg_n1to1_global((int)*fn, (int)*B, &i_n1to1_global);
    *n1to1_global = i_n1to1_global;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_read_global_f, CG_1TO1_READ_GLOBAL_F) (cgsize_t *fn,
	cgsize_t *B, STR_PSTR(connectname), STR_PSTR(zonename), STR_PSTR(donorname),
	cgsize_t *range, cgsize_t *donor_range, cgsize_t *transform,
	cgsize_t *ier STR_PLEN(connectname) STR_PLEN(zonename) STR_PLEN(donorname))
{
    int n, i, step, len;
    int cell_dim, phys_dim;     /* number of dimension for model    */
    int Ndim;           /* indexDimension           */
    int Nglobal;            /* number of 1to1 interface in base     */
    char **c_connectname, **c_zonename, **c_donorname;
    char basename[CGIO_MAX_NAME_LENGTH+1];
    cgsize_t **c_range, **c_donor_range;
    int **c_transform;

     /* get number of dimension for model: Ndim */
    *ier = cg_base_read((int)*fn, (int)*B, basename, &cell_dim, &phys_dim);
    if (*ier) return;

     /* For structured grid: */
    Ndim = cell_dim;

     /* get number of 1to1 interface in base:  Nglobal */
    *ier = cg_n1to1_global((int)*fn, (int)*B, &Nglobal);
    if (*ier) return;
    if (Nglobal < 1) {
        cgi_error("Number of interface must equal 1 or more");
        *ier = 1;
        return;
    }
     /* allocate memory for C-arrays (ptr-to-ptr) */
    if ((c_connectname = (char **)malloc(Nglobal*sizeof(char *)))==NULL ||
        (c_zonename    = (char **)malloc(Nglobal*sizeof(char *)))==NULL ||
        (c_donorname   = (char **)malloc(Nglobal*sizeof(char *)))==NULL ||
        (c_range       = (cgsize_t **)malloc(Nglobal*sizeof(cgsize_t *)))==NULL ||
        (c_donor_range = (cgsize_t **)malloc(Nglobal*sizeof(cgsize_t *)))==NULL ||
        (c_transform   = (int **)malloc(Nglobal*sizeof(int *)))==NULL) {
        cgi_error("Error allocating memory...");
        *ier = 1;
        return;
    }
    len = CGIO_MAX_NAME_LENGTH+1;
    for (n = 0; n < Nglobal; n++) {
        if ((c_connectname[n] = (char *)malloc(len*sizeof(char)))==NULL ||
            (c_zonename[n]    = (char *)malloc(len*sizeof(char)))==NULL ||
            (c_donorname[n]   = (char *)malloc(len*sizeof(char)))==NULL ||
            (c_range[n]       = (cgsize_t *)malloc(6*sizeof(cgsize_t)))==NULL ||
            (c_donor_range[n] = (cgsize_t *)malloc(6*sizeof(cgsize_t)))==NULL ||
            (c_transform[n]   = (int *)malloc(3*sizeof(int)))==NULL) {
            cgi_error("Error allocating memory...");
            *ier = 1;
            return;
        }
    }
     /* get all 1to1 interfaces */
    *ier = cg_1to1_read_global((int)*fn, (int)*B, c_connectname, c_zonename,
               c_donorname, c_range, c_donor_range, c_transform);
     /* copy C-arrays in Fortran arrays */
    if (*ier == 0) {
        for (n = 0; n < Nglobal; n++) {
            step = n*CGIO_MAX_NAME_LENGTH;
            string_2_F_string(c_connectname[n], STR_PTR(connectname)+step,
                              CGIO_MAX_NAME_LENGTH, ier);
            if (*ier) break;
            string_2_F_string(c_zonename[n],    STR_PTR(zonename)   +step,
                              CGIO_MAX_NAME_LENGTH, ier);
            if (*ier) break;
            string_2_F_string(c_donorname[n],   STR_PTR(donorname)  +step,
                              CGIO_MAX_NAME_LENGTH, ier);
            if (*ier) break;

            for (i = 0; i < Ndim; i++) {
                step = Ndim*2*n;
                range[step+i] = c_range[n][i];
                range[step+i+Ndim] = c_range[n][i+Ndim];
                donor_range[step+i] = c_donor_range[n][i];
                donor_range[step+i+Ndim] = c_donor_range[n][i+Ndim];
                transform[Ndim*n+i] = c_transform[n][i];
            }
#if DEBUG_FTOC
            printf("c_connectname[n]='%s'\n",c_connectname[n]);
            printf("c_zonename   [n]='%s'\n",c_zonename[n]);
            printf("c_donorname  [n]='%s'\n",c_donorname[n]);
            printf("..................12345678901234567890123456789012\n");
#endif
        }
#if DEBUG_FTOC
        printf("connectname='%s'\n",connectname);
        printf("zonename   ='%s'\n",zonename);
        printf("donorname  ='%s'\n",donorname);
        printf(".............+2345678901234567890123456789012+2345678901234567890123456789012\n");
#endif
    }
    for (n = 0; n < Nglobal; n++) {
        free(c_connectname[n]);
        free(c_zonename[n]);
        free(c_donorname[n]);
        free(c_range[n]);
        free(c_donor_range[n]);
        free(c_transform[n]);
    }
    free(c_connectname);
    free(c_zonename);
    free(c_donorname);
    free(c_range);
    free(c_donor_range);
    free(c_transform);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write BC_t Nodes                                        *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_nbocos_f, CG_NBOCOS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *nbocos, cgsize_t *ier)
{
    int i_nbocos;

    *ier = cg_nbocos(*fn, *B, *Z, &i_nbocos);
    *nbocos = i_nbocos;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_info_f, CG_BOCO_INFO_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *BC, STR_PSTR(boconame), cgsize_t *bocotype,
	cgsize_t *ptset_type, cgsize_t *npnts, cgsize_t *NormalIndex,
	cgsize_t *NormalListFlag, cgsize_t *NormalDataType, cgsize_t *ndataset,
	cgsize_t *ier STR_PLEN(boconame))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(BCType_t) i_bocotype;
    CGNS_ENUMT(PointSetType_t) i_ptset_type;
    CGNS_ENUMT(DataType_t) i_NormalDataType;
    int n, i_ndataset, index_dim, i_NormalIndex[3];

    *ier = cg_index_dim((int)*fn, (int)*B, (int)*Z, &index_dim);
    if (*ier) return;
    *ier = cg_boco_info((int)*fn, (int)*B, (int)*Z, (int)*BC, c_name,
               &i_bocotype, &i_ptset_type, npnts, i_NormalIndex,
               NormalListFlag, &i_NormalDataType, &i_ndataset);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(boconame), STR_LEN(boconame), ier);
    *bocotype = i_bocotype;
    *ptset_type = i_ptset_type;
    *NormalDataType = i_NormalDataType;
    *ndataset = i_ndataset;
    for (n = 0; n < index_dim; n++)
        NormalIndex[n] = i_NormalIndex[n];
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_read_f, CG_BOCO_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *BC, cgsize_t *pnts, void *NormalList, cgsize_t *ier)
{
    *ier = cg_boco_read((int)*fn, (int)*B, (int)*Z, (int)*BC, pnts, NormalList);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_id_f, CG_BOCO_ID_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *BC, double *boco_id, cgsize_t *ier)
{
    *ier = cg_boco_id((int)*fn, (int)*B, (int)*Z, (int)*BC, boco_id);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_write_f, CG_BOCO_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(boconame), cgsize_t *bocotype,
	cgsize_t *ptset_type, cgsize_t *npnts, cgsize_t *pnts,
	cgsize_t *BC, cgsize_t *ier STR_PLEN(boconame))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_BC;

    string_2_C_string(STR_PTR(boconame), STR_LEN(boconame),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("boconame='%s'\n", c_name);
#endif
    *ier = cg_boco_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(BCType_t))*bocotype,
               (CGNS_ENUMT(PointSetType_t))*ptset_type,
               *npnts, pnts, &i_BC);
    *BC = i_BC;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_normal_write_f, CG_BOCO_NORMAL_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *BC,
	cgsize_t *NormalIndex, cgsize_t *NormalListFlag,
	cgsize_t *NormalDataType, void *NormalList, cgsize_t *ier)
{
    int n, index_dim, i_NormalIndex[3];

    *ier = cg_index_dim((int)*fn, (int)*B, (int)*Z, &index_dim);
    if (*ier) return;
    for (n = 0; n < index_dim; n++)
        i_NormalIndex[n] = (int)NormalIndex[n];
    *ier = cg_boco_normal_write((int)*fn, (int)*B, (int)*Z, (int)*BC,
               i_NormalIndex, (int)*NormalListFlag,
               (CGNS_ENUMT(DataType_t))*NormalDataType, NormalList);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_gridlocation_read_f, CG_BOCO_GRIDLOCATION_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *BC,
	cgsize_t *location, cgsize_t *ier)
{
    CGNS_ENUMT(GridLocation_t) i_location;

    *ier = cg_boco_gridlocation_read((int)*fn, (int)*B, (int)*Z,
               (int)*BC, &i_location);
    *location = i_location;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_boco_gridlocation_write_f, CG_BOCO_GRIDLOCATION_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *BC,
	cgsize_t *location, cgsize_t *ier)
{
    *ier = cg_boco_gridlocation_write((int)*fn, (int)*B, (int)*Z,
               (int)*BC, (CGNS_ENUMT(GridLocation_t))*location);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write BCProperty_t/WallFunction_t Nodes                 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_bc_wallfunction_read_f, CG_BC_WALLFUNCTION_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *BC,
	cgsize_t *WallFunctionType, cgsize_t *ier)
{
    CGNS_ENUMT(WallFunctionType_t) i_WallFunctionType;

    *ier = cg_bc_wallfunction_read((int)*fn, (int)*B, (int)*Z, (int)*BC,
               &i_WallFunctionType);
    *WallFunctionType = i_WallFunctionType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_bc_wallfunction_write_f, CG_BC_WALLFUNCTION_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *BC,
	cgsize_t *WallFunctionType, cgsize_t *ier)
{
    *ier = cg_bc_wallfunction_write((int)*fn, (int)*B, (int)*Z, (int)*BC,
               (CGNS_ENUMT(WallFunctionType_t))*WallFunctionType);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write BCProperty_t/Area_t Nodes                         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_bc_area_read_f, CG_BC_AREA_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *BC, cgsize_t *AreaType,
	float *SurfaceArea, STR_PSTR(RegionName),
	cgsize_t *ier STR_PLEN(RegionName))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(AreaType_t) i_AreaType;

    *ier = cg_bc_area_read((int)*fn, (int)*B, (int)*Z, (int)*BC, &i_AreaType,
               SurfaceArea, c_name);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(RegionName), STR_LEN(RegionName), ier);
    *AreaType = i_AreaType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_bc_area_write_f, CG_BC_AREA_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *BC, cgsize_t *AreaType,
	float *SurfaceArea, STR_PSTR(RegionName),
	cgsize_t *ier STR_PLEN(RegionName))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
#if DEBUG_FTOC
    int n;
#endif

    string_2_C_string(STR_PTR(RegionName), STR_LEN(RegionName),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("RegionName='");
    for (n=0; n<32; n++) printf("%c",*((STR_PTR(RegionName))+n));
        printf("', c_name='%s'\n", c_name);
#endif
    *ier = cg_bc_area_write((int)*fn, (int)*B, (int)*Z, (int)*BC,
               (CGNS_ENUMT(AreaType_t))*AreaType,
               *SurfaceArea, c_name);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GridConnectivityProperty_t/Periodic_t Nodes       *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_conn_periodic_read_f, CG_CONN_PERIODIC_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	float *RotationCenter, float *RotationAngle, float *Translation,
	cgsize_t *ier)
{
    *ier = cg_conn_periodic_read((int)*fn, (int)*B, (int)*Z, (int)*I,
               RotationCenter, RotationAngle, Translation);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_periodic_write_f, CG_CONN_PERIODIC_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	float *RotationCenter, float *RotationAngle, float *Translation,
	cgsize_t *ier)
{
    *ier = cg_conn_periodic_write((int)*fn, (int)*B, (int)*Z, (int)*I,
               RotationCenter, RotationAngle, Translation);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_periodic_read_f, CG_1TO1_PERIODIC_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	float *RotationCenter, float *RotationAngle, float *Translation,
	cgsize_t *ier)
{
    *ier = cg_1to1_periodic_read((int)*fn, (int)*B, (int)*Z, (int)*I,
               RotationCenter, RotationAngle, Translation);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_periodic_write_f, CG_1TO1_PERIODIC_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	float *RotationCenter, float *RotationAngle, float *Translation,
	cgsize_t *ier)
{
    *ier = cg_1to1_periodic_write((int)*fn, (int)*B, (int)*Z, (int)*I,
               RotationCenter, RotationAngle, Translation);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *   Read and write GridConnectivityProperty_t/AverageInterface_t Nodes  *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_conn_average_read_f, CG_CONN_AVERAGE_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	cgsize_t *AverageInterfaceType, cgsize_t *ier)
{
    CGNS_ENUMT(AverageInterfaceType_t) i_AverageInterfaceType;

    *ier = cg_conn_average_read((int)*fn, (int)*B, (int)*Z, (int)*I,
               &i_AverageInterfaceType);
    *AverageInterfaceType = i_AverageInterfaceType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conn_average_write_f, CG_CONN_AVERAGE_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	cgsize_t *AverageInterfaceType, cgsize_t *ier)
{
    *ier = cg_conn_average_write((int)*fn, (int)*B, (int)*Z, (int)*I,
               (CGNS_ENUMT(AverageInterfaceType_t))*AverageInterfaceType);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_average_read_f, CG_1TO1_AVERAGE_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	cgsize_t *AverageInterfaceType, cgsize_t *ier)
{
    CGNS_ENUMT(AverageInterfaceType_t) i_AverageInterfaceType;

    *ier = cg_1to1_average_read((int)*fn, (int)*B, (int)*Z, (int)*I,
               &i_AverageInterfaceType);
    *AverageInterfaceType = i_AverageInterfaceType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_1to1_average_write_f, CG_1TO1_AVERAGE_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *I,
	cgsize_t *AverageInterfaceType, cgsize_t *ier)
{
    *ier = cg_1to1_average_write((int)*fn, (int)*B, (int)*Z, (int)*I,
               (CGNS_ENUMT(AverageInterfaceType_t))*AverageInterfaceType);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write BCDataSet_t Nodes                                 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_dataset_read_f, CG_DATASET_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *BC, cgsize_t *DSet,
	STR_PSTR(Dataset_name), cgsize_t *BCType, cgsize_t *DirichletFlag,
	cgsize_t *NeumannFlag, cgsize_t *ier STR_PLEN(Dataset_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(BCType_t) i_BCType;
    int i_DirichletFlag, i_NeumannFlag;

    *ier = cg_dataset_read((int)*fn, (int)*B, (int)*Z, (int)*BC, (int)*DSet,
               c_name, &i_BCType, &i_DirichletFlag, &i_NeumannFlag);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(Dataset_name), STR_LEN(Dataset_name), ier);
    *BCType = i_BCType;
    *DirichletFlag = i_DirichletFlag;
    *NeumannFlag = i_NeumannFlag;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_dataset_write_f, CG_DATASET_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *BC, STR_PSTR(Dataset_name),
	cgsize_t *BCType, cgsize_t *Dset, cgsize_t *ier STR_PLEN(Dataset_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_Dset;

    string_2_C_string(STR_PTR(Dataset_name), STR_LEN(Dataset_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("Dataset_name='%s'\n", c_name);
#endif
    *ier = cg_dataset_write((int)*fn, (int)*B, (int)*Z, (int)*BC, c_name,
               (CGNS_ENUMT(BCType_t))*BCType, &i_Dset);
    *Dset = i_Dset;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_bcdataset_write_f, CG_BCDATASET_WRITE_F) (
	STR_PSTR(Dataset_name), cgsize_t *BCType, cgsize_t *BCDataType,
	cgsize_t *ier STR_PLEN(Dataset_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(Dataset_name), STR_LEN(Dataset_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
#if DEBUG_FTOC
    printf("Dataset_name='%s'\n", c_name);
#endif
    *ier = cg_bcdataset_write(c_name, (CGNS_ENUMT(BCType_t))*BCType,
               (CGNS_ENUMT(BCDataType_t))*BCDataType);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_bcdataset_info_f, CG_BCDATASET_INFO_F) (
	cgsize_t *ndataset, cgsize_t *ier STR_PLEN(Dataset_name))
{
    int i_ndataset;

    *ier = cg_bcdataset_info(&i_ndataset);
    *ndataset = i_ndataset;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_bcdataset_read_f, CG_BCDATASET_READ_F) (
	cgsize_t *index, STR_PSTR(Dataset_name), cgsize_t *BCType,
	cgsize_t *DirichletFlag, cgsize_t *NeumannFlag,
	cgsize_t *ier STR_PLEN(Dataset_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(BCType_t) i_BCType;
    int i_DirichletFlag, i_NeumannFlag;

    *ier = cg_bcdataset_read((int)*index, c_name, &i_BCType,
               &i_DirichletFlag, &i_NeumannFlag);
    if (*ier) return;
    *BCType = i_BCType;
    *DirichletFlag = i_DirichletFlag;
    *NeumannFlag = i_NeumannFlag;
    string_2_F_string(c_name, STR_PTR(Dataset_name), STR_LEN(Dataset_name), ier);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write BCData_t Nodes                                    *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_bcdata_write_f, CG_BCDATA_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *BC, cgsize_t *Dset,
	cgsize_t *BCDataType, cgsize_t *ier)
{
    *ier = cg_bcdata_write((int)*fn, (int)*B, (int)*Z, (int)*BC,
               (int)*Dset, (CGNS_ENUMT(BCDataType_t))*BCDataType);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write RigidGridMotion_t Nodes                           *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_n_rigid_motions_f, CG_N_RIGID_MOTIONS_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *n_rigid_motions, cgsize_t *ier)
{
    int i_n_rigid_motions;

    *ier = cg_n_rigid_motions((int)*fn, (int)*B, (int)*Z, &i_n_rigid_motions);
    *n_rigid_motions = i_n_rigid_motions;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_rigid_motion_read_f, CG_RIGID_MOTION_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *R, STR_PSTR(rmotion_name),
	cgsize_t *type, cgsize_t *ier STR_PLEN(rmotion_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(RigidGridMotionType_t) i_type;

    *ier = cg_rigid_motion_read((int)*fn, (int)*B, (int)*Z, (int)*R,
               c_name, &i_type);
    if (*ier) return;
    *type = i_type;
    string_2_F_string(c_name, STR_PTR(rmotion_name), STR_LEN(rmotion_name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_rigid_motion_write_f, CG_RIGID_MOTION_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(rmotion_name),
	cgsize_t *type, cgsize_t *R, cgsize_t *ier STR_PLEN(rmotion_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_R;

    string_2_C_string(STR_PTR(rmotion_name), STR_LEN(rmotion_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_rigid_motion_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(RigidGridMotionType_t))*type, &i_R);
    *R = i_R;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write ArbitraryGridMotion_t Nodes                       *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_n_arbitrary_motions_f, CG_N_ARBITRARY_MOTIONS_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *n_arbitrary_motions,
	cgsize_t *ier)
{
    int i_n_arbitrary_motions;

    *ier = cg_n_arbitrary_motions((int)*fn, (int)*B, (int)*Z,
               &i_n_arbitrary_motions);
    *n_arbitrary_motions = i_n_arbitrary_motions;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_arbitrary_motion_read_f, CG_ARBITRARY_MOTION_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, cgsize_t *A,
	STR_PSTR(amotion_name), cgsize_t *type,
	cgsize_t *ier STR_PLEN(amotion_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(ArbitraryGridMotionType_t) i_type;

    *ier = cg_arbitrary_motion_read((int)*fn, (int)*B, (int)*Z, (int)*A,
               c_name, &i_type);
    if (*ier) return;
    *type = i_type;
    string_2_F_string(c_name, STR_PTR(amotion_name), STR_LEN(amotion_name), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_arbitrary_motion_write_f, CG_ARBITRARY_MOTION_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *Z, STR_PSTR(amotion_name),
	cgsize_t *type, cgsize_t *A, cgsize_t *ier STR_PLEN(amotion_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_A;

    string_2_C_string(STR_PTR(amotion_name), STR_LEN(amotion_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_arbitrary_motion_write((int)*fn, (int)*B, (int)*Z, c_name,
               (CGNS_ENUMT(ArbitraryGridMotionType_t))*type, &i_A);
    *A = i_A;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write GridCoordinates_t Nodes                           *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_ngrids_f, CG_NGRIDS_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, cgsize_t *ngrids, cgsize_t *ier)
{
    int i_ngrids;

    *ier = cg_ngrids((int)*fn, (int)*B, (int)*Z, &i_ngrids);
    *ngrids = i_ngrids;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_grid_read_f, CG_GRID_READ_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, cgsize_t *G, STR_PSTR(gridname),
	cgsize_t *ier STR_PLEN(gridname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_grid_read((int)*fn, (int)*B, (int)*Z, (int)*G, c_name);
    if (!*ier)
        string_2_F_string(c_name, STR_PTR(gridname), STR_LEN(gridname), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_grid_write_f, CG_GRID_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, cgsize_t *Z, STR_PSTR(gridname), cgsize_t *G,
	cgsize_t *ier STR_PLEN(gridname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_G;

    string_2_C_string(STR_PTR(gridname), STR_LEN(gridname),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_grid_write((int)*fn, (int)*B, (int)*Z, c_name, &i_G);
    *G = i_G;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write SimulationType_t Node                             *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_simulation_type_read_f, CG_SIMULATION_TYPE_READ_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *type, cgsize_t *ier)
{
    CGNS_ENUMT(SimulationType_t) i_type;

    *ier = cg_simulation_type_read((int)*fn, (int)*B, &i_type);
    *type = i_type;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_simulation_type_write_f, CG_SIMULATION_TYPE_WRITE_F) (
	cgsize_t *fn, cgsize_t *B, cgsize_t *type, cgsize_t *ier)
{
    *ier = cg_simulation_type_write((int)*fn, (int)*B,
               (CGNS_ENUMT(SimulationType_t))*type);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write BaseIterativeData_t Node                          *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_biter_read_f, CG_BITER_READ_F) (cgsize_t *fn,
	cgsize_t *B, STR_PSTR(bitername), cgsize_t *nsteps,
	cgsize_t *ier STR_PLEN(bitername))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_nsteps;

    *ier = cg_biter_read((int)*fn, (int)*B, c_name, &i_nsteps);
    if (*ier) return;
    *nsteps = i_nsteps;
    string_2_F_string(c_name, STR_PTR(bitername), STR_LEN(bitername), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_biter_write_f, CG_BITER_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, STR_PSTR(bitername), cgsize_t *nsteps,
	cgsize_t *ier STR_PLEN(bitername))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(bitername), STR_LEN(bitername),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_biter_write((int)*fn, (int)*B, c_name, (int)*nsteps);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write ZoneIterativeData_t Nodes                         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_ziter_read_f, CG_ZITER_READ_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(zitername), cgsize_t *ier STR_PLEN(zitername))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_ziter_read((int)*fn, (int)*B, (int)*Z, c_name);
    if (*ier == 0)
        string_2_F_string(c_name, STR_PTR(zitername), STR_LEN(zitername), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_ziter_write_f, CG_ZITER_WRITE_F) (cgsize_t *fn, cgsize_t *B,
	cgsize_t *Z, STR_PSTR(zitername), cgsize_t *ier STR_PLEN(zitername))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(zitername), STR_LEN(zitername),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_ziter_write((int)*fn, (int)*B, (int)*Z, c_name);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write Gravity_t Node                                    *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_gravity_read_f, CG_GRAVITY_READ_F) (cgsize_t *fn,
	cgsize_t *B, float *gravity_vector, cgsize_t *ier)
{
    *ier = cg_gravity_read((int)*fn, (int)*B, gravity_vector);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_gravity_write_f, CG_GRAVITY_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, float *gravity_vector, cgsize_t *ier)
{
   *ier = cg_gravity_write((int)*fn, (int)*B, gravity_vector);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write Axisymmetry_t Node                                *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_axisym_read_f, CG_AXISYM_READ_F) (cgsize_t *fn,
	cgsize_t *B, float *ref_point, float *axis, cgsize_t *ier)
{
    *ier = cg_axisym_read((int)*fn, (int)*B, ref_point, axis);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_axisym_write_f, CG_AXISYM_WRITE_F) (cgsize_t *fn,
	cgsize_t *B, float *ref_point, float *axis, cgsize_t *ier)
{
    *ier = cg_axisym_write((int)*fn, (int)*B, ref_point, axis);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write RotatingCoordinates_t Node                        *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_rotating_read_f, CG_ROTATING_READ_F) (
	float *rot_rate, float *rot_center, cgsize_t *ier)
{
    *ier = cg_rotating_read(rot_rate, rot_center);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_rotating_write_f, CG_ROTATING_WRITE_F) (
	float *rot_rate, float *rot_center, cgsize_t *ier)
{
        *ier = cg_rotating_write(rot_rate, rot_center);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Read and write  IndexArray/Range_t Nodes                         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_ptset_info_f, CG_PTSET_INFO_F) (
	cgsize_t *ptset_type, cgsize_t *npnts, cgsize_t *ier)
{
    CGNS_ENUMT(PointSetType_t) i_ptset_type;

    *ier = cg_ptset_info(&i_ptset_type, npnts);
    *ptset_type = i_ptset_type;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_ptset_read_f, CG_PTSET_READ_F) (
	cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_ptset_read(pnts);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_ptset_write_f, CG_PTSET_WRITE_F) (
	cgsize_t *ptset_type, cgsize_t *npnts, cgsize_t *pnts, cgsize_t *ier)
{
    *ier = cg_ptset_write((CGNS_ENUMT(PointSetType_t))*ptset_type, *npnts, pnts);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Go - To Function                                                 *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef WIN32_FORTRAN
CGNSDLL void __stdcall cg_goto_f(cgsize_t *fn, cgsize_t *B, cgsize_t *ier, ...)
#else
CGNSDLL void FMNAME(cg_goto_f, CG_GOTO_F)(cgsize_t *fn, cgsize_t *B, cgsize_t *ier, ...)
#endif
{
#ifdef _CRAY
    _fcd cray_string;
#endif
    char *f_label[CG_MAX_GOTO_DEPTH], *label[CG_MAX_GOTO_DEPTH];
    int index[CG_MAX_GOTO_DEPTH], n, i, len[CG_MAX_GOTO_DEPTH];
    va_list ap;

     /* initialize ap to the last parameter before the variable argument list */
     /* Note:  On HP, print statements btw va_start and va_end create major problems */

    va_start(ap, ier);

     /* read arguments */
    for (n = 0; n < CG_MAX_GOTO_DEPTH; n++)  {
#ifdef _CRAY
        cray_string = va_arg(ap, _fcd);
        f_label[n] = _fcdtocp(cray_string);
        len[n] = _fcdlen(cray_string);
#else
        f_label[n] = va_arg(ap, char *);
# ifdef WIN32_FORTRAN
     /* In Windows, the arguments appear in a different order: char*, len, index,...*/
        len[n] = va_arg(ap, int);
# endif
#endif
        if (f_label[n][0] == ' ' || 0 == strncmp(f_label[n],"end",3) ||
            0 == strncmp(f_label[n],"END",3)) break;

        index[n] = (int)*(va_arg(ap, cgsize_t *));
        if (index[n] < 0) {
            cgi_error("Incorrect input to function cg_goto_f");
            *ier = 1;
            return;
        }
    }

#if !defined(_CRAY) && !defined(WIN32_FORTRAN)
    for (i=0; i<n; i++) {
        len[i] = va_arg(ap, int);
    }
#endif
    va_end(ap);

     /* convert strings to C-strings */
    for (i=0; i < n; i++) {
        label[i] = CGNS_NEW(char,len[i]+1);
        string_2_C_string(f_label[i], len[i], label[i], len[i], ier);
    }

#if DEBUG_GOTO
    printf("\nIn cg_ftoc.c: narguments=%d\n",n);
    for (i=0; i<n; i++) printf("\targ %d: '%s' #%d\n",i,label[i], index[i]);
#endif

    *ier = cgi_set_posit((int)*fn, (int)*B, n, index, label);

    for (i=0; i<n; i++) CGNS_FREE(label[i]);
    return;
}

/*-----------------------------------------------------------------------*/

#ifdef WIN32_FORTRAN
CGNSDLL void __stdcall cg_gorel_f(cgsize_t *fn, cgsize_t *ier, ...)
#else
CGNSDLL void FMNAME(cg_gorel_f, CG_GOREL_F)(cgsize_t *fn, cgsize_t *ier, ...)
#endif
{
#ifdef _CRAY
    _fcd cray_string;
#endif
    char *f_label[CG_MAX_GOTO_DEPTH], *label[CG_MAX_GOTO_DEPTH];
    int index[CG_MAX_GOTO_DEPTH], n, i, len[CG_MAX_GOTO_DEPTH];
    va_list ap;

    if (posit == 0) {
        cgi_error ("position not set with cg_goto");
        *ier = CG_ERROR;
        return;
    }
    if ((int)*fn != posit_file) {
        cgi_error("current position is in the wrong file");
        *ier = CG_ERROR;
        return;
    }

     /* initialize ap to the last parameter before the variable argument list */
     /* Note:  On HP, print statements btw va_start and va_end create major problems */

    va_start(ap, ier);

     /* read arguments */
    for (n = 0; n < CG_MAX_GOTO_DEPTH; n++)  {
#ifdef _CRAY
        cray_string = va_arg(ap, _fcd);
        f_label[n] = _fcdtocp(cray_string);
        len[n] = _fcdlen(cray_string);
#else
        f_label[n] = va_arg(ap, char *);
# ifdef WIN32_FORTRAN
     /* In Windows, the arguments appear in a different order: char*, len, index,...*/
        len[n] = va_arg(ap, int);
# endif
#endif
        if (f_label[n][0] == ' ' || 0 == strncmp(f_label[n],"end",3) ||
            0 == strncmp(f_label[n],"END",3)) break;

        index[n] = (int)*(va_arg(ap, cgsize_t *));
        if (index[n] < 0) {
            cgi_error("Incorrect input to function cg_goto_f");
            *ier = 1;
            return;
        }
    }

#if !defined(_CRAY) && !defined(WIN32_FORTRAN)
    for (i=0; i<n; i++) {
        len[i] = va_arg(ap, int);
    }
#endif
    va_end(ap);

     /* convert strings to C-strings */
    for (i=0; i < n; i++) {
        label[i] = CGNS_NEW(char,len[i]+1);
        string_2_C_string(f_label[i], len[i], label[i], len[i], ier);
    }

#if DEBUG_GOTO
    printf("\nIn cg_ftoc.c: narguments=%d\n",n);
    for (i=0; i<n; i++) printf("\targ %d: '%s' #%d\n",i,label[i], index[i]);
#endif

    *ier = cgi_update_posit(n, index, label);

    for (i=0; i<n; i++) CGNS_FREE(label[i]);
    return;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_gopath_f, CG_GOPATH_F) (cgsize_t *fn,
	STR_PSTR(path), cgsize_t *ier STR_PLEN(path))
{
    int length;
    char *c_path;

    length = (int) STR_LEN(path);
    c_path = CGNS_NEW(char, length+1);

    string_2_C_string(STR_PTR(path), STR_LEN(path), c_path, length, ier);
    if (*ier == 0)
        *ier = cg_gopath((int)*fn, c_path);
    CGNS_FREE(c_path);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *              Read Multiple path nodes                         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_famname_read_f, CG_FAMNAME_READ_F) (
	STR_PSTR(famname), cgsize_t *ier STR_PLEN(famname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_famname_read(c_name);
    if (*ier == 0)
        string_2_F_string(c_name, STR_PTR(famname), STR_LEN(famname), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_convergence_read_f, CG_CONVERGENCE_READ_F) (
	cgsize_t *iterations, STR_PSTR(NormDefinitions),
	cgsize_t *ier STR_PLEN(NormDefinitions))
{
    char *c_descr_text;
    int i_iterations;

    *ier = cg_convergence_read(&i_iterations, &c_descr_text);
    if (*ier) return;
    string_2_F_string(c_descr_text, STR_PTR(NormDefinitions),
        STR_LEN(NormDefinitions), ier);
    *iterations = i_iterations;
    free(c_descr_text);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_state_size_f, CG_STATE_SIZE_F) (
	cgsize_t *size, cgsize_t *ier)
{
    char *c_descr_text;

    *ier = cg_state_read(&c_descr_text);
    if (*ier) return;
    *size = (cgsize_t)strlen(c_descr_text);
    free(c_descr_text);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_state_read_f, CG_STATE_READ_F) (
	STR_PSTR(StateDescription), cgsize_t *ier STR_PLEN(StateDescription))
{
    char *c_descr_text;

    *ier = cg_state_read(&c_descr_text);
    if (*ier) return;
    string_2_F_string(c_descr_text, STR_PTR(StateDescription),
        STR_LEN(StateDescription), ier);
    free(c_descr_text);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_equationset_read_f, CG_EQUATIONSET_READ_F) (
	cgsize_t *EquationDimension, cgsize_t *GoverningEquationsFlag,
	cgsize_t *GasModelFlag, cgsize_t *ViscosityModelFlag,
	cgsize_t *ThermalConductivityModelFlag,
	cgsize_t *TurbulenceClosureFlag, cgsize_t *TurbulenceModelFlag,
	cgsize_t *ier)
{
    int i_EquationDimension, i_GoverningEquationsFlag, i_GasModelFlag;
    int i_ViscosityModelFlag, i_ThermalConductivityModelFlag;
    int i_TurbulenceClosureFlag, i_TurbulenceModelFlag;

    *ier = cg_equationset_read(&i_EquationDimension,
               &i_GoverningEquationsFlag,
               &i_GasModelFlag, &i_ViscosityModelFlag,
               &i_ThermalConductivityModelFlag,
               &i_TurbulenceClosureFlag, &i_TurbulenceModelFlag);
#if DEBUG_FTOC
    printf("in cg_ftoc, EquationDimension=%d\n",*EquationDimension);
#endif
    *EquationDimension = i_EquationDimension;
    *GoverningEquationsFlag = i_GoverningEquationsFlag;
    *GasModelFlag = i_GasModelFlag;
    *ViscosityModelFlag = i_ViscosityModelFlag;
    *ThermalConductivityModelFlag = i_ThermalConductivityModelFlag;
    *TurbulenceClosureFlag = i_TurbulenceClosureFlag;
    *TurbulenceModelFlag = i_TurbulenceModelFlag;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_equationset_chemistry_read_f, CG_EQUATIONSET_CHEMISTRY_READ_F) (
	cgsize_t *ThermalRelaxationFlag, cgsize_t *ChemicalKineticsFlag, cgsize_t *ier)
{
    int i_ThermalRelaxationFlag, i_ChemicalKineticsFlag;

    *ier = cg_equationset_chemistry_read(&i_ThermalRelaxationFlag,
               &i_ChemicalKineticsFlag);
    *ThermalRelaxationFlag = i_ThermalRelaxationFlag;
    *ChemicalKineticsFlag = i_ChemicalKineticsFlag;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_equationset_elecmagn_read_f, CG_EQUATIONSET_ELECMAGN_READ_F) (
	cgsize_t *ElecFldModelFlag, cgsize_t *MagnFldModelFlag,
	cgsize_t *ConductivityModelFlag, cgsize_t *ier)
{
    int i_ElecFldModelFlag, i_MagnFldModelFlag, i_ConductivityModelFlag;

    *ier = cg_equationset_elecmagn_read(&i_ElecFldModelFlag,
               &i_MagnFldModelFlag, &i_ConductivityModelFlag);
    *ElecFldModelFlag = i_ElecFldModelFlag;
    *MagnFldModelFlag = i_MagnFldModelFlag;
    *ConductivityModelFlag = i_ConductivityModelFlag;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_governing_read_f, CG_GOVERNING_READ_F) (
	cgsize_t *EquationsType, cgsize_t *ier)
{
    CGNS_ENUMT(GoverningEquationsType_t) i_EquationsType;

    *ier = cg_governing_read(&i_EquationsType);
    *EquationsType = i_EquationsType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_diffusion_read_f, CG_DIFFUSION_READ_F) (
	cgsize_t *diffusion_model, cgsize_t *ier)
{
    int n, index_dim, ndata, i_diffusion_model[6];

    index_dim = cgi_posit_index_dim();
    if (index_dim == 1) ndata = 1;
    else if (index_dim == 2) ndata = 3;
    else if (index_dim == 3) ndata = 6;
    else {
        cgi_error("invalid value for IndexDimension");
        *ier = CG_ERROR;
        return;
    }
    *ier = cg_diffusion_read(i_diffusion_model);
    if (*ier) return;
    for (n = 0; n < ndata; n++)
        diffusion_model[n] = i_diffusion_model[n];
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_model_read_f, CG_MODEL_READ_F) (STR_PSTR(ModelLabel),
	cgsize_t *ModelType, cgsize_t *ier STR_PLEN(ModelLabel))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(ModelType_t) i_ModelType;

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(ModelLabel), STR_LEN(ModelLabel),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_model_read(c_name, &i_ModelType);
    *ModelType = i_ModelType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_narrays_f, CG_NARRAYS_F) (cgsize_t *narrays, cgsize_t *ier)
{
    int i_narrays;

    *ier = cg_narrays(&i_narrays);
    *narrays = i_narrays;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_array_info_f, CG_ARRAY_INFO_F) (cgsize_t *A,
	STR_PSTR(ArrayName), cgsize_t *DataType, cgsize_t *DataDimension,
	cgsize_t *DimensionVector, cgsize_t *ier STR_PLEN(ArrayName))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
    int i_DataDimension;
    CGNS_ENUMT(DataType_t) i_DataType;

    *ier = cg_array_info((int)*A, c_name, &i_DataType, &i_DataDimension,
                         DimensionVector);
    if (*ier) return;
    string_2_F_string(c_name, STR_PTR(ArrayName), STR_LEN(ArrayName), ier);
    *DataType = i_DataType;
    *DataDimension = i_DataDimension;
}

/*-----------------------------------------------------------------------*/

#ifdef WIN32_FORTRAN
CGNSDLL void __stdcall cg_array_read_f(cgsize_t *A, void *Data, ...)
{
    va_list ap;
    cgsize_t *ier;
    int DataDimension;
    cgsize_t DimensionVector[CGIO_MAX_DIMENSIONS];
    char ArrayName[CGIO_MAX_NAME_LENGTH+1];
    CGNS_ENUMT(DataType_t) DataType;

    cg_array_info((int)*A, ArrayName, &DataType, &DataDimension, DimensionVector);

    va_start(ap, Data);
    if (DataType == CGNS_ENUMV(Character)) (void) va_arg(ap, int);
    ier = va_arg(ap, cgsize_t *);
    va_end(ap);
#else
CGNSDLL void FMNAME(cg_array_read_f, CG_ARRAY_READ_F) (cgsize_t *A,
	void *Data, cgsize_t *ier)
{
#endif
    *ier = cg_array_read((int)*A, Data);
}

/*-----------------------------------------------------------------------*/

#ifdef WIN32_FORTRAN
CGNSDLL void __stdcall cg_array_read_as_f(cgsize_t *A, cgsize_t *type,
	void *Data, ...)
{
    va_list ap;
    cgsize_t *ier;
    va_start(ap, Data);
    if ((CGNS_ENUMT(DataType_t))*type == CGNS_ENUMV(Character))
        (void) va_arg(ap, int);
    ier = va_arg(ap, cgsize_t *);
    va_end(ap);
#else
CGNSDLL void FMNAME(cg_array_read_as_f, CG_ARRAY_READ_AS_F) (cgsize_t *A,
	cgsize_t *type, void *Data, cgsize_t *ier)
{
#endif
    *ier = cg_array_read_as((int)*A, (CGNS_ENUMT(DataType_t))*type, Data);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_nintegrals_f, CG_NINTEGRALS_F) (
	cgsize_t *nintegrals, cgsize_t *ier)
{
    int i_nintegrals;

    *ier = cg_nintegrals(&i_nintegrals);
    *nintegrals = i_nintegrals;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_integral_read_f, CG_INTEGRAL_READ_F) (
	cgsize_t *IntegralDataIndex, STR_PSTR(IntegralDataName),
	cgsize_t *ier STR_PLEN(IntegralDataName))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_integral_read((int)*IntegralDataIndex, c_name);
    if (!*ier)
        string_2_F_string(c_name, STR_PTR(IntegralDataName),
            STR_LEN(IntegralDataName), ier);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_rind_read_f, CG_RIND_READ_F) (
	cgsize_t *RindData, cgsize_t *ier)
{
    int n, index_dim, i_RindData[6];

    index_dim = cgi_posit_index_dim();
    *ier = cg_rind_read(i_RindData);
    if (*ier) return;
    for (n = 0; n < 2*index_dim; n++)
        RindData[n] = i_RindData[n];
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_ndescriptors_f, CG_NDESCRIPTORS_F) (
	cgsize_t *ndescriptors, cgsize_t *ier)
{
    int i_ndescriptors;

    *ier = cg_ndescriptors(&i_ndescriptors);
    *ndescriptors = i_ndescriptors;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_descriptor_size_f, CG_DESCRIPTOR_SIZE_F) (
	cgsize_t *descr_no, cgsize_t *descr_size, cgsize_t *ier)
{
    char *c_descr_text;
    char descr_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_descriptor_read((int)*descr_no, descr_name, &c_descr_text);
    if (!*ier) {
        *descr_size = (cgsize_t)strlen(c_descr_text);
        free(c_descr_text);
    }
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_descriptor_read_f, CG_DESCRIPTOR_READ_F) (
	cgsize_t *descr_no, STR_PSTR(descr_name), STR_PSTR(descr_text),
	cgsize_t *ier STR_PLEN(descr_name)  STR_PLEN(descr_text))
{
    char *c_descr_text;
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_descriptor_read((int)*descr_no, c_name, &c_descr_text);
    if (*ier) return;
#if DEBUG_FTOC
    printf("In cg_descriptor_read_f, descr_no=%d, descr_name='%s', c_descr_text='%s'\n",
        *descr_no, c_name, c_descr_text);
#endif
    string_2_F_string(c_name, STR_PTR(descr_name), STR_LEN(descr_name), ier);
    if (!*ier)
        string_2_F_string(c_descr_text, STR_PTR(descr_text),
            STR_LEN(descr_text), ier);
    free(c_descr_text);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_nunits_f, CG_NUNITS_F) (cgsize_t *nunits, cgsize_t *ier)
{
    int i_nunits;

    *ier = cg_nunits(&i_nunits);
    *nunits = i_nunits;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_units_read_f, CG_UNITS_READ_F) (
	cgsize_t *mass, cgsize_t *length, cgsize_t *time,
	cgsize_t *temperature, cgsize_t *angle, cgsize_t *ier)
{
    CGNS_ENUMT(MassUnits_t) i_mass;
    CGNS_ENUMT(LengthUnits_t) i_length;
    CGNS_ENUMT(TimeUnits_t) i_time;
    CGNS_ENUMT(TemperatureUnits_t) i_temperature;
    CGNS_ENUMT(AngleUnits_t) i_angle;

    *ier = cg_units_read(&i_mass, &i_length, &i_time, &i_temperature, &i_angle);
    *mass = i_mass;
    *length = i_length;
    *time = i_time;
    *temperature = i_temperature;
    *angle = i_angle;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_unitsfull_read_f, CG_UNITSFULL_READ_F) (
	cgsize_t *mass, cgsize_t *length, cgsize_t *time,
	cgsize_t *temperature, cgsize_t *angle, cgsize_t *current,
	cgsize_t *amount, cgsize_t *intensity, cgsize_t *ier)
{
    CGNS_ENUMT(MassUnits_t) i_mass;
    CGNS_ENUMT(LengthUnits_t) i_length;
    CGNS_ENUMT(TimeUnits_t) i_time;
    CGNS_ENUMT(TemperatureUnits_t) i_temperature;
    CGNS_ENUMT(AngleUnits_t) i_angle;
    CGNS_ENUMT(ElectricCurrentUnits_t) i_current;
    CGNS_ENUMT(SubstanceAmountUnits_t) i_amount;
    CGNS_ENUMT(LuminousIntensityUnits_t) i_intensity;

    *ier = cg_unitsfull_read(&i_mass, &i_length, &i_time, &i_temperature,
                             &i_angle, &i_current, &i_amount, &i_intensity);
    *mass = i_mass;
    *length = i_length;
    *time = i_time;
    *temperature = i_temperature;
    *angle = i_angle;
    *current = i_current;
    *amount = i_amount;
    *intensity = i_intensity;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_exponents_info_f, CG_EXPONENTS_INFO_F) (
	cgsize_t *DataType, cgsize_t *ier)
{
    CGNS_ENUMT(DataType_t) i_DataType;

    *ier = cg_exponents_info(&i_DataType);
    *DataType = i_DataType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_nexponents_f, CG_NEXPONENTS_F) (
	cgsize_t *nexps, cgsize_t *ier)
{
    int i_nexps;

    *ier = cg_nexponents(&i_nexps);
    *nexps = i_nexps;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_exponents_read_f, CG_EXPONENTS_READ_F) (
	void *exponents, cgsize_t *ier)
{
    *ier = cg_exponents_read(exponents);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_expfull_read_f, CG_EXPFULL_READ_F) (
	void *exponents, cgsize_t *ier)
{
    *ier = cg_expfull_read(exponents);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conversion_info_f, CG_CONVERSION_INFO_F) (
	cgsize_t *DataType, cgsize_t *ier)
{
    CGNS_ENUMT(DataType_t) i_DataType;

    *ier = cg_conversion_info(&i_DataType);
    *DataType = i_DataType;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conversion_read_f, CG_CONVERSION_READ_F) (
	void *ConversionFactors, cgsize_t *ier)
{
    *ier = cg_conversion_read(ConversionFactors);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_dataclass_read_f, CG_DATACLASS_READ_F) (
	cgsize_t *dataclass, cgsize_t *ier)
{
    CGNS_ENUMT(DataClass_t) i_dataclass;

    *ier = cg_dataclass_read(&i_dataclass);
    *dataclass = i_dataclass;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_gridlocation_read_f, CG_GRIDLOCATION_READ_F) (
	cgsize_t *GridLocation, cgsize_t *ier)
{
    CGNS_ENUMT(GridLocation_t) i_GridLocation;

    *ier = cg_gridlocation_read(&i_GridLocation);
    *GridLocation = i_GridLocation;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_ordinal_read_f, CG_ORDINAL_READ_F) (
	cgsize_t *Ordinal, cgsize_t *ier)
{
    int i_Ordinal;

    *ier = cg_ordinal_read(&i_Ordinal);
    *Ordinal = i_Ordinal;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_npe_f, CG_NPE_F) (cgsize_t *type,
	cgsize_t *npe, cgsize_t *ier)
{
    int i_npe;

    *ier = cg_npe((CGNS_ENUMT(ElementType_t))*type, &i_npe);
    *npe = i_npe;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_is_link_f, CG_IS_LINK_F) (
	cgsize_t *path_length, cgsize_t *ier)
{
    int i_path_length;

    *ier = cg_is_link(&i_path_length);
    *path_length = i_path_length;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_link_read_f, CG_LINK_READ_F) (
	STR_PSTR(filename), STR_PSTR(link_path), cgsize_t *ier
	STR_PLEN(filename)  STR_PLEN(link_path))
{
    char *f_name, *l_name;

    *ier = cg_link_read(&f_name, &l_name);
    if (*ier) return;
    string_2_F_string(f_name, STR_PTR(filename), STR_LEN(filename), ier);
    if (*ier == 0)
        string_2_F_string(l_name, STR_PTR(link_path), STR_LEN(link_path), ier);
    free(f_name);
    free(l_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_nuser_data_f, CG_NUSER_DATA_F) (
	cgsize_t *nuser_data, cgsize_t *ier)
{
    int i_nuser_data;

    *ier = cg_nuser_data(&i_nuser_data);
    *nuser_data = i_nuser_data;
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_user_data_read_f, CG_USER_DATA_READ_F) (cgsize_t *index,
	STR_PSTR(dataname), cgsize_t *ier STR_PLEN(dataname))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    *ier = cg_user_data_read((int)*index, c_name);
    if (*ier == 0)
        string_2_F_string(c_name, STR_PTR(dataname), STR_LEN(dataname), ier);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *                   Write Multiple path nodes                           *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_famname_write_f, CG_FAMNAME_WRITE_F) (
	STR_PSTR(family_name), cgsize_t *ier STR_PLEN(family_name))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(family_name), STR_LEN(family_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_famname_write(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_convergence_write_f, CG_CONVERGENCE_WRITE_F) (
	cgsize_t *iterations, STR_PSTR(NormDefinitions),
	cgsize_t *ier STR_PLEN(NormDefinitions))
{
    char *c_string;
    int len;

    len = STR_LEN(NormDefinitions);
     /* convert Fortran-text-string to a C-string */
    c_string = CGNS_NEW(char, len+1);
    string_2_C_string(STR_PTR(NormDefinitions), len, c_string, len, ier);
    if (*ier == 0) {
#if DEBUG_FTOC
        printf("In cg_ftoc: c_NormDefinitions = '%s'",c_string);
#endif
        *ier = cg_convergence_write((int)*iterations, c_string);
    }
    free(c_string);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_state_write_f, CG_STATE_WRITE_F) (STR_PSTR(StateDescription),
	cgsize_t *ier STR_PLEN(StateDescription))
{
    char *c_string;
    int len;

    len = STR_LEN(StateDescription);
     /* convert Fortran-text-string to a C-string */
    c_string = CGNS_NEW(char, len+1);
    string_2_C_string(STR_PTR(StateDescription), len, c_string, len, ier);
    if (*ier == 0) {
#if DEBUG_FTOC
        printf("In cg_ftoc: C_StateDescription = '%s'",c_string);
#endif
        *ier = cg_state_write(c_string);
    }
    free(c_string);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_equationset_write_f, CG_EQUATIONSET_WRITE_F) (
	cgsize_t *EquationDimension, cgsize_t *ier)
{
#if DEBUG_FTOC
    printf("In cg_ftoc: EquationDimension=%d\n",*EquationDimension);
#endif
    *ier = cg_equationset_write((int)*EquationDimension);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_governing_write_f, CG_GOVERNING_WRITE_F) (
	cgsize_t *Equationstype, cgsize_t *ier)
{
    *ier = cg_governing_write((CGNS_ENUMT(GoverningEquationsType_t))*Equationstype);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_diffusion_write_f, CG_DIFFUSION_WRITE_F) (
	cgsize_t *diffusion_model, cgsize_t *ier)
{
    int n, index_dim, ndata, i_diffusion_model[6];

    index_dim = cgi_posit_index_dim();
    if (index_dim == 1) ndata = 1;
    else if (index_dim == 2) ndata = 3;
    else if (index_dim == 3) ndata = 6;
    else {
        cgi_error("invalid value for IndexDimension");
        *ier = CG_ERROR;
        return;
    }
    for (n = 0; n < ndata; n++)
        i_diffusion_model[n] = (int)diffusion_model[n];
    *ier = cg_diffusion_write(i_diffusion_model);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_model_write_f, CG_MODEL_WRITE_F) (STR_PSTR(ModelLabel),
	cgsize_t *ModelType, cgsize_t *ier STR_PLEN(ModelLabel))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(ModelLabel), STR_LEN(ModelLabel),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_model_write(c_name, (CGNS_ENUMT(ModelType_t))*ModelType);
}

/*-----------------------------------------------------------------------*/

#ifdef WIN32_FORTRAN
CGNSDLL void __stdcall cg_array_write_f(STR_PSTR(ArrayName),
	cgsize_t *DataType, cgsize_t *DataDimension,
	cgsize_t *DimensionVector, void *Data, ...)
{
    va_list ap;
    cgsize_t *ier;
    char c_name[CGIO_MAX_NAME_LENGTH+1];

    va_start(ap, Data);
    if ((CGNS_ENUMT(DataType_t))*DataType == CGNS_ENUMV(Character))
        (void) va_arg(ap, int);
    ier = va_arg(ap, cgsize_t *);
    va_end(ap);
#else
CGNSDLL void FMNAME(cg_array_write_f, CG_ARRAY_WRITE_F) (STR_PSTR(ArrayName),
	cgsize_t *DataType, cgsize_t *DataDimension, cgsize_t *DimensionVector,
	void *Data, cgsize_t *ier STR_PLEN(ArrayName))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];
#endif

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(ArrayName), STR_LEN(ArrayName),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_array_write(c_name, (CGNS_ENUMT(DataType_t))*DataType,
                              (int)*DataDimension, DimensionVector, Data);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_integral_write_f, CG_INTEGRAL_WRITE_F) (
	STR_PSTR(IntegralDataName), cgsize_t *ier STR_PLEN(IntegralDataName))
{
    char c_name[CGIO_MAX_NAME_LENGTH+1];

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(IntegralDataName), STR_LEN(IntegralDataName),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_integral_write(c_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_rind_write_f, CG_RIND_WRITE_F) (
	cgsize_t *RindData, cgsize_t *ier)
{
    int n, index_dim, i_RindData[6];

    index_dim = cgi_posit_index_dim();
    for (n = 0; n < 2*index_dim; n++)
        i_RindData[n] = (int)RindData[n];
    *ier = cg_rind_write(i_RindData);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_descriptor_write_f, CG_DESCRIPTOR_WRITE_F) (
	STR_PSTR(descr_name), STR_PSTR(descr_text),
	cgsize_t *ier STR_PLEN(descr_name) STR_PLEN(descr_text))
{
    char *c_descr_text, c_descr_name[CGIO_MAX_NAME_LENGTH+1];
    int len;

/*  On Linux, the string found in STR_PTR(descr_text) is not terminated.
    Therefore it can be much longer and can't be used as is.  The value returned
    by STR_LEN(descr_text) is correct.  So the string can be terminatted properly:
    char *terminated_descr_text;
    terminated_descr_text=(char*)malloc(strlen(STR_PTR(descr_text))+1);
    strcpy(terminated_descr_text,STR_PTR(descr_text));
    terminated_descr_text[STR_LEN(descr_text)]='\0';
    It's not necessary to do this here since we're always calling get_adf_c_name()
    which truncates the string to STR_LEN(descr_name) or shorter.
*/

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(descr_name), STR_LEN(descr_name),
        c_descr_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;

    len = STR_LEN(descr_text);
    c_descr_text = CGNS_NEW(char, len+1);
    string_2_C_string(STR_PTR(descr_text), len, c_descr_text, len, ier);
    if (*ier == 0) {
#if DEBUG_FTOC
        printf("c_descr_name='%s', c_descr_text='%s'\n",c_descr_name, c_descr_text);
#endif

         /* Call C-routine */
        *ier = cg_descriptor_write(c_descr_name, c_descr_text);
    }
    free(c_descr_text);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_units_write_f, CG_UNITS_WRITE_F) (
	cgsize_t *mass, cgsize_t *length, cgsize_t *time,
	cgsize_t *temperature, cgsize_t *angle, cgsize_t *ier)
{
    *ier = cg_units_write((CGNS_ENUMT(MassUnits_t))*mass,
               (CGNS_ENUMT(LengthUnits_t))*length,
               (CGNS_ENUMT(TimeUnits_t))*time,
               (CGNS_ENUMT(TemperatureUnits_t))*temperature,
               (CGNS_ENUMT(AngleUnits_t))*angle);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_unitsfull_write_f, CG_UNITSFULL_WRITE_F) (
	cgsize_t *mass, cgsize_t *length, cgsize_t *time,
	cgsize_t *temperature, cgsize_t *angle, cgsize_t *current,
	cgsize_t *amount, cgsize_t *intensity, cgsize_t *ier)
{
    *ier = cg_unitsfull_write((CGNS_ENUMT(MassUnits_t))*mass,
               (CGNS_ENUMT(LengthUnits_t))*length,
               (CGNS_ENUMT(TimeUnits_t))*time,
               (CGNS_ENUMT(TemperatureUnits_t))*temperature,
               (CGNS_ENUMT(AngleUnits_t))*angle,
               (CGNS_ENUMT(ElectricCurrentUnits_t))*current,
               (CGNS_ENUMT(SubstanceAmountUnits_t))*amount,
               (CGNS_ENUMT(LuminousIntensityUnits_t))*intensity);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_exponents_write_f, CG_EXPONENTS_WRITE_F) (
	cgsize_t *DataType, void *exponents, cgsize_t *ier)
{
    *ier = cg_exponents_write((CGNS_ENUMT(DataType_t))*DataType, exponents);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_expfull_write_f, CG_EXPFULL_WRITE_F) (
	cgsize_t *DataType, void *exponents, cgsize_t *ier)
{
    *ier = cg_expfull_write((CGNS_ENUMT(DataType_t))*DataType, exponents);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_conversion_write_f, CG_CONVERSION_WRITE_F) (
	cgsize_t *DataType, void *ConversionFactors, cgsize_t *ier)
{
    *ier = cg_conversion_write((CGNS_ENUMT(DataType_t))*DataType,
               ConversionFactors);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_dataclass_write_f, CG_DATACLASS_WRITE_F) (
	cgsize_t *dataclass, cgsize_t *ier)
{
    *ier = cg_dataclass_write((CGNS_ENUMT(DataClass_t))*dataclass);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_gridlocation_write_f, CG_GRIDLOCATION_WRITE_F) (
	cgsize_t *GridLocation, cgsize_t *ier)
{
    *ier = cg_gridlocation_write((CGNS_ENUMT(GridLocation_t))*GridLocation);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_ordinal_write_f, CG_ORDINAL_WRITE_F) (
	cgsize_t *Ordinal, cgsize_t *ier)
{
    *ier = cg_ordinal_write((int)*Ordinal);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_link_write_f, CG_LINK_WRITE_F) (
	STR_PSTR(nodename), STR_PSTR(filename), STR_PSTR(name_in_file), cgsize_t *ier
	STR_PLEN(nodename)  STR_PLEN(filename)  STR_PLEN(name_in_file))
{
    char n_name[CGIO_MAX_NAME_LENGTH+1];
    char f_name[CGIO_MAX_FILE_LENGTH+1];
    char i_name[CGIO_MAX_LINK_LENGTH+1];

    string_2_C_string(STR_PTR(nodename), STR_LEN(nodename),
        n_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    string_2_C_string(STR_PTR(filename), STR_LEN(filename),
        f_name, CGIO_MAX_FILE_LENGTH, ier);
    if (*ier) return;
    string_2_C_string(STR_PTR(name_in_file), STR_LEN(name_in_file),
        i_name, CGIO_MAX_LINK_LENGTH, ier);
    if (*ier) return;

    *ier = cg_link_write(n_name, f_name, i_name);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_user_data_write_f, CG_USER_DATA_WRITE_F) (
	STR_PSTR(dataname), cgsize_t *ier STR_PLEN(dataname))
{
    char d_name[CGIO_MAX_NAME_LENGTH+1];

    string_2_C_string(STR_PTR(dataname), STR_LEN(dataname),
        d_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier == 0)
        *ier = cg_user_data_write(d_name);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      General Delete Function                      *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_delete_node_f, CG_DELETE_NODE_F) (STR_PSTR(node_name),
	cgsize_t *ier STR_PLEN(node_name))
{
/* ici */
    char c_name[CGIO_MAX_NAME_LENGTH+1];

     /* convert Fortran-text-string to a C-string */
    string_2_C_string(STR_PTR(node_name), STR_LEN(node_name),
        c_name, CGIO_MAX_NAME_LENGTH, ier);
    if (*ier) return;
    *ier = cg_delete_node(c_name);
#if DEBUG_FTOC
    printf("\n  Deleting node ='%s'\n", c_name);
#endif
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *\
 *      Error Handling Functions                                         *
\* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

CGNSDLL void FMNAME(cg_get_error_f, CG_GET_ERROR_F) (
	STR_PSTR(errmsg) STR_PLEN(errmsg))
{
    cgsize_t ierr;

    string_2_F_string ((char *)cg_get_error(), STR_PTR(errmsg),
        STR_LEN(errmsg), &ierr);
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_error_exit_f, CG_ERROR_EXIT_F) ()
{
    cg_error_exit();
}

/*-----------------------------------------------------------------------*/

CGNSDLL void FMNAME(cg_error_print_f, CG_ERROR_PRINT_F) ()
{
    cg_error_print();
}

/*-----------------------------------------------------------------------*/

static void exit_on_error(int is_fatal, char *errmsg)
{
    if (is_fatal) {
        fprintf(stderr, "FATAL ERROR:%s\n", errmsg);
        exit(1);
    }
}

CGNSDLL void FMNAME(cg_exit_on_error_f, CG_EXIT_ON_ERROR_F) (cgsize_t *flag)
{
    cg_error_handler(*flag ? exit_on_error : NULL);
}

