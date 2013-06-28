/*
 * cgnsutil.h - CGNS utility header
 */

#ifndef _CGNSUTIL_H_
#define _CGNSUTIL_H_

#include "cgnslib.h"

#ifndef CGNSTYPES_H
# define cgsize_t int
# define CG_MAX_INT32 0x7FFFFFFF
#endif
#ifndef CGNS_ENUMT
# define CGNS_ENUMT(e) e
# define CGNS_ENUMV(e) e
#endif

typedef struct _DESC {
    int id;
    char name[33];
    char *desc;
} DESC;

typedef struct _VERTEX {
    cgsize_t id;
    double x, y, z, w;
} VERTEX;

typedef struct _ELEMSET {
    int id;
    char name[33];
    int type;
    cgsize_t start;
    cgsize_t end;
    int nbndry;
    cgsize_t *conn;
    cgsize_t *parent;
} ELEMSET;

typedef struct _INTERFACE {
    int id;
    char name[33];
    cgsize_t range[3][2];
    char d_name[33];
    cgsize_t d_range[3][2];
    int transform[3];
    int d_zone;
} INTERFACE;

typedef struct _CONNECT {
    int id;
    char name[33];
    int type;
    int location;
    int ptype;
    cgsize_t npnts;
    cgsize_t *pnts;
    char d_name[33];
    int d_ztype;
    int d_ptype;
    cgsize_t d_npnts;
    cgsize_t *d_pnts;
    int d_zone;
} CONNECT;

typedef struct _BOCO {
    int id;
    char name[33];
    int type;
    int ptype;
    cgsize_t npnts;
    cgsize_t *pnts;
    int n_index[3];
    cgsize_t n_cnt;
    int n_type;
    double *n_list;
} BOCO;

typedef struct _FIELD {
    int id;
    char name[33];
    int datatype;
    int units[5];
    int dataclass;
    int convtype;
    double dataconv[2];
    int exptype;
    double exponent[5];
    double *data;
} FIELD;

typedef struct _SOLUTION {
    int id;
    char name[33];
    int location;
    int rind[3][2];
    cgsize_t size;
    int units[5];
    int dataclass;
    int nflds;
    FIELD *flds;
    int ndesc;
    DESC *desc;
} SOLUTION;

typedef struct _ZONE {
    int id;
    char name[33];
    int type;
    int idim;
    cgsize_t dim[3];
    int units[5];
    int dataclass;
    int datatype;
    int vertflags;
    cgsize_t nverts;
    VERTEX *verts;
    int nesets;
    ELEMSET *esets;
    int nints;
    INTERFACE *ints;
    int nconns;
    CONNECT *conns;
    int nbocos;
    BOCO *bocos;
    int nsols;
    SOLUTION *sols;
    void *user;
    int ndesc;
    DESC *desc;
} ZONE;

extern int nZones;
extern ZONE *Zones;

extern int baseunits[5];
extern int baseclass;

extern int cgnsfn;
extern int cgnsbase;

extern int element_node_counts[];

#ifdef __cplusplus
extern "C" {
#endif

void FATAL (char *procname, char *errmsg);

ZONE *new_zone (int count);
VERTEX *new_vertex (cgsize_t nverts);
ELEMSET *new_elemset (int nsets);
INTERFACE *new_interface (int nints);
CONNECT *new_connect (int nconns);
BOCO *new_boco (int nbocos);
DESC *new_desc (int ndesc);
SOLUTION *new_solution (int nsols);
FIELD *new_field (int nfields, cgsize_t size);

cgsize_t vertex_index (ZONE *z, cgsize_t i, cgsize_t j, cgsize_t k);
cgsize_t cell_index (ZONE *z, cgsize_t i, cgsize_t j, cgsize_t k);
cgsize_t solution_index (ZONE *z, SOLUTION *s, cgsize_t i, cgsize_t j, cgsize_t k);

int file_exists (char *file);
int is_executable (char *file);
char *find_executable (char *exename);
char *find_file (char *filename, char *exename);
int same_file (char *file1, char *file2);
char *temporary_file (char *basename);
void copy_file (char *oldfile, char *newfile);

int open_cgns (char *cgnsfile, int read_only);
int find_base (char *basename);

void read_cgns (void);
int read_zones (void);
void read_zone_data (int izone);
cgsize_t read_zone_grid (int izone);
int read_zone_element (int izone);
int structured_elements (int izone);
int read_zone_interface (int izone);
int read_zone_connect (int izone);
int read_zone_boco (int izone);
int read_zone_solution (int izone);
cgsize_t read_solution_field (int izone, int isol, int ifld);
int read_units (int units[5]);

void write_cgns (void);
void write_zones (void);
void write_zone_data (int izone);
void write_zone_grid (int izone);
void write_zone_element (int izone);
void write_zone_interface (int izone);
void write_zone_connect (int izone);
void write_zone_boco (int izone);
void write_zone_solution (int izone, int isol);
void write_solution_field (int izone, int isol, int ifld);

double volume_tet (VERTEX *v1, VERTEX *v2, VERTEX *v3, VERTEX *v4);
double volume_pyr (VERTEX *v1, VERTEX *v2, VERTEX *v3,
                   VERTEX *v4, VERTEX *v5);
double volume_wdg (VERTEX *v1, VERTEX *v2, VERTEX *v3,
                   VERTEX *v4, VERTEX *v5, VERTEX *v6);
double volume_hex (VERTEX *v1, VERTEX *v2, VERTEX *v3, VERTEX *v4,
                   VERTEX *v5, VERTEX *v6, VERTEX *v7, VERTEX *v8);
double volume_element (int nnodes, VERTEX *v[]);
double cell_volume (ZONE *z, cgsize_t i, cgsize_t j, cgsize_t k);
void vertex_volumes (int izone);

void cell_vertex_solution (int izone, int isol, int weight);
void cell_center_solution (int izone, int isol, int weight);

#ifdef __cplusplus
}
#endif

#endif

