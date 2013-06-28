/*
 * cgnsversion - this changes the version number of a CGNS file.
 *
 * This only works for versions 1.2 and later, since changes
 * were made in the ADF routines between 1.1 and 1.2 that
 * prevent the 1.1 CGNS library being able to read a file
 * written using the later ADF routines (CGNS version 1.2+)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
# include <io.h>
# define access _access
# define unlink _unlink
#else
# include <unistd.h>
#endif
#include "getargs.h"
#include "cgnslib.h"
#include "cgns_header.h"
#include "cgns_io.h"

#ifndef CGNSTYPES_H
# define cgsize_t int
#endif
#ifndef CGNS_ENUMT
# define CGNS_ENUMT(e) e
# define CGNS_ENUMV(e) e
#endif

static cgns_base *CurrentBase;
static cgns_zone *CurrentZone;

static char *inpfile = NULL;
static char *outfile = NULL;

static cgns_file *cgfile;

static int inpcgio, outcgio;
static double inproot = -1.0;
static double outroot = -1.0;

static float FloatVersion;
static int LibraryVersion = CGNS_VERSION;
static int FileVersion;
static int FromVersion;

static int VersionList[] = {
    1200, 1270, 2000, 2100, 2200, 2300, 2400, 2420, 2460,
    2500, 2510, 2520, 2530, 2540, 2550, 3000, 3100
};
#define nVersions (sizeof(VersionList)/sizeof(int))

static char datatype[CGIO_MAX_NAME_LENGTH+1];
static char label[CGIO_MAX_NAME_LENGTH+1];
static char linkfile[CGIO_MAX_FILE_LENGTH+1];
static char linkpath[CGIO_MAX_LINK_LENGTH+1];

static int numdims;
static cgsize_t dimvals[CGIO_MAX_DIMENSIONS];

/* command line options */

static int verbose = 0;
static int keep_links = 1;
static int keep_nodes = 0;

static char options[] = "vrk";

static char *usgmsg[] = {
    "usage  : cgnsversion [options] Version CGNSfile [CGNSoutfile]",
    "options:",
    "   -v : verbose output",
    "   -r : remove links",
    "   -k : keep all nodes",
    NULL
};


/*--------------------------------------------------------------------*/

static void error_exit (char *msg, int errcode)
{
    if (errcode > 0) {
        char errmsg[128];
        cgio_error_message (errmsg);
        fprintf (stderr, "%s:%s\n", msg, errmsg);
    }
    else if (errcode < 0)
        fprintf (stderr, "%s:%s\n", msg, cg_get_error());
    else
        fprintf (stderr, "%s\n", msg);

    /* close files, delete temporary file, and exit */

    cgio_cleanup ();
    if (outfile != NULL) unlink (outfile);
    exit (1);
}

/*--------------------------------------------------------------------*/

static size_t get_size (char *type, int nd, cgsize_t *dims)
{
    int n;
    size_t size = 1;

    if (nd < 1) return 0;
    for (n = 0; n < nd; n++)
        size *= (size_t)dims[n];
    switch (*type) {
        case 'B':
        case 'C':
            if (type[1] == '1') return size;
            break;
        case 'I':
        case 'U':
            if (type[1] == '4') return size * sizeof(int);
            if (type[1] == '8') return size * sizeof(cglong_t);
            break;
        case 'R':
            if (type[1] == '4') return size * sizeof(float);
            if (type[1] == '8') return size * sizeof(double);
            break;
        case 'X':
            if (type[1] == '4') return 2 * size * sizeof(float);
            if (type[1] == '8') return 2 * size * sizeof(double);
            break;
    }
    return 0;
}

/*--------------------------------------------------------------------*/

static void copy_1200 (double inpid, double outid)
{
    size_t size;
    void *data;

    if (cgio_get_label (inpcgio, inpid, label))
        error_exit ("cgio_get_label", 1);
    if (cgio_set_label (outcgio, outid, label))
        error_exit ("cgio_set_label", 1);

    if (cgio_get_data_type (inpcgio, inpid, datatype))
        error_exit ("cgio_get_data_type", 1);
    if (cgio_get_dimensions (inpcgio, inpid, &numdims, dimvals))
        error_exit ("cgio_get_dimensions", 1);
    if (numdims == 0) return;
    if (FromVersion == 1200 && numdims == 1) {
        numdims = 2;
        dimvals[1] = dimvals[0];
        dimvals[0] = 1;
    }
    if (FileVersion == 1200 && numdims == 2) {
        numdims = 1;
        dimvals[0] = dimvals[1];
    }
    if (cgio_set_dimensions (outcgio, outid, datatype, numdims, dimvals))
        error_exit ("cgio_set_dimensions", 1);

    size = get_size (datatype, numdims, dimvals);
    if (size == 0) return;
    if ((data = malloc (size)) == NULL)
        error_exit ("malloc failed for data", 0);
    if (cgio_read_all_data (inpcgio, inpid, data))
        error_exit ("cgio_read_all_data", 1);
    if (cgio_write_all_data (outcgio, outid, data))
        error_exit ("cgio_write_all_data", 1);
    free (data);
}

/*--------------------------------------------------------------------*/

static double get_child_id (int cgio, double parid, int nchild)
{
    int cnt;
    double childid;

    if (cgio_children_ids (cgio, parid, nchild, 1, &cnt, &childid))
        error_exit ("cgio_children_ids", 1);
    return childid;
}

/*--------------------------------------------------------------------*/

static double create_node (double inpid, double parid, int recurse)
{
    int nchild, nc, len_ret;
    char name[CGIO_MAX_NAME_LENGTH+1],type[CGIO_MAX_NAME_LENGTH+1];
    double outid, inpchild;

    if (cgio_get_name (inpcgio, inpid, name))
        error_exit ("cgio_get_name", 1);
    if (cgio_get_label (inpcgio, inpid, type))
        error_exit ("cgio_get_label", 1);

    if (!keep_nodes) {
        if (FileVersion < 2100 &&
            0 == strcmp (type, "UserDefinedData_t"))
            return 0.0;
        if (FileVersion < 2400 &&
           (0 == strcmp (type, "AdditionalUnits_t") ||
            0 == strcmp (type, "AdditionalExponents_t")))
            return 0.0;
    }

    if (cgio_is_link (inpcgio, inpid, &len_ret))
        error_exit ("cgio_is_link", 1);

    /* create link */

    if (len_ret > 0 && keep_links) {
        if (cgio_get_link (inpcgio, inpid, linkfile, linkpath))
            error_exit ("cgio_get_link", 1);
        if (cgio_create_link (outcgio, parid, name, linkfile,
                linkpath, &outid))
            error_exit ("cgio_create_link", 1);
        return 0.0;
    }

    /* create node */

    if (cgio_create_node (outcgio, parid, name, &outid))
        error_exit ("cgio_create_node", 1);

/* check on Celcius <-> Celsius */

    if (cgio_copy_node (inpcgio, inpid, outcgio, outid))
        error_exit ("cgio_copy_node", 1);

    if (!recurse) return outid;

    /* recurse on children */

    if (cgio_number_children (inpcgio, inpid, &nchild))
        error_exit ("cgio_number_children", 1);
    if (nchild < 1) return outid;

    if (FileVersion < 2100 &&
        0 == strcmp (type, "FlowEquationSet_t")) {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            if (cgio_get_label (inpcgio, inpchild, label))
                error_exit ("cgio_get_label", 1);
            if (strcmp (label, "ThermalRelaxationModel_t") &&
                strcmp (label, "ChemicalKineticsModel_t") &&
                strcmp (label, "EMElectricFieldModel_t") &&
                strcmp (label, "EMMagneticFieldModel_t") &&
                strcmp (label, "EMConductivityModel_t"))
                create_node (inpchild, outid, 1);
        }
    }
    else if (FileVersion < 2400 &&
        0 == strcmp (type, "FlowEquationSet_t")) {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            if (cgio_get_label (inpcgio, inpchild, label))
                error_exit ("cgio_get_label", 1);
            if (strcmp (label, "EMElectricFieldModel_t") &&
                strcmp (label, "EMMagneticFieldModel_t") &&
                strcmp (label, "EMConductivityModel_t"))
                create_node (inpchild, outid, 1);
        }
    }
    else if (FileVersion < 2400 &&
        0 == strcmp (type, "UserDefinedData_t")) {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            if (cgio_get_label (inpcgio, inpchild, label))
                error_exit ("cgio_get_label", 1);
            if (strcmp (label, "UserDefinedData_t") &&
                strcmp (label, "FamilyName_t") &&
                strcmp (label, "Ordinal_t") &&
                strcmp (label, "GridLocation_t") &&
                strcmp (label, "IndexArray_t") &&
                strcmp (label, "IndexRange_t"))
                create_node (inpchild, outid, 1);
        }
    }
    else if (FileVersion < 2400 &&
        0 == strcmp (type, "GridConnectivity1to1_t")) {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            if (cgio_get_label (inpcgio, inpchild, label))
                error_exit ("cgio_get_label", 1);
            if (strcmp (label, "GridConnectivityProperty_t"))
                create_node (inpchild, outid, 1);
        }
    }
    else if (FileVersion < 2400 &&
        0 == strcmp (type, "Family_t")) {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            if (cgio_get_label (inpcgio, inpchild, label))
                error_exit ("cgio_get_label", 1);
            if (strcmp (label, "RotatingCoordinates_t"))
                create_node (inpchild, outid, 1);
        }
    }
    else if (FileVersion < 2400 &&
        0 == strcmp (type, "FamilyBC_t")) {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            if (cgio_get_label (inpcgio, inpchild, label))
                error_exit ("cgio_get_label", 1);
            if (strcmp (label, "BCDataSet_t"))
                create_node (inpchild, outid, 1);
        }
    }
    else {
        for (nc = 1; nc <= nchild; nc++) {
            inpchild = get_child_id (inpcgio, inpid, nc);
            create_node (inpchild, outid, 1);
        }
    }
    return outid;
}

/*--------------------------------------------------------------------*/

void fix_name (int cgio, double pid, double id, char *newname)
{
    char oldname[CGIO_MAX_NAME_LENGTH+1];

    if (cgio_get_name (cgio, id, oldname))
        error_exit ("cgio_get_name", 1);
    if (strcmp (oldname, newname) &&
        cgio_set_name (cgio, pid, id, newname))
        error_exit ("cgio_set_name", 1);
}

/*--------------------------------------------------------------------*/

void CGI_write_dataset(double parent_id, cgns_dataset *dataset,
                       CGNS_ENUMT(GridLocation_t) location)
{
    int n;
    cgsize_t dim_vals;
    const char *type_name;

    if (dataset->link && keep_links) {
        create_node (dataset->id, parent_id, 0);
        return;
    }

    /* BCDataSet_t */
    type_name = cg_BCTypeName(dataset->type);
    dim_vals = (int)strlen(type_name);
    if (cgio_new_node (outcgio, parent_id, dataset->name, "BCDataSet_t",
            "C1", 1, &dim_vals, type_name, &dataset->id))
        error_exit ("cgio_new_node", 1);

    /* GridLocation */
    if (FileVersion >= 2400) location = dataset->location;
    if (location != CGNS_ENUMV(Vertex) &&
        (FileVersion == 1200 || FileVersion >= 2400)) {
        double dummy_id;
        type_name = cg_GridLocationName(location);
        dim_vals = (cgsize_t)strlen(type_name);
        if (cgio_new_node (outcgio, dataset->id, "GridLocation",
                "GridLocation_t", "C1", 1, &dim_vals,
                type_name, &dummy_id))
            error_exit ("cgio_new_node", 1);
    }

    /* DirichletData */
    if (dataset->dirichlet)
        create_node (dataset->dirichlet->id, dataset->id, 1);

    /* NeumannData */
    if (dataset->neumann)
        create_node (dataset->neumann->id, dataset->id, 1);

     /* Descriptor_t */
    for (n=0; n<dataset->ndescr; n++)
        create_node (dataset->descr[n].id, dataset->id, 1);

    /* ReferenceState_t */
    if (dataset->state) create_node (dataset->state->id, dataset->id, 1);

    /* DataClass_t */
    if (dataset->data_class &&
        cgi_write_dataclass (dataset->id, dataset->data_class))
        error_exit ("cgi_write_dataclass", -1);

    /* DimensionalUnits_t */
    if (dataset->units) create_node (dataset->units->id, dataset->id, 1);

    if (keep_nodes || FileVersion >= 2100) {
        /* UserDefinedData_t */
        for (n=0; n<dataset->nuser_data; n++)
            create_node (dataset->user_data[n].id, dataset->id, 1);
    }
    if (dataset->ptset && (keep_nodes || FileVersion >= 2400))
        create_node (dataset->ptset->id, dataset->id, 1);
}

/*--------------------------------------------------------------------*/

void CGI_write_boco (double parent_id, cgns_boco *boco)
{
    int n;
    cgsize_t dim_vals;
    double dummy_id;
    const char *type_name;
    CGNS_ENUMT(BCType_t) type;
    CGNS_ENUMT(GridLocation_t) location = boco->location;

    if (boco->link && keep_links) {
        create_node (boco->id, parent_id, 0);
        return;
    }

    /* BC_t */
    type = boco->type;
    if (type == CGNS_ENUMV(FamilySpecified) && FileVersion < 2000) {
        printf ("WARNING:BC type FamilySpecified changed to UserDefined\n");
        printf ("        for boundary condition \"%s\"\n", boco->name);
        type = CGNS_ENUMV(BCTypeUserDefined);
    }
    type_name = cg_BCTypeName(type);
    dim_vals = (int)strlen(type_name);
    if (cgio_new_node (outcgio, parent_id, boco->name, "BC_t",
            "C1", 1, &dim_vals, type_name, &boco->id))
        error_exit ("cgio_new_node", 1);

    if (boco->ptset) {
        CGNS_ENUMT(PointSetType_t) ptype = boco->ptset->type;

        /* handle the various combinations of PointSetType and GridLocation */

        if (ptype == CGNS_ENUMV(PointList) ||
            ptype == CGNS_ENUMV(PointRange)) {
            if (location != CGNS_ENUMV(Vertex)) {
                if (FileVersion == 1200 || FileVersion >= 2300) {
                    if (ptype == CGNS_ENUMV(PointList))
                      ptype = CGNS_ENUMV(ElementList);
                    else
                      ptype = CGNS_ENUMV(ElementRange);
                    location = CGNS_ENUMV(Vertex);
                }
            }
        }
        else if (ptype == CGNS_ENUMV(ElementList) ||
                 ptype == CGNS_ENUMV(ElementRange)) {
            if (FileVersion > 1200 && FileVersion < 2300) {
                if (ptype == CGNS_ENUMV(ElementList))
                    ptype = CGNS_ENUMV(PointList);
                else
                    ptype = CGNS_ENUMV(PointRange);
                location = CGNS_ENUMV(FaceCenter);
            }
            else
                location = CGNS_ENUMV(Vertex);
        }
        else
            ptype = 0;

        if (ptype) {
            type_name = cg_PointSetTypeName(ptype);
            if (cgio_create_node (outcgio, boco->id, type_name, &dummy_id))
                error_exit ("cgio_create_node", 1);
            if (CurrentZone->type == CGNS_ENUMV(Unstructured) &&
                ((FileVersion == 1200 &&
                  (ptype == CGNS_ENUMV(ElementList) ||
                   ptype == CGNS_ENUMV(ElementRange))) ||
                 (FromVersion == 1200 &&
                  (boco->ptset->type == CGNS_ENUMV(ElementList) ||
                   boco->ptset->type == CGNS_ENUMV(ElementRange)))))
                copy_1200 (boco->ptset->id, dummy_id);
            else {
                if (cgio_copy_node (inpcgio, boco->ptset->id,
                        outcgio, dummy_id))
                    error_exit ("cgio_copy_node", 1);
            }
        }
    }

    /* GridLocation */
    if (location != CGNS_ENUMV(Vertex) && FileVersion > 1200) {
        type_name = cg_GridLocationName(location);
        dim_vals = (cgsize_t)strlen(type_name);
        if (cgio_new_node (outcgio, boco->id, "GridLocation",
               "GridLocation_t", "C1", 1, &dim_vals, type_name, &dummy_id))
            error_exit ("cgio_new_node", 1);
    }

     /* FamilyName_t */
    if (boco->family_name[0]!='\0') {
        if (FileVersion == 1200) {
            if (cgio_new_node (outcgio, boco->id, boco->family_name,
                    "FamilyName_t", "MT",0, &dim_vals, NULL, &dummy_id))
                error_exit ("cgio_new_node", 1);
        }
        else {
            dim_vals = (cgsize_t)strlen(boco->family_name);
            if (cgio_new_node (outcgio, boco->id, "FamilyName",
                    "FamilyName_t", "C1",1, &dim_vals,
                    boco->family_name, &dummy_id))
                error_exit ("cgio_new_node", 1);
        }
    }

    /* BCDataSet_t */
    for (n=0; n<boco->ndataset; n++)
        CGI_write_dataset (boco->id, &boco->dataset[n], boco->location);

     /* InwardNormalIndex */
    if (boco->Nindex) {
        dim_vals = CurrentZone->index_dim;
        if (cgio_new_node (outcgio, boco->id, "InwardNormalIndex",
                "\"int[IndexDimension]\"", "I4", 1, &dim_vals,
                boco->Nindex, &boco->index_id))
            error_exit ("cgio_new_node", 1);
    }

     /* InwardNormalList */
    if (boco->normal) create_node (boco->normal->id, boco->id, 0);

    /* Descriptor_t */
    for (n=0; n<boco->ndescr; n++)
        create_node (boco->descr[n].id, boco->id, 1);

    /* ReferenceState_t */
    if (boco->state) create_node (boco->state->id, boco->id, 1);

    /* DataClass_t */
    if (boco->data_class &&
        cgi_write_dataclass(boco->id, boco->data_class))
        error_exit (NULL, 0);

    /* DimensionalUnits_t */
    if (boco->units) create_node (boco->units->id, boco->id, 1);

     /* Ordinal_t */
    if (boco->ordinal && cgi_write_ordinal(boco->id, boco->ordinal))
        error_exit ("cgi_write_ordinal", -1);

    if (keep_nodes || FileVersion >= 2100) {
        /* BCProperty_t */
        if (boco->bprop && (keep_nodes || FileVersion >= 2200))
            create_node (boco->bprop->id, boco->id, 1);
        /* UserDefinedData_t */
        for (n=0; n<boco->nuser_data; n++)
            create_node (boco->user_data[n].id, boco->id, 1);
    }
}

/*--------------------------------------------------------------------*/

void CGI_write_zboco(double parent_id, cgns_zboco *zboco) {
    int n;

    if (zboco->link && keep_links) {
        create_node (zboco->id, parent_id, 0);
        return;
    }

    /* ZoneBC_t */
    if (cgio_new_node (outcgio, parent_id, "ZoneBC", "ZoneBC_t",
            "MT", 0, 0, 0, &zboco->id))
        error_exit ("cgio_new_node", 1);

    /* BC_t */
    for (n=0; n<zboco->nbocos; n++)
        CGI_write_boco (zboco->id, &zboco->boco[n]);

    /* Descriptor_t */
    for (n=0; n<zboco->ndescr; n++)
        create_node (zboco->descr[n].id, zboco->id, 1);

    /* ReferenceState_t */
    if (zboco->state) create_node (zboco->state->id, zboco->id, 1);

    /* DataClass_t */
    if (zboco->data_class &&
        cgi_write_dataclass (zboco->id, zboco->data_class))
        error_exit (NULL, 0);

    /* DimensionalUnits_t */
    if (zboco->units) create_node (zboco->units->id, zboco->id, 1);

    if (keep_nodes || FileVersion >= 2100) {
        /* UserDefinedData_t */
        for (n=0; n<zboco->nuser_data; n++)
            create_node (zboco->user_data[n].id, zboco->id, 1);
    }
}

/*--------------------------------------------------------------------*/

void CGI_write_conns (double parent_id, cgns_conn *conn)
{
    int n;
    cgsize_t dim_vals;
    double dummy_id;
    const char *type_name;
    cgns_ptset *ptset;
    CGNS_ENUMT(GridLocation_t) location = conn->location;

    if (conn->link && keep_links) {
        create_node (conn->id, parent_id, 0);
        return;
    }

    dim_vals = (cgsize_t)strlen(conn->donor);
    if (cgio_new_node (outcgio, parent_id, conn->name, "GridConnectivity_t",
            "C1", 1, &dim_vals, conn->donor, &conn->id))
        error_exit ("cgio_new_node", 1);

     /* GridConnectivityType_t */
    type_name = cg_GridConnectivityTypeName(conn->type);
    dim_vals = (cgsize_t)strlen(type_name);
    if (cgio_new_node (outcgio, conn->id, "GridConnectivityType",
            "GridConnectivityType_t", "C1", 1, &dim_vals,
            type_name, &dummy_id))
        error_exit ("cgio_new_node", 1);

     /* write GridLocation */
    if (location != CGNS_ENUMV(Vertex)) {
        if (location != CGNS_ENUMV(CellCenter) && FileVersion < 2300) {
            if (FileVersion == 2200)
                location = CGNS_ENUMV(FaceCenter);
            else                                
                location = CGNS_ENUMV(CellCenter);
        }
        type_name = cg_GridLocationName(location);
        dim_vals = (cgsize_t)strlen(type_name);
        if (cgio_new_node (outcgio, conn->id, "GridLocation",
                "GridLocation_t", "C1", 1, &dim_vals,
                type_name, &dummy_id))
            error_exit ("cgio_new_node", 1);
    }
    if (location != conn->location) {
        printf ("WARNING:GridLocation changed from %s to %s\n",
            cg_GridLocationName(conn->location), type_name);
        printf ("        for connectivity \"%s\"\n", conn->name);
    }

    /* PointRange or PointList */

    ptset = &(conn->ptset);
    create_node (ptset->id, conn->id, 0);

    /* Cell or Point ListDonor */
    ptset = &(conn->dptset);
    if (FileVersion > 1200) {
        create_node (ptset->id, conn->id, 0);
        if (conn->interpolants)
            create_node (conn->interpolants->id, conn->id, 1);
    }
    else {
        double donorid;
        CGNS_ENUMT(ZoneType_t) dtype = CGNS_ENUMV(ZoneTypeNull);
        CGNS_ENUMT(PointSetType_t) ptype = ptset->type;
        for (n = 0; n < CurrentBase->nzones; n++) {
            if (0 == strcmp (conn->donor, CurrentBase->zone[n].name)) {
                dtype = CurrentBase->zone[n].type;
                break;
            }
        }
        if (dtype == CGNS_ENUMV(Structured)) {
            if (cgio_new_node (outcgio, conn->id, "StructuredDonor",
                    "StructuredDonor_t", "MT", 0,
                    &dim_vals, NULL, &dummy_id))
                error_exit ("cgio_new_node", 1);
            donorid = create_node (ptset->id, dummy_id, 0);
            ptype = CGNS_ENUMV(PointListDonor);
        }
        else if (dtype == CGNS_ENUMV(Unstructured)) {
            if (cgio_new_node (outcgio, conn->id, "UnstructuredDonor",
                    "UnstructuredDonor_t", "MT", 0,
                    &dim_vals, NULL, &dummy_id))
                error_exit ("cgio_new_node", 1);
            donorid = create_node (ptset->id, dummy_id, 0);
            if (ptype != CGNS_ENUMV(PointListDonor) &&
                ptype != CGNS_ENUMV(CellListDonor))
                ptype = (location == CGNS_ENUMV(Vertex) ||
                         conn->interpolants == NULL) ?
                         CGNS_ENUMV(PointListDonor) :
                         CGNS_ENUMV(CellListDonor);
            if (ptype == CGNS_ENUMV(CellListDonor) && conn->interpolants)
                create_node (conn->interpolants->id, dummy_id, 0);
        }
        else {
            printf ("ERROR:couldn't find donor zone \"%s\" for \"%s\"\n",
                conn->donor, conn->name);
            printf ("     therefore donor data was not written\n");
            donorid = 0.0;
        }
        if (ptype != ptset->type && donorid > 0.0) {
            type_name = cg_PointSetTypeName(ptype);
            fix_name (outcgio, dummy_id, donorid, (char *)type_name);
            printf ("WARNING:%s changed to %s\n",
                cg_PointSetTypeName(ptset->type), type_name);
            printf ("        for connectivity \"%s\"\n", conn->name);
        }
    }

    /* Descriptor_t */
    for (n=0; n<conn->ndescr; n++)
        create_node (conn->descr[n].id, conn->id, 1);

     /* Ordinal_t */
    if (conn->ordinal && cgi_write_ordinal(conn->id, conn->ordinal))
        cg_error_exit();

    if (keep_nodes || FileVersion >= 2100) {
        /* GridConnectivityProperty_t */
        if (conn->cprop && (keep_nodes || FileVersion >= 2200))
            create_node (conn->cprop->id, conn->id, 1);
        /* UserDefinedData_t */
        for (n=0; n<conn->nuser_data; n++)
            create_node (conn->user_data[n].id, conn->id, 1);
    }
}

/*--------------------------------------------------------------------*/

void CGI_write_zconn(double parent_id, cgns_zconn *zconn) {
    int n;

    if (zconn->link && keep_links) {
        create_node (zconn->id, parent_id, 0);
        return;
    }

    /* ZoneGridConnectivity_t */
    if (cgio_new_node (outcgio, parent_id, "ZoneGridConnectivity",
            "ZoneGridConnectivity_t", "MT", 0, 0, 0, &zconn->id))
        error_exit ("cgio_new_node", 1);

    /* GridConnectivity1to1_t */
    for (n=0; n<zconn->n1to1; n++)
        create_node (zconn->one21[n].id, zconn->id, 1);

    /* GridConnectivity_t */
    for (n=0; n<zconn->nconns; n++)
        CGI_write_conns (zconn->id, &zconn->conn[n]);

    /* OversetHoles_t */
    for (n=0; n<zconn->nholes; n++)
        create_node (zconn->hole[n].id, zconn->id, 1);

    /* Descriptor_t */
    for (n=0; n<zconn->ndescr; n++)
        create_node (zconn->descr[n].id, zconn->id, 1);

    if (keep_nodes || FileVersion >= 2100) {
        /* UserDefinedData_t */
        for (n=0; n<zconn->nuser_data; n++)
            create_node (zconn->user_data[n].id, zconn->id, 1);
    }
}

/*--------------------------------------------------------------------*/

void CGI_write_section(double parent_id, cgns_section *section)
{
    int n, data[2];
    cgsize_t dim_vals;
    double dummy_id;
    cgsize_t *elems;

    if ((section->link && keep_links) ||
        ((FromVersion < 3000 || FromVersion >= 3100) &&
         (FileVersion < 3000 || FileVersion >= 3100) &&
         section->el_type < CGNS_ENUMV(MIXED)) ||
        (FromVersion >= 3100 && FileVersion >= 3100)) {
        create_node (section->id, parent_id, 1);
        return;
    }

    if (FileVersion < 3000 && section->el_type > CGNS_ENUMV(MIXED)) {
        printf ("can't convert %s element set - skipping\n",
            cg_ElementTypeName(section->el_type));
        return;
    }

    if (section->connect->data) free(section->connect->data);
    elems = CGNS_NEW(cgsize_t, section->connect->dim_vals[0]);
    if (cgi_read_int_data(section->connect->id, section->connect->data_type,
            section->connect->dim_vals[0], elems))
        error_exit ("cgi_read_int_data", 1);

    if (section->el_type == CGNS_ENUMV(MIXED)) {
        cgsize_t ne, nelems, cnt = 0;
        CGNS_ENUMT(ElementType_t) type;
        int npe;

        nelems = section->range[1] - section->range[0] + 1;
        for (ne = 0; ne < nelems; ne++) {
            type = (CGNS_ENUMT(ElementType_t))elems[cnt];
            if (FromVersion < 3000) {
                if (type > CGNS_ENUMV(MIXED)) type++;
            }
            else if (FromVersion < 3100) {
                if (type > CGNS_ENUMV(PYRA_5) && type < CGNS_ENUMV(NGON_n)) {
                    if (type == CGNS_ENUMV(PYRA_14))
                        type = CGNS_ENUMV(PYRA_13);
                    else
                        type--;
                }
            }
            if (type > CGNS_ENUMV(NGON_n))
                npe = type - CGNS_ENUMV(NGON_n);
            else
                cg_npe (type, &npe);
            if (npe <= 0)
                error_exit("invalid element found in MIXED elements", 0);
            if (FileVersion < 3000) {
                if (type > CGNS_ENUMV(MIXED)) type--;
            }
            else if (FileVersion < 3100) {
                if (type > CGNS_ENUMV(PYRA_5) && type < CGNS_ENUMV(NGON_n)) {
                    if (type == CGNS_ENUMV(PYRA_13))
                        type = CGNS_ENUMV(PYRA_14);
                    else
                        type++;
                }
            }
            elems[cnt++] = (cgsize_t)type;
            cnt += npe;
        }
    }
    section->connect->data = (void *)elems;

    if (FileVersion >= 3000 && FileVersion < 3100 &&
        section->el_type > CGNS_ENUMV(PYRA_5) &&
        section->el_type < CGNS_ENUMV(NGON_n)) {
        if (section->el_type == CGNS_ENUMV(PYRA_13))
            section->el_type = CGNS_ENUMV(PYRA_14);
        else
            section->el_type++;
    }

    dim_vals = 2;
    data[0] = section->el_type;
    data[1] = section->el_bound;
    if (cgio_new_node (outcgio, parent_id, section->name, "Elements_t",
            "I4", 1, &dim_vals, data, &section->id))
        error_exit ("cgio_new_node", 1);

     /* ElementRange */
    if (cgio_new_node (outcgio, section->id, "ElementRange", "IndexRange_t",
            "I4", 1, &dim_vals, section->range, &dummy_id))
        error_exit ("cgio_new_node", 1);

     /* ElementConnectivity */
    if (cgi_write_array(section->id, section->connect))
        cg_error_exit();

     /* ParentData */
#ifdef CG_SPLIT_PARENT_DATA
    if (section->parelem && cgi_write_array(section->id, section->parelem))
        cg_error_exit();
    if (section->parface && cgi_write_array(section->id, section->parface))
        cg_error_exit();
#else
    if (section->parent && cgi_write_array(section->id, section->parent))
        cg_error_exit();
#endif

     /* Descriptor_t */
    for (n=0; n<section->ndescr; n++)
        create_node (section->descr[n].id, section->id, 1);

    if (keep_nodes || FileVersion >= 2100) {
        /* UserDefinedData_t */
        for (n=0; n<section->nuser_data; n++)
            create_node (section->user_data[n].id, section->id, 1);
    }
    CGNS_FREE(elems);
}

/*--------------------------------------------------------------------*/

void CGI_write_zone(double parent_id, cgns_zone *zone)
{
    int n;
    cgsize_t dim_vals[2];
    double dummy_id;
    const char *type_name;

    CurrentZone = zone;
    if (zone->link && keep_links) {
        create_node (zone->id, parent_id, 0);
        return;
    }

     /* Create the Zone_t nodes */
    dim_vals[0]= zone->index_dim;
    dim_vals[1]= 3;
    if (cgio_new_node (outcgio, parent_id, zone->name, "Zone_t",
            "I4", 2, dim_vals, zone->nijk, &zone->id))
        error_exit ("cgio_new_node", 1);

     /* write ZoneType */
    type_name = cg_ZoneTypeName(zone->type);
    dim_vals[0] = (cgsize_t)strlen(type_name);
    if (cgio_new_node (outcgio, zone->id, "ZoneType", "ZoneType_t",
            "C1", 1, dim_vals, type_name, &dummy_id))
        error_exit ("cgio_new_node", 1);

    /* GridCoordinates_t */
    if (zone->nzcoor >= 1) {
        if (FileVersion >= 2000) {
            for (n=0; n<zone->nzcoor; n++)
                create_node (zone->zcoor[n].id, zone->id, 1);
        }
        else {
            dummy_id = create_node (zone->zcoor->id, zone->id, 1);
            fix_name (outcgio, zone->id, dummy_id, "GridCoordinates");
        }
    }

     /* FamilyName_t */
    if (zone->family_name[0]!='\0') {
        if (FileVersion == 1200) {
            if (cgio_new_node (outcgio, zone->id, zone->family_name,
                    "FamilyName_t", "MT",0, dim_vals, NULL, &dummy_id))
                error_exit ("cgio_new_node", 1);
        }
        else {
            dim_vals[0] = (cgsize_t)strlen(zone->family_name);
            if (cgio_new_node (outcgio, zone->id, "FamilyName",
                    "FamilyName_t", "C1",1, dim_vals,
                    zone->family_name, &dummy_id))
                error_exit ("cgio_new_node", 1);
        }
    }

    /* Elements_t */
    for (n=0; n<zone->nsections; n++)
        CGI_write_section (zone->id, &zone->section[n]);

    /* FlowSolution_t */
    for (n=0; n<zone->nsols; n++)
        create_node (zone->sol[n].id, zone->id, 1);

    /* ZoneGridConnectivity_t */
    if (zone->zconn) CGI_write_zconn(zone->id, zone->zconn);

    /* ZoneBC_t */
    if (zone->zboco) CGI_write_zboco(zone->id, zone->zboco);

    /* DescreteData_t */
    for (n=0; n<zone->ndiscrete; n++)
        create_node (zone->discrete[n].id, zone->id, 1);

    /* Descriptor_t */
    for (n=0; n<zone->ndescr; n++)
        create_node (zone->descr[n].id, zone->id, 1);

    /* ReferenceState_t */
    if (zone->state) create_node (zone->state->id, zone->id, 1);

    /* DataClass_t */
    if (zone->data_class && cgi_write_dataclass (zone->id, zone->data_class))
        error_exit ("cgi_write_dataclass", -1);

    /* DimensionalUnits_t */
    if (zone->units) create_node (zone->units->id, zone->id, 1);

    /* ConvergenceHistory_t */
    if (zone->converg) create_node (zone->converg->id, zone->id, 1);

    /* FlowEquationSet_t */
    if (zone->equations) create_node (zone->equations->id, zone->id, 1);

    /* IntegralData_t */
    for (n=0; n<zone->nintegrals; n++)
        create_node (zone->integral[n].id, zone->id, 1);

     /* Ordinal_t */
    if (zone->ordinal && cgi_write_ordinal (zone->id, zone->ordinal))
        cg_error_exit();

    if (keep_nodes || FileVersion >= 2000) {
        /* RigidGridMotion_t */
        for (n=0; n<zone->nrmotions; n++)
            create_node (zone->rmotion[n].id, zone->id, 1);

        /* ArbitraryGridMotion_t */
        for (n=0; n<zone->namotions; n++)
            create_node (zone->amotion[n].id, zone->id, 1);

        /* ZoneIterativeData_t */
        if (zone->ziter) create_node (zone->ziter->id, zone->id, 1);

        if (keep_nodes || FileVersion >= 2100) {
            /* UserDefinedData_t */
            for (n=0; n<zone->nuser_data; n++)
                create_node (zone->user_data[n].id, zone->id, 1);
             /* RotatingCoordinates_t */
            if (zone->rotating && (keep_nodes || FileVersion >= 2200))
                create_node (zone->rotating->id, zone->id, 1);
        }
    }
}

/*--------------------------------------------------------------------*/

void CGI_write ()
{
    cgns_base *base;
    int n, b;
    cgsize_t dim_vals;
    int data[2];
    double dummy_id;

    /* write version number */
    dim_vals = 1;
    if (cgio_new_node (outcgio, cgfile->rootid, "CGNSLibraryVersion",
            "CGNSLibraryVersion_t", "R4", 1, &dim_vals,
            &FloatVersion, &dummy_id))
        error_exit ("cgio_new_node", 1);

    /* write all CGNSBase_t nodes in ADF file */
    for (b=0; b < cgfile->nbases; b++) {
        CurrentBase = base = &(cgfile->base[b]);

        data[0]=base->cell_dim;
        data[1]=base->phys_dim;

        /* Create the CGNSBase_t nodes */
        dim_vals=2;
        if (cgio_new_node (outcgio, cgfile->rootid, base->name, "CGNSBase_t",
                "I4", 1, &dim_vals, data, &base->id))
            error_exit ("cgio_new_node", 1);

        /* Descriptor_t */
        for (n=0; n<base->ndescr; n++)
            create_node (base->descr[n].id, base->id, 1);

        /* ReferenceState_t */
        if (base->state) create_node (base->state->id, base->id, 1);

        if (keep_nodes || FileVersion >= 2200) {
            /* Gravity_t */
            if (base->gravity) create_node (base->gravity->id, base->id, 1);

            /* Axisymmetry_t */
            if (base->axisym) create_node (base->axisym->id, base->id, 1);

            /* RotatingCoordinates_t */
            if (base->rotating) create_node (base->rotating->id, base->id, 1);
        }

        /* Zone_t */
        for (n=0; n<base->nzones; n++)
            CGI_write_zone(base->id, &base->zone[n]);

       /* Family_t */
        for (n=0; n<base->nfamilies; n++)
            create_node (base->family[n].id, base->id, 1);

       /* DataClass_t */
        if (base->data_class &&
            cgi_write_dataclass (base->id, base->data_class))
            error_exit ("cgi_write_dataclass", -1);

       /* DimensionalUnits_t */
        if (base->units) create_node (base->units->id, base->id, 1);

       /* ConvergenceHistory_t */
        if (base->converg) create_node (base->converg->id, base->id, 1);

       /* FlowEquationSet_t */
        if (base->equations) create_node (base->equations->id, base->id, 1);

       /* IntegralData_t */
        for (n=0; n<base->nintegrals; n++)
            create_node (base->integral[n].id, base->id, 1);

        if (keep_nodes || FileVersion >= 2000) {
            /* SimulationType_t */
            if (base->type) {
                const char *sim_name;
                sim_name = cg_SimulationTypeName(base->type);
                dim_vals = (cgsize_t)strlen(sim_name);
                if (cgio_new_node (outcgio, base->id, "SimulationType",
                        "SimulationType_t","C1", 1,
                        &dim_vals, sim_name, &base->type_id))
                    error_exit ("cgio_new_node", 1);
            }

            /* BaseIterativeData_t */
            if (base->biter) create_node (base->biter->id, base->id, 1);

            /* UserDefinedData_t */
            if (keep_nodes || FileVersion >= 2100) {
                for (n=0; n<base->nuser_data; n++)
                    create_node (base->user_data[n].id, base->id, 1);
            }
        }
    }
}

/*--------------------------------------------------------------------*/

static char *temporary_file (char *basename)
{
    char *p, *temp;
    int n;

    if (basename == NULL || !*basename)
        basename = "cgnstmpfile";
    temp = (char *) malloc (strlen(basename) + 10);
    if (temp == NULL) {
        fprintf (stderr, "malloc failed for temp filename\n");
        exit (1);
    }
    sprintf (temp, "%s.tmp", basename);
    p = temp + strlen(temp);
    for (n = 0; n < 1000; n++) {
        sprintf (p, "%3.3d~", n);
        if (access (temp, 0)) return temp;
    }
    fprintf (stderr, "failed to create temporary filename\n");
    exit (1);
}

/*--------------------------------------------------------------------*/

int main (int argc, char **argv)
{
    int n, inpfn, filetype;
    float file_version;

#ifdef CGNSTYPES_H
# if CG_BUILD_64BIT
    fprintf(stderr, "doesn't work in 64-bit mode yet\n");
    return 1;
# endif
#endif

    if (argc < 3)
        print_usage (usgmsg, NULL);

    while ((n = getargs (argc, argv, options)) > 0) {
        switch (n) {
            case 'v':
                verbose = 1;
                break;
            case 'r':
                keep_links = 0;
                break;
            case 'k':
                keep_nodes = 1;
                break;
        }
    }

    if (argind >= argc - 1)
        print_usage (usgmsg, "version and/or CGNSfile not given");

    /* get the file version */

    FloatVersion = (float)atof(argv[argind++]);
    FileVersion = (int)(1000.0 * FloatVersion + 0.5);

    if (FileVersion < 1200) {
        fprintf (stderr,
            "ADF incompatibilities do not allow versions prior to 1.2\n");
        exit (1);
    }

    for (n = 0; n < nVersions; n++) {
        if (FileVersion == VersionList[n]) break;
    }
    if (n >= nVersions || VersionList[n] > LibraryVersion) {
        fprintf (stderr, "Version %g is not valid\n", FloatVersion);
        fprintf (stderr, "Version is one of 1.2");
        for (n = 1; n < nVersions; n++) {
            if (VersionList[n] < LibraryVersion)
                fprintf (stderr, ", %g", 0.001 * (float)VersionList[n]);
        }
        fprintf (stderr, " or %g\n", 0.001 * (float)LibraryVersion);
        exit (1);
    }

    inpfile = argv[argind++];
    if (access (inpfile, 0)) {
        fprintf (stderr, "input file \"%s\" not found\n", inpfile);
        exit (1);
    }

    /* read the input file using the CGNS routines in order to
       fill in the cgns_* structures. */

    if (cg_open (inpfile, CG_MODE_READ, &inpfn)) cg_error_exit();

    /* get version and check if valid */

    if (cg_version (inpfn, &file_version)) cg_error_exit();
    FromVersion = (int)(1000.0 * file_version + 0.5);
    if (LibraryVersion < FromVersion) {
        cg_close (inpfn);
        fprintf (stderr,
            "file version is more recent than then CGNS library version\n");
        exit (1);
    }
    if (FileVersion == FromVersion) {
        cg_close (inpfn);
        fprintf (stderr, "file version is already at version %g\n",
            FloatVersion);
        exit (1);
    }

    printf ("converting \"%s\" from version %g to %g\n",
        inpfile, file_version, FloatVersion);

    if (FileVersion < 2100) keep_links = 0;

    /* set the default cgns_file structure pointer */

    cgfile = cgi_get_file(inpfn);
    if (cgio_get_root_id(cgfile->cgio, &inproot))
        cgio_error_exit("cgio_get_root_id");
    if (cgio_get_file_type(cgfile->cgio, &filetype))
        cgio_error_exit("cgio_get_file_type");

    /* write output to temporary file */

    outfile = temporary_file (inpfile);
    printf ("writing output to temporary file \"%s\"\n", outfile);

    /* use cgio routines to write the file */

    if (cgio_open_file(outfile, 'w', filetype, &outcgio))
        cgio_error_exit("cgio_file_open");
    if (cgio_get_root_id(outcgio, &outroot))
        cgio_error_exit("cgio_get_root_id");

    inpcgio = cgfile->cgio;
    cgfile->cgio = outcgio;
    cgfile->rootid = outroot;
    CGI_write ();

    /* close using cgio routines */

    cgio_close_file(inpcgio);
    cgio_close_file(outcgio);

    /* rename the temporary file to the output name */

    if (argind < argc)
        inpfile = argv[argind];
    printf ("renaming \"%s\" to \"%s\"\n", outfile, inpfile);

    unlink (inpfile);
    if (rename (outfile, inpfile)) {
        fprintf (stderr, "rename %s -> %s failed", outfile, inpfile);
        exit (1);
    }
    free (outfile);

    return 0;
}

