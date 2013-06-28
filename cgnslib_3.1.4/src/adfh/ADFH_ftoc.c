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
#include <string.h>
#include "fortran_macros.h"
#include "ADFH.h"

#if defined(_WIN32) && defined(BUILD_DLL)
# define CGNSDLL _declspec(dllexport)
#else
# define CGNSDLL
#endif

/*-------------------------------------------------------------------*/

static void StringFtoC(char *string, int string_length,
    char *c_string, int max_len, int *ierr) {
    int i, iend;

    if (string == NULL || c_string == NULL) {
        *ierr = NULL_STRING_POINTER;
        return;
    }

    /** Skip and trailing blanks **/
    for (iend = string_length-1; iend >= 0; iend--) {
        if (string[iend] != ' ') break;
    }
    if (iend >= max_len) {
        *ierr = STRING_LENGTH_TOO_BIG;
        return;
    }

    /** Copy the non-trailing blank portion of the string **/
    for (i = 0; i <= iend; i++)
        c_string[i] = string[i];

    /** NULL terminate the C string **/
    c_string[i] = '\0';
    *ierr = NO_ERROR;
}

/*-------------------------------------------------------------------*/

static void StringCtoF(char *c_string, char *string,
    int string_length, int *ierr) {
    int i, len;

    if (c_string == NULL || string == NULL) {
        *ierr = NULL_STRING_POINTER;
        return;
    }
    len = strlen(c_string);
    /* string truncated - should this be an error ? */
    if (len > string_length) len = string_length;

    for (i = 0; i < len; i++)
        string[i] = c_string[i];
    while (i < string_length)
        string[i++] = ' ';
    *ierr = NO_ERROR;
}

/*----- ADFH Fortran to C interface -----*/

CGNSDLL void FMNAME(adfcnam, ADFCNAM) (double *ID, int *istart, int *inum,
    int *inamlen, int *inumret, STR_PSTR(names), int *ier STR_PLEN(names)) {
    char *tmp, *p = STR_PTR(names);
    int n, len = *inamlen;

    *inumret = 0;
    if (*inum < 1) {
        *ier = NO_ERROR;
        return;
    }
    tmp = (char *) malloc (*inum * (len + 1));
    if (tmp == NULL) {
        *ier = MEMORY_ALLOCATION_FAILED;
        return;
    }
    ADFH_Children_Names(*ID, *istart, *inum, len+1, inumret, tmp, ier);
    if (*ier == NO_ERROR) {
        for (n = 0; n < *inumret; n++)
            StringCtoF(&tmp[n*(len+1)], &p[n*len], len, ier);
    }
    free(tmp);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfcids, ADFCIDS) (double *ID, int *istart, int *inum,
    int *inumret, double *cIDs, int *ier) {

    ADFH_Children_IDs(*ID, *istart, *inum, inumret, cIDs, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfrid, ADFRID) (double *ID) {

    ADFH_Release_ID(*ID);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfcre, ADFCRE) (double *PID, STR_PSTR(name), double *ID,
    int *ier STR_PLEN(name)) {
    char c_name[ADF_NAME_LENGTH+1];

    StringFtoC(STR_PTR(name), STR_LEN(name), c_name, ADF_NAME_LENGTH, ier);
    if (*ier == NO_ERROR)
        ADFH_Create(*PID, c_name, ID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdclo, ADFDCLO) (double *RootID, int *ier) {

    ADFH_Database_Close(*RootID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfddel, ADFDDEL) (STR_PSTR(filename),
    int *ier STR_PLEN(filename)) {
    int len = STR_LEN(filename);
    char *c_name;

    if ((c_name = (char *) malloc (len + 1)) == NULL) {
        *ier = MEMORY_ALLOCATION_FAILED;
        return;
    }
    StringFtoC(STR_PTR(filename), len, c_name, len, ier);
    if (*ier == NO_ERROR)
        ADFH_Database_Delete(c_name, ier);
    free (c_name);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdel, ADFDEL) (double *PID, double *ID, int *ier) {

    ADFH_Delete(*PID, *ID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdgc, ADFDGC) (double *ID, int *ier) {

    ADFH_Database_Garbage_Collection(*ID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdgf, ADFDGF) (double *RootID, STR_PSTR(format),
    int *ier STR_PLEN(format)) {
    char c_format[ADF_FORMAT_LENGTH+1];

    ADFH_Database_Get_Format(*RootID, c_format, ier);
    if (*ier == NO_ERROR)
        StringCtoF(c_format, STR_PTR(format), STR_LEN(format), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdopn, ADFDOPN) (STR_PSTR(filename), STR_PSTR(status),
    STR_PSTR(format), double *RootID, int *ier STR_PLEN(filename)
    STR_PLEN(status) STR_PLEN(format)) {
    int len = STR_LEN(filename);
    char *c_filename;
    char c_status[ADF_STATUS_LENGTH+1];
    char c_format[ADF_FORMAT_LENGTH+1];

    StringFtoC(STR_PTR(status), STR_LEN(status), c_status,
        ADF_STATUS_LENGTH, ier);
    if (*ier != NO_ERROR) return;
    StringFtoC(STR_PTR(format), STR_LEN(format), c_format,
        ADF_FORMAT_LENGTH, ier);
    if (*ier != NO_ERROR) return;
    if ((c_filename = (char *) malloc (len + 1)) == NULL) {
        *ier = MEMORY_ALLOCATION_FAILED;
        return;
    }
    StringFtoC(STR_PTR(filename), len, c_filename, len, ier);
    if (*ier == NO_ERROR)
        ADFH_Database_Open(c_filename, c_status, c_format, RootID, ier);
    free (c_filename);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdsf, ADFDSF) (double *RootID, STR_PSTR(format),
    int *ier STR_PLEN(format)) {
    char c_format[ADF_FORMAT_LENGTH+1];

    StringFtoC(STR_PTR(format), STR_LEN(format), c_format,
        ADF_FORMAT_LENGTH, ier);
    if (*ier == NO_ERROR)
        ADFH_Database_Set_Format(*RootID, c_format, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfdver, ADFDVER) (double *RootID, STR_PSTR(version),
    STR_PSTR(cdate), STR_PSTR(mdate), int *ier STR_PLEN(version)
    STR_PLEN(cdate) STR_PLEN(mdate)) {
    char c_version[ADF_VERSION_LENGTH+1];
    char c_cdate[ADF_DATE_LENGTH+1];
    char c_mdate[ADF_DATE_LENGTH+1];

    ADFH_Database_Version(*RootID, c_version, c_cdate, c_mdate, ier);
    if (*ier == NO_ERROR) {
        StringCtoF(c_version, STR_PTR(version), STR_LEN(version), ier);
        StringCtoF(c_cdate, STR_PTR(cdate), STR_LEN(cdate), ier);
        StringCtoF(c_mdate, STR_PTR(mdate), STR_LEN(mdate), ier);
    }
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adferr, ADFERR) (int *ier, STR_PSTR(errstr) STR_PLEN(errstr)) {
    int err;
    char errmsg[ADF_MAX_ERROR_STR_LENGTH+1];

    ADFH_Error_Message(*ier, errmsg);
    StringCtoF(errmsg, STR_PTR(errstr), STR_LEN(errstr), &err);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfftd, ADFFTD) (double *ID, int *ier) {

    ADFH_Flush_to_Disk(*ID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfgdt, ADFGDT) (double *ID, STR_PSTR(dtype),
    int *ier STR_PLEN(dtype)) {
    char c_type[ADF_DATA_TYPE_LENGTH+1];

    ADFH_Get_Data_Type(*ID, c_type, ier);
    if (*ier == NO_ERROR)
        StringCtoF(c_type, STR_PTR(dtype), STR_LEN(dtype), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfgdv, ADFGDV) (double *ID, int *dvals, int *ier) {

    ADFH_Get_Dimension_Values(*ID, dvals, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfges, ADFGES) (int *estate, int *ier) {

    ADFH_Get_Error_State(estate, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfglb, ADFGLB) (double *ID, STR_PSTR(label),
    int *ier STR_PLEN(label)) {
    char c_label[ADF_LABEL_LENGTH+1];

    ADFH_Get_Label(*ID, c_label, ier);
    if (*ier == NO_ERROR)
        StringCtoF(c_label, STR_PTR(label), STR_LEN(label), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfglkp, ADFGLKP) (double *ID, STR_PSTR(file), STR_PSTR(name),
    int *ier STR_PLEN(file) STR_PLEN(name)) {
    char c_file[ADF_FILENAME_LENGTH+1];
    char c_name[ADF_MAX_LINK_DATA_SIZE+1];

    ADFH_Get_Link_Path(*ID, c_file, c_name, ier);
    if (*ier == NO_ERROR) {
        StringCtoF(c_file, STR_PTR(file), STR_LEN(file), ier);
        StringCtoF(c_name, STR_PTR(name), STR_LEN(name), ier);
    }
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfgnam, ADFGNAM) (double *ID, STR_PSTR(name),
    int *ier STR_PLEN(name)) {
    char c_name[ADF_NAME_LENGTH+1];

    ADFH_Get_Name(*ID, c_name, ier);
    if (*ier == NO_ERROR)
        StringCtoF(c_name, STR_PTR(name), STR_LEN(name), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfgnd, ADFGND) (double *ID, int *ndims, int *ier) {

    ADFH_Get_Number_of_Dimensions(*ID, ndims, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfgnid, ADFGNID) (double *PID, STR_PSTR(name), double *ID,
    int *ier STR_PLEN(name)) {
    int len = STR_LEN(name);
    char *c_name;

    if ((c_name = (char *) malloc (len+1)) == NULL) {
        *ier = MEMORY_ALLOCATION_FAILED;
        return;
    }
    StringFtoC(STR_PTR(name), len, c_name, len, ier);
    if (*ier == NO_ERROR)
        ADFH_Get_Node_ID(*PID, c_name, ID, ier);
    free(c_name);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfgrid, ADFGRID) (double *ID, double *RootID, int *ier) {

    ADFH_Get_Root_ID(*ID, RootID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfislk, ADFISLK) (double *ID, int *lplen, int *ier) {

    ADFH_Is_Link(*ID, lplen, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adflink, ADFLINK) (double *PID, STR_PSTR(name), STR_PSTR(file),
    STR_PSTR(nfile), double *ID, int *ier STR_PLEN(name)
    STR_PLEN(file) STR_PLEN(nfile)) {
    char c_name[ADF_NAME_LENGTH+1];
    char c_file[ADF_FILENAME_LENGTH+1];
    char c_nfile[ADF_MAX_LINK_DATA_SIZE+1];

    StringFtoC(STR_PTR(name), STR_LEN(name), c_name, ADF_NAME_LENGTH, ier);
    if (*ier != NO_ERROR) return;
    StringFtoC(STR_PTR(file), STR_LEN(file), c_file, ADF_FILENAME_LENGTH, ier);
    if (*ier != NO_ERROR) return;
    StringFtoC(STR_PTR(nfile), STR_LEN(nfile), c_nfile,
        ADF_MAX_LINK_DATA_SIZE, ier);
    if (*ier == NO_ERROR)
        ADFH_Link(*PID, c_name, c_file, c_nfile, ID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adflver, ADFLVER) (STR_PSTR(version), int *ier STR_PLEN(version)) {
    char c_version[ADF_VERSION_LENGTH+1];

    ADFH_Library_Version(c_version, ier);
    if (*ier == NO_ERROR)
        StringCtoF(c_version, STR_PTR(version), STR_LEN(version), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfmove, ADFMOVE) (double *PID, double *ID, double *NPID, int *ier) {

    ADFH_Move_Child(*PID, *ID, *NPID, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfncld, ADFNCLD) (double *ID, int *numcld, int *ier) {

    ADFH_Number_of_Children(*ID, numcld, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfpdim, ADFPDIM) (double *ID, STR_PSTR(dtype), int *dims,
    int *dvals, int *ier STR_PLEN(dtype)) {
    char c_type[ADF_DATA_TYPE_LENGTH+1];

    StringFtoC(STR_PTR(dtype), STR_LEN(dtype), c_type,
        ADF_DATA_TYPE_LENGTH, ier);
    if (*ier == NO_ERROR)
        ADFH_Put_Dimension_Information(*ID, c_type, *dims, dvals, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfpnam, ADFPNAM) (double *PID, double *ID, STR_PSTR(name),
    int *ier STR_PLEN(name)) {
    char c_name[ADF_NAME_LENGTH+1];

    StringFtoC(STR_PTR(name), STR_LEN(name), c_name, ADF_NAME_LENGTH, ier);
    if (*ier == NO_ERROR)
        ADFH_Put_Name(*PID, *ID, c_name, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfses, ADFSES) (int *estate, int *ier) {

    ADFH_Set_Error_State(*estate, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfslb, ADFSLB) (double *ID, STR_PSTR(label),
    int *ier STR_PLEN(label)) {
    char c_label[ADF_LABEL_LENGTH+1];

    StringFtoC(STR_PTR(label), STR_LEN(label), c_label, ADF_LABEL_LENGTH, ier);
    if (*ier == NO_ERROR)
        ADFH_Set_Label(*ID, c_label, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfread, ADFREAD) (double *ID, int *s_start, int *s_end,
    int *s_stride, int *m_num_dims, int *m_dims, int *m_start, int *m_end,
    int *m_stride, char *data, int *ier) {

    ADFH_Read_Data(*ID, s_start, s_end, s_stride, *m_num_dims,
        m_dims, m_start, m_end, m_stride, data, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfrall, ADFRALL) (double *ID, char *data, int *ier) {

    ADFH_Read_All_Data(*ID, data, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfrblk, ADFRBLK) (double *ID, int *b_start, int *b_end,
    char *data, int *ier) {

    ADFH_Read_Block_Data(*ID, *b_start, *b_end, data, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfwrit, ADFWRIT) (double *ID, int *s_start, int *s_end,
    int *s_stride, int *m_num_dims, int *m_dims, int *m_start, int *m_end,
    int *m_stride, char *data, int *ier) {

    ADFH_Write_Data(*ID, s_start, s_end, s_stride, *m_num_dims,
        m_dims, m_start, m_end, m_stride, data, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfwall, ADFWALL) (double *ID, char *data, int *ier) {

    ADFH_Write_All_Data(*ID, data, ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfwblk, ADFRWBLK) (double *ID, int *b_start, int *b_end,
    char *data, int *ier) {

    ADFH_Write_Block_Data(*ID, *b_start, *b_end, data, ier);
}

/*----- added mainly to handle character output under windows -----*/

CGNSDLL void FMNAME(adfreadc, ADFREADC) (double *ID, int *s_start, int *s_end,
    int *s_stride, int *m_num_dims, int *m_dims, int *m_start, int *m_end,
    int *m_stride, STR_PSTR(data), int *ier STR_PLEN(data)) {

    ADFH_Read_Data(*ID, s_start, s_end, s_stride, *m_num_dims,
        m_dims, m_start, m_end, m_stride, STR_PTR(data), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfrallc, ADFRALLC) (double *ID, STR_PSTR(data),
    int *ier STR_PLEN(data)) {

    ADFH_Read_All_Data(*ID, STR_PTR(data), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfrblkc, ADFRBLKC) (double *ID, int *b_start, int *b_end,
    STR_PSTR(data), int *ier STR_PLEN(data)) {

    ADFH_Read_Block_Data(*ID, *b_start, *b_end, STR_PTR(data), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfwritc, ADFWRITC) (double *ID, int *s_start, int *s_end,
    int *s_stride, int *m_num_dims, int *m_dims, int *m_start, int *m_end,
    int *m_stride, STR_PSTR(data), int *ier STR_PLEN(data)) {

    ADFH_Write_Data(*ID, s_start, s_end, s_stride, *m_num_dims,
        m_dims, m_start, m_end, m_stride, STR_PTR(data), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfwallc, ADFWALLC) (double *ID, STR_PSTR(data),
    int *ier STR_PLEN(data)) {

    ADFH_Write_All_Data(*ID, STR_PTR(data), ier);
}

/*-------------------------------------------------------------------*/

CGNSDLL void FMNAME(adfwblkc, ADFRBWLKC) (double *ID, int *b_start, int *b_end,
    STR_PSTR(data), int *ier STR_PLEN(data)) {

    ADFH_Write_Block_Data(*ID, *b_start, *b_end, STR_PTR(data), ier);
}

