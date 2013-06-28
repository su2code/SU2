/* created by combine 2.0 */
/* file ADFAAA2.c */
/***********************************************************************
   ADF Core:
   Glue routines between the FORTRAN interface and the C interface.

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ADF.h"
#include "ADF_internals.h"
#ifdef MEM_DEBUG
#include "cg_malloc.h"
#endif

/* end of file ADFAAA2.c */
/* file ADFCNA2.c */
/***********************************************************************
ADFCNAM   ADF_Children_Names:
***********************************************************************/
void    FNAME(adfcna2,ADFCNA2)(
                const Fdouble *PID,
                const Fint *istart,
                const Fint *imaxnum,
                const Fint *idim,
                const Fint *name_length,
                Fint *inum_ret,
                Fchar names,
                Fint *error_return )
{
int i;
char *pstr, *tmp_names;

if( *idim != *name_length ) {  /** inconsistency **/
   *error_return = BAD_DIMENSION_VALUE ;  /** not really what this error code
                                              was meant for but it might do **/
   return ;
   } /* end if */

pstr = F2CP(names);

tmp_names = (char *) malloc( (*imaxnum) * (*name_length + 1) * (sizeof(char)) ) ;
if( tmp_names == NULL ) {
   *error_return = MEMORY_ALLOCATION_FAILED ;
   return ;
   } /* end if */

ADF_Children_Names( *PID, *istart, *imaxnum, *name_length, inum_ret,
   tmp_names, error_return ) ;
if( *error_return != NO_ERROR ) {
   free( tmp_names ) ;
   return ;
   } /* end if */

for( i=0; i<*inum_ret; i++ )   { /* Blank-fill the names */
   if( strlen( &tmp_names[ i * (*name_length+1) ]) == *name_length ) {
/** string is maximum length, do not null terminate or blank fill **/
      strncpy( &pstr[ i * (*name_length) ], &tmp_names[ i * (*name_length+1) ],
               *name_length ) ;
      }
   else {
/** string is short enough, blank fill remainder **/
      strcpy( &pstr[ i * (*name_length) ], &tmp_names[ i * (*name_length+1) ] ) ;
      ADFI_blank_fill_string( &pstr[ i * (*name_length) ], *name_length ) ;
      } /* end if */
   } /* end for */

free( tmp_names ) ;
}
/* end of file ADFCNA2.c */
/* file ADFCID2.c */
/***********************************************************************
ADFCNAM   ADF_Children_IDs:
***********************************************************************/
void    FNAME(adfcid2,ADFCID2)(
                const Fdouble *PID,
                const Fint *istart,
                const Fint *imaxnum,
                Fint *inum_ret,
                Fdouble *CIDs,
                Fint *error_return )
{
ADF_Children_IDs( *PID, *istart, *imaxnum, inum_ret, CIDs, error_return ) ;
}
/* end of file ADFCID2.c */
/* file ADFCRE2.c */
/***********************************************************************
ADFCRE    ADF_Create:
***********************************************************************/
void    FNAME(adfcre2,ADFCRE2)(
                const Fdouble *PID,
                const Fchar name,
		const Fint *name_length,
                Fdouble *ID,
                Fint *error_return )
{
char	c_name[ ADF_NAME_LENGTH + 1 ] ;

ADFI_string_2_C_string( F2CP(name), MIN(ADF_NAME_LENGTH, *name_length), c_name,
	error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADF_Create( *PID, c_name, ID, error_return ) ;
}
/* end of file ADFCRE2.c */
/* file ADFDCL2.c */
/***********************************************************************
ADFDCLO   ADF_Database_Close:
***********************************************************************/
void    FNAME(adfdcl2,ADFDCL2)(
                const Fdouble *Root_ID,
                Fint *error_return )
{
ADF_Database_Close( *Root_ID, error_return ) ;
}
/* end of file ADFDCL2.c */
/* file ADFDDE2.c */
/***********************************************************************
ADFDDEL   ADF_Database_Delete:
***********************************************************************/
void    FNAME(adfdde2,ADFDDE2)(
                const Fchar filename,
		const Fint *name_length,
                Fint *error_return )
{
char	c_name[ ADF_FILENAME_LENGTH + 1 ] ;

ADFI_string_2_C_string( F2CP(filename), MIN(ADF_FILENAME_LENGTH, *name_length),
	c_name, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
ADF_Database_Delete( c_name, error_return ) ;
}
/* end of file ADFDDE2.c */
/* file ADFDEL2.c */
/***********************************************************************
ADFDEL    ADF_Delete:
***********************************************************************/
void    FNAME(adfdel2,ADFDEL2)(
                const Fdouble *PID,
                const Fdouble *ID,
                Fint *error_return )
{
ADF_Delete( *PID, *ID, error_return ) ;
}
/* end of file ADFDEL2.c */
/* file ADFDGC2.c */
/***********************************************************************
ADFDGC    ADF_Database_Garbage_Collection:
***********************************************************************/
void    FNAME(adfdgc2,ADFDGC2)(
                const Fdouble *ID,
                Fint *error_return )
{
ADF_Database_Garbage_Collection( *ID, error_return ) ;
}
/* end of file ADFDGC2.c */
/* file ADFDGF2.c */
/***********************************************************************
ADFDGF    ADF_Database_Get_Format:
***********************************************************************/
void    FNAME(adfdgf2,ADFDGF2)(
                const Fdouble *Root_ID,
                Fchar format,
                const Fint *format_length,
                Fint *error_return )
{
ADF_Database_Get_Format( *Root_ID, F2CP(format), error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADFI_blank_fill_string( F2CP(format), *format_length ) ;
}
/* end of file ADFDGF2.c */
/* file ADFDOP2.c */
/***********************************************************************
ADFDOPN   ADF_Database_Open:
***********************************************************************/
void    FNAME(adfdop2,ADFDOP2)(
                const Fchar filename,
		const Fint *filename_length,
                Fchar status_in,
		const Fint *status_length,
                const Fchar format,
		const Fint *format_length,
                Fdouble *Root_ID,
                Fint *error_return )
{
char	c_filename[ ADF_FILENAME_LENGTH + 1 ],
	c_status[ ADF_NAME_LENGTH+1 ],
	c_format[ ADF_NAME_LENGTH+1 ] ;

ADFI_string_2_C_string( F2CP(filename),
        MIN(ADF_FILENAME_LENGTH, *filename_length),
	c_filename, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
ADFI_string_2_C_string( F2CP(status_in), 
        MIN(ADF_NAME_LENGTH, *status_length),
	c_status, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
ADFI_string_2_C_string( F2CP(format), 
        MIN(ADF_NAME_LENGTH, *format_length), c_format,
	error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADF_Database_Open( c_filename, c_status, c_format,
                Root_ID, error_return ) ;
}
/* end of file ADFDOP2.c */
/* file ADFDSF2.c */
/***********************************************************************
ADFDSF    ADF_Database_Set_Format:
***********************************************************************/
void    FNAME(adfdsf2,ADFDSF2)(
                const Fdouble *Root_ID,
                const Fchar format,
		const Fint *format_length,
                Fint *error_return )
{
char	c_format[ ADF_NAME_LENGTH + 1 ] ;

ADFI_string_2_C_string( F2CP(format), 
        MIN(ADF_NAME_LENGTH, *format_length), c_format,
	error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADF_Database_Set_Format( *Root_ID, c_format, error_return ) ;
}
/* end of file ADFDSF2.c */
/* file ADFDVE2.c */
/***********************************************************************
ADVDVER   ADF_Database_Version:
***********************************************************************/
void    FNAME(adfdve2,ADFDVE2)(
                const Fdouble *Root_ID,
                Fchar version,
                Fchar creation_date,
                Fchar modification_date,
		const Fint *v_length,
		const Fint *c_length,
		const Fint *m_length,
                Fint *error_return )
{
ADF_Database_Version( *Root_ID, F2CP(version), F2CP(creation_date), 
                F2CP(modification_date), error_return ) ;

ADFI_blank_fill_string( F2CP(version), *v_length ) ;
ADFI_blank_fill_string( F2CP(creation_date), *c_length ) ;
ADFI_blank_fill_string( F2CP(modification_date), *m_length ) ;
}
/* end of file ADFDVE2.c */
/* file ADFERR2.c */
/***********************************************************************
ADFERR    ADF_Error_Message:
***********************************************************************/
void    FNAME(adferr2,ADFERR2)(
                const Fint *error_return_input,
                Fchar error_string,
                const Fint *str_length )
{
char  msg_buf[ADF_MAX_ERROR_STR_LENGTH+1] ;

ADF_Error_Message( *error_return_input, msg_buf ) ;
strncpy( F2CP(error_string), msg_buf, *str_length ) ;
ADFI_blank_fill_string( F2CP(error_string), *str_length ) ;
}
/* end of file ADFERR2.c */
/* file ADFFTD2.c */
/***********************************************************************
ADFFTD    ADF_Flush_to_Disk:
***********************************************************************/
void    FNAME(adfftd2,ADFFTD2)(
                const Fdouble *ID,
                Fint *error_return )
{
ADF_Flush_to_Disk( *ID, error_return ) ;
}
/* end of file ADFFTD2.c */
/* file ADFGDT2.c */
/***********************************************************************
ADFGDT    ADF_Get_Data_Type:
***********************************************************************/
void    FNAME(adfgdt2,ADFGDT2)(
                const Fdouble *ID,
                Fchar data_type,
                const Fint *data_type_length,
                Fint *error_return )
{
char  ctype[ ADF_DATA_TYPE_LENGTH + 1 ] ;

ADF_Get_Data_Type( *ID, ctype, error_return ) ;
if( *error_return == NO_ERROR ) {
   if( strlen( ctype ) < *data_type_length ) {
      strcpy( F2CP(data_type), ctype ) ;
      ADFI_blank_fill_string( F2CP(data_type), *data_type_length ) ;
      }
   else {
      strncpy( F2CP(data_type), ctype, *data_type_length ) ;
      } /* end if */
   } /* end if */
}
/* end of file ADFGDT2.c */
/* file ADFGDV2.c */
/***********************************************************************
ADFGDV    ADF_Get_Dimension_Values:
***********************************************************************/
void    FNAME(adfgdv2,ADFGDV2)(
                const Fdouble *ID,
                Fint dim_vals[],
                Fint *error_return )
{
ADF_Get_Dimension_Values( *ID, dim_vals, error_return ) ;
}
/* end of file ADFGDV2.c */
/* file ADFGES2.c */
/***********************************************************************
ADFGES    ADF_Get_Error_State:
***********************************************************************/
void    FNAME(adfges2,ADFGES2)(
                Fint *error_state,
                Fint *error_return )
{
ADF_Get_Error_State( error_state, error_return ) ;
}
/* end of file ADFGES2.c */
/* file ADFGLB2.c */
/***********************************************************************
ADFGLB    ADF_Get_Label:
***********************************************************************/
void    FNAME(adfglb2,ADFGLB2)(
                const Fdouble *ID,
                Fchar label,
                const Fint *label_length,
                Fint *error_return )
{
char  clabel[ ADF_LABEL_LENGTH + 1 ] ;

ADF_Get_Label( *ID, clabel, error_return ) ;
if( *error_return == NO_ERROR ) {
   if( strlen( clabel ) < *label_length ) {
      strcpy( F2CP(label), clabel ) ;
      ADFI_blank_fill_string( F2CP(label), *label_length ) ;
      }
   else {
      strncpy( F2CP(label), clabel, *label_length ) ;
      } /* end if */
   } /* end if */
}
/* end of file ADFGLB2.c */
/* file ADFGLK2.c */
/***********************************************************************
ADFGLKP   ADF_Get_Link_Path:
***********************************************************************/
void    FNAME(adfglk2,ADFGLK2)(
                const Fdouble *ID,
                Fchar filename,
                const Fint *filename_length,
                Fchar link_path,
                const Fint *link_path_length,
                Fint *error_return )
{
char  cpath[ ADF_MAX_LINK_DATA_SIZE + 1 ],
      cfilename[ ADF_FILENAME_LENGTH + 1 ] ;

ADF_Get_Link_Path( *ID, cfilename, cpath, error_return ) ;
if( *error_return == NO_ERROR ) {
   if( strlen(cfilename) < *filename_length ) {
      strcpy( F2CP(filename), cfilename ) ;
      ADFI_blank_fill_string( F2CP(filename), *filename_length ) ;
      }
   else {
      strncpy( F2CP(filename), cfilename, *filename_length ) ;
      } /* end if */

   if( strlen(cpath) < *link_path_length ) {
      strcpy( F2CP(link_path), cpath ) ;
      ADFI_blank_fill_string( F2CP(link_path), *link_path_length ) ;
      }
   else {
      strncpy( F2CP(link_path), cpath, *link_path_length ) ;
      } /* end if */
   } /* end if */
}
/* end of file ADFGLK2.c */
/* file ADFGNA2.c */
/***********************************************************************
ADFGNAM   ADF_Get_Name:
***********************************************************************/
void    FNAME(adfgna2,ADFGNA2)(
                const Fdouble *ID,
                Fchar name,
                const Fint *name_length,
                Fint *error_return )
{
char  cname[ ADF_NAME_LENGTH + 1 ] ;

ADF_Get_Name( *ID, cname, error_return ) ;
if( *error_return == NO_ERROR ) {
   if( strlen( cname ) < *name_length ) {
      strcpy( F2CP(name), cname ) ;
      ADFI_blank_fill_string( F2CP(name), *name_length ) ;
      }
   else {
      strncpy( F2CP(name), cname, *name_length ) ;
      } /* end if */
   } /* end if */
}
/* end of file ADFGNA2.c */
/* file ADFGND2.c */
/***********************************************************************
ADFGND    ADF_Get_Number_of_Dimensions:
***********************************************************************/
void    FNAME(adfgnd2,ADFGND2)(
                const Fdouble *ID,
                Fint *num_dims,
                Fint *error_return )
{
ADF_Get_Number_of_Dimensions( *ID, num_dims, error_return ) ;
}
/* end of file ADFGND2.c */
/* file ADFGNI2.c */
/***********************************************************************
ADFGNID   ADF_Get_Node_ID:
***********************************************************************/
void    FNAME(adfgni2,ADFGNI2)(
                const Fdouble *PID,
                const Fchar name,
		const Fint *name_length,
                Fdouble *ID,
                Fint *error_return )
{
char	c_name[ ADF_FILENAME_LENGTH + 1 ] ;

ADFI_string_2_C_string( F2CP(name), MIN(ADF_FILENAME_LENGTH, *name_length),
	c_name, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
ADF_Get_Node_ID( *PID, c_name, ID, error_return ) ;
}
/* end of file ADFGNI2.c */
/* file ADFGRI2.c */
/***********************************************************************
ADFGRID   ADF_Get_Root_ID:
***********************************************************************/
void    FNAME(adfgri2,ADFGRI2)(
                const Fdouble *ID,
                Fdouble *Root_ID,
                Fint *error_return )
{
ADF_Get_Root_ID( *ID, Root_ID, error_return ) ;
}
/* end of file ADFGRI2.c */
/* file ADFISL2.c */
/***********************************************************************
ADFISLK   ADF_Is_Link:
***********************************************************************/
void    FNAME(adfisl2,ADFISL2)(
                const Fdouble *ID,
                Fint *link_path_length,
                Fint *error_return )
{
ADF_Is_Link( *ID, link_path_length, error_return ) ;
}
/* end of file ADFISL2.c */
/* file ADFLIN2.c */
/***********************************************************************
ADFLINK   ADF_Link:
***********************************************************************/
void    FNAME(adflin2,ADFLIN2)(
                const Fdouble *PID,
                const Fchar name,
                const Fchar file,
                const Fchar name_in_file,
		const Fint *name_length,
		const Fint *file_length,
		const Fint *nfile_length,
                Fdouble *ID,
                Fint *error_return )
{
char	c_name[ ADF_FILENAME_LENGTH + 1 ],
	c_file[ ADF_FILENAME_LENGTH + 1 ],
	c_nfile[ ADF_MAX_LINK_DATA_SIZE + 1 ] ;

ADFI_string_2_C_string( F2CP(name), 
                MIN(ADF_FILENAME_LENGTH, *name_length),
		c_name, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
ADFI_string_2_C_string( F2CP(file), 
                MIN(ADF_FILENAME_LENGTH, *file_length),
		c_file, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
ADFI_string_2_C_string( F2CP(name_in_file), 
                MIN(ADF_MAX_LINK_DATA_SIZE, *nfile_length),
		c_nfile,
	error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADF_Link( *PID, c_name, c_file, c_nfile, ID, error_return ) ;
}
/* end of file ADFLIN2.c */
/* file ADFLVE2.c */
/***********************************************************************
ADFLVER   ADF_Library_Version:
***********************************************************************/
void    FNAME(adflve2,ADFLVE2)(
                Fchar version,
		const Fint *version_length,
                Fint *error_return )
{
ADF_Library_Version( F2CP(version), error_return ) ;
ADFI_blank_fill_string ( F2CP(version), *version_length );
}
/* end of file ADFLVE2.c */
/* file ADFMOV2.c */
/***********************************************************************
ADFMOVE   ADF_Move_Child:
***********************************************************************/
void    FNAME(adfmov2,ADFMOV2)(
                const Fdouble *PID,
                const Fdouble *ID,
                const Fdouble *NPID,
                Fint *error_return )
{
ADF_Move_Child( *PID, *ID, *NPID, error_return ) ;
}
/* end of file ADFMOV2.c */
/* file ADFNCL2.c */
/***********************************************************************
ADFNCLD   ADF_Number_of_Children:
***********************************************************************/
void    FNAME(adfncl2,ADFNCL2)(
                const Fdouble *ID,
                Fint *num_children,
                Fint *error_return )
{
ADF_Number_of_Children( *ID, num_children, error_return ) ;
}
/* end of file ADFNCL2.c */
/* file ADFPDI2.c */
/***********************************************************************
ADFPDIM   ADF_Put_Dimension_Information:
***********************************************************************/
void    FNAME(adfpdi2,ADFPDI2)(
                const Fdouble *ID,
                const Fchar data_type,
                const Fint *data_type_length,
                const Fint *dims,
                const Fint dim_vals[],
                Fint *error_return )
{
char    c_data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;

/* CMLU */
int i;
for (i=0;i<ADF_DATA_TYPE_LENGTH;i++)
   c_data_type[i] = ' ';
c_data_type[i] = '\0';


ADFI_string_2_C_string( F2CP(data_type), 
        MIN(ADF_DATA_TYPE_LENGTH, *data_type_length),
        c_data_type, error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADF_Put_Dimension_Information( *ID, c_data_type, *dims, dim_vals,
                error_return ) ;
}
/* end of file ADFPDI2.c */
/* file ADFPNA2.c */
/***********************************************************************
ADFPNAM   ADF_Put_Name:
***********************************************************************/
void    FNAME(adfpna2,ADFPNA2)(
                const Fdouble *PID,
                const Fdouble *ID,
                const Fchar name,
                const Fint *name_length,
                Fint *error_return )
{
char    c_name[ ADF_NAME_LENGTH + 1 ] ;

ADFI_string_2_C_string( F2CP(name), 
                MIN(ADF_NAME_LENGTH, *name_length), c_name,
                error_return ) ;
if( *error_return != NO_ERROR )
   return ;

ADF_Put_Name( *PID, *ID, c_name, error_return ) ;
}
/* end of file ADFPNA2.c */
/* file ADFRAL2.c */
/* end of file ADFRAL2.c */
/* file ADFRALL.c */
/***********************************************************************
ADFRALL   ADF_Read_All_Data:
   Read all data from a Node.  Reads all the node's data and returns
   it into a contiguous memory space.

   input:   real*8          ID             The ID of the node to use.
   output:  character *(*)  data           The start of the data in memory.
   output:  integer         error_return   Error flag.
***********************************************************************/
void    FNAME(adfrall,ADFRALL)(
                const Fdouble *ID,
                Fchar data,
                Fint *error_return )
{

#if  defined(cray) && defined(_ADDR64)

int  local_error ;
char data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;
char errmsg[ ADF_MAX_ERROR_STR_LENGTH + 1 ] ;

/** see ADFWALL() for more details **/

   ADF_Get_Data_Type( *ID, data_type, &local_error ) ;
   if( local_error != NO_ERROR ) {
      ADF_Error_Message( local_error, errmsg ) ;
      printf( "%s\n", errmsg ) ;
      printf( "Unrecoverable ADF error.  ADFRALL\n" ) ;
      printf( "Cannot determine data type, so cannot determine function\n" ) ;
      printf( "argument list (character arrays are different than other\n" ) ;
      printf( "types in this environemnt), so cannot set error_return.\n" ) ;
      abort () ;
      } /* end if */

   if ( ADFI_stridx_c( data_type, "C1" ) >= 0 ) {
/** expecting character data type **/

      ADF_Read_All_Data( *ID, F2CP(data), error_return ) ;
      }
   else {
/** expecting non-character data type **/

      ADF_Read_All_Data( *ID, data.c_pointer, (int *)data.fcd_len ) ;
      }

#else

   ADF_Read_All_Data( *ID, F2CP(data), error_return ) ;

#endif
}
/* end of file ADFRALL.c */
/* file ADFRBLK.c */
/***********************************************************************
ADFRBLK   ADF_Read_Block_Data:
   Read block of data from a Node.  Reads a block of the node's data and
   returns it into a contiguous memory space.

   input:   real*8  ID               The ID of the node to use.
   input:   const int b_start	     The starting point in block in token space
   input:   const int b_end          The ending point in block in token space
   output:  character *(*) data      The start of the data in memory.
   output:  integer   error_return   Error flag.
***********************************************************************/
void    FNAME(adfrblk,ADFRBLK)(
                const Fdouble *ID,
                const Fint *b_start,
                const Fint *b_end,
                Fchar data,
                Fint *error_return )
{

#if  defined(cray) && defined(_ADDR64)

int  local_error ;
char data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;
char errmsg[ ADF_MAX_ERROR_STR_LENGTH + 1 ] ;

/** see ADFWALL() for more details **/

   ADF_Get_Data_Type( *ID, data_type, &local_error ) ;
   if( local_error != NO_ERROR ) {
      ADF_Error_Message( local_error, errmsg ) ;
      fprintf(stderr,"%s\n", errmsg ) ;
      fprintf(stderr,"Unrecoverable ADF error.  ADFRBLK\n" ) ;
      fprintf(stderr,"Cannot determine data type, so cannot determine function\n" ) ;
      fprintf(stderr,"argument list (character arrays are different than other\n" ) ;
      fprintf(stderr,"types in this environemnt), so cannot set error_return.\n" ) ;
      abort () ;
      } /* end if */

   if ( ADFI_stridx_c( data_type, "C1" ) >= 0 ) {
/** expecting character data type **/

      ADF_Read_Block_Data( *ID, (long) *b_start, (long) *b_end,
               F2CP(data), error_return ) ;
      }
   else {
/** expecting non-character data type **/

      ADF_Read_Block_Data( *ID, (long) *b_start, (long) *b_end,
               data.c_pointer, (int *)data.fcd_len ) ;
      }

#else

   ADF_Read_Block_Data( *ID, (long) *b_start, (long) *b_end,
               F2CP(data), error_return ) ;

#endif
}
/* end of file ADFRBLK.c */
/* file ADFREA2.c */
/* end of file ADFREA2.c */
/* file ADFREAD.c */
/***********************************************************************
ADFREAD   ADF_Read_Data:

   A 1-based system is used with all index values  (the first element has
   an index of 1, not 0).
   R1: Read data from a node, with partial capabilities.  The partial 
   capabilities are both in the node's data and also in memory.  
   Vectors of integers are used to indicate the data to be accessed 
   from the node, and another set of integer vectors is used to 
   describe the memory location for the data.  
    Note:  If the data-type of the node is a compound data-type ("I4[3],R8") 
   for example, the partial capabilities will access one or more of 
   these 20-byte data entities.  You cannot access a subset of an 
   occurrence of the data-type.

   f77: ADFREAD( ID, sstart[], send[], sstrid[], mnumd, 
       mdims[], mstart[], mend[], mstrid[], data, ierr )
   input:  real*8  ID            The ID of the node to use.
   input:  integer sstart(12)    The starting dimension values to use 
                                 in the database (node).
   input:  integer send(12)      The ending dimension values to use in 
                                 the database (node).
   input:  integer sstrid(12)    The stride values to use in the 
                                 database (node).
   input:  integer mnumd         The number of dimensions to use in memory.
   input:  integer mdims(mnumd)  The dimensionality to use in memory.
   input:  integer mstart(mnumd) The starting dimension values
                                 to use in memory.
   input:  integer mend(mnumd)   The ending dimension values to
                                 use in memory.
   input:  integer mstrid(mnumd) The stride values to use 
                                 in memory.
   output:  character*(*) data   The start of the data in memory.
   output:  integer ierr

***********************************************************************/
void    FNAME(adfread,ADFREAD)(
                const Fdouble *ID,
                const Fint s_start[],
                const Fint s_end[],
                const Fint s_stride[],
                const Fint *m_num_dims,
                const Fint m_dims[],
                const Fint m_start[],
                const Fint m_end[],
                const Fint m_stride[],
                Fchar data,
                Fint *error_return )
{

#if  defined(cray) && defined(_ADDR64)

int  local_error ;
char data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;
char errmsg[ ADF_MAX_ERROR_STR_LENGTH + 1 ] ;

/** see ADFWALL() for more details **/

   ADF_Get_Data_Type( *ID, data_type, &local_error ) ;
   if( local_error != NO_ERROR ) {
      ADF_Error_Message( local_error, errmsg ) ;
      fprintf(stderr,"%s\n", errmsg ) ;
      fprintf(stderr,"Unrecoverable ADF error.  ADFREAD\n" ) ;
      fprintf(stderr,"Cannot determine data type, so cannot determine function\n" ) ;
      fprintf(stderr,"argument list (character arrays are different than other\n" ) ;
      fprintf(stderr,"types in this environemnt), so cannot set error_return.\n" ) ;
      abort () ;
      } /* end if */

   if ( ADFI_stridx_c( data_type, "C1" ) >= 0 ) {
/** expecting character data type **/

      ADF_Read_Data( *ID, s_start, s_end, s_stride, *m_num_dims, m_dims,
                     m_start, m_end, m_stride, F2CP(data), error_return ) ;
   }
   else {
/** expecting non-character data type **/

      ADF_Read_Data( *ID, s_start, s_end, s_stride, *m_num_dims, m_dims,
                     m_start, m_end, m_stride,
                     data.c_pointer, (int *)data.fcd_len ) ;
   }

#else

   ADF_Read_Data( *ID, s_start, s_end, s_stride, *m_num_dims, m_dims,
                  m_start, m_end, m_stride, F2CP(data), error_return ) ;

#endif
}
/* end of file ADFREAD.c */
/* file ADFSES2.c */
/***********************************************************************
ADFSES    ADF_Set_Error_State:
***********************************************************************/
void    FNAME(adfses2,ADFSES2)(
                const Fint *error_state,
                Fint *error_return )
{
ADF_Set_Error_State( *error_state, error_return ) ;
}
/* end of file ADFSES2.c */
/* file ADFSLB2.c */
/***********************************************************************
ADFSLB    ADF_Set_Label:
***********************************************************************/
void    FNAME(adfslb2,ADFSLB2)(
                const Fdouble *ID,
                const Fchar label,
                const Fint *label_length,
                Fint *error_return )
{
char    c_label[ ADF_LABEL_LENGTH + 1 ] ;

ADFI_string_2_C_string( F2CP(label), 
                MIN(ADF_LABEL_LENGTH, *label_length),
                c_label, error_return ) ;
if( *error_return != NO_ERROR )
   return ;
   
ADF_Set_Label( *ID, c_label, error_return ) ;
}
/* end of file ADFSLB2.c */
/* file ADFWAL2.c */
/* end of file ADFWAL2.c */
/* file ADFWALL.c */
/***********************************************************************
ADFWALL   ADF_Write_All_Data:
   Write all data to a Node.  Writes all the node's data from a
   contiguous memory space.

   input:    real*8         ID             The node's id.
   input:    character *(*) data           The start of data in memory.
   output:   integer        error_return   Error flag.
***********************************************************************/
void    FNAME(adfwall,ADFWALL)(
                const Fdouble *ID,
                const Fchar data,
                Fint *error_return )
{

#if  defined(cray) && defined(_ADDR64)

int  local_error ;
char data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;
char errmsg[ ADF_MAX_ERROR_STR_LENGTH + 1 ] ;

/** On Crays with 64 bit addresses (like Tritons), Fortran
    character arrays occupy two 64 bit words on argument stacks
    (see <fortran.h>).  Other data (INTEGER, REAL, etc) occupy one word.
    To accommodate both situations with one function, some tricks
    need to be employed... **/

/** First, find the data type that ADF is expecting and assume the
    function is being called with that kind of argument.  This
    function is defined with character data in mind and deviates
    from that if necessary. **/

   ADF_Get_Data_Type( *ID, data_type, &local_error ) ;
   if( local_error != NO_ERROR ) {
      ADF_Error_Message( local_error, errmsg ) ;
      fprintf(stderr,"%s\n", errmsg ) ;
      fprintf(stderr,"Unrecoverable ADF error.  ADFWALL\n" ) ;
      fprintf(stderr,"Cannot determine data type, so cannot determine function\n" ) ;
      fprintf(stderr,"argument list (character arrays are different than other\n" ) ;
      fprintf(stderr,"types in this environemnt), so cannot set error_return.\n" ) ;
      abort () ;
      } /* end if */

   if ( ADFI_stridx_c( data_type, "C1" ) >= 0 ) {
/** expecting character data type **/

      ADF_Write_All_Data( *ID, F2CP(data), error_return ) ;
   }
   else {

/** expecting data type other than character - the stack is not as
    long as is with character data and so "error_return" does not
    correspond to anything valid.  Since "data" is declared Fchar but
    the actual argument is not, the second half of "data" contains
    the error flag argument (the stack is assumed to be contiguous). **/

      ADF_Write_All_Data( *ID, data.c_pointer, (int *)data.fcd_len ) ;
   }

#else

/** In other CRAY environments the declared length is encoded in
    unused portions of pointers and the F2CP macro handles its
    conversion.
    In other environments, the character array declared length
    mystery argument is assumed to be on the end of the stack and
    is ignored here (F2CP macro does nothing). **/

   ADF_Write_All_Data( *ID, F2CP(data), error_return ) ;

#endif
}
/* end of file ADFWALL.c */
/* file ADFWBLK.c */
/***********************************************************************
ADFWBLK   ADF_Write_Block_Data:
   Write block of data from a Node.  Writes a block of the node's data and
   returns it into a contiguous memory space.

   input:   real*8    ID               The ID of the node to use.
   input:   const int b_start          The starting point in block in token space
   input:   const int b_end            The ending point in block in token space
   output:  character *(*)  data       The start of the data in memory.
   output:  integer   error_return     Error flag.
***********************************************************************/
void    FNAME(adfwblk,ADFWBLK)(
                const Fdouble *ID,
                const Fint *b_start,
                const Fint *b_end,
                Fchar data,
                Fint *error_return )
{

#if  defined(cray) && defined(_ADDR64)

int  local_error ;
char data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;
char errmsg[ ADF_MAX_ERROR_STR_LENGTH + 1 ] ;

/** see ADFWALL() for more details **/

   ADF_Get_Data_Type( *ID, data_type, &local_error ) ;
   if( local_error != NO_ERROR ) {
      ADF_Error_Message( local_error, errmsg ) ;
      fprintf(stderr,"%s\n", errmsg ) ;
      fprintf(stderr,"Unrecoverable ADF error.  ADFWBLK\n" ) ;
      fprintf(stderr,"Cannot determine data type, so cannot determine function\n" ) ;
      fprintf(stderr,"argument list (character arrays are different than other\n" ) ;
      fprintf(stderr,"types in this environemnt), so cannot set error_return.\n" ) ;
      abort () ;
      } /* end if */

   if ( ADFI_stridx_c( data_type, "C1" ) >= 0 ) {
/** expecting character data type **/

      ADF_Write_Block_Data( *ID, (long) *b_start, (long) *b_end,
		 	    F2CP(data), error_return ) ;
      }
   else {
/** expecting non-character data type **/

      ADF_Write_Block_Data( *ID, (long) *b_start, (long) *b_end,
                data.c_pointer, (int *)data.fcd_len ) ;
      }

#else

   ADF_Write_Block_Data( *ID, (long) *b_start, (long) *b_end,
                F2CP(data), error_return ) ;

#endif
}
/* end of file ADFWBLK.c */
/* file ADFWRI2.c */
/* end of file ADFWRI2.c */
/* file ADFWRIT.c */
/***********************************************************************
ADFWRIT   ADF_Write_Data:
   Write data to a Node, with partial capabilities.
   See ADF_Read_Data for description.

   f77: ADFWRIT( ID, sstart[], send[], sstrid[], mnumd, 
                 mdims[], mstart[], mend[], mstrid[], data, ierr )
   input:  real*8 ID             The ID of the node to use.
   input:  integer sstart(12)    The starting dimension values to use 
                                 in the database (node).
   input:  integer send(12)      The ending dimension values to use in 
                                 the database (node).
   input:  integer sstrid(12)    The stride values to use in the 
                                 database (node).
   input:  integer mnumd         The number of dimensions to use in memory.
   input:  integer mdims(mnumd)  The dimensionality to use in memory.
   input:  integer mstart(mnumd) The starting dimension values to use
                                 in memory.
   input:  integer mend(mnumd)   The ending dimension values to use in
                                  memory.
   input:  integer mstrid(mnumd) The stride values to use in memory.
   input:  character*(*) data    The start of the data in memory.
   output: integer ierr


***********************************************************************/
void    FNAME(adfwrit,ADFWRIT)(
                const Fdouble *ID,
                const Fint s_start[],
                const Fint s_end[],
                const Fint s_stride[],
                const Fint *m_num_dims,
                const Fint m_dims[],
                const Fint m_start[],
                const Fint m_end[],
                const Fint m_stride[],
                const Fchar data,
                Fint *error_return )
{

#if  defined(cray) && defined(_ADDR64)

int  local_error ;
char data_type[ ADF_DATA_TYPE_LENGTH + 1 ] ;
char errmsg[ ADF_MAX_ERROR_STR_LENGTH + 1 ] ;

/** see ADFWALL() for more details **/

   ADF_Get_Data_Type( *ID, data_type, &local_error ) ;
   if( local_error != NO_ERROR ) {
      ADF_Error_Message( local_error, errmsg ) ;
      fprintf(stderr,"%s\n", errmsg ) ;
      fprintf(stderr,"Unrecoverable ADF error.  ADFWRIT\n" ) ;
      fprintf(stderr,"Cannot determine data type, so cannot determine function\n" ) ;
      fprintf(stderr,"argument list (character arrays are different than other\n" ) ;
      fprintf(stderr,"types in this environemnt), so cannot set error_return.\n" ) ;
      abort () ;
      } /* end if */

   if ( ADFI_stridx_c( data_type, "C1" ) >= 0 ) {
/** expecting character data type **/

      ADF_Write_Data( *ID, s_start, s_end, s_stride, *m_num_dims, m_dims,
                      m_start, m_end, m_stride, F2CP(data), error_return ) ;
   }
   else {
/** expecting non-character data type **/

      ADF_Write_Data( *ID, s_start, s_end, s_stride, *m_num_dims, m_dims,
                      m_start, m_end, m_stride,
                      data.c_pointer, (int *)data.fcd_len ) ;
   }

#else
   ADF_Write_Data( *ID, s_start, s_end, s_stride, *m_num_dims, m_dims,
                   m_start, m_end, m_stride, F2CP(data), error_return ) ;

#endif
}
/* end of file ADFWRIT.c */
/* end of combine 2.0 */
