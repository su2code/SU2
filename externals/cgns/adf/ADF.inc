C
C      File:	ADF.inc
C      ----------------------------------------------------------------------
C                                BOEING
C      ----------------------------------------------------------------------
C        Project: CGNS
C        Authors: Tom Dickens   865-6122    tpd6908@yak.ca.boeing.com
C                 Chuck Keagle  234-0990    keagle@cfdd38.cfd.ca.boeing.com
C        Date: 2/7/96
C        Purpose: Provide FORTRAN definitions for the ADF-Core routines.
C     ----------------------------------------------------------------------
C     ----------------------------------------------------------------------

      implicit none

C***********************************************************************
C	Defines:  These defines are used within the ADF core routines
C	to specify the size in bytes of various items.
C  Caution:  Simply changing a define here may not correctly adjust the
C	ADF core code.  These sizes are provided for reference only!
C***********************************************************************/
      integer ADF_DATA_TYPE_LENGTH
      integer ADF_DATE_LENGTH
      integer ADF_FILENAME_LENGTH
      integer ADF_FORMAT_LENGTH
      integer ADF_LABEL_LENGTH
      integer ADF_MAXIMUM_LINK_DEPTH
      integer ADF_MAX_DIMENSIONS
      integer ADF_MAX_ERROR_STR_LENGTH
      integer ADF_MAX_LINK_DATA_SIZE
      integer ADF_NAME_LENGTH
      integer ADF_STATUS_LENGTH
      integer ADF_VERSION_LENGTH

      parameter (ADF_DATA_TYPE_LENGTH = 32)
      parameter (ADF_DATE_LENGTH = 32 )
      parameter (ADF_FILENAME_LENGTH = 1024)
      parameter (ADF_FORMAT_LENGTH = 20)
      parameter (ADF_LABEL_LENGTH = 32)
      parameter (ADF_MAXIMUM_LINK_DEPTH = 100)
      parameter (ADF_MAX_DIMENSIONS = 12)
      parameter (ADF_MAX_ERROR_STR_LENGTH = 80)
      parameter (ADF_MAX_LINK_DATA_SIZE = 4096)
      parameter (ADF_NAME_LENGTH = 32)
      parameter (ADF_STATUS_LENGTH = 32)
      parameter (ADF_VERSION_LENGTH = 32 )

!***********************************************************************
!   Error-return values
!   These values need to be kept in sync with the error strings in
!   file ADF_interface.c
!***********************************************************************
!      Don't use zero since you can assign zero to a pointer

      integer NO_ERROR                      
      integer NUMBER_LESS_THAN_MINIMUM      
      integer NUMBER_GREATER_THAN_MAXIMUM   
      integer STRING_LENGTH_ZERO            
      integer STRING_LENGTH_TOO_BIG         
      integer STRING_NOT_A_HEX_STRING       
      integer TOO_MANY_ADF_FILES_OPENED     
      integer ADF_FILE_STATUS_NOT_RECOGNIZED
      integer FILE_OPEN_ERROR               
      integer ADF_FILE_NOT_OPENED           
      integer FILE_INDEX_OUT_OF_RANGE     
      integer BLOCK_OFFSET_OUT_OF_RANGE     
      integer NULL_STRING_POINTER           
      integer FSEEK_ERROR                
      integer FWRITE_ERROR               
      integer FREAD_ERROR                
      integer ADF_MEMORY_TAG_ERROR       
      integer ADF_DISK_TAG_ERROR         
      integer REQUESTED_NEW_FILE_EXISTS  
      integer ADF_FILE_FORMAT_NOT_RECOGNIZED
      integer FREE_OF_ROOT_NODE             
      integer FREE_OF_FREE_CHUNK_TABLE      
      integer REQUESTED_OLD_FILE_NOT_FOUND  
      integer UNIMPLEMENTED_CODE            
      integer SUB_NODE_TABLE_ENTRIES_BAD    
      integer MEMORY_ALLOCATION_FAILED      
      integer DUPLICATE_CHILD_NAME          
      integer ZERO_DIMENSIONS               
      integer BAD_NUMBER_OF_DIMENSIONS      
      integer CHILD_NOT_OF_GIVEN_PARENT     
      integer DATA_TYPE_TOO_LONG            
      integer INVALID_DATA_TYPE             
      integer NULL_POINTER                  
      integer NO_DATA                       
      integer ERROR_ZEROING_OUT_MEMORY      
      integer REQUESTED_DATA_TOO_LONG       
      integer END_OUT_OF_DEFINED_RANGE      
      integer BAD_STRIDE_VALUE              
      integer MINIMUM_GT_MAXIMUM            
      integer MACHINE_FORMAT_NOT_RECOGNIZED 
      integer CANNOT_CONVERT_NATIVE_FORMAT  
      integer CONVERSION_FORMATS_EQUAL      
      integer DATA_TYPE_NOT_SUPPORTED       
      integer FILE_CLOSE_ERROR              
      integer NUMERIC_OVERFLOW              
      integer START_OUT_OF_DEFINED_RANGE    
      integer ZERO_LENGTH_VALUE             
      integer BAD_DIMENSION_VALUE           
      integer BAD_ERROR_STATE               
      integer UNEQUAL_MEMORY_AND_DISK_DIMS  
      integer LINKS_TOO_DEEP                
      integer NODE_IS_NOT_A_LINK            
      integer LINK_TARGET_NOT_THERE         
      integer LINKED_TO_FILE_NOT_THERE      
      integer NODE_ID_ZERO                  
      integer INCOMPLETE_DATA               
      integer INVALID_NODE_NAME             
      integer INVALID_VERSION               
      integer NODES_NOT_IN_SAME_FILE        
      integer PRISTK_NOT_FOUND              
      integer MACHINE_FILE_INCOMPATABLE     
      integer FFLUSH_ERROR
      integer NULL_NODEID_POINTER
      integer MAX_FILE_SIZE_EXCEEDED

      parameter (NO_ERROR                      = -1)
      parameter (NUMBER_LESS_THAN_MINIMUM      =  1)
      parameter (NUMBER_GREATER_THAN_MAXIMUM   =  2)
      parameter (STRING_LENGTH_ZERO            =  3)
      parameter (STRING_LENGTH_TOO_BIG         =  4)
      parameter (STRING_NOT_A_HEX_STRING       =  5)
      parameter (TOO_MANY_ADF_FILES_OPENED     =  6)
      parameter (ADF_FILE_STATUS_NOT_RECOGNIZED=  7)
      parameter (FILE_OPEN_ERROR               =  8)
      parameter (ADF_FILE_NOT_OPENED           =  9)
      parameter (FILE_INDEX_OUT_OF_RANGE       = 10)
      parameter (BLOCK_OFFSET_OUT_OF_RANGE     = 11)
      parameter (NULL_STRING_POINTER           = 12)
      parameter (FSEEK_ERROR                   = 13)
      parameter (FWRITE_ERROR                  = 14)
      parameter (FREAD_ERROR                   = 15)
      parameter (ADF_MEMORY_TAG_ERROR          = 16)
      parameter (ADF_DISK_TAG_ERROR            = 17)
      parameter (REQUESTED_NEW_FILE_EXISTS     = 18)
      parameter (ADF_FILE_FORMAT_NOT_RECOGNIZED= 19)
      parameter (FREE_OF_ROOT_NODE             = 20)
      parameter (FREE_OF_FREE_CHUNK_TABLE      = 21)
      parameter (REQUESTED_OLD_FILE_NOT_FOUND  = 22)
      parameter (UNIMPLEMENTED_CODE            = 23)
      parameter (SUB_NODE_TABLE_ENTRIES_BAD    = 24)
      parameter (MEMORY_ALLOCATION_FAILED      = 25)
      parameter (DUPLICATE_CHILD_NAME          = 26)
      parameter (ZERO_DIMENSIONS               = 27)
      parameter (BAD_NUMBER_OF_DIMENSIONS      = 28)
      parameter (CHILD_NOT_OF_GIVEN_PARENT     = 29)
      parameter (DATA_TYPE_TOO_LONG            = 30)
      parameter (INVALID_DATA_TYPE             = 31)
      parameter (NULL_POINTER                  = 32)
      parameter (NO_DATA                       = 33)
      parameter (ERROR_ZEROING_OUT_MEMORY      = 34)
      parameter (REQUESTED_DATA_TOO_LONG       = 35)
      parameter (END_OUT_OF_DEFINED_RANGE      = 36)
      parameter (BAD_STRIDE_VALUE              = 37)
      parameter (MINIMUM_GT_MAXIMUM            = 38)
      parameter (MACHINE_FORMAT_NOT_RECOGNIZED = 39)
      parameter (CANNOT_CONVERT_NATIVE_FORMAT  = 40)
      parameter (CONVERSION_FORMATS_EQUAL      = 41)
      parameter (DATA_TYPE_NOT_SUPPORTED       = 42)
      parameter (FILE_CLOSE_ERROR              = 43)
      parameter (NUMERIC_OVERFLOW              = 44)
      parameter (START_OUT_OF_DEFINED_RANGE    = 45)
      parameter (ZERO_LENGTH_VALUE             = 46)
      parameter (BAD_DIMENSION_VALUE           = 47)
      parameter (BAD_ERROR_STATE               = 48)
      parameter (UNEQUAL_MEMORY_AND_DISK_DIMS  = 49)
      parameter (LINKS_TOO_DEEP                = 50)
      parameter (NODE_IS_NOT_A_LINK            = 51)
      parameter (LINK_TARGET_NOT_THERE         = 52)
      parameter (LINKED_TO_FILE_NOT_THERE      = 53)
      parameter (NODE_ID_ZERO                  = 54)
      parameter (INCOMPLETE_DATA               = 55)
      parameter (INVALID_NODE_NAME             = 56)
      parameter (INVALID_VERSION               = 57)
      parameter (NODES_NOT_IN_SAME_FILE        = 58)
      parameter (PRISTK_NOT_FOUND              = 59)
      parameter (MACHINE_FILE_INCOMPATABLE     = 60)
      parameter (FFLUSH_ERROR                  = 61)
      parameter (NULL_NODEID_POINTER           = 62)
      parameter (MAX_FILE_SIZE_EXCEEDED        = 63)
