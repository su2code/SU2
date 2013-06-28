c   
c***********************************************************************
c   R1: Get Children names of a Node.  Return the name of children nodes 
c   directly associated with a parent node.  The names of the children 
c   are NOT guaranteed to be returned in any particular order.  If a 
c   new child is added, it is NOT guaranteed to be returned as the last 
c   child.  Note:  To start with the first child, use an istart value of 1.
c
c   f77: ADFCNAM( ID, istart, inum, inamlen, inumret, names, ierr )
c   input:  real*8 ID.         The ID of the Node to use.
c   input:  integer istart.    The Nth child's name to start with (first=1).
c   input:  integer inum.      Maximum number of names to return.
c   input:  integer inamlen.   Length of names.
c   output:  integer inumret.  The number of names returned.
c   output:  character*(*) names   Array of names.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFCNAM( ID, istart, inum, inamlen, inumret, names,
     1                    ierr )
      IMPLICIT NONE
      real*8 ID
      integer istart
      integer inum
      integer inamlen
      integer inumret
      character*(*) names
      integer ierr

      call ADFCNA2( ID, istart, inum, inamlen, len(names), inumret,
     1              names, ierr )
      return
      end

c***********************************************************************
c   R1: Create a Node.  Create a new node (not a link-node) as a child
c   of a given parent.  Default values in this new node are: label=blank, 
c   number of sub-nodes = 0, data-type = "MT", number of dimensions 
c   = 0, data = NULL.
c
c   f77: ADFCRE( PID, name, ID, ierr )
c   input:  real*8 PID.  The ID of the parent node, to whom we are 
c                        creating a new child node.
c   input:  character*(*) name.
c   output:  real*8 ID.  The ID of the newly created node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFCRE( PID, name, ID, ierr )
      IMPLICIT NONE
      real*8 PID
      character*(*) name
      real*8 ID
      integer ierr

      call ADFCRE2( PID, name, LEN( name ), ID, ierr )
      return
      end

c***********************************************************************
c   R1: Close an opened database.  If the ADF database spans multiple
c   files then all files used will also be closed.  If an ADF file which
c   is linked to by this database is also opened through another 
c   database, only the opened file stream associated with this database 
c   will be closed. 
c
c   f77: ADFDCLO( RootID, ierr )
c   input:  real*8 RootID  Root ID of the database.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDCLO( RootID, ierr )
      IMPLICIT NONE
      real*8 RootID
      integer ierr

      call ADFDCL2( RootID, ierr )
      return
      end

c***********************************************************************
c   Rn: Delete an existing database.  This will delete one or more ADF
c   files which are linked together under file top ADF file named 
c   "filename".
c
c   f77: ADFDDEL( filename, ierr )
c   input:  character*(*) filename
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDDEL( filename, ierr )
      IMPLICIT NONE
      character*(*) filename
      integer ierr

      call ADFDDE2( filename, LEN( filename ), ierr )
      return
      end

c***********************************************************************
c   R2: Delete a Node.   If the node is NOT a link, then the specified
c   node and all 
c   sub-nodes anywhere under it are also deleted.  For a link, and also 
c   for links farther down in the tree, the link-node will be deleted, 
c   but the node which the link is linked to is not affected.  When a 
c   node is deleted, other link-nodes which point to it are left 
c   dangling.  For example, if N13 is deleted, then L1 and L2 point to a 
c   non-existing node.  This is OK until L1 and L2 are used.
c
c   f77: ADFDEL( PID, ID, ierr )
c   input:  real*8 PID.  The ID of the node's parent.
c   input:  real*8 ID.  The ID of the node to use.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDEL( PID, ID, ierr )
      IMPLICIT NONE
      real*8 PID
      real*8 ID
      integer ierr

      call ADFDEL2( PID, ID, ierr )
      return
      end

c***********************************************************************
c   Rn: Garbage Collection.  Redistribute data in the file to  use free-
c   space which is not located at the end of the file.  Neighboring free 
c   spaces will be merged.  Note:  For better file compaction a utility 
c   could be written to copy an ADF file, creating a new ADF file 
c   without any wasted space.
c
c   f77: ADFDGC( ID, ierr )
c   input:  real*8 ID.  The ID of a node in the ADF file in which to do 
c                       garbage collection.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDGC( ID, ierr )
      IMPLICIT NONE
      real*8 ID
      integer ierr

      call ADFDGC2( ID, ierr )
      return
      end

c***********************************************************************
c   R1: Get the data format used in an existing database.
c
c   f77: ADFDGF( RootID, format, ierr )
c   input:  real*8 RootID The rootID of the ADF file.
c   output:  character*20 format.  See format for ADFDOPN.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDGF( RootID, format, ierr )
      IMPLICIT NONE
      real*8 RootID
      character*(*) format
      integer ierr

      call ADFDGF2( RootID, format, LEN( format ), ierr )
      return
      end

c***********************************************************************
c   R1: Open a database.  Open either a new or an existing ADF file.
c   If links to other ADF files are used, these additional files will be
c   opened automatically as required.
c   
c   f77: ADFDOPN( filename, status, format, rootID, ierr)
c   input:  character*(*) filename.  Not used if status SCRATCH is 
c               used.  Filename must be a legal name and may 
c               include a relative or absolute path.  It must be 
c               directly usable by the C fopen() system 
c               routine (no environment expansion is done).
c   input:  character*(*) status.  Like FORTRAN OPEN() status.  
c                This field is required.  Allowable values are:
c           READ_ONLY - File must exist.  Writing NOT allowed.
c           OLD - File must exist.  Reading and writing allowed.
c           NEW - File must not exist.
c           SCRATCH - New file.  Filename is ignored.
c           UNKNOWN - OLD if file exists, else NEW is used.
c   input:  character*(*) format.  Specifies the numeric format for 
c               the file.  If blank, the machine's native 
c               format is used.  This field is only used when a 
c               file is created.
c           NATIVE - Determine the format on the machine.  If the 
c               native format is not one of the formats 
c               supported, the created file cannot be used on 
c               other machines.
c           IEEE_BIG - Use the IEEE big ENDIAN format.
c           IEEE_LITTLE - Use the IEEE little ENDIAN format.
c           CRAY - Use the native Cray format.
c   output:  real*8 rootID
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDOPN( filename, status, format, rootID, ierr)
      IMPLICIT NONE
      character*(*) filename
      character*(*) status
      character*(*) format
      real*8 rootID
      integer ierr

      call ADFDOP2( filename, len( filename ), status, len( status ),
     1      format, len( format ), rootID, ierr )
      return
      end

c***********************************************************************
c   R1: Set the data format used in an existing database.
c    Note:  Use with extreme caution.  Needed only 
c   for data conversion utilities and NOT intended 
c   for the general user!!!
c
c   f77: ADFDSF( RootID, format, ierr )
c   input:  real*8 RootID The rootID if the ADF file.
c   input:  character*(*) format.  See format for ADFDOPN.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDSF( RootID, format, ierr )
      IMPLICIT NONE
      real*8 RootID
      character*(*) format
      integer ierr

      call ADFDSF2( RootID, format, len( format ), ierr )
      return
      end

c***********************************************************************
c   R1: Get ADF File Version ID.  This is the version number of the ADF
c   library routines which created an ADF database.  Modified ADF databases 
c   will take on the version ID of the current ADF library version if 
c   it is higher than the version indicated in the file.
c    The format of the version ID is:  "ADF Database Version 000.01"
c
c   f77: ADFDVER( RootID, version, cdate, mdate, ierr )
c   input:  real*8 RootID.  The ID of the root node in the ADF file.
c   output:  character(32) version.  A 32-byte character string 
c               containing the version ID.
c   output:  character(32) cdate.  A 32-byte character string 
c               containing the creation date of the file.
c   output:  character(32) mdate.  A 32-byte character 
c               string containing the last modification date of the file.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFDVER( RootID, version, cdate, mdate, ierr )
      IMPLICIT NONE
      real*8 RootID
      character*(*) version
      character*(*) cdate
      character*(*) mdate
      integer ierr

      call ADFDVE2( RootID, version, cdate, mdate,
     1  len( version ), len( cdate ), len( mdate ), ierr )
      return
      end

c***********************************************************************
c   R1: Return Error Message.  Given an ierr from an ADF routine, 
c   get a textual description of the error.
c
c   f77: ADFERR( ierr, errstr )
c   input:  integer ierr.
c   output:  character(80) errstr.  An 80-byte description of 
c             the specified error.  If the number is bad, the 
c             string "Unknown error #nnn" is returned.
c***********************************************************************
      subroutine ADFERR( ierr, errstr )
      IMPLICIT NONE
      integer ierr
      character*(*) errstr

      call ADFERR2( ierr, errstr, len( errstr ) )
      return
      end

c***********************************************************************
c   R1: Flush data to disk.  This routine is used to force any modified 
c   information to be flushed to the physical disk.  This ensures that 
c   data will not be lost if a program aborts.  This control of when to 
c   flush all data to disk is provided to the user rather than to flush 
c   the data every time it is modified, which would result in reduced 
c   performance.
c
c   f77: ADFFTD( ID, ierr )
c   input:  real*8 ID.  The ID of a node in the ADF file in which to flush.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFFTD( ID, ierr )
      IMPLICIT NONE
      real*8 ID
      integer ierr

      call ADFFTD2( ID, ierr )
      return
      end

c***********************************************************************
c   R1: Get Data Type.  Return the 32 character string in a node's data-
c   type field.  In C, the data-type will be null terminated after the last 
c   non-blank character.  A maximum of 33 characters may be used 
c   (32 for the data-type plus 1 for the null).
c
c   f77: ADFGDT( ID, dtype, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  character*(*) dtype.  The 32-character data-type of the node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGDT( ID, dtype, ierr )
      IMPLICIT NONE
      real*8 ID
      character*(*) dtype
      integer ierr

      call ADFGDT2( ID, dtype, len( dtype ), ierr )
      return
      end

c***********************************************************************
c   R1: Get Dimension Values.  Return the dimension values for a node.
c   Values will only be returned for the number of dimensions defined in 
c   the node.  If the number of dimensions for the node is zero, an 
c   error is returned.
c
c   f77: ADFGDV( ID, dvals, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  integer dvals(12).
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGDV( ID, dvals, ierr )
      IMPLICIT NONE
      real*8 ID
      integer dvals(12)
      integer ierr

      call ADFGDV2( ID, dvals, ierr )
      return
      end

c***********************************************************************
c   R1      Get Error State.  Return the active error state.
c
c   f77: ADFGES( estate, ierr )
c   output:  integer estate.  Flag for ABORT on error (1) or return 
c                 error status (0).  Set on a per database basis.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGES( estate, ierr )
      IMPLICIT NONE
      integer estate
      integer ierr

      call ADFGES2( estate, ierr )
      return
      end

c***********************************************************************
c   R1: Get Label.  Return the 32 character string in a node's label field.
c   In C, the label will be null terminated after the last non-blank 
c   character.  A maximum of 33 characters may be used (32 for the 
c   label plus 1 for the null).
c
c   f77: ADFGLB( ID, label, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  character*(*) label.  The 32-character label of the node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGLB( ID, label, ierr )
      IMPLICIT NONE
      real*8 ID
      character*(*) label
      integer ierr

      call ADFGLB2( ID, label, len( label ), ierr )
      return
      end

c***********************************************************************
c   R1: Get path information from a link.  If the node is a link-node,
c   return the path information.  Else, return an error.
c   If the link is in the same file, then a blank filename is returned.
c
c   f77: ADFGLKP( ID, file, name, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  character*(*) file.  The filename. 
c   output:  character*(*) name.  The name of node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGLKP( ID, file, name, ierr )
      IMPLICIT NONE
      real*8 ID
      character*(*) file
      character*(*) name 
      integer ierr

      call ADFGLK2( ID, file, len( file ), name, len( name ), ierr )
      return
      end

c***********************************************************************
c   R1: Get Name of a Node.  Given a node's ID, return the 32 character 
c   name of that node.  In C, the name will be null terminated after the
c   last non-blank character.  A maximum of 33 characters may be used 
c   (32 for the name plus 1 for the null).
c
c   f77: ADFGNAM( ID, name, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  character*(*) name.  The simple name of the node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGNAM( ID, name, ierr )
      IMPLICIT NONE
      real*8 ID
      character*(*) name
      integer ierr

      call ADFGNA2( ID, name, len( name ), ierr )
      return
      end

c***********************************************************************
c   R1: Get Number of Dimensions.  Return the number of data dimensions 
c   used in a node.  Valid values are from 0 to 12.
c
c   f77: ADFGND( ID, ndims, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  integer ndims.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGND( ID, ndims, ierr )
      IMPLICIT NONE
      real*8 ID
      integer ndims
      integer ierr

      call ADFGND2( ID, ndims, ierr )
      return
      end

c***********************************************************************
c   R1: Get Unique-Identifier of a Node.  Given a parent node ID and a
c   a name of child node, return the ID of the child.  To return the ID
c   of the root-node in an ADF file, use any known ID in the ADF file
c   and a name of "/".
c
c   f77: ADFGNID( PID, name, ID, ierr )
c   input:  real*8 PID.  The ID of the Node's parent.
c   input:  character*(*) name.  The name of the node.  Compound 
c                names including path information use a slash 
c                "/" notation between node names.  If a 
c                leading slash is used, then PID can be any 
c                valid node ID in the ADF database.
c   output:  real*8 ID.  The ID of the named node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGNID( PID, name, ID, ierr )
      IMPLICIT NONE
      real*8 PID
      character*(*) name
      real*8 ID
      integer ierr

      call ADFGNI2( PID, name, len( name ), ID, ierr )
      return
      end

c***********************************************************************
c   R1: Get root-ID for an ADF system from any ID in the system.
c
c   f77: ADFGRID( ID, RootID, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output: real*8 RootID.  The ID of the root node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFGRID( ID, RootID, ierr )
      IMPLICIT NONE
      real*8 ID
      real*8 RootID
      integer ierr

      call ADFGRI2( ID, RootID, ierr )
      return
      end

c***********************************************************************
c   R1: Test if a Node is a link.  If the actual data-type of the node
c   is "LK" (created with ADF_Link), return the link path length.  
c   Otherwise, return 0.
c
c   f77: ADFISLK( ID, lplen, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   input:  integer lplen.  0 if the node is NOT a link.  If 
c             the node is a link, the length of the path string is returned.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFISLK( ID, lplen, ierr )
      IMPLICIT NONE
      real*8 ID
      integer lplen
      integer ierr

      call ADFISL2( ID, lplen, ierr )
      return
      end

c***********************************************************************
c   R1: Create a link.  Note:  The Node linked to does not have to exist
c   when the link is created (but it may exist and that is OK).  However,
c   when the link is used, an error will occur if the linked to node does
c   not exist.
c
c   f77: ADFLINK( PID, name, file, nfile, ID, ierr )
c   input:  real*8 PID.  The ID of the Node's parent.
c   input:  character*(*) name.  The name of the link node.
c   input:  character*(*) file.  The filename to use for the link 
c             (directly usable by a C open() routine).  If 
c             blank (null), the link will be within the same file.
c   input:  character*(*) nfile.  The name of the node which 
c             the link will point to.  This can be a simple or 
c             compound name.
c   output:  real*8 ID.  The ID of the link-node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFLINK( PID, name, file, nfile, ID, ierr )
      IMPLICIT NONE
      real*8 PID
      character*(*) name
      character*(*) file
      character*(*) nfile
      real*8 ID
      integer ierr

      call ADFLIN2( PID, name, file, nfile, len( name ),
     1      len( file ), len( nfile ), ID, ierr )
      return
      end

c***********************************************************************
c   R1: Get ADF Library Version ID.  This is the version number of the
c   ADF library routines which your program is currently using.
c    The format of the version ID is:  "ADF Library  Version 000.01"
c    Note:  There is a double space between Library and Version.  This
c   lines up the number with the Database Version string.
c
c   f77: ADFLVER( version, ierr )
c   output:  character(32) version.  A 32-byte character string 
c                containing the ADF Library version ID information.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFLVER( version, ierr )
      IMPLICIT NONE
      character*(*) version
      integer ierr

      call ADFLVE2( version, len( version ), ierr )
      return
      end

c***********************************************************************
c   R2: Change Parent (move a Child Node).  The node and the 2 parents
c   must all exist within a single ADF file.  If the node is pointed to
c   by a link-node, changing the node's parent will break the link.
c
c   f77: ADFMOVE( PID, ID, NPID, ierr )
c   input:  real*8 PID.  The ID of the Node's parent.
c   input:  real*8 ID.  The ID of the node to use.
c   input:  real*8 NPID.  The ID of the Node's New Parent 
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFMOVE( PID, ID, NPID, ierr )
      IMPLICIT NONE
      real*8 PID
      real*8 ID
      real*8 NPID
      integer ierr

      call ADFMOV2( PID, ID, NPID, ierr )
      return
      end

c***********************************************************************
c   R1: Get Number of Children of a Node.  Return the number of children 
c   nodes directly associated with a parent node.
c
c   f77: ADFNCLD( ID, numcld, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   output:  integer numcld.  The number of children directly 
c                associated with this node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFNCLD( ID, numcld, ierr )
      IMPLICIT NONE
      real*8 ID
      integer numcld
      integer ierr

      call ADFNCL2( ID, numcld, ierr )
      return
      end

c***********************************************************************
c   R1: Set/change the data-type and Dimension Information of a Node.
c   Valid user-definable data-types are:
c   
c   No data    MT
c   Integer 32    I4
c   Integer 64    I8
c   Unsigned Int 32   U4
c   Unsigned Int 64   U8
c   Real 32    R4
c   Real 64    R8
c   Complex 64    X4
c   Complex 128    X8
c   Character (unsigned byte) C1
c   Byte (unsigned byte)  B1
c   Compound data-types can be used which combine types 
c   ("I4,I4,R8"), define an array ("I4[25]"), or a combination of these 
c   ("I4,C1[20],R8[3]").
c   dims can be a number from 0 to 12.
c   dim_vals is an array of integers.  The number of integers used is 
c   determined by the dims argument.  If dims is zero, the dim_values 
c   are not used.  Valid range for dim_values are from 1 to 
c   2,147,483,648.  The total data size, calculated by the data-type-size 
c   times the dimension value(s), cannot exceed 2,147,483,648.
c   Note:  When this routine is called and the data-type or the 
c   number of dimensions changes, any data currently associated 
c   with the node is lost!!   The dimension values can be changed and 
c   the data space will be extended as needed.
c-
c   f77: ADFPDIM( ID, dtype, dims, dvals, ierr)
c   input:  real*8 ID.  The ID of the node.
c   input:  character*(*) dtype.
c   input:  integer dims.  The number of dimensions this node has.
c   input:  integer dvals(12).  The dimension values for this node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFPDIM( ID, dtype, dims, dvals, ierr)
      IMPLICIT NONE
      real*8 ID
      character*(*) dtype
      integer dims
      integer dvals(12)
      integer ierr

      call ADFPDI2( ID, dtype, len( dtype ), dims, dvals, ierr)
      return
      end

c***********************************************************************
c   R2: Put (change) Name of a Node.  Warning:  If the node is pointed
c   to by a link-node, changing the node's name will break the link.
c
c   f77: ADFPNAM( PID, ID, name, ierr )
c   input:  real*8 PID.  The ID of the Node's parent.
c   input:  real*8 ID.  The ID of the node to use.
c   input:  character*(*) name.  The new name of the node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFPNAM( PID, ID, name, ierr )
      IMPLICIT NONE
      real*8 PID
      real*8 ID
      character*(*) name
      integer ierr

      call ADFPNA2( PID, ID, name, len( name ), ierr )
      return
      end

c***********************************************************************
c   R1      Set Error State.  For all ADF calls, set the error handling
c   convention; either return error codes, or abort the program on an error.
c   The default state for the ADF interface is to return error codes and NOT
c   abort.
c
c   f77: ADFSES( estate, ierr )
c   input:  integer estate.  Flag for ABORT on error (1) or return 
c                 error status (0).  Set on a per database basis.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFSES( estate, ierr )
      IMPLICIT NONE
      integer estate
      integer ierr

      call ADFSES2( estate, ierr )
      return
      end

c***********************************************************************
c   R1: Set Label.  Set the 32 character string in a node's label field.
c
c   f77: ADFSLB( ID, label, ierr )
c   input:  real*8 ID.  The ID of the node to use.
c   input:  character*(*) label.  The 32-character label of the node.
c   output:  integer ierr.
c***********************************************************************
      subroutine ADFSLB( ID, label, ierr )
      IMPLICIT NONE
      real*8 ID
      character*(*) label
      integer ierr

      call ADFSLB2( ID, label, len( label ), ierr )
      return
      end
