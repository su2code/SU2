	program write_cgns_1
	implicit none

!       author: Diane Poirier (diane@icemcfd.com)
!       last revised on March 8 2000

!       This example test the complete SIDS for multi-block data.
!	It creates a dummy mesh composed of 2 structured blocks in 3D.
#include "cgnstypes_f.h"
#ifdef WINNT
	include 'cgnswin_f.h'
#endif
	include 'cgnslib_f.h'

	integer Ndim
	parameter (Ndim = 3)
	cgsize_t one
	parameter (one = 1)
	
	cgsize_t index_dim, cell_dim, phys_dim
	cgsize_t base_no, zone_no, coord_no, sol_no, discr_no, conn_no
	cgsize_t hole_no, boco_no, field_no, dset_no
	cgsize_t num, NormalIndex(Ndim)
	cgsize_t ndims, npnts,size(Ndim*3)
	cgsize_t cg, ier, zone, coord, i, j, k, n, pos, sol, field
	cgsize_t pnts(Ndim,120), donor_pnts(Ndim,120)
	cgsize_t transform(Ndim)
	cgsize_t nptsets, nrmlistflag
	real*4 data(120), normals(360)
	double precision Dxyz(120), values(120)
	character*32 zonename, solname, fieldname
	character*32 coordname(Ndim)
	character*32 donorname

	coordname(1) = 'CoordinateX'
	coordname(2) = 'CoordinateY'
	coordname(3) = 'CoordinateZ'

! *** initialize
	ier = 0
	index_dim=Ndim
	cell_dim=Ndim
	phys_dim=Ndim

! *** open CGNS file for writing
 	call cg_open_f('cgtest.cgns', CG_MODE_WRITE, cg, ier)
 	if (ier .eq. ERROR) call cg_error_exit_f

! *** base
 	call cg_base_write_f(cg, 'Basename', cell_dim, phys_dim,
     &                       base_no, ier)
 	if (ier .eq. ERROR) call cg_error_exit_f

! *** zone
	do zone=1, 2
    	    write(zonename,'(a5,i1)') 'zone#',zone
    	    num = 1
	    do i=1,index_dim          		! zone#1: 3*4*5, zone#2: 4*5*6
                size(i) = i+zone+1		! nr of nodes in i,j,k
		size(i+Ndim) = size(i)-1	! nr of elements in i,j,k
	 	size(i+2*Ndim) = 0		! nr of bnd nodes if ordered
       		num = num * size(i)		! nr of nodes
    	    enddo
!234567890!234567890!234567890!234567890!234567890!234567890!23456789012
    	    call cg_zone_write_f(cg, base_no, zonename, size,
     &                           Structured, zone_no, ier)
    	    if (ier .eq. ERROR) call cg_error_exit_f

! *** coordinate
     	    do coord=1, phys_dim
 		do k=1, size(3)
 		do j=1, size(2)
 		do i=1, size(1)
		    pos = i + (j-1)*size(1) + (k-1)*size(1)*size(2)
		! * make up some dummy coordinates just for the test:
 	    	    if (coord.eq.1) Dxyz(pos) = i
 	    	    if (coord.eq.2) Dxyz(pos) = j
 	    	    if (coord.eq.3) Dxyz(pos) = k
 		enddo
 		enddo
 		enddo

         	call cg_coord_write_f(cg, base_no, zone_no, RealDouble,
     &                            coordname(coord), Dxyz, coord_no, ier)
         	if (ier .eq. ERROR) call cg_error_exit_f

     	    enddo

! *** solution
     	    do sol=1, 2
 		write(solname,'(a5,i1,a5,i1)') 'Zone#',zone,' sol#',sol
 		call cg_sol_write_f(cg, base_no, zone_no, solname, 
     &                              Vertex, sol_no, ier)
 		if (ier .eq. ERROR) call cg_error_exit_f

! *** solution field
 		do field=1, 2
		    ! make up some dummy solution values
 	    	    do i=1, num
 		  	values(i) = i*field*sol
 	    	    enddo
 	    	    write(fieldname,'(a6,i1)') 'Field#',field
 	    	    call cg_field_write_f(cg, base_no, zone_no, sol_no,
     &                  RealDouble, fieldname, values, field_no, ier)
 	    	    if (ier .eq. ERROR) call cg_error_exit_f

 		enddo				! field loop
            enddo				! solution loop

! *** discrete data
	    call cg_discrete_write_f(cg, base_no, zone_no, 'discrete#1',
     &                               discr_no, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

  	! *** discrete data arrays, defined on vertices:
	    call cg_goto_f(cg, base_no, ier, 'Zone_t', zone,
     &                     'DiscreteData_t', discr_no, 'end')
            if (ier .eq. ERROR) call cg_error_exit_f

	    do 123 k=1, size(3)
            do 123 j=1, size(2)
            do 123 i=1, size(1)
                pos = i + (j-1)*size(1) + (k-1)*size(1)*size(2)
	        data(pos) = pos	! * make up some dummy data
 123	    continue
	    call cg_array_write_f('arrayname', RealSingle, index_dim,
     &                             size, data, ier)
            if (ier .eq. ERROR) call cg_error_exit_f

  	! *** discrete data arrays attribute: GOTO DataArray node
	    call cg_goto_f(cg, base_no, ier, 'Zone_t', zone,
     &	            'DiscreteData_t', discr_no, 'DataArray_t', 1, 'end')
	    if (ier .eq. ERROR) call cg_error_exit_f

	    call cg_units_write_f(Kilogram, Meter, Second, Kelvin, 
     &                            Radian, ier)
            if (ier .eq. ERROR) call cg_error_exit_f
		
! *** overset holes
	 !  create dummy data
	    do i=1,3
	      ! Define 2 separate PointRange, for 2 patches in the hole
		pnts(i,1)=1
		pnts(i,2)=size(i)
	      ! second PointRange of hole
		pnts(i,3)=2
		pnts(i,4)=size(i)
	    enddo
	  ! Hole defined with 2 point set type PointRange, so 4 points:
	    nptsets = 2
	    npnts = 4
	    call cg_hole_write_f(cg, base_no, zone_no, 'hole#1', Vertex,
     &                           PointRange, nptsets, npnts, pnts,
     &                           hole_no, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! *** general connectivity
	    do 100 n=1, 5
	    do 100 i=1,3
	        pnts(i,n)=i		! * dummy data
		donor_pnts(i,n)=i*2
 100        continue
	  ! create a point matching connectivity 
	    npnts = 5
	    call cg_conn_write_f(cg, base_no, zone_no, 'Connect#1', 
     &          Vertex, Abutting1to1, PointList, npnts, pnts, 'zone#2',
     &          Structured, PointListDonor, Integer, npnts, donor_pnts,
     &          conn_no, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! *** connectivity 1to1
	 !  generate data
	    do i=1,3
	     !**make up some dummy data:
                pnts(i,1)=1
                pnts(i,2)=size(i)
		donor_pnts(i,1)=1
		donor_pnts(i,2)=size(i)
		transform(i)=i*(-1)
            enddo
	    if (zone .eq. 1) then
		donorname='zone#2'
	    else if (zone .eq. 2) then
		donorname='zone#1'
	    endif

	    call cg_1to1_write_f(cg, base_no, zone_no, '1to1_#1',
     &		donorname, pnts, donor_pnts, transform, conn_no, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! *** ZoneGridConnectivity attributes:  GOTO ZoneGridConnectivity_t node
	    call cg_goto_f(cg, base_no, ier, 'Zone_t', zone,
     &                     'ZoneGridConnectivity_t', one, 'end')
	    if (ier .eq. ERROR) call cg_error_exit_f

! *** ZoneGridConnectivity attributes: Descriptor_t
!234567890!234567890!234567890!234567890!234567890!234567890!23456789012
	    call cg_descriptor_write_f('DescriptorName',
     &                                 'Zone Connectivity', ier)

! *** bocos
            npnts = 2
	    call cg_boco_write_f(cg, base_no, zone_no, 'boco#1',
     &           BCInflow, PointRange, npnts, pnts, boco_no, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! *** boco normal
	    npnts = 1
	    do i=1,Ndim
	        NormalIndex(i)=0
		! compute nr of points on bc patch:
		npnts = npnts * (pnts(i,2)-pnts(i,1)+1)
	    enddo
	    NormalIndex(1)=1
	    do i=1,phys_dim*npnts
		normals(i)=i 
	    enddo

            nrmlistflag = 1
	    call cg_boco_normal_write_f(cg, base_no, zone_no, boco_no,
     &         NormalIndex, nrmlistflag, RealSingle, normals, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! ** boundary condition attributes: GOTO BC_t node
	    call cg_goto_f(cg, base_no, ier, 'Zone_t', zone, 'ZoneBC_t',
     &          one, 'BC_t', boco_no, 'end')
	    if (ier .eq. ERROR) call cg_error_exit_f

! ** boundary condition attributes:  GridLocation_t
	    call cg_gridlocation_write_f(Vertex, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! ** boundary condition dataset
            call cg_dataset_write_f(cg, base_no, zone,
     &         boco_no, 'DataSetName', BCInflow, dset_no, ier)
            if (ier .eq. ERROR) call cg_error_exit_f

! ** boundary condition data:
	    call cg_bcdata_write_f(cg, base_no, zone,
     &         boco_no, dset_no, Neumann, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	     
! ** boundary condition data arrays: GOTO BCData_t node
	    call cg_goto_f(cg, base_no, ier, 'Zone_t', zone_no, 
     &          'ZoneBC_t', one, 'BC_t', boco_no, 'BCDataSet_t',
     &          dset_no, 'BCData_t', Neumann, 'end') 
	    if (ier .eq. ERROR) call cg_error_exit_f

	    do i=1, npnts
		data(i) = i
	    enddo
	    ndims = 1
	    call cg_array_write_f('dataset_arrayname', RealSingle, 
     &           ndims, npnts, data, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

! ** boundary condition data attributes: 
	    call cg_dataclass_write_f(NormalizedByDimensional, ier)
            if (ier .eq. ERROR) call cg_error_exit_f

 	enddo					! zone loop

! *** close CGNS file
	call cg_close_f(cg, ier)
	if (ier .eq. ERROR) call cg_error_exit_f

	end
