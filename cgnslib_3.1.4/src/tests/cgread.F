
	program read_cgns_1
	implicit none

! This program reads a 3D mesh, structured or unstructured.

#include "cgnstypes_f.h"
#ifdef WINNT
	include 'cgnswin_f.h'
#endif
	include 'cgnslib_f.h'

	integer Ndim, Nglobal
	parameter (Ndim = 3)
	parameter (Nglobal = 500)

	cgsize_t i, narrays, iarray, nintegrals, integral
	cgsize_t ndescriptors, idescr, nptsets, nzonetype
	cgsize_t ndonor_ptset_type, ndonor_data_type
	cgsize_t idataset, dirichletflag, neumannflag
!234567890!234567890!234567890!234567890!234567890!234567890!23456789012
	cgsize_t IndexDim, CellDim, PhysDim
	cgsize_t ier, n, zonetype
	cgsize_t nbases, nzones
	cgsize_t rmin(3), DataSize(Ndim)
	cgsize_t size(Ndim*3)
	cgsize_t ncoords, type, nsols, nfields, location
	cgsize_t nholes, nconns, n1to1, n1to1_global, nbocos
	cgsize_t ptset_type
	cgsize_t npnts, pnts(100000), donor_pnts(100000)
        cgsize_t npnts_donor
        cgsize_t bocotype, datatype
	character*32 basename, zonename, solname, fieldname 
	character*32 coordname, holename, connectname
	character*32 boconame, donorname
	cgsize_t cg, base, zone, coord, sol, field, discr
	cgsize_t hole, conn, one21, boco
	cgsize_t range(Ndim, 2), donor_range(Ndim, 2)
 	cgsize_t transform(Ndim)
  	cgsize_t G_range(Ndim*2, Nglobal) 
	cgsize_t G_donor_range(Ndim*2, Nglobal)
	cgsize_t G_transform(Ndim, Nglobal)
	character*32 G_connectname(Nglobal), G_zonename(Nglobal)
	character*32 G_donorname(Nglobal)
	character*32 name, filename
	character*40 text, NormDefinitions, StateDescription
	cgsize_t equation_dimension, GoverningEquationsFlag
	cgsize_t GasModelFlag, ViscosityModelFlag
	cgsize_t ThermalConductivityModelFlag 
	cgsize_t TurbulenceClosureFlag, TurbulenceModelFlag
	cgsize_t diffusion_model(6), niterations
	cgsize_t nndim, rind(6), ndiscrete, num
	cgsize_t dim_vals(12)
	cgsize_t mass, length, time, temp, deg
	cgsize_t NormalIndex(3), ndataset
	cgsize_t NormalListFlag
	real*4 data_single(100000)
	double precision data_double(100000)
	real*4 version

        cgsize_t one
        parameter (one = 1)

! *** open file
!	write(6,*) 'Input filename'
!	read(5,600) filename
 	write(filename,'(a)')'cgtest.cgns'
	call cg_open_f(filename, CG_MODE_READ, cg, ier)
 	if (ier .eq. ERROR) call cg_error_exit_f
	write(6,600)'READING FILE ',filename

! *** CGNS Library Version used for file creation:
	call cg_version_f(cg, version, ier)
	if (ier .eq. ERROR) call cg_error_exit_f
	write(6,102)
     &     'Library Version used for file creation:',version

! *** base
	call cg_nbases_f(cg, nbases, ier)
	if (ier .eq. ERROR) call cg_error_exit_f
	write(6,200)'nbases=',nbases
	
	do base=1, nbases

!234567890!234567890!234567890!234567890!234567890!234567890!23456789012
	  call cg_base_read_f(cg, base, basename, CellDim, PhysDim, ier)
	  if (ier .eq. ERROR) call cg_error_exit_f
	  write(6,300)'BaseName = "',basename,'"', 
     &                'cell_dimension=',CellDim

! *** base attribute:  GOTO base node
	  call cg_goto_f(cg, base, ier, 'end')
          if (ier .eq. ERROR) call cg_error_exit_f

! ***     base attribute:  Descriptor
          call cg_descriptor_read_f(one, name, text, ier)
          if (ier .eq. ERROR) call cg_error_exit_f
          if (ier.eq.ALL_OK) then
            write(6,400)'Base Descriptor_t Information:'
            write(6,500)' DescriptorName="',name,'"',
     &                  ' DescriptorText="',text,'"'
          endif

! ***     base attribute: flow equation set:
	  call cg_equationset_read_f(equation_dimension, 
     &      GoverningEquationsFlag,  GasModelFlag, 
     &      ViscosityModelFlag, ThermalConductivityModelFlag, 
     &      TurbulenceClosureFlag,  TurbulenceModelFlag, ier)
	  if (ier .eq. ERROR) then
	    call cg_error_exit_f
	  elseif (ier .eq. NODE_NOT_FOUND) then
	    write(6,200)
     &      'FlowEquationSet_t not defined under CGNSBase_t #',base
	  elseif (ier .eq. INCORRECT_PATH) then
  	    write(6,400)'Incorrect path input to cg_goto_f'
	  else
	    write(6,400) 'FlowEquationSet_t Information:' 
	    write(6,100)' equation_dimension=',equation_dimension

! ***       flow equation set attributes:  GOTO FlowEquationSet_t node
            call cg_goto_f(cg,base,ier,'FlowEquationSet_t',one,'end')
            if (ier .eq. ERROR) call cg_error_exit_f

! ***       flow equation set attribute: Descriptor
            call cg_descriptor_read_f(one, name,text,ier)
            if (ier .eq. ERROR) call cg_error_exit_f
            if (ier .eq. ALL_OK) write(6,500)
     &        ' DescriptorName="',name,'"',' DescriptorText="',text,'"'

! ***       flow equation set attribute: Gas Model Type
	    if (GasModelFlag.eq.1) then
                call cg_model_read_f('GasModel_t', type, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                if (ier .eq. ALL_OK) write(6,600)
     &            ' GasModelType="',ModelTypeName(type),'"'
            endif

! ***       flow equation set attribute: ViscosityModel Type
            if (ViscosityModelFlag.eq.1) then
              call cg_model_read_f('ViscosityModel_t', type, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              if (ier .eq. ALL_OK) write(6,600)
     &          ' ViscosityModelType="',ModelTypeName(type),'"'
            endif

! ***       flow equation set attribute:  TypmlConductivityModel Type
            if (ThermalConductivityModelFlag.eq.1) then
              call cg_model_read_f('ThermalConductivityModel_t',
     &                  type, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              if (ier .eq. ALL_OK) write(6,600)
     &          ' ThermalConductivityModelType=',
     &            ModelTypeName(type),'"'
            endif

! ***   flow equation set attribute: TurbulenceClosureType
            if (TurbulenceClosureFlag.eq.1) then
              call cg_model_read_f('TurbulenceClosure_t', type, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              if (ier .eq. ALL_OK) write(6,600)
     &          ' TurbulenceClosureType="', ModelTypeName(type),'"'
            endif

! ***   flow equation set attribute: TurbulenceModelType
            if (TurbulenceModelFlag.eq.1) then
              call cg_model_read_f('TurbulenceModel_t', type, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              if (ier .eq. ALL_OK) write(6,600)
     &          ' TurbulenceModelType="',ModelTypeName(type),'"'
            endif

! ***   flow equation set attribute: Governing Equations Type
            if (GoverningEquationsFlag .eq. 1) then
              call cg_governing_read_f(type, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              if (ier.eq.ALL_OK) 
     &          write(6,600)' GoverningEquationsType="',
     &                          GoverningEquationsTypeName(type),'"'

    	! *** Governing Equations attribute:  GOTO GoverningEquations_t node
	      call cg_goto_f(cg,base,ier, 'FlowEquationSet_t', one,
     &            'GoverningEquations_t', one ,'end')
	      if (ier .eq. ERROR) call cg_error_exit_f
	      

        ! *** Governing Equations attribute:  Diffusion model
	      call cg_diffusion_read_f(diffusion_model, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      if (ier.eq.ALL_OK)write(6,103)'   Diffusion model=',
     &                           (diffusion_model(i), i=1,6)
            endif       ! If Governing Equations are defined
	  endif		! If FlowEquationSet_t exists under CGNSBase_t


	  write(6,400)'                              *     *     *'

	  call cg_nzones_f(cg, base, nzones, ier)
	  if (ier .eq. ERROR) call cg_error_exit_f
	  write(6,200)'nzones=',nzones

! *** zone
	  do zone=1, nzones
	    call cg_zone_read_f(cg, base, zone, zonename, size, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	    write(6,104)'Name of Zone',zone,' is "',zonename,'"'

	    call cg_zone_type_f(cg, base, zone, zonetype, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	    write(6,600)'  Zone type is ', ZoneTypeName(zonetype)


	    if (zonetype.eq.Structured) then
		IndexDim=CellDim
	    else 
		IndexDim=1
	    endif

	
	    write(6,104)'  IndexDimension=',IndexDim

      ! *** zone attribute:  GOTO zone node
	    call cg_goto_f(cg, base, ier, 'Zone_t', zone, 'end')	
	    if (ier .eq. ERROR) call cg_error_exit_f

      ! *** zone attribute:  ordinal
	    call cg_ordinal_read_f(num, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	    if (ier .eq. ALL_OK) 
     &          write(6,200)' Zone ordinal=',num


      ! *** zone attribute:  convergence history
	    call cg_convergence_read_f(niterations, 
     &            NormDefinitions, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f

	    if (ier .eq. ALL_OK) then
	      write(6,600)'Convergence History of ',zonename
	      write(6,104) ' niterations=',niterations,
     &                     ' NormDefinitions="',NormDefinitions,'"'
	
         ! ** ConvergenceHistory_t attributes:
	      call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &                      'ConvergenceHistory_t', one, 'end')
	      if (ier .eq. ERROR) call cg_error_exit_f

         ! ** ConvergenceHistory_t attributes: DataArray_t
	      call cg_narrays_f(narrays, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,105) 'ConvergenceHistory_t contains ',
     &                      narrays,' array(s)'
	      do iarray=1, narrays  
                call cg_array_info_f(iarray, name, datatype, 
     &                             nndim, dim_vals, ier)
                if (ier .eq. ERROR) call cg_error_exit_f 

	        write(6,600) ' DataArrayName="',name,'"'
	        write(6,600) ' DataType="',DataTypeName(datatype),'"'
	        write(6,200) ' DataNdim=',nndim
	        write(6,200) ' DataDim=',dim_vals(1)

	        write(6,105) ' Data:'
	        if (datatype .eq. RealSingle) then
	          call cg_array_read_f(iarray, data_single, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f 
	          write(6,106) (data_single(n),n=1,dim_vals(1))
	        elseif (datatype .eq. RealDouble) then
		  call cg_array_read_f(iarray, data_double, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f 
	 	  write(6,106) (data_double(n),n=1,dim_vals(1))
	        endif
	      enddo

         ! ** ConvergenceHistory_t attributes: DataClass_t
	      call cg_dataclass_read_f(type,ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,600)'DataClassName=',DataClassName(type)

         ! ** ConvergenceHistory_t attributes: DimensionalUnits_t
	      call cg_units_read_f(mass, length, time, temp, deg, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      if (ier .eq. ALL_OK) then
	        write(6,100) 
     &	 	  'Dimensional Units:',
     &            MassUnitsName(mass), LengthUnitsName(length),
     &            TemperatureUnitsName(temp), TimeUnitsName(time),
     &            AngleUnitsName(deg)
	      endif
	    endif
            write(6,400)'                             *     *     *'

! *** zone attribute:  return to Zone_t node
	    call cg_goto_f(cg, base, ier, 'Zone_t', zone, 'end')
	    if (ier .eq. ERROR) call cg_error_exit_f
	    write(6,401)'Integral Data Information of ',zonename

	    call cg_nintegrals_f(nintegrals, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	    write(6,107) nintegrals, ' IntegralData_t node in ',
     &                    zonename

! *** zone attribute:  IntegralData_t
	    do integral=1, nintegrals
	      call cg_integral_read_f(integral, name, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,104) 'IntegralData_t #',integral,
     &                             ' is named "', name,'"'

        ! *** IntegralData_t attribute:  GOTO IntegralData_t node
	      call cg_goto_f(cg, base, ier, 'Zone_t', zone, 
     &                       'IntegralData_t', integral, 'end')
              if (ier .eq. ERROR) call cg_error_exit_f

	      call cg_narrays_f(narrays, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
              write(6,108) 'IntegralData_t #',integral,
     &                     ' contains ', narrays,' data'

	      do iarray=1, narrays

          ! *** IntegralData_t attribute: DataArray_t
                call cg_array_info_f(iarray, name, datatype,
     &                               nndim, dim_vals, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                write(6,600) ' DataArrayName="',name,'"'
                write(6,600) ' DataType=',DataTypeName(datatype)
                write(6,108) ' DataNdim=',nndim,
     &                       ', DataDim=',dim_vals(1)

                if (datatype .eq. RealSingle) then
                  call cg_array_read_f(iarray, data_single, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  write(6,109) ' integraldata=',data_single(1)
                elseif (datatype .eq. RealDouble) then
                  call cg_array_read_f(iarray, data_double, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  write(6,109) 'integraldata=',data_double(1)
                endif

	  ! *** DattaArray_t attribute: GOTO DataArray_t
		call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &                         'IntegralData_t', integral, 
     &                         'DataArray_t', iarray, 'end')
                if (ier .eq. ERROR) call cg_error_exit_f


          ! *** DattaArray_t attribute: DimensionalExponents_t
	        call cg_exponents_info_f(datatype, ier)
                if (ier .eq. ERROR) then
	          call cg_error_exit_f
	        elseif (ier .eq. ALL_OK) then
		  write(6,600)' Datatype for exponents is ',
     &               DataTypeName(datatype)
                  if (datatype .eq. RealSingle) then
	            call cg_exponents_read_f(data_single, ier)
		    if (ier .eq. ERROR) call cg_error_exit_f
		    write(6,110)' Exponents:',(data_single(n),n=1,5)
                  elseif (datatype .eq. RealDouble) then
                    call cg_exponents_read_f(data_double, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,110)' Exponents:',(data_double(n),n=1,5)
                  endif
	        endif

          ! *** DattaArray_t attribute: DataConversion_t
	        call cg_conversion_info_f(datatype, ier)
                if (ier .eq. ERROR) then
		  call cg_error_exit_f
                elseif (ier .eq. ALL_OK) then
		  write(6,600)' Datatype for conversion is ',
     &               DataTypeName(datatype)
                  if (datatype .eq. RealSingle) then
                    call cg_conversion_read_f(data_single, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,110)' Conversion:',(data_single(n),n=1,2)
                  elseif (datatype .eq. RealDouble) then
                    call cg_conversion_read_f(data_double, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,110)' Conversion:',(data_double(n),n=1,2)
                  endif
		endif

	      enddo	! loop through DataArray_t
	    enddo	! loop through IntegralData_t

            write(6,400)'                             *     *     *'

      ! *** zone coordinate attribute:  GOTO GridCoordinates_t node
            call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &            'GridCoordinates_t', one, 'end')
            if (ier .eq. ERROR) call cg_error_exit_f
	    if (ier .eq. ALL_OK) then

        ! *** GridCoordinates_t attribute: dimensional units
              call cg_units_read_f(mass, length, time, temp, deg, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
	      if (ier .eq. ALL_OK) write(6,400) 
     &            'Dimensional Units for GridCoordinates_t: ',
     &            LengthUnitsName(length)

	! *** GridCoordinates_t attribute:  Rind
              call cg_rind_read_f(rind, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              write(6,103)'GC Rind Data is ',(rind(i),i=1,6)

        ! *** coordinate array
	      call cg_narrays_f(narrays, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
              write(6,105) 'GridCoordinates_t contains ',
     &                            narrays,' arrays'
	      do iarray=1,narrays

		call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &            'GridCoordinates_t', one, 'end')
                if (ier .eq. ERROR) call cg_error_exit_f

	  ! *** GridCoordinates_t attribute: DataArray_t
                call cg_array_info_f(iarray, name, datatype,
     &                               nndim, dim_vals, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                write(6,600)' DataArrayName="',name,'"'
                write(6,600)' DataType=',DataTypeName(datatype)
                write(6,104)' DataNdim=',nndim
		do i=1,nndim
                  write(6,111)' DataDim(',i,')=',dim_vals(i)
		enddo

	  ! *** Compute nr of data in data array:
	   	num = 1
		do i=1,nndim
		  num = num*dim_vals(i)
		enddo

                if (datatype .eq. RealSingle) then
                  call cg_array_read_f(iarray, data_single, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  write(6,106) (data_single(i),i=1,2)
                  write(6,106) (data_single(i),i=num-1,num)
                elseif (datatype .eq. RealDouble) then
                  call cg_array_read_f(iarray, data_double, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  write(6,106) (data_double(i),i=1,2)
                  write(6,106) (data_double(i),i=num-1,num)
                endif

          ! *** coordinate attribute:  GOTO coordinate array node
                call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &          'GridCoordinates_t', one, 'DataArray_t', iarray, 'end')
                if (ier .eq. ERROR) call cg_error_exit_f

                call cg_ndescriptors_f(ndescriptors, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                write(6,105) 'No. of descriptors=',ndescriptors
                do idescr=1, ndescriptors
                  call cg_descriptor_read_f(idescr, name, text, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  write(6,500) ' DescriptorName="',name,'"',
     &                               ' DescriptorText="',text,'"'
                enddo

	      enddo	! loop through data arrays

        ! *** read coordinates using coordinate arrays' specific functions:

	      write(6,400)'Specific functions to read coordinates arrays'
	      call cg_ncoords_f(cg, base, zone, ncoords, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,103)'no. of coordinates=',ncoords

	 ! ** Compute the nr of data to be read
              do i=1,IndexDim
                rmin(i)=1
                DataSize(i)=size(i) + rind(2*i-1) + rind(2*i)
              enddo

	      do coord=1, ncoords
	  	call cg_coord_info_f(cg, base, zone, coord, datatype,
     &                               coordname, ier)
		if (ier .eq. ERROR) call cg_error_exit_f
 		write(6,112)'coord #',coord, 
     &            '   datatype=',DataTypeName(datatype), 
     &            '   name="',coordname,'"'

                if (datatype .eq. RealSingle) then
        	    call cg_coord_read_f(cg, base, zone, coordname,
     &                  RealSingle, rmin, DataSize, data_single, ier)
        	    if (ier .eq. ERROR) call cg_error_exit_f

                elseif (datatype .eq. RealDouble) then
                    call cg_coord_read_f(cg, base, zone, coordname,
     &                  RealDouble, rmin, DataSize, data_double, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                endif
	      enddo
	    endif 	! if GridCoordinates_t exists

            write(6,400)'                             *     *     *'

! *** solution
	   
	    call cg_nsols_f(cg, base, zone, nsols, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	    write(6,113) nsols,' FlowSolution_t node(s)',
     &              'found for ',zonename

      ! *** Read solution with general cg_array_read function
	    do sol=1, nsols
	      call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &            'FlowSolution_t', sol, 'end')
	      if (ier .eq. ERROR) call cg_error_exit_f

	! *** FlowSolution_t attribute:  DataArray_t
	      call cg_narrays_f(narrays, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
              write(6,108) ' FlowSolution_t #',sol,
     &              ' contains ',narrays,' solution arrays'

	! *** FlowSolution_t attribute:  GridLocation
	      call cg_gridlocation_read_f(location, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,600)'  The solution data are recorded at the ',
     &                 GridLocationName(location)

	! *** FlowSolution_t attribute:  Rind
	      call cg_rind_read_f(rind, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,103)'  The Rind Data is ',(rind(i),i=1,6)

              do iarray=1,narrays
	        call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &            'FlowSolution_t', sol, 'end')
                if (ier .eq. ERROR) call cg_error_exit_f

                call cg_array_info_f(iarray, name, datatype,
     &                               nndim, dim_vals, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
	        write(6,114) '  DataArray #',iarray
                write(6,600) '   Name="',name,'"'
                write(6,600) '   DataType=',DataTypeName(datatype)
                write(6,103) '   DataNdim=',nndim
                do i=1,nndim
                  write(6,111)'   DataDim(',i,')=',dim_vals(i)
                enddo

	 ! *** For dynamic memory allocation, compute the number of data to be read:
                num = 1
                do i=1,nndim
                    num = num*dim_vals(i)
                enddo
                write(6,200) 'Nr of data in solution vector=',num

                if (datatype .eq. RealSingle) then
                  call cg_array_read_f(iarray, data_single, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  !write(6,106) (data_single(i),i=1,num)
                elseif (datatype .eq. RealDouble) then
                  call cg_array_read_f(iarray, data_double, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  !write(6,106) (data_double(i),i=1,num)
                endif

          ! *** solution field attribute:  GOTO solution array node
                call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &          'FlowSolution_t',sol,'DataArray_t',iarray,'end')
                    if (ier .eq. ERROR) call cg_error_exit_f

          ! *** solution field attribute:  DimensionalUnits
                call cg_units_read_f(mass, length, time, temp,
     &                   deg, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                if (ier .eq. ALL_OK) then
                write(6,100)
     &		  '   Dimensional Units:',
     &            MassUnitsName(mass), LengthUnitsName(length),
     &            TemperatureUnitsName(temp), TimeUnitsName(time),
     &            AngleUnitsName(deg)
                endif

	      enddo	! loop through DataArray_t
	      write(6,103)' '

      ! *** Reading solution data with solution specific functions:
              call cg_sol_info_f(cg, base, zone, sol, solname,
     &                           location, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              write(6,115)'sol #',sol,':',
     &           '   solname="',solname,'"',
     &           '   location=',GridLocationName(location)

        ! *** Compute the nr of data to be read
	
              if (zonetype.eq.Structured) then
                do i=1,3
                  DataSize(i)=size(i) + rind(2*i-1) + rind(2*i)
                  if (location.eq.CellCenter) DataSize(i)=DataSize(i)-1
                enddo
              else
                DataSize(1)=size(2)
              endif

        ! *** solution field
              call cg_nfields_f(cg, base, zone, sol, nfields, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              write(6,105)'  nfields=',nfields

              do field=1, nfields
                call cg_field_info_f(cg, base, zone, sol, field,
     &                               type, fieldname, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                write(6,115)'  field #',field,':',
     &                    '   fieldname="',fieldname,'"',
     &                    '   datatype=',DataTypeName(type)

          ! *** read entire range of solution data and record in double precision
                call cg_field_read_f(cg, base, zone, sol, fieldname,
     &              RealDouble, rmin, DataSize, data_double, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
              enddo                             ! field loop

            enddo	! loop through FlowSolution_t

            write(6,400)'                             *     *     *'

! *** discrete data under zone
	    call cg_ndiscrete_f(cg, base, zone, ndiscrete, ier)
	    if (ier .eq. ERROR) call cg_error_exit_f
	    if (ier .eq. ALL_OK) write(6,113)ndiscrete, 
     &        ' DiscreteData_t node(s) found under ',zonename

	    do discr=1, ndiscrete
	      call cg_discrete_read_f(cg, base,zone, discr, name, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,600)' name=',name

        ! *** discrete data attribute:  GOTO DiscreteData_t node
	      call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &                'DiscreteData_t',  discr, 'end')
	      if (ier .eq. ERROR) call cg_error_exit_f

        ! *** discrete data attribute:  GridLocation_t
	      call cg_gridlocation_read_f(location, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      if (ier .eq. ALL_OK) write(6,600)
     &           ' The location of the DiscreteData vector is ',
     &                 GridLocationName(location)

  	! *** discrete data arrays:
              call cg_narrays_f(narrays, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              write(6,116) ' DiscreteData #', discr,
     &                         ' contains ', narrays,' arrays'
	      do iarray=1, narrays
	        call cg_array_info_f(iarray, name, datatype,
     &                              nndim, dim_vals, ier)
                if (ier .eq. ERROR) call cg_error_exit_f

		write(6,116) 'DataArray #',iarray,':'
	        write(6,600)'  Name =',name
		write(6,600)'  Datatype=',
     &            DataTypeName(datatype)

	  ! *** compute nr of data to be read
		num=1
		do n=1, nndim
		  num=num*dim_vals(n)
		enddo

                if (datatype .eq. RealSingle) then
                  call cg_array_read_f(iarray, data_single, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  !write(6,*) (data_single(n),n=1,num)
                elseif (datatype .eq. RealDouble) then
                  call cg_array_read_f(iarray, data_double, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  !write(6,*) (data_double(n),n=1,num)
                endif

          ! *** discrete data arrays attribute: GOTO DataArray node
                call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &            'DiscreteData_t', discr, 'DataArray_t', iarray, 'end')
	        if (ier .eq. ERROR) call cg_error_exit_f

	        call cg_units_read_f(mass, length, time, temp, deg, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
		if (ier .eq. ALL_OK) then
                  write(6,100)
     &		    '  Dimensional Units for DiscreteData_t:',
     &              MassUnitsName(mass), LengthUnitsName(length),
     &              TemperatureUnitsName(temp), TimeUnitsName(time),
     &              AngleUnitsName(deg)
                endif
	      enddo		! loop through DataArray_t
	    enddo

            write(6,400)'                             *     *     *'

! *** Interblock Connectivity:
	    write(6,401)'Interblock Connectivity for ',zonename

! *** ZoneGridConnectivity attributes:  GOTO ZoneGridConnectivity_t node
            call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &         'ZoneGridConnectivity_t', one, 'end')
            if (ier .eq. ERROR) call cg_error_exit_f

	    if (ier.eq. ALL_OK) then
! *** ZoneGridConnectivity attributes: Descriptor_t
              call cg_ndescriptors_f(ndescriptors, ier)
              if (ier .ne. 0) call cg_error_exit_f
              write(6,117)
     &        ndescriptors, ' descriptors for ZoneGridConnectivity_t'
              do idescr=1, ndescriptors
                call cg_descriptor_read_f(idescr, name, text, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                write(6,500) '     DescriptorName="',name,'"',
     &                             '     DescriptorText="',text,'"'
              enddo


        ! *** overset holes
	      call cg_nholes_f(cg, base, zone, nholes, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
	      write(6,107) nholes, ' holes found'

	      do hole=1, nholes
	        call cg_hole_info_f(cg, base, zone, hole, holename,
     &             location, ptset_type, nptsets, npnts, ier)
	        if (ier .eq. ERROR) call cg_error_exit_f
	        write(6,118)
     &	      	  '  hole #',hole,':', '   holename="',holename,'"',
     &            '   data location=',GridLocationName(location),
     &            '   nptsets = ',nptsets, 
     &            ', total no. of points =',npnts

	        if (npnts .lt. 30000) then
	 	  call cg_hole_read_f(cg, base, zone, hole, pnts, ier)
	  	  if (ier .eq. ERROR) call cg_error_exit_f
	        endif

          ! *** overset holes attributes:  GOTO OversetHoles_t node
	        call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &             'ZoneGridConnectivity_t', one,
     &             'OversetHoles_t', hole, 'end')
                if (ier .ne. 0) call cg_error_exit_f

	        call cg_ndescriptors_f(ndescriptors, ier)
	        if (ier .ne. 0) call cg_error_exit_f
	        write(6,117)
     &            ndescriptors, ' descriptors for ',holename
                do idescr=1, ndescriptors
                  call cg_descriptor_read_f(idescr, name, text, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f
                  write(6,500) '     DescriptorName="',name,'"',
     &                               '     DescriptorText="',text,'"'
                enddo
	      enddo	!hole loop
		


! *** general connectivity
	      call cg_nconns_f(cg, base, zone, nconns, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              write(6,107) nconns,' GridConnectivity_t found'

	      do conn=1, nconns
                call cg_conn_info_f(cg, base, zone, conn, connectname, 
     &             location, type, ptset_type, npnts, donorname,
     &             nzonetype, ndonor_ptset_type, ndonor_data_type,
     &             npnts_donor, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
	        write(6, 101)
!234567890!234567890!234567890!234567890!234567890!234567890!23456789012
     &            '  GridConnectivity #',conn,':',
     &            '   connect name=',connectname,
     &            '   Grid location=',GridLocationName(location),
     &            '   Connect-type=',GridConnectivityTypeName(type),
     &            '   ptset type="',PointSetTypeName(ptset_type),'"',
     &            '   npnts=',npnts,'   donorname="',donorname,'"',
     &            '   donor zonetype=',ZoneTypeName(nzonetype),
     &            '   donor ptset type=',
     &                PointSetTypeName(ndonor_ptset_type),
     &            '   npnts_donor=',npnts_donor

                call cg_conn_read_f(cg, base, zone, conn, pnts,
     &                            Integer, donor_pnts, ier)
                if (ier .eq. ERROR) call cg_error_exit_f

	        write(6,119) '   Current zone:',
     &           '    first point:', pnts(1),pnts(2),pnts(3),
     &           '    last point :', pnts(3*npnts-2), pnts(3*npnts-1), 
     &                               pnts(3*npnts)
	        write(6,119) '   Donor zone:',
     &           '    first point:', donor_pnts(1),donor_pnts(2),
     &                               donor_pnts(3),
     &           '    last point :', donor_pnts(3*npnts-2), 
     &                               donor_pnts(3*npnts-1), 
     &                               donor_pnts(3*npnts)

! *** general connectivity attributes:  GOTO GridConnectivity_t node
	        call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &              'ZoneGridConnectivity_t', one,
     &              'GridConnectivity_t', conn, 'end')
	        if (ier .eq. ERROR) call cg_error_exit_f

	        call cg_ordinal_read_f(num, ier)
	        if (ier .eq. ERROR) call cg_error_exit_f
	        if (ier .eq. ALL_OK) write(6,200)'  Ordinal=',num
              enddo

! *** connectivity 1to1
	      call cg_n1to1_f(cg, base, zone, n1to1, ier)
	      if (ier .eq. ERROR) call cg_error_exit_f
              write(6,107) n1to1,' GridConnectivity1to1_t found'

              do one21=1, n1to1
                call cg_1to1_read_f(cg, base, zone, one21, connectname,
     &             donorname, range, donor_range, transform, ier)
                if (ier .eq. ERROR) call cg_error_exit_f

                write(6,105) 'GridConnectivity1to1 #',one21
                write(6,600) 'connectname="',connectname,'"'
                write(6,600) 'donorname  ="',donorname,'"'

                write(6,120) ' range: ',
     &           '(',range(1,1),',',range(2,1),',',range(3,1),
     &           ') to (',range(1,2),',',range(2,2),',',range(3,2),')'

                write(6,121)' donor_range: ',
     &          '(', donor_range(1,1), ',', donor_range(2,1), ',',
     &            donor_range(3,1), ') to (',
     &            donor_range(1,2), ',', donor_range(2,2), ',',
     &            donor_range(3,2), ')'

                write(6,122) ' Transform: ', '(',
     &            transform(1), ',',
     &            transform(2), ',', transform(3), ')'

	
! *** connectivity 1to1 attributes:  GOTO GridConnectivity1to1_t node
	        call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &              'ZoneGridConnectivity_t', one,
     &           'GridConnectivity1to1_t', one21, 'end')
		  if (ier .eq. ERROR) call cg_error_exit_f
		  if (ier .eq. ALL_OK) then

! *** connectivity 1to1 attributes:  Descriptor_t
	            call cg_ndescriptors_f(ndescriptors, ier)
                    if (ier .ne. 0) call cg_error_exit_f
                    write(6,117)
     &              ndescriptors, ' descriptors for ',connectname
                    do idescr=1, ndescriptors
                      call cg_descriptor_read_f(idescr, name, text, ier)
                      if (ier .eq. ERROR) call cg_error_exit_f
                      write(6,500) '   DescriptorName="',name,'"',
!234567890!234567890!234567890!234567890!234567890!234567890!23456789012
     &                               '   DescriptorText="',text,'"'
                    enddo
		endif
              enddo
	    endif	! if ZoneGridConnectivity exists

            write(6,400)'                             *     *     *'

! *** bocos
	    write(6,600)'Boundary Conditions for ',zonename


! *** Zone bound. condition attributes: GOTO ZoneBC_t node
            call cg_goto_f(cg, base,ier, 'Zone_t', zone,
     &                     'ZoneBC_t', one, 'end')
            if (ier .eq. ERROR) call cg_error_exit_f
	    if (ier .eq. ALL_OK) then

  	! *** Zone bound. condition attributes: ReferenceState_t
              call cg_state_read_f(StateDescription, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
	      if (ier.eq.ALL_OK) then
		write(6,600)' ReferenceState defined under ZoneBC_t'
                write(6,600)'  StateDescription=',StateDescription

         ! ** ReferenceState_t attributes:  GOTO ReferenceState_t
                call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &            'ZoneBC_t', one, 'ReferenceState_t', one, 'end')
                if (ier .eq. ERROR) call cg_error_exit_f

                call cg_narrays_f(narrays, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                write(6,105) '  ReferenceState_t contains ',
     &                          narrays,' array(s)'

                do iarray=1, narrays

                  call cg_array_info_f(iarray, name, datatype,
     &              nndim, dim_vals, ier)
                  if (ier .eq. ERROR) call cg_error_exit_f

		  write(6,105) '   DataArray #',iarray,':'
                  write(6,600)'    Name =',name
                  write(6,600)'    Datatype=',DataTypeName(datatype)

		  write(6,600)'    Data:'
                  if (datatype .eq. RealSingle) then
                    call cg_array_read_f(iarray, data_single, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,124) data_single(1)
                  elseif (datatype .eq. RealDouble) then
                    call cg_array_read_f(iarray, data_double, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,124) data_double(1)
                  endif
                enddo


           ! ** ReferenceState_t attributes: DimensionalUnits_t
                call cg_units_read_f(mass, length, time, temp,
     &                   deg, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                if (ier .eq. ALL_OK) then
                  write(6,100)
     &		    '  Dimensional Units:',
     &              MassUnitsName(mass), LengthUnitsName(length),
     &              TemperatureUnitsName(temp), TimeUnitsName(time),
     &              AngleUnitsName(deg)
                endif
	      endif	!if ReferenceState exists under ZoneBC_t

	      call cg_nbocos_f(cg, base, zone, nbocos, ier)
              if (ier .eq. ERROR) call cg_error_exit_f
              write(6,113)nbocos,' bound. conditions found for ',
     &                   zonename

	      do boco=1, nbocos
	        call cg_boco_info_f(cg, base, zone, boco, boconame, 
     &             bocotype, ptset_type, npnts, 
     &             NormalIndex, NormalListFlag, datatype,
     &             ndataset, ier)
	   	if (ier .eq. ERROR) call cg_error_exit_f
           	write(6,105) ' boundary condition #',boco
		write(6,600) '  boconame=',boconame
	  	write(6,600) '  bocotype=',BCTypeName(bocotype)
		write(6,600) '  ptset_type=',
     &            PointSetTypeName(ptset_type)
		write(6,103) '  NormalIndex=',
     &            NormalIndex(1),NormalIndex(2), NormalIndex(3)
		write(6,104) '  NormalListFlag=',NormalListFlag
		write(6,600) '  datatype for normals=',
     &            DataTypeName(datatype)

	       ! read patch points and InwardNormalList
	        if (datatype.eq.RealSingle) then
		   call cg_boco_read_f(cg, base, zone, boco, pnts, 
     &              data_single, ier)
		   if (ier .eq. ERROR) call cg_error_exit_f
	        elseif (datatype.eq.RealDouble) then
		   call cg_boco_read_f(cg, base, zone, boco, pnts,
     &              data_double, ier)
		   if (ier .eq. ERROR) call cg_error_exit_f
		endif

                write(6,119) '   Bound. Condition Patch:',
     &           '    first point:', pnts(1),pnts(2),pnts(3),
     &           '    last point :', pnts(3*npnts-2), pnts(3*npnts-1),
     &                               pnts(3*npnts)

	 	if (NormalListFlag .ne. 0) then
		  if (datatype.eq.RealSingle) 
     &              write(6,126) '   Normals:',
     &           '    first point:', data_single(1),data_single(2),
     &                               data_single(3),
     &           '    last point :', data_single(3*npnts-2), 
     &                               data_single(3*npnts-1),
     &                               data_single(3*npnts)
		  if (datatype.eq.RealDouble)
     &              write(6,126) '   Normals:',
     &           '    first point:', data_double(1),data_double(2),
     &                               data_double(3),
     &           '    last point :', data_double(3*npnts-2), 
     &                               data_double(3*npnts-1),
     &                               data_double(3*npnts)
		endif
         ! ***  bound. condition attributes: GOTO BC_t node
		call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &              'ZoneBC_t', one, 'BC_t', boco, 'end')
		if (ier .eq. ERROR) call cg_error_exit_f

         ! ***  bound. condition attributes: DataClass_t
	        call cg_dataclass_read_f(type,ier)
                if (ier .eq. ERROR) call cg_error_exit_f
                if (ier.eq.ALL_OK) 
     &            write(6,600)'  B.C. DataClass=',
     &                           DataClassName(type)

         ! ***  boundary condition attributes:  GridLocation_t
                call cg_gridlocation_read_f(location, ier)
                if (ier .eq. ERROR) call cg_error_exit_f
		if (ier.eq.ALL_OK)
     &            write(6,600)'    data location=',
     &                   GridLocationName(location)

           ! ** boundary condition dataset
	        write(6,103) '  ndataset=',ndataset
		do idataset=1, ndataset
                  call cg_dataset_read_f(cg, base, zone, boco,idataset,
     &            name, type, DirichletFlag, NeumannFlag, ier)
		  if (ier .eq. ERROR) call cg_error_exit_f
		    
		  write(6,103)'   Dataset #',idataset
		  write(6,600)'    Name=',name
		  write(6,600)'    BCType=',BCTypeName(type)

             ! ** boundary condition data:  GOTO BCData_t node
		  if (DirichletFlag.eq.1) then
		    call cg_goto_f(cg, base, ier, 'Zone_t', zone,
     &                'ZoneBC_t', one, 'BC_t', boco, 'BCDataSet_t',
     &                idataset,'BCData_t',Dirichlet,'end')
		      if (ier .eq. ERROR) call cg_error_exit_f

               ! ** boundary condition data attributes: DataClass_t
		    write(6,401)'   Dirichlet DataSet:'
		    call cg_dataclass_read_f(type,ier)
		    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,600)'    DataClass=',
     &                              DataClassName(type)

	       ! ** boundary condition data attributes: DataArray_t
		    call cg_narrays_f(narrays, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,127) '    DirichletData', 
     &                         ' contains ', narrays,' data arrays'
                    do iarray=1, narrays
                      call cg_array_info_f(iarray, name, datatype,
     &                              nndim, dim_vals, ier)
                      if (ier .eq. ERROR) call cg_error_exit_f

                      write(6,105) '    DataArray #',iarray,':'
                      write(6,600)'     Name =',name
                      write(6,600)'     Datatype=',
     &                  DataTypeName(datatype)

		      write(6,105)'    Dirichlet Data:'
                      if (datatype .eq. RealSingle) then
                        call cg_array_read_f(iarray, data_single, ier)
                        if (ier .eq. ERROR) call cg_error_exit_f
                        write(6,106) 
     &                      (data_single(n),n=1,dim_vals(1))

                      elseif (datatype .eq. RealDouble) then
                        call cg_array_read_f(iarray, data_double, ier)
                        if (ier .eq. ERROR) call cg_error_exit_f
                        write(6,106) 
     &                      (data_double(n),n=1,dim_vals(1))
                      endif
		    enddo
		  endif

		  if (NeumannFlag.eq.1) then
		    call cg_goto_f(cg, base, ier, 'Zone_t', zone, 
     &                  'ZoneBC_t', one, 'BC_t', boco, 'BCDataSet_t',
     &                  idataset, 'BCData_t', Neumann,'end')
                      if (ier .eq. ERROR) call cg_error_exit_f

               ! ** boundary condition data attributes: DataClass_t
                    call cg_dataclass_read_f(type,ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,600)'    DataClass=',
     &                              DataClassName(type)

               ! ** boundary condition data attributes: DataArray_t
                    call cg_narrays_f(narrays, ier)
                    if (ier .eq. ERROR) call cg_error_exit_f
                    write(6,105) 
     &              '    Neumann Data contains ', narrays,' data arrays'
                    do iarray=1, narrays
                      call cg_array_info_f(iarray, name, datatype,
     &                              nndim, dim_vals, ier)
                      if (ier .eq. ERROR) call cg_error_exit_f

                      write(6,105) '    DataArray #',iarray,':'
                      write(6,600)'     Name =',name
                      write(6,600)'     Datatype=',
     &                  DataTypeName(datatype)

                      write(6,400)'    Neumann Data:'
                      if (datatype .eq. RealSingle) then
                        call cg_array_read_f(iarray, data_single, ier)
                        if (ier .eq. ERROR) call cg_error_exit_f
                        write(6,106) 
     &                      (data_single(n),n=1,dim_vals(1))

                      elseif (datatype .eq. RealDouble) then
                        call cg_array_read_f(iarray, data_double, ier)
                        if (ier .eq. ERROR) call cg_error_exit_f
                        write(6,106) 
     &                      (data_double(n),n=1,num)
                      endif

                    enddo	! loop through DataArray
		  endif		! if Neumann
		enddo		! loop through dataset
              enddo		! loop through boco		
	    endif		! if ZoneBC_t exists
 	  enddo			! zone loop

          write(6,400)'                             *     *     *'

! *** connectivity 1to1 - Global
          write(6,600)' Reading 1to1 connectivity for entire Base'
          call cg_n1to1_global_f(cg, base, n1to1_global, ier)
          if (ier .eq. ERROR) call cg_error_exit_f
          write(6,200)'n1to1_global=',n1to1_global

	  if (n1to1_global .gt. 0) then
	  call cg_1to1_read_global_f(cg, base, 
     &       G_connectname, G_zonename, G_donorname,
     &       G_range, G_donor_range, G_transform, ier)
	  if (ier .eq. ERROR) call cg_error_exit_f

	  do i=1, n1to1_global
	    write(6,600) ' '
	    write(6,130) '*** interface #',i,' ***'
	    write(6,600) 'G_connectname="',G_connectname(i),'"'
	    write(6,600) 'G_zonename   ="',G_zonename(i),'"'
	    write(6,600) 'G_donorname  ="',G_donorname(i),'"'

	    write(6,131) 'G_range: ', 
     &           '(',G_range(1,i),',',G_range(2,i),',',G_range(3,i), 
     &      ') to (',G_range(4,i),',',G_range(5,i),',',G_range(6,i),')'

            write(6,132) 'G_donor_range: ',
     &     '(', G_donor_range(1,i), ',', G_donor_range(2,i), ',',
     &          G_donor_range(3,i), ') to (', 
     &          G_donor_range(4,i), ',', G_donor_range(5,i), ',',
     &          G_donor_range(6,i), ')'

	   write(6,133) 'Transform: ', '(', 
     &          G_transform(1,i), ',',
     &          G_transform(2,i), ',', G_transform(3,i), ')'
	  enddo
	  endif


 	enddo    				! loop through bases

        write(6,400)'                             *     *     *'

        call cg_close_f(cg, ier)
        if (ier .eq. ERROR) call cg_error_exit_f

 100 	format(a/,'    Mass units: ',a/,'    Length units: ',a/,
     &    '    Temperature units: ',a/,'    Time units: ',a/,
     &    '    Angle units:',a)
 101	format(a,i1,a,/2a,/2a,/2a,/3a,/a,i4,3a,/2a,/2a,/2a,/a,i4)
 102 	format(a,f5.3)
 103	format(a,6i2)
 104	format(a,i5,3a)
 105	format(a,i2,a)
 106    format(6f10.3)
 107	format(i2,2a)
 108    format(a,i2,a,i2,a)
 109	format(a,f5.1)
 110	format(a,5f5.1)
 111	format(a,i1,a,i8)
 112	format(a,i1/2a/3a)
 113	format(i1,3a)
 114	format(/a, i1)
 115	format(a,i1,a/3a/2a)
 116	format(a,i1,a,i1,a)
 117	format(/i4,2a)
 118	format(a,i1,a/3a/2a/a,i1,a,i5)
 119	format(a/a,3i2/a,3i2)
 120	format(a10, 3(a1,i1),a6,3(i1,a1))
 121 	format(a16,3(a1,i1),a6,3(i1,a1))
 122	format(a12,3(a1,i2),a1)
 124	format(4x, f7.2)
 126	format(a/a,3f5.2/a,3f5.2)
 127	format(2a,i1,a)
 130	format(a15, i2, a4)
 131	format(a10, 3(a1,i1),a6,3(i1,a1))
 132	format(a16,3(a1,i1),a6,3(i1,a1))
 133	format(a12,3(a1,i2),a1)
 200    format(a,i5)
 300	format(3a/a,i2)
 400	format(/a/)
 401	format(/2a/)
 500	format(3a/3a)
 600	format(3a)

 9999	end
