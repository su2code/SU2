      program test64f
      implicit none
#ifdef WINNT
      include 'cgnswin_f.h'
#endif
      include 'cgnslib_f.h'

      integer*8 nnodes,nelems
      integer dotest

      parameter (nnodes = 40000000)
      parameter (nelems = 40000000)
      parameter (dotest = 1)

      integer*8 cgfile,cgbase,cgzone,cgcoord,cgsect
      integer*8 cdim,pdim,i,ier,type,nb,pflag
      integer*8 n,ns,ne,nerrs
      integer*8 size(3),elems(4,nelems)
      real*4 nodes(nnodes)

      if (CG_BUILD_64BIT .eq. 0) then
        print *, 'must be built in 64 bit mode'
        stop
      endif
      
      call cg_open_f('test64.cgns',CG_MODE_WRITE,cgfile,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      cdim = 3
      pdim = 3
      call cg_base_write_f(cgfile,'Base',cdim,pdim,cgbase,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      size(1) = nnodes
      size(2) = nelems
      size(3) = 0
      call cg_zone_write_f(cgfile,cgbase,'Zone',size,
     &                     Unstructured, cgzone, ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      print *,'writing',nnodes,' coordinates'
      do n=1,nnodes
        nodes(n) = n
      enddo
      call cg_coord_write_f(cgfile,cgbase,cgzone,RealSingle,
     &                      'Coordinates',nodes,cgcoord,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      print *,'writing',nelems,' tetrahedra'
      do n=1,nelems
        do i=1,4
          elems(i,n) = n
        enddo
      enddo
      ns = 1
      ne = nelems
      nb = 0
      call cg_section_write_f(cgfile,cgbase,cgzone,'Elements',
     &                        TETRA_4,ns,ne,nb,elems,cgsect,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      call cg_close_f(cgfile,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      if (dotest.eq.0) goto 100

      call cg_open_f('test64.cgns',CG_MODE_READ,cgfile,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f
      
      cgbase = 1
      cgzone = 1

      print *,'reading/checking',nnodes,' coordinates'
      ns = 1
      ne = nnodes
      call cg_coord_read_f(cgfile,cgbase,cgzone,'Coordinates',
     &                     RealSingle,ns,ne,nodes,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f
      nerrs = 0
      do n=1,nnodes
        if (nodes(n).ne.n) nerrs = nerrs + 1
      enddo
      print *,nerrs,' errors'

      print *,'reading/checking',nelems,' tetrahedra'
      call cg_elements_read_f(cgfile,cgbase,cgzone,cgsect,
     &                        elems,0,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f
      nerrs = 0
      do n=1,nelems
        do i=1,4
          if (elems(i,n).ne.n) nerrs = nerrs + 1
        enddo
      enddo
      print *,nerrs,' errors'

      call cg_close_f(cgfile,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

  100 continue
      end
