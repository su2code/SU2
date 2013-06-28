      program cgsubreg

#include "cgnstypes_f.h"
#ifdef WINNT
      include 'cgnswin_f.h'
#endif
      include 'cgnslib_f.h'

      cgsize_t ierr, cgfile, cgbase, cgzone, cgcoord
      cgsize_t csub, nsub, dim
      cgsize_t i, j, k, n, size(9)
      cgsize_t ptrange(6), ptlist(125)
      cgsize_t loc, ptype, npts, bclen, gclen
      real*4 x(125), y(125), z(125)
      character*32 sname, cname

c open and write base
      dim = 3
      call cg_open_f('cgsubreg.cgns', CG_MODE_WRITE, cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f
      call cg_base_write_f(cgfile, 'Base', dim, dim, cgbase, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

c create some bogus data
      do n=1,3
        size(n)   = 5
        size(n+3) = 4
        size(n+6) = 0
      enddo
      n = 0
      do k=1,5
        do j=1,5
          do i=1,5
            n = n + 1
            x(n) = i
            y(n) = j
            z(n) = k
          enddo
        enddo
      enddo
      do n=1,3
        ptrange(n) = 1
        ptrange(n+3) = 5
      enddo
      ptrange(4) = 1
      n = 0
      do j=1,5
        do i=1,5
          n = n + 1
          ptlist(n) = i
          n = n + 1
          ptlist(n) = j
          n = n + 1
          ptlist(n) = 1
        enddo
      enddo

c write zone
      call cg_zone_write_f(cgfile, cgbase, 'Zone', size,
     &                     Structured, cgzone, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f
c write coordinates
      call cg_coord_write_f(cgfile, cgbase, cgzone, RealSingle,
     &                      'CoordinateX', x, cgcoord, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f
      call cg_coord_write_f(cgfile, cgbase, cgzone, RealSingle,
     &                      'CoordinateY', y, cgcoord, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f
      call cg_coord_write_f(cgfile, cgbase, cgzone, RealSingle,
     &                      'CoordinateZ', z, cgcoord, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

c write first ZoneSubRegion as BCRegionName
      dim = 2
      call cg_subreg_bcname_write_f(cgfile, cgbase, cgzone, 'reg1',
     &                              dim, 'bcname', csub, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

c write second ZoneSubRegion as GridConnectivityRegionName
      call cg_subreg_gcname_write_f(cgfile, cgbase, cgzone, 'reg2',
     &                              dim, 'gcname', csub, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

c write third ZoneSubRegion as PointRange
      npts = 2
      call cg_subreg_ptset_write_f(cgfile, cgbase, cgzone, 'reg3',
     &                             dim, Vertex, PointRange, npts,
     &                             ptrange, csub, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

c write fourth ZoneSubRegion as PointList
      npts = 25
      call cg_subreg_ptset_write_f(cgfile, cgbase, cgzone, 'reg4',
     &                             dim, EdgeCenter, PointList, npts,
     &                             ptlist, csub, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

c close the file and reopen in read mode
      call cg_close_f(cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

      call cg_open_f('cgsubreg.cgns', CG_MODE_READ, cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

      cgbase = 1
      cgzone = 1
      call cg_nsubregs_f(cgfile, cgbase, cgzone, nsub, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f
      if (nsub .ne. 4) then
        print *,'expecting nsub=4 - got',nsub
        stop
      endif

      do csub=1,4
        call cg_subreg_info_f(cgfile, cgbase, cgzone, csub, sname, dim,
     &                        loc, ptype, npts, bclen, gclen, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f

        if (sname .eq. 'reg1') then
          if (bclen .ne. 6 .or. gclen .ne. 0 .or.
     &        ptype .ne. CG_NULL .or. npts .ne. 0) then
            print *,'bad info for',sname
            stop
          endif
          call cg_subreg_bcname_read_f(cgfile, cgbase, cgzone, csub,
     &                                 cname, ierr)
          if (ierr .ne. CG_OK) call cg_error_exit_f
          if (cname .ne. 'bcname') then
            print *,'expecting bcname - got ',cname
            stop
          endif

        elseif (sname .eq. 'reg2') then
          if (bclen .ne. 0 .or. gclen .ne. 6 .or.
     &        ptype .ne. CG_NULL .or. npts .ne. 0) then
            print *,'bad info for',sname
            stop
          endif
          call cg_subreg_gcname_read_f(cgfile, cgbase, cgzone, csub,
     &                                 cname, ierr)
          if (ierr .ne. CG_OK) call cg_error_exit_f
          if (cname .ne. 'gcname') then
            print *,'expecting gcname - got ',cname
            stop
          endif

        elseif (sname .eq. 'reg3') then
          if (bclen .ne. 0 .or. gclen .ne. 0 .or.
     &        ptype .ne. PointRange .or. npts .ne. 2 .or.
     &        loc .ne. Vertex) then
            print *,'bad info for',sname
            stop
          endif
          call cg_subreg_ptset_read_f(cgfile, cgbase, cgzone, csub,
     &                                ptrange, ierr)
          if (ierr .ne. CG_OK) call cg_error_exit_f

        elseif (sname .eq. 'reg4') then
          if (bclen .ne. 0 .or. gclen .ne. 0 .or.
     &        ptype .ne. PointList .or. npts .ne. 25 .or.
     &        loc .ne. EdgeCenter) then
            print *,'bad info for',sname
            stop
          endif
          call cg_subreg_ptset_read_f(cgfile, cgbase, cgzone, csub,
     &                                ptlist, ierr)
          if (ierr .ne. CG_OK) call cg_error_exit_f

        else
          print *,'unknown subregion ',sname
          stop
        endif
      enddo

      call cg_close_f(cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

      end
