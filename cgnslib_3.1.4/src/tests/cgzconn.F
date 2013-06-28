      program cgzconn

#include "cgnstypes_f.h"
#ifdef WINNT
      include 'cgnswin_f.h'
#endif
      include 'cgnslib_f.h'

      cgsize_t ierr, cgfile, cgbase, cgzone, cgcoord
      cgsize_t cgz1, cgz2, cgconn, nzconn, nconn, n1to1
      cgsize_t i, j, k, n, size(9), dim
      cgsize_t ptrange(6), ptlist(125), transform(3)
      cgsize_t loc, type, ptype, npts, dztype, dptype, ddtype, dnpts
      real*4 x(125), y(125), z(125)
      character*32 zname, dname, cname

c open and write base
      dim = 3
      call cg_open_f('cgzconn.cgns', CG_MODE_WRITE, cgfile, ierr)
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
        transform(n) = n
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

c loop over zones
      do n=1,2
        if (n .eq. 1) then
          zname = 'Zone1'
          dname = 'Zone2'
        else
          zname = 'Zone2'
          dname = 'Zone1'
        endif
c write zone
        call cg_zone_write_f(cgfile, cgbase, zname, size,
     &                       Structured, cgzone, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
c write coordinates
        call cg_coord_write_f(cgfile, cgbase, cgzone, RealSingle,
     &                        'CoordinateX', x, cgcoord, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        call cg_coord_write_f(cgfile, cgbase, cgzone, RealSingle,
     &                        'CoordinateY', y, cgcoord, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        call cg_coord_write_f(cgfile, cgbase, cgzone, RealSingle,
     &                        'CoordinateZ', z, cgcoord, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f

c write first ZoneGridConnectivity - will be active
        call cg_zconn_write_f(cgfile, cgbase, cgzone,
     &                        'ZoneConnectivity1', cgz1, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f

c write second ZoneGridConnectivity - will be active
        call cg_zconn_write_f(cgfile, cgbase, cgzone,
     &                        'ZoneConnectivity2', cgz2, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
c write general connectivity
        npts = 2
        dnpts = 25
        call cg_conn_write_f(cgfile, cgbase, cgzone, 'conn', Vertex,
     &                       Abutting1to1, PointRange, npts, ptrange,
     &                       dname, Structured, PointListDonor,
     &                       Integer, dnpts, ptlist, cgconn, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
c set it back to previous ZoneGridConnectivity and write 1to1
        call cg_zconn_set_f(cgfile, cgbase, cgzone, cgz1, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        call cg_1to1_write_f(cgfile, cgbase, cgzone, '1to1', dname,
     &                       ptrange, ptrange, transform, cgconn, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
      enddo

c close the file and reopen in read mode
      call cg_close_f(cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

      call cg_open_f('cgzconn.cgns', CG_MODE_READ, cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f
      cgbase = 1
      cgz1 = 1
      cgz2 = 2
      cgconn = 1
      do cgzone=1,2
        call cg_nzconns_f(cgfile, cgbase, cgzone, nzconn, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (nzconn .ne. 2) then
          print *,'nzconn != 2'
          stop
        endif
c read should make ZoneGridConnectivity active
        call cg_zconn_read_f(cgfile, cgbase, cgzone, cgz2, zname, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (zname .ne. 'ZoneConnectivity2') then
          print *,'expecting Zoneconnectivity2 - got',zname
          stop
        endif
        call cg_nconns_f(cgfile, cgbase, cgzone, nconn, ierr);
        if (ierr .ne. CG_OK) call cg_error_exit_f
        call cg_n1to1_f(cgfile, cgbase, cgzone, n1to1, ierr);
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (nconn .ne. 1 .or. n1to1 .ne. 0) then
          print *,'expecting nconn=1,n1to1=0 - got',nconn,n1to1
          stop
        endif
        call cg_conn_info_f(cgfile, cgbase, cgzone, cgconn, cname,
     &                      loc, type, ptype, npts, dname, dztype,
     &                      dptype, ddtype, dnpts, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (cname .ne. 'conn' .or. loc .ne. Vertex .or.
     &      type .ne. Abutting1to1 .or. ptype .ne. PointRange .or.
     &      npts .ne. 2 .or. dztype .ne. Structured .or.
     &      dptype .ne. PointListDonor .or. dnpts .ne. 25) then
          print *,'invalid conn data'
          stop
        endif

c read should make ZoneGridConnectivity active
        call cg_zconn_read_f(cgfile, cgbase, cgzone, cgz1, zname, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (zname .ne. 'ZoneConnectivity1') then
          print *,'expecting Zoneconnectivity1 - got',zname
          stop
        endif
        call cg_nconns_f(cgfile, cgbase, cgzone, nconn, ierr);
        if (ierr .ne. CG_OK) call cg_error_exit_f
        call cg_n1to1_f(cgfile, cgbase, cgzone, n1to1, ierr);
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (nconn .ne. 0 .or. n1to1 .ne. 1) then
          print *,'expecting nconn=0,n1to1=1 - got',nconn,n1to1
          stop
        endif
        call cg_1to1_read_f(cgfile, cgbase, cgzone, cgconn, cname,
     &                      dname, ptrange, ptrange, transform, ierr)
        if (ierr .ne. CG_OK) call cg_error_exit_f
        if (cname .ne. '1to1') then
          print *,'invalid 1to1 data'
          stop
        endif
      enddo

      call cg_close_f(cgfile, ierr)
      if (ierr .ne. CG_OK) call cg_error_exit_f

      end
