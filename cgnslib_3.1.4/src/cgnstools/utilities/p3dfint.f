      subroutine openf (iread,fname)
      integer iread
      character *(*) fname
      if (iread .eq. 0) then
        open (2,file=fname,form='unformatted',status='new')
      else
        open (2,file=fname,form='unformatted',status='old')
      endif
      end
c
      subroutine closef
      close (2)
      end
c
      subroutine readif (icnt,idata,ierr)
      integer icnt,idata(*)
      read (2,end=10) (idata(i),i=1,icnt)
      ierr = 0
      return
   10 ierr = 1
      end
c
      subroutine readff (icnt,rdata,ierr)
      integer icnt
      real*4 rdata(*)
      read (2,end=10) (rdata(i),i=1,icnt)
      ierr = 0
      return
   10 ierr = 1
      end
c
      subroutine readdf (icnt,rdata,ierr)
      integer icnt
      real*8 rdata(*)
      read (2,end=10) (rdata(i),i=1,icnt)
      ierr = 0
      return
   10 ierr = 1
      end
c
      subroutine readgff (icnt,rdata,idata,ierr)
      integer icnt,idata(*),ierr
      real*4 rdata(*)
      read (2,end=10) (rdata(i),i=1,3*icnt),(idata(i),i=1,icnt)
      ierr = 0
      return
   10 ierr = 1
      end
c
      subroutine readgdf (icnt,rdata,idata,ierr)
      integer icnt,idata(*),ierr
      real*8 rdata(*)
      read (2,end=10) (rdata(i),i=1,3*icnt),(idata(i),i=1,icnt)
      ierr = 0
      return
   10 ierr = 1
      end
c
      subroutine writeif (icnt,idata,ierr)
      integer icnt,idata(*)
      write (2) (idata(i),i=1,icnt)
      ierr = 0
      end
c
      subroutine writeff (icnt,rdata,ierr)
      integer icnt
      real*4 rdata(*)
      write (2) (rdata(i),i=1,icnt)
      ierr = 0
      end
c
      subroutine writedf (icnt,rdata,ierr)
      integer icnt
      real*8 rdata(*)
      write (2) (rdata(i),i=1,icnt)
      ierr = 0
      end
c
      subroutine writegff (icnt,rdata,idata,ierr)
      integer icnt,idata(*),ierr
      real*4 rdata(*)
      write (2) (rdata(i),i=1,3*icnt),(idata(i),i=1,icnt)
      ierr = 0
      end
c
      subroutine writegdf (icnt,rdata,idata,ierr)
      integer icnt,idata(*),ierr
      real*8 rdata(*)
      write (2) (rdata(i),i=1,3*icnt),(idata(i),i=1,icnt)
      ierr = 0
      end

