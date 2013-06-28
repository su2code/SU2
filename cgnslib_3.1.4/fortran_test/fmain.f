      program fmain
      integer*4 i,j
      character*32 string
      i = 1
      j = 3
      string = 'string'
      call cg_sub(i,string,j)
      call adfsub(i,string,j)
      stop
      end
