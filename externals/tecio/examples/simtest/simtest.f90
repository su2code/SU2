!
! Simple example fortran program to write a
! binary datafile for tecplot.  This example
! does the following:
!
!   1.  Open a datafile called "t.plt"
!   2.  Assign values for X,Y, and P
!   3.  Write out a zone dimensioned 4x5
!   4.  Close the datafile.
!
!
      program test

      INCLUDE 'tecio.f90'

      character*1 NULLCHR
      Integer*4   Debug,III,NPts,NElm,isDouble

      Dimension X(4,5), Y(4,5), P(4,5)
      Real*8    SolTime
      Integer*4 VIsDouble, FileType
      Integer*4 FileFormat
      Integer*4 ZoneType,StrandID,ParentZn,IsBlock
      Integer*4 ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer*4 Null(*)

      NULLCHR = CHAR(0)
      NullPtr = 0
      Debug   = 1
      FileType = 0
      FileFormat = 0 ! 0 = PLT, 1 = SZPLT
      VIsDouble = 0
      IMax    = 4
      JMax    = 5
      KMax    = 1
      ZoneType = 0
      SolTime = 360.0
      StrandID = 0
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
!
!... Open the file and write the tecplot datafile 
!... header information.
!
      I = TecIni142('SIMPLE DATASET'//NULLCHR, &
                    'X Y P'//NULLCHR, &
                    'simtestf90-t.plt'//NULLCHR, &
                    '.'//NULLCHR, &
                    FileFormat, &
                    FileType, &
                    Debug, &
                    VIsDouble)

      Do 10 I = 1,4
      Do 10 J = 1,5
        X(I,J) = I
        Y(I,J) = J
        P(I,J) = I*J
   10 Continue
!
!... Write the zone header information.
!
      I = TecZne142('Simple Zone'//NULLCHR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    KMax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    ShrConn)
!
!... Write out the field data.
!
      III = IMax*JMax
      isDouble = 0
      I   = TecDat142(III,X,0)
      I   = TecDat142(III,Y,0)
      I   = TecDat142(III,P,0)

      I = TecEnd142()
      Stop
      End
