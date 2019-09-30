! 
!  Complex example FORTRAN program to write a
!  binary data file for Tecplot. This example
!  does the following:
! 
!    1.  Open a data file called "field.plt."
!    2.  Open a data file called "line.plt."
!    3.  Assign values for X, Y and P. These will be used
!        in both the ordered and FE data files.
!    4.  Write out an ordered zone dimensioned 4 x 5 to "field.plt."
!    5.  Assign values for XL and YL arrays.
!    6.  Write out data for line plot to "line.plt."  Make the data
!        use double precision.
!    7.  Write out a finite element zone to "field.plt."
!    8.  Write out a text record to "field.plt."
!    9.  Write out a geometry (circle) record to "field.plt."
!   10.  Close file 1.
!   11.  Close file 2.
!  
      Program ComplexTest

      Include "tecio.f90"

      REAL*4      X(4,5),  Y(4,5), P(4,5)  
      REAL*8      XL(50), YL(50)
      REAL*4      XLDummy(1), YLDummy(1)
      EQUIVALENCE (XLDummy(1), XL(1))
      EQUIVALENCE (YLDummy(1), YL(1))
      REAL*8      SolTime
      INTEGER*4   Debug,I,J,K,L,III,NPts,NElm,DIsDouble,VIsDouble,FileType
      INTEGER*4   FileFormat
      INTEGER*4   IMax,JMax,KMax,NM(4,12)   
      INTEGER*4   StrandID,ParentZn
      INTEGER*4   SharingZone(3)
      REAL*8      XP, YP, ZP, FH, LineSpacing, PatternLength
      REAL*8      BoxMargin, BoxLineThickness, TextAngle
      INTEGER*4   AttachToZone, Zone, Scope, PositionCoordSys
      INTEGER*4   Clipping
      INTEGER*4   FontType, HeightUnits, Anchor, BoxType
      INTEGER*4   IsFilled, GeomType, LinePattern, NumEllipsePts
      INTEGER*4   BoxColor, BoxFillColor, TextColor, Color, FillColor
      INTEGER*4   ArrowheadStyle, ArrowheadAttachment, NumSegments
      INTEGER*4   NumSegPts(1)
      REAL*8      LineThickness, ArrowheadSize, ArrowheadAngle
      REAL*4      XGeomData(1), YGeomData(1), ZGeomData(1)
      CHARACTER*1 NULCHAR
      INTEGER*4   Zero
      POINTER     (NullPtr,Null)
      INTEGER*4   Null(*)

      Debug      = 2
      VIsDouble  = 0
      FileType   = 0
      FileFormat = 0 ! 0 = PLT, 1 = SZPLT
      DIsDouble  = 0
      NULCHAR    = CHAR(0)
      Zero       = 0
      NullPtr    = 0
!
! Open field.plt and write the header information.
! 
      I = TECINI142('DATASET WITH 1 ORDERED ZONE, '// &
                    '1 QUAD ZONE OVER 2 TIME STEPS'//NULCHAR, &
                    'X Y P'//NULCHAR, &
                    'comtestf90-field.plt'//NULCHAR, &
                    '.'//NULCHAR, &
                     FileFormat, &
                     FileType, &
                     Debug, &
                     VIsDouble)
!  
!  Open line.plt and write the header information.
!  
      VIsDouble = 1
      I = TECINI142('DATASET WITH ONE I-ORDERED ZONE'//NULCHAR, &
                    'X Y'//NULCHAR, &
                    'comtestf90-line.plt'//NULCHAR, &
                    '.'//NULCHAR, &
                     FileFormat, &
                     FileType, &
                     Debug, &
                     VIsDouble)

!  
!  Calculate values for the field variables.
!  
      Do 10 J = 1,5
      Do 10 I = 1,4
          X(I,J) = I
          Y(I,J) = J
          P(I,J) = I*J
   10 Continue

!  
!  Make sure writing to file #1.
!  
      III = 1
      I = TECFIL142(III)

!  
!  Write the zone header information for the ordered zone.
!  
      IMax = 4
      JMax = 5
      KMax = 1
      SolTime = 10.0
      StrandID = 1
      ParentZn = 0
      I = TECZNE142('Ordered Zone 1'//NULCHAR, &
                     0,  & ! ZONETYPE
                     IMax, &
                     JMax, &
                     KMax, &
                     0, &
                     0, &
                     0, &
                     SolTime, &
                     StrandID, &
                     ParentZn, &
                     1, &     ! ISBLOCK
                     0, &     ! NumFaceConnections
                     0, &     ! FaceNeighborMode
                     0, &     ! TotalNumFaceNodes
                     0, &     ! NumConnectedBoundaryFaces
                     0, &     ! TotalNumBoundaryConnections
                     Null, &  ! PassiveVarList
                     Null, &  ! ValueLocation
                     Null, &  ! ShareVarFromZone
                     0)       ! ShareConnectivityFromZone)

!  
!  Write out the field data for the ordered zone.
!  
      III = IMax*JMax
      I   = TECDAT142(III,X,DIsDouble)
      I   = TECDAT142(III,Y,DIsDouble)
      I   = TECDAT142(III,P,DIsDouble)

!   
!  Calculate values for the I-ordered zone.
!  

      Do 20 I = 1,50
         XL(I) = I
         YL(I) = sin(I/20.0)
   20 Continue
!  
!  Switch to the 'line.plt' file (file number 2)
!  and write out the line plot data.
!  
      III = 2
      I = TECFIL142(III)
!  
!  Write the zone header information for the XY-data.
!  
      IMax = 50
      JMax = 1
      KMax = 1
      SolTime = 0.0
      StrandID = 0
      I = TECZNE142('XY Line plot'//NULCHAR, &
                    0, &
                    IMax, &
                    JMax, &
                    KMax, &
                    0, &
                    0, &
                    0, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    1, &
                    0, &
                    0, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    0)
!  
!  Write out the line plot.
!  
      DIsDouble = 1
      III = IMax
      I   = TECDAT142(III,XLDummy,DIsDouble)
      I   = TECDAT142(III,YLDummy,DIsDouble)

!  
!  Switch back to the field plot file and write out
!  the finite-element zone.
!  
      III = 1
      I = TECFIL142(III)
!
!  Move the coordinates so this zone's not on top of the other
!
      Do 30 J = 1,5
      Do 30 I = 1,4
          X(I,J) = I+5
          Y(I,J) = J
          P(I,J) = I*J
   30 Continue
!  
!  Write the zone header information for the finite-element zone.
!  
      NPts      = 20 
      NElm      = 12 
      KMax      = 1  
      SolTime   = 10.0
      StrandID  = 2
      I = TECZNE142('Finite Zone 1'//NULCHAR, &
                    3, &  ! FEQUADRILATERAL
                    NPts, &
                    NElm, &
                    KMax, &
                    0, &
                    0, &
                    0, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    1, &
                    0, &
                    0, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    0)
!  
!  Write out the field data for the finite-element zone.
!  
      IMax      = 4
      JMax      = 5
      III       = IMax*JMax
      DIsDouble = 0
      I    = TECDAT142(III,X,DIsDouble)
      I    = TECDAT142(III,Y,DIsDouble)
      I    = TECDAT142(III,P,DIsDouble)

!  
!  Calculate and then write out the connectivity list.
!  Note: The NM array references cells starting with
!        offset of 1.
!  

      Do 40 I = 1,IMax-1
      Do 40 J = 1,JMax-1
          K = I+(J-1)*(IMax-1)
          L = I+(J-1)*IMax
          NM(1,K) = L
          NM(2,K) = L+1
          NM(3,K) = L+IMax+1
          NM(4,K) = L+IMax
   40 Continue

      I = TECNOD142(NM)
!
!  Calculate vlues for the new solution variable.
!
      Do 50 J = 1,5
      Do 50 I = 1,4
          P(I,J) = 2*I*J
   50 Continue
!
! Write the zone header information for time step 2
!
      IMax           = 4
      JMax           = 5
      KMax           = 1
      SolTime        = 20.0
      StrandID       = 1
      SharingZone(1) = 1
      SharingZone(2) = 1
      SharingZone(3) = 0
      I = TECZNE142('Ordered Zone 2'//NULCHAR, &
                    0, &  ! ORDERED
                    IMax, &
                    JMax, &
                    KMax, &
                    0, &
                    0, &
                    0, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    1, &
                    0, &
                    0, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    SharingZone, &
                    0)
!
!  Write out the solution variable the grid variables are shared.
!
      IMax      = 4
      JMax      = 5
      III       = IMax*JMax
      DIsDouble = 0
      I   = TECDAT142(III,P,DIsDouble)
!
!  Calculate values for the new solution variable.
!
      Do 60 J = 1,5
      Do 60 I = 1,4
          P(I,J) = 3*I*J
   60 Continue
!
!  Write another time step for the FEZone and share from the first
!
      SolTime = 20.0
      StrandID = 2
      KMax = 0
      SharingZone(1) = 2
      SharingZone(2) = 2
      SharingZone(3) = 0
      I = TECZNE142('Finite Zone 2'//NULCHAR, &
                    3, &  ! FEQUADRILATERAL
                    NPts, &
                    NElm, &
                    KMax, &
                    0, &
                    0, &
                    0, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    1, &
                    0, &
                    0, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    SharingZone, &
                    2)
!
!  Write out the solution variable the grid variables are shared.
!
      IMax      = 4
      JMax      = 5
      III       = IMax*JMax
      DIsDouble = 0
      I   = TECDAT142(III,P,DIsDouble)

!  
!  Prepare to write out text record. Text is positioned
!  at 50, 50 in frame units and has a height 5 frame units.
!  
      XP               = 50 
      YP               = 50 
      FH               = 5
      Scope            = 1 
      Clipping         = 0
      PositionCoordSys = 1 
      FontType         = 1 
      HeightUnits      = 1 
      AttachToZone     = 0
      Zone             = 0
      BoxType          = 0 
      BoxMargin        = 5.0
      BoxLineThickness = 0.5
      BoxColor         = 3
      BoxFillColor     = 7
      TextAngle        = 0.0
      Anchor           = 0 
      LineSpacing      = 1.5
      TextColor        = 0 
    
      III =  TECTXT142(XP, &
                       YP, &
                       0.0d0, &
                       PositionCoordSys, &
                       AttachToZone, &
                       Zone, &
                       FontType, &
                       HeightUnits, &
                       FH, &
                       BoxType, &
                       BoxMargin, &
                       BoxLineThickness, &
                       BoxColor, &
                       BoxFillColor, &
                       TextAngle, &
                       Anchor, &
                       LineSpacing, &
                       TextColor, &
                       Scope, &
                       Clipping, &
                      'Hi Mom'//NULCHAR, &
                      NULCHAR)

!  
!  Prepare to write out geometry record (circle). Circle is 
!  positioned at 25, 25 in frame units and has a radius of 20.
!  Circle is drawn using a dashed line pattern.
!  


      XP                  = 25
      YP                  = 25
      ZP                  = 0.0
      IsFilled            = 0
      Color               = 0
      FillColor           = 7
      GeomType            = 3 
      LinePattern         = 1 
      LineThickness       = 0.3
      PatternLength       = 1.5
      NumEllipsePts       = 72
      ArrowheadStyle      = 0
      ArrowheadAttachment = 0
      ArrowheadSize       = 0.0
      ArrowheadAngle      = 15.0
      NumSegments         = 1
      NumSegPts(1)        = 1
    
      XGeomData(1) = 20
      YGeomData(1) = 0.0
      ZGeomData(1) = 0.0
    
    
      III =  TECGEO142(XP, &
                       YP, &
                       ZP, &
                       PositionCoordSys, &
                       AttachToZone, &
                       Zone, &
                       Color, &
                       FillColor, &
                       IsFilled, &
                       GeomType, &
                       LinePattern, &
                       PatternLength, &
                       LineThickness, &
                       NumEllipsePts, &
                       ArrowheadStyle, &
                       ArrowheadAttachment, &
                       ArrowheadSize, &
                       ArrowheadAngle, &
                       Scope, &
                       Clipping, &
                       NumSegments, &
                       NumSegPts, &
                       XGeomData, &
                       YGeomData, &
                       ZGeomData, &
                       NULCHAR)
          
! 
!  Close out file 1.
! 
      I = TECEND142() 

!  
!  Close out file 2.
!  
      III = 2
      I = TECFIL142(III)
      I = TECEND142() 
      STOP
      END
