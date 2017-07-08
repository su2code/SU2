#!MC 1200
# Created by Tecplot 360 build 12.0.0.4231
$!VarSet |MFBD| = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Bakhtar_Nozzle_S'
$!EXTRACTFROMPOLYLINE 
  EXTRACTTHROUGHVOLUME = NO
  EXTRACTLINEPOINTSONLY = NO
  INCLUDEDISTANCEVAR = NO
  NUMPTS = 1000
  EXTRACTTOFILE = YES
  FNAME = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Bakhtar_Nozzle_S/streamline.dat'
  RAWDATA
2
-0.056 0.00 0
0.0762 0.00 0
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!RemoveVar |MFBD|
