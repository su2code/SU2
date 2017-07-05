#!MC 1200
# Created by Tecplot 360 build 12.0.0.4231
$!VarSet |MFBD| = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Moore_Nozzle/Moore_Nozzle_verification'
$!EXTRACTFROMPOLYLINE 
  EXTRACTTHROUGHVOLUME = NO
  EXTRACTLINEPOINTSONLY = NO
  INCLUDEDISTANCEVAR = NO
  NUMPTS = 1000
  EXTRACTTOFILE = YES
  FNAME = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Moore_Nozzle/Moore_Nozzle_verification/streamline.dat'
  RAWDATA
2
-0.10 0.00 0
0.5 0.00 0
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!RemoveVar |MFBD|
