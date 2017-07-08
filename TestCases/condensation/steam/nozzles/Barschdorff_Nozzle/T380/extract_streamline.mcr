#!MC 1200
# Created by Tecplot 360 build 12.0.0.4231
$!VarSet |MFBD| = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Barschdorff_Nozzle/T380'
$!EXTRACTFROMPOLYLINE 
  EXTRACTTHROUGHVOLUME = NO
  EXTRACTLINEPOINTSONLY = NO
  INCLUDEDISTANCEVAR = NO
  NUMPTS = 10000
  EXTRACTTOFILE = YES
  FNAME = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Barschdorff_Nozzle/T380/streamline.dat'
  RAWDATA
2
-0.20 0.00 0
0.12 0.00 0
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!RemoveVar |MFBD|
