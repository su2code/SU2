#!MC 1200
# Created by Tecplot 360 build 12.0.0.4231
$!VarSet |MFBD| = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Gyarmathy_Nozzle_18C'
$!EXTRACTFROMPOLYLINE 
  EXTRACTTHROUGHVOLUME = NO
  EXTRACTLINEPOINTSONLY = NO
  INCLUDEDISTANCEVAR = NO
  NUMPTS = 1000
  EXTRACTTOFILE = YES
  FNAME = '/home/lazzini/SOURCE/SU2/TestCases/condensation/nozzles/Gyarmathy_Nozzle_18C/streamline.dat'
  RAWDATA
2
-0.02 0.0005 0
 0.02 0.0005 0
$!PICK SETMOUSEMODE
  MOUSEMODE = SELECT
$!RemoveVar |MFBD|
