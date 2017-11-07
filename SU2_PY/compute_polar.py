#!/usr/bin/env python

## \file Compute_polar.py
#  \brief Python script for performing polar sweep.
#  \author E Arad (based on T. Lukaczyk and  F. Palacios script)
#  \version 5.0.0 "Raven"
#
# SU2 Original Developers: Dr. Francisco D. Palacios.
#                          Dr. Thomas D. Economon.
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# imports
import numpy as np
from optparse import OptionParser
import os, sys, shutil, copy, os.path
sys.path.append(os.environ['SU2_RUN'])
import SU2
from SU2.util.polarSweepLib import *

def main():

# Command Line Options
   parser = OptionParser()
   parser.add_option("-c", "--ctrl", dest="ctrlFile",
                     help="reads polar control parameters from FILE (default:polarCtrl.in) ",
                     metavar="FILE", default="polarCtrl.in")
   parser.add_option("-n", "--partitions", dest="partitions", default=2,
                     help="number of PARTITIONS", metavar="PARTITIONS")
   parser.add_option("-i", "--iterations", dest="iterations", default=99999,
                     help="number of ITERATIONS", metavar="ITERATIONS")
   parser.add_option("-d", "--dimmension", dest="geomDim", default=2,
                     help="Geometry dimension (2 or 3)", metavar="geomDim")
   parser.add_option("-w","--Wind", action="store_true", dest="Wind",default=False,
                     help=" Wind system (default is body system" )
   parser.add_option("-v","--Verbose", action="store_true", dest="verbose",default=False,
                     help=" Verbose printout (if activated)" )

   (options, args)=parser.parse_args()
   options.partitions = int( options.partitions )
   options.iterations = int( options.iterations )
   options.geomDim = int(options.geomDim )

   d2r=np.pi/180

   #--------------- now read the parameters control file and parse it
    
   fc=open(options.ctrlFile,'r')
   ctrl=fc.readlines()
   nc=np.size(ctrl)
   fc.close()

   print str(nc)+" lines read from control file: "+options.ctrlFile

   PA,polarSweepType,velDirOption,nAlpha,nBeta,nPhi,nMach,alpha,beta,phi,MachList,polarVar=setPolaraType(ctrl,nc,options.verbose)

   if options.verbose:
      velDirOptionLegend=['V(alpha,phi)','V(alpha,beta)']
      print '>>>  Control file details: Pitch axis is '+PA+'. Polar sweep type is '+str(polarSweepType)+\
         '; polarVar = '+polarVar
      print '>>>  Velocity definiton: '+velDirOptionLegend[velDirOption-1]
      print '>>>  nAalpha = '+str(nAlpha)+'; nBeta = '+str(nBeta)+'; nPhi = '+str(nPhi)+'; nMach = '+str(nMach)
   if polarSweepType < 4 :
      nPolara = max(nAlpha,nPhi)
   else:
      nPolara = nMach

   #-------------Configuration base file ----------------------- 
   inputbaseFileString='input base file'
   keyWordInputbaseFile=inputbaseFileString.lower()
   iBaseInputF = parLocator(keyWordInputbaseFile,ctrl,nc,-1,options.verbose)
   bIFLine=ctrl[iBaseInputF]
   icol=bIFLine.index(':')
   sBIF=bIFLine[icol+1:]
   inputbaseFile=sBIF.strip(' ')
   inputbaseFile=inputbaseFile.strip('\n')

   print ' '
   print '--------------------------------------------------------------------------------------------------------------------'
   print ' '
   print 'Configuration file: '+inputbaseFile
   print 'PolarSweepType = '+str(polarSweepType)+' Polar sweep in '+polarVar+' using '+str(nPolara)+' angles/Mach No '
   print ' '
   print '--------------------------------------------------------------------------------------------------------------------'
   print ' '
      
   if polarSweepType == 4:
      nPolara=1 # prevent angles inner loop
   if options.geomDim not in [2,3]:
      raise SystemExit('ERROR: dimension can be either 2 or 3 (-d parameter)  ')

   if options.Wind:
      outSystem='Wind'
   else:
      outSystem='Body'

   print " "
   print "==============================================================================="
   print "   Polar sweep in "+str(options.geomDim)+"D ; output in "+outSystem+" system"
   print "==============================================================================="
   print " "
   

   # load config, start state
   config = SU2.io.Config(inputbaseFile)
   state  = SU2.io.State()

   # prepare config
   config.NUMBER_PART = options.partitions
   config.EXT_ITER    = options.iterations
   config.NZONES      = 1

   # find solution files if they exist
   state.find_files(config)

   # start results data
   results = SU2.util.bunch()

   if nMach==0:
      if 'MACH_NUMBER' in config:
         MachList.append(config.MACH_NUMBER)
      else:
         MachList.append(0.5)
      nMach=1

   if nAlpha == 0:
      if 'AOA' in config:
         alpha.append(config.AOA)
      else:
         alpha.append(0.0)
      nAlpha=1

   if nPhi == 0:
      phi.append(0.0)
      nPhi=1
      noPhi_in_CTRL=True
   else:
      noPhi_in_CTRL=False

   if nBeta ==0:
      if noPhi_in_CTRL:
         if 'SIDESLIP_ANGLE' in config:
            beta.append(config.SIDESLIP_ANGLE)
         else:
             beta.append(0.0)
         nBeta=1
      else:
         if polarSweepType < 4:  # alpha sweep with phi set
            tAlpha= [np.tan(d2r*x) for x in alpha]
            tPhi= [np.tan(d2r*x) for x in phi]
            tb=[x*y  for y in tAlpha for x in tPhi]
            beta=[np.arctan(x)/d2r for x in tb]
            nBeta=np.size(beta)
         else:   # Mach ramp
            if 'SIDESLIP_ANGLE' in config:
               beta.append(config.SIDESLIP_ANGLE)
            else:
               beta.append(0.0)
            nBeta=1

   if options.verbose:
      print '>>> alpha: '+str(alpha)
      print '>>> beta:  '+str(beta)
      print '>>> phi:   '+str(phi)
      print '>>> Mach   '+str(MachList)

      
   results.AOA = alpha
   results.MACH = MachList
   results.SIDESLIP_ANGLE=beta
   
   if options.geomDim == 3:
      results.MOMENT_X = []
      results.MOMENT_Y = []
   if options.Wind:
      results.DRAG = []
      results.LIFT = []
      if options.geomDim == 3:
         results.SIDEFORCE = []
   else:
      results.FORCE_X = []
      results.FORCE_Y = []
      if options.geomDim == 3:
         results.FORCE_Z = []
    
   results.MOMENT_Z = []

   if polarSweepType ==4:
      outFile='machRamp_aoa' + str(alpha[0]) + '.dat'
   else:
      outFile='Polar_M' + str(MachList[0]) + '.dat'
   f = open(outFile, 'w')
   if options.verbose:
      print 'Opening polar sweep file: '+outFile
   if options.Wind:
      f.write('%  AOA, Mach, CL, CD,  ')
      if options.geomDim == 3:
         f.write('CY,  ')
   else:
      f.write('%  AOA, Mach, CX, CY,  ')
      if options.geomDim == 3: 
         f.write('CZ,  ')
        
   if options.geomDim == 3:
      f.write('Cmx,  Cmy,   ')

   f.write('Cmz \n')
  
   # iterate mach
   for MachNumber in MachList:

      # iterate angles
      for j in range(0,nPolara):
         if polarSweepType < 3:
            AngleAttack=alpha[j]
            SIDESLIP_ANGLE = beta[0]
         elif polarSweepType == 3:
            AngleAttack=alpha[0]
            SIDESLIP_ANGLE = beta[j]
         else:
             AngleAttack=alpha[0]
             SIDESLIP_ANGLE = beta[0]

         if options.verbose:
            print 'Sweep step '+str(j)+': Mach = '+str(MachNumber)+', aoa = ',str(AngleAttack)+', beta = '\
               +str(SIDESLIP_ANGLE)
             
         # local config and state
         konfig = copy.deepcopy(config)
         ztate  = copy.deepcopy(state)
         #
         # The eval functions below requires definition of various optimization
         # variables, though we are handling here only a direct solution.
         # So, if they are missing in the cfg file (and only then), some dummy values are
         # introduced here
         if not 'OBJECTIVE_FUNCTION' in konfig:
            konfig.OBJECTIVE_FUNCTION='DRAG'
         if not 'DV_KIND' in konfig:
            konfig.DV_KIND=['FFD_SETTING']
         if not 'DV_PARAM' in konfig:
            konfig.DV_PARAM={'FFDTAG': ['1'], 'PARAM': [[0.0, 0.5]], 'SIZE': [1]}
         if not 'DEFINITION_DV' in konfig:
            konfig.DEFINITION_DV={'FFDTAG': [[]],
                                  'KIND': ['HICKS_HENNE'],
                                  'MARKER': [['WING']],
                                  'PARAM': [[0.0, 0.05]],
                                  'SCALE': [1.0],
                                  'SIZE': [1]}
         if not 'OPT_OBJECTIVE' in konfig:
            obj = {}
            obj['DRAG'] = {'SCALE':1.e-2,'OBJTYPE':'DEFAULT'}
            konfig.OPT_OBJECTIVE  = obj
         #
         # --------- end of dummy optimization variables definition section -----------------------
         #
         
    
         # set angle of attack and side-slip angle
         konfig.AOA = AngleAttack
         konfig.SIDESLIP_ANGLE = SIDESLIP_ANGLE
         konfig.MACH_NUMBER = MachNumber
         caseName='DIRECT_M_'+str(MachNumber)+'_AOA_'+str(AngleAttack)
         print 'Mach = ' , konfig.MACH_NUMBER , 'AOA = ' , konfig.AOA
         print 'case :'+caseName
    
         # run su2
         if options.Wind:
            drag = SU2.eval.func('DRAG',konfig,ztate)
            lift = SU2.eval.func('LIFT',konfig,ztate)
            if options.geomDim == 3:
               sideforce=SU2.eval.func('SIDEFORCE',konfig,ztate)
         else:
            force_x=SU2.eval.func('FORCE_X',konfig,ztate)
            force_y=SU2.eval.func('FORCE_Y',konfig,ztate)
            if options.geomDim == 3:
               force_z=SU2.eval.func('FORCE_Z',konfig,ztate)
          
         momentz = SU2.eval.func('MOMENT_Z',konfig,ztate)
         if options.geomDim == 3:
            momentx = SU2.eval.func('MOMENT_X',konfig,ztate)
            momenty = SU2.eval.func('MOMENT_Y',konfig,ztate)

         # append results

         if options.Wind:
            results.DRAG.append(drag)
            results.LIFT.append(lift)
            if options.geomDim == 3:
               results.SIDEFORCE.append(sideforce)
         else:
            results.FORCE_X.append(force_x)
            results.FORCE_Y.append(force_y)
            if options.geomDim == 3:
               results.FORCE_Z.append(force_z)

         results.MOMENT_Z.append(momentz)
         if options.geomDim == 3:
            results.MOMENT_X.append(momentz)
            results.MOMENT_Y.append(momenty)

         output = str(AngleAttack) + ", "+str(MachNumber)+", "

         if options.Wind:
            output = output+ str(lift) + ", " + str(drag)
            if options.geomDim == 3:
               output = output+", "+str(sideforce)
         else:
            output = output+ str(force_x) + ", " + str(force_y)
            if options.geomDim == 3:
               output = output+", "+str(force_z)
         if options.geomDim == 3:
            output = output+", "+str(momentx)+", "+str(momenty)
         output = output+", "+str(momentz)+" \n"

         f.write(output)
         if os.path.isdir(caseName):
            os.system('rm -R '+caseName)
         command='mv DIRECT '+caseName
         if options.verbose:
            print command
         os.system(command)

   f.close()

   # Close open file



#         sys.exit(0)

      #----------------------------------------------------------#
      
   #: for each angle

   # plotting
   #plt.figure()
   #plt.plot( results.MACH_NUMBER, results.AOA , results.LIFT , results.DRAG )
   #plt.show()

   # save data
   SU2.io.save_data('results.pkl',results)
   
if __name__ == "__main__":
    main()
