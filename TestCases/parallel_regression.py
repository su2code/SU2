#!/usr/bin/env python 

## \file autotest.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author Aniket C. Aranake, Alejandro Campos, Thomas D. Economon
#  \version 2.0.2.
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys,time, os, subprocess, datetime, signal, os.path

class testcase:

  def __init__(self,tag_in):

    datestamp = time.strftime("%Y%m%d", time.gmtime())
    self.tag  = "%s_%s"%(tag_in,datestamp)  # Input, string tag that identifies this run

    # Configuration file path/filename
    self.cfg_dir  = "/home/ale11"
    self.cfg_file = "default.cfg"

    # The test condition. These must be set after initialization
    self.test_iter = 1
    self.test_vals = []  

    # These can be optionally varied 
    self.su2_dir     = "/home/ale11"
    self.su2_exec    = "default" 
    self.timeout     = 300
    self.tol         = 0.001
    self.outputdir   = "/home/ale11"

  def run_test(self):

    passed       = True
    exceed_tol   = False
    timed_out    = False
    iter_missing = True
    start_solver = True

    # Adjust the number of iterations in the config file   
    self.do_adjust_iter()

    # Assemble the shell command to run SU2
    self.su2_exec = os.path.join("$SU2_RUN", self.su2_exec)
    command_base = "%s %s > outputfile"%(self.su2_exec, self.cfg_file)
    command      = "%s"%(command_base)

    # Run SU2
    os.chdir(os.path.join('./',self.cfg_dir)) 
    start   = datetime.datetime.now()
    process = subprocess.Popen(command, shell=True)  # This line launches SU2

    while process.poll() is None:
      time.sleep(0.1)
      now = datetime.datetime.now()
      if (now - start).seconds> self.timeout:
        try:
          process.kill()
          os.system('killall %s' % self.su2_exec)   # In case of parallel execution
        except AttributeError: # popen.kill apparently fails on some versions of subprocess... the killall command should take care of things!
          pass
        timed_out = True
        passed    = False

    # Examine the output
    f = open('outputfile','r')
    output = f.readlines()
    delta_vals = []
    sim_vals = []
    if not timed_out:
      start_solver = False
      for line in output:
        if not start_solver: # Don't bother parsing anything before --Start solver ---
          if line.find('Begin solver') > -1:
            start_solver=True
        else:   # Found the --Begin solver --- line; parse the input
          raw_data = line.split()
          try:
            iter_number = int(raw_data[0])
            data        = raw_data[2:]    # Take the last 4 columns for comparison
          except ValueError:
            continue
          except IndexError:
            continue
         
          if iter_number == self.test_iter:  # Found the iteration number we're checking for
            iter_missing = False
            if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
              print "Error in test_vals!"
              passed = False
              break
            for j in range(len(data)):
              sim_vals.append( float(data[j]) )
              delta_vals.append( abs(float(data[j])-self.test_vals[j]) )
              if delta_vals[j] > self.tol:
                exceed_tol = True
                passed     = False
            break
          else:
            iter_missing = True

      if not start_solver:
        passed = False
        
      if iter_missing:
        passed = False

    print '=========================================================\n'
      
    # Write the test results 
    #for j in output:
    #  print j

    if passed:
      print "%s: PASSED"%self.tag
    else:
      print "%s: FAILED"%self.tag

    print 'execution command: %s'%command

    if timed_out:
      print 'ERROR: Execution timed out. timeout=%d'%self.timeout

    if exceed_tol:
      print 'ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f'%self.tol

    if not start_solver:
      print 'ERROR: The code was not able to get to the "Begin solver" section.'

    if iter_missing:
      print 'ERROR: The iteration number %d could not be found.'%self.test_iter

    print 'test_iter=%d, test_vals: '%self.test_iter,
    for j in self.test_vals:
      print '%f '%j,
    print '\n',

    print 'sim_vals: ',
    for j in sim_vals:
      print '%f '%j,
    print '\n',
  
    print 'delta_vals: ',
    for j in delta_vals:
      print '%f '%j,
    print '\n'
    
    os.chdir('../../../')
    return passed

  def do_adjust_iter(self):
  
    # Read the cfg file
    self.cfg_file = os.path.join(os.environ['SU2_HOME'], self.cfg_dir, self.cfg_file)
    file_in = open(self.cfg_file, 'r')
    lines   = file_in.readlines()
    file_in.close()
  
    # Rewrite the file with a .autotest extension
    self.cfg_file = "%s.autotest"%self.cfg_file
    file_out = open(self.cfg_file,'w')
    file_out.write('%% This file automatically generated by autotest.py\n')
    file_out.write('%% Number of iterations changed to %d\n'%(self.test_iter+1))
    for line in lines:
      if line.find("EXT_ITER")==-1:
        file_out.write(line)
      else:
        file_out.write("EXT_ITER=%d\n"%(self.test_iter+1))
    file_out.close()

    

if __name__=="__main__":
  '''This program runs SU^2 and ensures that the output matches specified values. This will be used to do nightly checks to make sure nothing is broken. '''

  # Build SU2_CFD in parallel using autoconf
  os.system('./configure --prefix=$SU2_HOME --with-MPI=mpicxx --with-Metis-lib=/home/ale11/tools/metis-5.0.1/lib --with-Metis-include=/home/ale11/tools/metis-5.0.1/include --with-Metis-version=5 CXXFLAGS="-O3"')
  os.system('make clean')
  os.system('make install')

  os.chdir(os.environ['SU2_RUN'])
  if not os.path.exists("./SU2_CFD"):
    print 'Could not build SU2_CFD'
    sys.exit(1)
  
  if not os.path.exists("./SU2_DDC"):
    print 'Could not build SU2_DDC'
    sys.exit(1)
    
  os.chdir(os.environ['SU2_HOME'])  
  os.chdir('../')

  #############
  ### EULER ###
  #############

  # Inviscid Channel
  channel           = testcase('channel')
  channel.cfg_dir   = "TestCases/euler/channel"
  channel.cfg_file  = "inv_channel_RK.cfg"
  channel.test_iter = 100
  channel.test_vals = [-2.361460,3.021777,0.007201,0.050694]
  channel.su2_exec  = "parallel_computation.py -f"
  channel.timeout   = 1600
  channel.tol       = 0.00001
  passed1           = channel.run_test()
  
  # Inviscid NACA0012 
  naca0012           = testcase('naca0012')
  naca0012.cfg_dir   = "TestCases/euler/naca0012"
  naca0012.cfg_file  = "inv_NACA0012_Roe.cfg"
  naca0012.test_iter = 100
  naca0012.test_vals = [-3.669700,-3.181916,0.129571,0.068786]
  naca0012.su2_exec  = "parallel_computation.py -f"
  naca0012.timeout   = 1600
  naca0012.tol       = 0.00001
  passed2            = naca0012.run_test()

  # Supersonic wedge 
  wedge           = testcase('wedge')
  wedge.cfg_dir   = "TestCases/euler/wedge"
  wedge.cfg_file  = "inv_wedge_HLLC.cfg"
  wedge.test_iter = 100
  wedge.test_vals = [-2.676757,3.186233,-0.252186,0.044415 ]
  wedge.su2_exec  = "parallel_computation.py -f"
  wedge.timeout   = 1600
  wedge.tol       = 0.00001
  passed3         = wedge.run_test()

  # Inviscid ONERA M6 Wing
  oneram6           = testcase('oneram6')
  oneram6.cfg_dir   = "TestCases/euler/oneram6"
  oneram6.cfg_file  = "inv_ONERAM6_JST.cfg"
  oneram6.test_iter = 10
  oneram6.test_vals = [-2.623942,-2.088731,0.283443,0.017103]
  oneram6.su2_exec  = "parallel_computation.py -f"
  oneram6.timeout   = 3200
  oneram6.tol       = 0.00001
  passed4           = oneram6.run_test()

  #############
  ###  N-S  ###
  #############

  # Laminar flat plate
  flatplate           = testcase('flatplate')
  flatplate.cfg_dir   = "TestCases/navierstokes/flatplate"
  flatplate.cfg_file  = "lam_flatplate_Roe.cfg"
  flatplate.test_iter = 100
  flatplate.test_vals = [-4.930464,0.551458,0.085523,0.015966]
  flatplate.su2_exec  = "parallel_computation.py -f"
  flatplate.timeout   = 1600
  flatplate.tol       = 0.00001
  passed5             = flatplate.run_test()


  # Laminar cylinder (steady)
  cylinder           = testcase('cylinder')
  cylinder.cfg_dir   = "TestCases/navierstokes/cylinder"
  cylinder.cfg_file  = "lam_cylinder_JST.cfg"
  cylinder.test_iter = 25
  cylinder.test_vals = [-9.663630,-8.686690,0.007463,4.333015]
  cylinder.su2_exec  = "parallel_computation.py -f"
  cylinder.timeout   = 1600
  cylinder.tol       = 0.00001
  passed6            = cylinder.run_test()

  ##########################
  ### Compressible RANS  ###
  ##########################

  # rae2822
  rae2822           = testcase('rae2822')
  rae2822.cfg_dir   = "TestCases/rans/rae2822"
  rae2822.cfg_file  = "turb_SA_RAE2822.cfg"
  rae2822.test_iter = 100
  rae2822.test_vals = [-3.352130,-5.432892,0.881425,0.024306] #last 4 columns
  rae2822.su2_exec  = "parallel_computation.py -f"
  rae2822.timeout   = 1600
  rae2822.tol       = 0.00001
  passed7           = rae2822.run_test()

  # flatplate
  turb_flatplate           = testcase('turb_flatplate')
  turb_flatplate.cfg_dir   = "TestCases/rans/flatplate"
  turb_flatplate.cfg_file  = "turb_SA_flatplate_Roe.cfg"
  turb_flatplate.test_iter = 100
  turb_flatplate.test_vals = [-5.091512,-6.319214,0.001744,0.013873] #last 4 columns
  turb_flatplate.su2_exec  = "parallel_computation.py -f"
  turb_flatplate.timeout   = 1600
  turb_flatplate.tol       = 0.00001
  passed8                  = turb_flatplate.run_test()

  # oneram6
  turb_oneram6           = testcase('turb_oneram6')
  turb_oneram6.cfg_dir   = "TestCases/rans/oneram6"
  turb_oneram6.cfg_file  = "turb_ONERAM6.cfg"
  turb_oneram6.test_iter = 20
  turb_oneram6.test_vals = [-4.706643,-11.453973,0.228641,0.114652] #last 4 columns
  turb_oneram6.su2_exec  = "parallel_computation.py -f"
  turb_oneram6.timeout   = 3200
  turb_oneram6.tol       = 0.00001
  passed9                = turb_oneram6.run_test()
  
  # naca0012
  turb_naca0012           = testcase('turb_naca0012')
  turb_naca0012.cfg_dir   = "TestCases/rans/naca0012"
  turb_naca0012.cfg_file  = "naca0012.cfg"
  turb_naca0012.test_iter = 20
  turb_naca0012.test_vals = [-8.694570,-9.210104,-0.000034,0.008124] #last 4 columns
  turb_naca0012.su2_exec  = "parallel_computation.py -f"
  turb_naca0012.timeout   = 3200
  turb_naca0012.tol       = 0.00001
  passed10                = turb_naca0012.run_test()
  
  ############################
  ### Incompressible RANS  ###
  ############################
  
  # Incompressible NACA0012
  inc_turb_naca0012           = testcase('inc_turb_naca0012')
  inc_turb_naca0012.cfg_dir   = "TestCases/incomp_rans/naca0012"
  inc_turb_naca0012.cfg_file  = "naca0012.cfg"
  inc_turb_naca0012.test_iter = 20
  inc_turb_naca0012.test_vals = [-9.066931,-8.386696,-0.000003,0.008181] #last 4 columns
  inc_turb_naca0012.su2_exec  = "parallel_computation.py -f"
  inc_turb_naca0012.timeout   = 1600
  inc_turb_naca0012.tol       = 0.00001
  passed11                    = inc_turb_naca0012.run_test()
  
  ########################
  ### Cont. Adj. Euler ###
  ########################

  # Inviscid NACA0012 (To be validated with finite differences)
  contadj_naca0012           = testcase('contadj_naca0012')
  contadj_naca0012.cfg_dir   = "TestCases/cont_adj_euler/naca0012"
  contadj_naca0012.cfg_file  = "inv_NACA0012.cfg"
  contadj_naca0012.test_iter = 100
  contadj_naca0012.test_vals = [-4.756880,-10.259729,0.005384,0.506290] #last 4 columns
  contadj_naca0012.su2_exec  = "parallel_computation.py -f"
  contadj_naca0012.timeout   = 1600
  contadj_naca0012.tol       = 0.00001
  passed12                   = contadj_naca0012.run_test()

  # Inviscid RAM-C (To be validated with finite differences)
  contadj_ram_c           = testcase('contadj_ram_c')
  contadj_ram_c.cfg_dir   = "TestCases/cont_adj_euler/ram_c"
  contadj_ram_c.cfg_file  = "inv_RAMC.cfg"
  contadj_ram_c.test_iter = 100
  contadj_ram_c.test_vals = [0.777333,-7.308160,-0.001880,0.080418] #last 4 columns
  contadj_ram_c.su2_exec  = "parallel_computation.py -f"
  contadj_ram_c.timeout   = 1600
  contadj_ram_c.tol       = 0.00001
  passed13                = contadj_ram_c.run_test()

  ######################
  ### Cont. Adj. N-S ###
  ######################

  # Adjoint laminar cylinder
  contadj_ns_cylinder           = testcase('contadj_ns_cylinder')
  contadj_ns_cylinder.cfg_dir   = "TestCases/cont_adj_navierstokes/cylinder"
  contadj_ns_cylinder.cfg_file  = "lam_cylinder.cfg"
  contadj_ns_cylinder.test_iter = 100
  contadj_ns_cylinder.test_vals = [2.533005,-2.050461,0.533940,136.580000] #last 4 columns
  contadj_ns_cylinder.su2_exec  = "parallel_computation.py -f"
  contadj_ns_cylinder.timeout   = 1600
  contadj_ns_cylinder.tol       = 0.00001
  passed14                      = contadj_ns_cylinder.run_test()

  # Adjoint laminar naca0012 (To be fixed)
  contadj_ns_naca0012           = testcase('contadj_ns_naca0012')
  contadj_ns_naca0012.cfg_dir   = "TestCases/cont_adj_navierstokes/naca0012"
  contadj_ns_naca0012.cfg_file  = "lam_NACA0012.cfg"
  contadj_ns_naca0012.test_iter = 100
  contadj_ns_naca0012.test_vals = [1.515646,-3.816185,8.9793e-01,1.7569e-01] #last 4 columns
  contadj_ns_naca0012.su2_exec  = "parallel_computation.py -f"
  contadj_ns_naca0012.timeout   = 1600
  contadj_ns_naca0012.tol       = 0.00001
  passed15                      = True #contadj_ns_naca0012.run_test() 

  ################################
  ### Cont. Adj. RANS (Frozen) ###
  ################################
  
  # Adjoint turbulent NACA0012 (To be validated with finite differences)
  contadj_rans_naca0012           = testcase('contadj_rans_naca0012')
  contadj_rans_naca0012.cfg_dir   = "TestCases/cont_adj_rans/naca0012"
  contadj_rans_naca0012.cfg_file  = "turb_nasa.cfg"
  contadj_rans_naca0012.test_iter = 100
  contadj_rans_naca0012.test_vals = [-5.329331,-8.633984,18.310000,-0.000000] #last 4 columns
  contadj_rans_naca0012.su2_exec  = "parallel_computation.py -f"
  contadj_rans_naca0012.timeout   = 1600
  contadj_rans_naca0012.tol       = 0.00001
  passed16                        = contadj_rans_naca0012.run_test()
  
  #######################################
  ### Cont. Adj. Incompressible Euler ###
  #######################################
  
  # Adjoint Incompressible Inviscid NACA0012
  contadj_incomp_NACA0012           = testcase('contadj_incomp_NACA0012')
  contadj_incomp_NACA0012.cfg_dir   = "TestCases/cont_adj_incomp_euler/naca0012"
  contadj_incomp_NACA0012.cfg_file  = "incomp_NACA0012.cfg"
  contadj_incomp_NACA0012.test_iter = 140
  contadj_incomp_NACA0012.test_vals = [-7.448212,-7.006721,0.003406,0.000000] #last 4 columns
  contadj_incomp_NACA0012.su2_exec  = "parallel_computation.py -f"
  contadj_incomp_NACA0012.timeout   = 1600
  contadj_incomp_NACA0012.tol       = 0.00001
  passed17                          = contadj_incomp_NACA0012.run_test()
  
  #####################################
  ### Cont. Adj. Incompressible N-S ###
  #####################################
  
  # Adjoint Incompressible Viscous Cylinder
  contadj_incomp_cylinder           = testcase('contadj_incomp_cylinder')
  contadj_incomp_cylinder.cfg_dir   = "TestCases/cont_adj_incomp_navierstokes/cylinder"
  contadj_incomp_cylinder.cfg_file  = "lam_incomp_cylinder.cfg"
  contadj_incomp_cylinder.test_iter = 25
  contadj_incomp_cylinder.test_vals = [-9.002409,-10.068685,0.048091,0.000000] #last 4 columns
  contadj_incomp_cylinder.su2_exec  = "parallel_computation.py -f"
  contadj_incomp_cylinder.timeout   = 1600
  contadj_incomp_cylinder.tol       = 0.00001
  passed18                          = contadj_incomp_cylinder.run_test()
  
  if (passed1 and passed2 and passed3 and passed4 and passed5 and passed6 and
      passed7 and passed8 and passed9 and passed10 and passed11 and passed12
      and passed13 and passed14 and passed15 and passed16 and passed17 and passed18):
    sys.exit(0)
  else:
    sys.exit(1)

