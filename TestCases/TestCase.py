#!/usr/bin/env python 

## \file TestCase.py
#  \brief Python class for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 4.0.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

import time, os, subprocess, datetime

class TestCase:

    def __init__(self,tag_in):

        self.tag  = tag_in  # Input, string tag that identifies this run

        # Configuration file path/filename
        self.cfg_dir  = "."
        self.cfg_file = "default.cfg"

        # Indicate if the test is unsteady
        self.unsteady = False
        
        # Indicates if the test uses discrete adjoint
	self.disc_adj = False
        
        # The test condition. These must be set after initialization
        self.test_iter = 1
        self.test_vals = []  

        # These can be optionally varied 
        self.su2_exec    = "SU2_CFD" 
        self.timeout     = 300
        self.tol         = 0.001

    def run_test(self):

        print '==================== Start Test: %s ===================='%self.tag
        passed       = True
        exceed_tol   = False
        timed_out    = False
        iter_missing = True
        start_solver = True

        # Adjust the number of iterations in the config file   
        self.adjust_iter()

        # Assemble the shell command to run SU2
        logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        command = "%s %s > %s" % (self.su2_exec, self.cfg_file,logfilename)

        # Run SU2
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        print os.getcwd()
        start   = datetime.datetime.now()
        process = subprocess.Popen(command, shell=True)  # This line launches SU2

        # check for timeout
        while process.poll() is None:
            time.sleep(0.1)
            now = datetime.datetime.now()
            running_time = (now - start).seconds
            if running_time > self.timeout:
                try:
                    process.kill()
                    os.system('killall %s' % self.su2_exec)   # In case of parallel execution
                except AttributeError: # popen.kill apparently fails on some versions of subprocess... the killall command should take care of things!
                    pass
                timed_out = True
                passed    = False

        # Examine the output
        f = open(logfilename,'r')
        output = f.readlines()
        delta_vals = []
        sim_vals = []
        if not timed_out:
            start_solver = False
            for line in output:
                if not start_solver: # Don't bother parsing anything before --Start solver ---
                    if line.find('Begin Solver') > -1:
                        start_solver=True
                else:   # Found the --Begin solver --- line; parse the input
                    raw_data = line.split()
                    try:
                        iter_number = int(raw_data[0])
                        if self.unsteady:
                            iter_number = int(raw_data[1])
                        data        = raw_data[len(raw_data)-4:]    # Take the last 4 columns for comparison
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

        # Write the test results 
        #for j in output:
        #  print j

        if passed:
            print "%s: PASSED"%self.tag
        else:
            print "%s: FAILED"%self.tag
            print 'Output for the failed case'
            subprocess.call(['cat', logfilename])      

        print 'execution command: %s'%command

        if timed_out:
            print 'ERROR: Execution timed out. timeout=%d'%self.timeout

        if exceed_tol:
            print 'ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f'%self.tol

        if not start_solver:
            print 'ERROR: The code was not able to get to the "Begin solver" section.'

        if iter_missing:
            print 'ERROR: The iteration number %d could not be found.'%self.test_iter

        print 'test_iter=%d \n'%self.test_iter,

        print 'test_vals (stored): ',
        for j in self.test_vals:
            print '%f,'%j,
        print '\n',

        print 'sim_vals (computed): ',
        for j in sim_vals:
            print '%f,'%j,
        print '\n',

        print 'delta_vals: ',
        for j in delta_vals:
            print '%f,'%j,
        print '\n',
        
        print 'test duration: %.2f min'%(running_time/60.0) 
        print '==================== End Test: %s ====================\n'%self.tag
        
        os.chdir(workdir)
        return passed

    def adjust_iter(self):

        # Read the cfg file
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        file_in = open(self.cfg_file, 'r')
        lines   = file_in.readlines()
        file_in.close()

        # Rewrite the file with a .autotest extension
        self.cfg_file = "%s.autotest"%self.cfg_file
        if self.disc_adj:
            self.cfg_file = self.cfg_file + '_disc_adj'
        file_out = open(self.cfg_file,'w')
        file_out.write('%% This file automatically generated by the regression script\n')
        file_out.write('%% Number of iterations changed to %d\n'%(self.test_iter+1))
        if self.disc_adj:
            file_out.write('%% MATH_PROBLEM changed to %s\n'%("DISCRETE_ADJOINT"))
            file_out.write('%% RESTART_SOL changed to %s\n'%("NO"))
        for line in lines:
            if all([line.startswith("MATH_PROBLEM"), self.disc_adj]):
                file_out.write("MATH_PROBLEM=DISCRETE_ADJOINT\n")
            elif all([line.startswith("RESTART_SOL"), self.disc_adj]):
                file_out.write("RESTART_SOL=NO\n")
            elif line.startswith("EXT_ITER"):
                file_out.write("EXT_ITER=%d\n"%(self.test_iter+1))
            else:
                file_out.write(line)
        file_out.close()
        os.chdir(workdir)

        return
