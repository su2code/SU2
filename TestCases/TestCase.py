#!/usr/bin/env python 

## \file TestCase.py
#  \brief Python class for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 4.2.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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

import time, os, subprocess, datetime, sys
import difflib

class TestCase:

    def __init__(self,tag_in):

        self.tag  = tag_in  # Input, string tag that identifies this run

        # Configuration file path/filename
        self.cfg_dir  = "."
        self.cfg_file = "default.cfg"

        # Indicate if the test is unsteady
        self.unsteady = False

        # The test condition. These must be set after initialization
        self.test_iter = 1
        self.test_vals = []  

        # These can be optionally varied 
        self.su2_exec    = "SU2_CFD" 
        self.timeout     = 300
        self.tol         = 0.001

        # Options for file-comparison tests
        self.reference_file = "of_grad.dat.ref"
        self.test_file      = "of_grad.dat"

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

    def run_filediff(self):
        print '==================== Start Test: %s ===================='%self.tag
        passed       = True
        timed_out    = False

        # Adjust the number of iterations in the config file
        self.adjust_iter()

        # Assemble the shell command to run
        logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        command = "%s -f %s > %s" % (self.su2_exec, self.cfg_file, logfilename)

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


        if not timed_out:
            # Compare files
            fromfile = self.reference_file
            tofile = self.test_file 

            fromdate = time.ctime(os.stat(fromfile).st_mtime)
            todate = time.ctime(os.stat(tofile).st_mtime)
            fromlines = open(fromfile, 'U').readlines()
            tolines = open(tofile, 'U').readlines()

            diff = list(difflib.unified_diff(fromlines, tolines, fromfile, tofile, fromdate, todate))


            if (diff==[]):
                passed=True
            else:
                for line in diff:
                    print line[:-1]
                passed=False

        else:
            passed = False

        print 'test duration: %.2f min'%(running_time/60.0)
        print '==================== End Test: %s ====================\n'%self.tag

        os.chdir(workdir)
        return passed

    def run_opt(self):

        print '==================== Start Test: %s ===================='%self.tag
        passed       = True
        exceed_tol   = False
        timed_out    = False
        iter_missing = True
        start_solver = True

        # Adjust the number of iterations in the config file   
        self.adjust_opt_iter()

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
                if not start_solver: # Don't bother parsing anything before optimizer starts
                    if line.find('OBJFUN') > -1:
                        start_solver=True
                else:   # Found the OBJFUN line; parse the input
                    raw_data = line.split()
                    try:
                        iter_number = int(raw_data[0])
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
            print 'ERROR: The code was not able to get to the "OBJFUN" section.'

        if iter_missing:
            print 'ERROR: The optimizer iteration number %d could not be found.'%self.test_iter

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

    def run_geo(self):

        print '==================== Start Test: %s ===================='%self.tag
        passed       = True
        exceed_tol   = False
        timed_out    = False
        start_solver = True
        iter_missing = True

        found_thick  = False
        found_area   = False
        found_twist  = False
        found_chord  = False

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
        data = []
        if not timed_out:
            start_solver = False
            for line in output:
                if not start_solver: # Don't bother parsing anything before SU2_GEO starts
                    if line.find('Section 1') > -1:
                        start_solver=True
                elif line.find('Section 2') > -1: # jump out of loop if we hit the next section
                    break
                else:   # Found the lines; parse the input

                    if line.find('Maximum thickness') > -1:
                        raw_data = line.split()
                        data.append(raw_data[-1][:-1])
                        found_thick = True

                    elif line.find('Area') > -1:
                        raw_data = line.split()
                        data.append(raw_data[-1][:-1])
                        found_area = True

                    elif line.find('Twist angle') > -1:
                        raw_data = line.split()
                        data.append(raw_data[-1][:-1])
                        found_twist = True

                    elif line.find('Chord') > -1:
                        raw_data = line.split()
                        data.append(raw_data[-1][:-1])
                        found_chord = True


            if found_thick and found_area and found_twist and found_chord:  # Found what we're checking for
                iter_missing = False
                if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
                    print "Error in test_vals!"
                    passed = False
                for j in range(len(data)):
                    sim_vals.append( float(data[j]) )
                    delta_vals.append( abs(float(data[j])-self.test_vals[j]) )
                    if delta_vals[j] > self.tol:
                        exceed_tol = True
                        passed     = False
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
            print 'ERROR: The code was not able to get to the "OBJFUN" section.'

        if iter_missing:
            print 'ERROR: The SU2_GEO values could not be found.'

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

    def run_def(self):
    
        print '==================== Start Test: %s ===================='%self.tag
        passed       = True
        exceed_tol   = False
        timed_out    = False
        iter_missing = True
        start_solver = True
    
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
                if not start_solver: # Don't bother parsing anything before -- Volumetric grid deformation ---
                    if line.find('Volumetric grid deformation') > -1:
                        start_solver=True
                else:   # Found the -- Volumetric grid deformation --- line; parse the input
                    raw_data = line.split()
                    try:
                        iter_number = int(raw_data[0])
                        data        = raw_data[len(raw_data)-1:]    # Take the last column for comparison
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
            print 'ERROR: Execution timed out. timeout=%d sec'%self.timeout
    
        if exceed_tol:
            print 'ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%e'%self.tol
    
        if not start_solver:
            print 'ERROR: The code was not able to get to the "Begin solver" section.'
    
        if iter_missing:
            print 'ERROR: The iteration number %d could not be found.'%self.test_iter
    
        print 'test_iter=%d \n'%self.test_iter,
    
        print 'test_vals (stored): ',
        for j in self.test_vals:
            print '%e'%j,
        print '\n',
    
        print 'sim_vals (computed): ',
        for j in sim_vals:
            print '%e'%j,
        print '\n',
    
        print 'delta_vals: ',
        for j in delta_vals:
            print '%e'%j,
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
        file_out = open(self.cfg_file,'w')
        file_out.write('%% This file automatically generated by the regression script\n')
        file_out.write('%% Number of iterations changed to %d\n'%(self.test_iter+1))
        for line in lines:
            if not line.startswith("EXT_ITER"):
                file_out.write(line)
            else:
                file_out.write("EXT_ITER=%d\n"%(self.test_iter+1))
        file_out.close()
        os.chdir(workdir)

        return

    def adjust_opt_iter(self):

        # Read the cfg file
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        file_in = open(self.cfg_file, 'r')
        lines   = file_in.readlines()
        file_in.close()

        # Rewrite the file with a .autotest extension
        self.cfg_file = "%s.autotest"%self.cfg_file
        file_out = open(self.cfg_file,'w')
        file_out.write('%% This file automatically generated by the regression script\n')
        file_out.write('%% Number of optimizer iterations changed to %d\n'%(self.test_iter))
        for line in lines:
            if not line.startswith("OPT_ITERATIONS"):
                file_out.write(line)
            else:
                file_out.write("OPT_ITERATIONS= %d\n"%(self.test_iter))
        file_out.close()
        os.chdir(workdir)

        return
