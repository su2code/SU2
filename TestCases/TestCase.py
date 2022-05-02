#!/usr/bin/env python 

## \file TestCase.py
#  \brief Python class for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 7.3.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
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
from __future__ import print_function, division, absolute_import
import time, os, subprocess, datetime, sys
import difflib


def print_vals(vals, name="Values"):
    """Print an array of floats."""
    print(name + ': ' + ', '.join('{:f}'.format(v) for v in vals))


class TestCase:

    def __init__(self,tag_in):

        self.tag  = tag_in  # Input, string tag that identifies this run

        # Configuration file path/filename
        self.cfg_dir  = "."
        self.cfg_file = "default.cfg"

        # Indicate if the test is unsteady
        self.unsteady = False

        # Indicate if the test is a polar run
        self.polar = False

        # Indicate whether to disable restart
        self.no_restart = False

        # Indicate whether the new output is used
        self.new_output = True   

        # multizone problem
        self.multizone = False

        # The test condition. These must be set after initialization
        self.test_iter = 1
        self.ntest_vals = 4
        self.test_vals = []  

        # These can be optionally varied 
        self.su2_exec    = "SU2_CFD" 
        self.timeout     = 300
        self.tol         = 0.001

        # Options for file-comparison tests
        self.reference_file = "of_grad.dat.ref"
        self.test_file      = "of_grad.dat"

    def run_test(self):

        print('==================== Start Test: %s ===================='%self.tag)
        passed       = True
        exceed_tol   = False
        timed_out    = False
        iter_missing = True
        start_solver = True

        # if root, add flag to mpirun
        if os.geteuid()==0:
            if self.su2_exec.startswith('mpirun'):
                self.su2_exec = self.su2_exec.replace('mpirun', 'mpirun --allow-run-as-root')

        # Adjust the number of iterations in the config file
        if len(self.test_vals) != 0:
            self.adjust_iter() 

        # Check for disabling the restart
        if self.no_restart:
            self.disable_restart()

        # Assemble the shell command to run SU2
        if len(self.test_vals) != 0:
            logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        else:
            logfilename = '%s_check.log' % os.path.splitext(self.cfg_file)[0]

        # Check for polar calls
        if self.polar:
             command = "%s > %s" % (self.su2_exec, logfilename)
        else:
            command = "%s %s > %s 2>&1" % (self.su2_exec, 
                                           self.cfg_file, 
                                           logfilename)

        # Run SU2
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        print(os.getcwd())
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
        if not timed_out and len(self.test_vals) != 0:
            start_solver = False
            for line in output:
                if not start_solver: # Don't bother parsing anything before --Start solver ---
                    if line.find('Begin Solver') > -1:
                        start_solver=True
                else:   # Found the --Begin solver --- line; parse the input
                    if self.new_output or self.multizone:
                        raw_data = line.strip() # Strip removes whitespaces head-tail
                        raw_data = raw_data[1:-1].split('|') # Remove heat-tail bars before splitting
                    else:
                        raw_data = line.split()
                    try:
                        iter_number = int(raw_data[0])
                        if self.unsteady and not self.multizone and not self.new_output:
                            iter_number = int(raw_data[1])
                        data = raw_data[len(raw_data) - len(self.test_vals):]
                    except ValueError:
                        continue
                    except IndexError:
                        continue

                    if iter_number == self.test_iter:  # Found the iteration number we're checking for
                        iter_missing = False
                        if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
                            print("Error in test_vals!")
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
        #  print(j)


        process.communicate()
        if process.returncode != 0:
            passed = False
        if passed:
            print("%s: PASSED"%self.tag)
        else:
            print("%s: FAILED"%self.tag)
            print('Output for the failed case')
            subprocess.call(['cat', logfilename])      

        print('execution command: %s'%command)

        if timed_out:
            print('ERROR: Execution timed out. timeout=%d'%self.timeout)

        if exceed_tol:
            print('ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f'%self.tol)

        if not start_solver:
            print('ERROR: The code was not able to get to the "Begin solver" section.')

        if iter_missing:
            print('ERROR: The iteration number %d could not be found.'%self.test_iter)

        if len(self.test_vals) != 0:
            print('test_iter=%d' % self.test_iter)

            print_vals(self.test_vals, name="test_vals (stored)")

            print_vals(sim_vals, name="sim_vals (computed)")

            print_vals(delta_vals, name="delta_vals")

        print('test duration: %.2f min'%(running_time/60.0))
        print('==================== End Test: %s ====================\n'%self.tag)

        sys.stdout.flush()
        os.chdir(workdir)
        return passed

    def run_filediff(self):
        print('==================== Start Test: %s ===================='%self.tag)
        passed       = True
        timed_out    = False

        # Adjust the number of iterations in the config file
        self.adjust_iter()

        # if root, add flag to mpirun
        if os.geteuid()==0:
            if self.su2_exec.startswith('mpirun'):
                self.su2_exec = self.su2_exec.replace('mpirun', 'mpirun --allow-run-as-root')

        # Assemble the shell command to run
        logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        command = "%s %s > %s 2>&1" % (self.su2_exec, self.cfg_file, logfilename)

        # Run SU2
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        print(os.getcwd())
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

        # Check for error output from that process
        if process.poll() != 0:
            passed = False
            print("ERROR")
            print("Output from the failed case:")
            subprocess.call(["cat", logfilename])

        if not timed_out and passed:
            # Compare files
            fromfile = self.reference_file
            tofile = self.test_file 
            # Initial value s.t. will fail if it does not get to diff step
            diff = ''
            try:
                fromdate = time.ctime(os.stat(fromfile).st_mtime)
                fromlines = open(fromfile, 'U').readlines()
                try: 
                    todate = time.ctime(os.stat(tofile).st_mtime)
                    tolines = open(tofile, 'U').readlines()
                    diff = list(difflib.unified_diff(fromlines, tolines, fromfile, tofile, fromdate, todate))
                except OSError:
                    print("OS error, most likely from missing reference file:", fromfile)
                    print("Current working directory contents:")
                    print(os.listdir("."))
            except OSError:
                print("OS error, most likely from missing reference file:", fromfile)
                print("Current working directory contents:")
                print(os.listdir("."))





            if (diff==[]):
                passed=True
            else:
                for line in diff:
                    print(line[:-1])
                passed=False

        else:
            passed = False

        print('test duration: %.2f min'%(running_time/60.0))
        print('==================== End Test: %s ====================\n'%self.tag)

        sys.stdout.flush()
        os.chdir(workdir)
        return passed

    def run_opt(self):

        print('==================== Start Test: %s ===================='%self.tag)
        passed       = True
        exceed_tol   = False
        timed_out    = False
        iter_missing = True
        start_solver = True

        # Adjust the number of iterations in the config file   
        self.adjust_opt_iter()

        # Assemble the shell command to run SU2
        logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        command = "%s %s > %s 2>&1" % (self.su2_exec, self.cfg_file, logfilename)

        # Run SU2
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        print(os.getcwd())
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
                        data = raw_data[len(raw_data) - len(self.test_vals):]
                    except ValueError:
                        continue
                    except IndexError:
                        continue

                    if iter_number == self.test_iter:  # Found the iteration number we're checking for
                        iter_missing = False
                        if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
                            print("Error in test_vals!")
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
        #  print(j)

        if passed:
            print("%s: PASSED"%self.tag)
        else:
            print("%s: FAILED"%self.tag)
            print('Output for the failed case')
            subprocess.call(['cat', logfilename])      

        print('execution command: %s'%command)

        if timed_out:
            print('ERROR: Execution timed out. timeout=%d'%self.timeout)

        if exceed_tol:
            print('ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f'%self.tol)

        if not start_solver:
            print('ERROR: The code was not able to get to the "OBJFUN" section.')

        if iter_missing:
            print('ERROR: The optimizer iteration number %d could not be found.'%self.test_iter)

        print('test_iter=%d' % self.test_iter)

        print_vals(self.test_vals, name="test_vals (stored)")

        print_vals(sim_vals, name="sim_vals (computed)")

        print_vals(delta_vals, name="delta_vals")

        print('test duration: %.2f min'%(running_time/60.0))
        print('==================== End Test: %s ====================\n'%self.tag)

        sys.stdout.flush()
        os.chdir(workdir)
        return passed

    def run_geo(self):

        print('==================== Start Test: %s ===================='%self.tag)
        passed       = True
        exceed_tol   = False
        timed_out    = False
        start_solver = True
        iter_missing = True

        found_thick  = False
        found_area   = False
        found_twist  = False
        found_chord  = False

        # if root, add flag to mpirun
        if os.geteuid()==0:
            if self.su2_exec.startswith('mpirun'):
                self.su2_exec = self.su2_exec.replace('mpirun', 'mpirun --allow-run-as-root')
                
        # Assemble the shell command to run SU2
        logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        command = "%s %s > %s 2>&1" % (self.su2_exec, self.cfg_file, logfilename)

        # Run SU2
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        print(os.getcwd())
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
                    if line.find('Station 1') > -1:
                        start_solver=True
                elif line.find('Station 2') > -1: # jump out of loop if we hit the next station
                    break
                else:   # Found the lines; parse the input

                    if line.find('Chord') > -1:
                        raw_data = line.replace(",", "").split()
                        data.append(raw_data[1])
                        found_chord = True
                        data.append(raw_data[5])
                        found_radius = True
                        data.append(raw_data[8])
                        found_toc = True
                        data.append(raw_data[10])
                        found_aoa = True

            if found_chord and found_radius and found_toc and found_aoa:  # Found what we're checking for
                iter_missing = False
                if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
                    print("Error in test_vals!")
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
        #  print(j)

        if passed:
            print("%s: PASSED"%self.tag)
        else:
            print("%s: FAILED"%self.tag)
            print('Output for the failed case')
            subprocess.call(['cat', logfilename])

        print('execution command: %s'%command)

        if timed_out:
            print('ERROR: Execution timed out. timeout=%d'%self.timeout)

        if exceed_tol:
            print('ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f'%self.tol)

        if not start_solver:
            print('ERROR: The code was not able to get to the "OBJFUN" section.')

        if iter_missing:
            print('ERROR: The SU2_GEO values could not be found.')

        print_vals(self.test_vals, name="test_vals (stored)")

        print_vals(sim_vals, name="sim_vals (computed)")

        print_vals(delta_vals, name="delta_vals")

        print('test duration: %.2f min'%(running_time/60.0))
        print('==================== End Test: %s ====================\n'%self.tag)

        sys.stdout.flush()
        os.chdir(workdir)
        return passed

    def run_def(self):
    
        print('==================== Start Test: %s ===================='%self.tag)
        passed       = True
        exceed_tol   = False
        timed_out    = False
        iter_missing = True
        start_solver = True
    
        # if root, add flag to mpirun
        if os.geteuid()==0:
            if self.su2_exec.startswith('mpirun'):
                self.su2_exec = self.su2_exec.replace('mpirun', 'mpirun --allow-run-as-root')

        # Assemble the shell command to run SU2
        logfilename = '%s.log' % os.path.splitext(self.cfg_file)[0]
        command = "%s %s > %s 2>&1" % (self.su2_exec, self.cfg_file, logfilename)
    
        # Run SU2
        workdir = os.getcwd()
        os.chdir(self.cfg_dir)
        print(os.getcwd())
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
                            print("Error in test_vals!")
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
        #  print(j)
    
        if passed:
            print("%s: PASSED"%self.tag)
        else:
            print("%s: FAILED"%self.tag)
            print('Output for the failed case')
            subprocess.call(['cat', logfilename])
    
        print('execution command: %s'%command)
    
        if timed_out:
            print('ERROR: Execution timed out. timeout=%d sec'%self.timeout)
    
        if exceed_tol:
            print('ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%e'%self.tol)
    
        if not start_solver:
            print('ERROR: The code was not able to get to the "Begin solver" section.')
    
        if iter_missing:
            print('ERROR: The iteration number %d could not be found.'%self.test_iter)
    
        print('test_iter=%d' % self.test_iter)

        print_vals(self.test_vals, name="test_vals (stored)")

        print_vals(sim_vals, name="sim_vals (computed)")

        print_vals(delta_vals, name="delta_vals")
 
        print('test duration: %.2f min'%(running_time/60.0))
        #print('==================== End Test: %s ====================\n'%self.tag)
    
        sys.stdout.flush()
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
        if (self.multizone or self.new_output) and self.unsteady:
            adjust_string = "TIME_ITER"
        elif self.multizone:
            adjust_string = "OUTER_ITER"
        else:
            adjust_string = "ITER"
        for line in lines:
            if not line.strip().split("=")[0].strip() == adjust_string:
                file_out.write(line)
            else:
                file_out.write(adjust_string+"=%d\n"%(self.test_iter+1))
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

    def disable_restart(self):

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
            if not line.startswith("RESTART_SOL"):
                file_out.write(line)
            else:
                file_out.write("RESTART_SOL= NO\n")
        file_out.close()
        os.chdir(workdir)

        return
