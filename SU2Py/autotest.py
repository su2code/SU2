#!/usr/bin/python 

## \file autotest.py
#  \brief Python script for automated testing of SU2 examples
#  \author Aniket C. Aranake
#  \version 1.1.
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

import time, os, subprocess, datetime, signal, os.path
import smtplib
from email.mime.text import MIMEText

class testcase:

  def __init__(self,tag_in):

    datestamp = time.strftime("%Y%m%d", time.gmtime())
    self.tag  = "%s_%s"%(tag_in,datestamp)  # Input, string tag that identifies this run

    # Configuration file path/filename
    self.cfg_dir  = "/home/su2_test/su2/trunk/SU2_CFD/bin/"
    self.cfg_file = "default.cfg"

    # The test condition. These must be set after initialization
    self.test_iter = 1
    self.test_vals = []  

    # These can be optionally varied 
    self.su2_dir     = "/home/su2_test/su2/trunk/SU2_CFD/bin/"
    self.su2_exec    = "SU2_CFD"
    self.timeout     = 300
    self.tol         = 0.001
    self.outputdir   = "/home/su2_test/test_results"
    self.adjust_iter = False

    # Who to e-mail if the test fails
    # self.recipients = ['susquared-dev@lists.stanford.edu']
    self.recipients = ['aniket@stanford.edu']

  def add_recipient(self,email):
    self.recipients.append(email)

  def run_test(self):

    passed     = True
    exceed_tol = False
    timed_out  = False

    # Optionally adjust the number of iterations in the config file
    if self.adjust_iter == True:
      self.do_adjust_iter()

    # Assemble the shell command to run SU2
    command_base = "%s %s"%(os.path.join(self.su2_dir,self.su2_exec), self.cfg_file)
    command      = "%s"%(command_base)

    # Run SU2
    os.chdir(self.cfg_dir) 
    start   = datetime.datetime.now()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)  # This line launches SU2
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
    output = process.stdout.readlines()
    start_solver = False
    delta_vals   = []
    for line in output:
      if not start_solver: # Don't bother parsing anything before --Start solver ---
        if line.find('Start solver') > -1:
          start_solver=True
      else:   # Found the --Start solver --- line; parse the input
        raw_data = line.split()
        try:
          iter_number = int(raw_data[0])
          data        = raw_data[2:]    # Take the last 4 columns for comparison
        except ValueError:
          continue
        except IndexError:
          continue
        
        if iter_number == self.test_iter:  # Found the iteration number we're checking for
          if not len(self.test_vals)==len(data):   # something went wrong... probably bad input
            print "Error in test_vals!"
            passed = False
            break
          for j in range(len(data)):
            delta_vals.append( abs(float(data[j])-self.test_vals[j]) )
            if delta_vals[j] > self.tol:
              exceed_tol = True
              passed     = False
          break

    # Write the test results and output to a file
    os.chdir(self.outputdir)
    file_out = open(self.tag, 'w')
    file_out.write('autotest.py output file, status: ')
    if passed:
      file_out.write('PASSED\n')
      print "%s: PASSED"%self.tag
    else:
      file_out.write('FAILED\n')
      print "%s: FAILED"%self.tag
    file_out.write('                                 ^^^^^^\n')

    file_out.write('execution command: %s\n'%command)

    if timed_out:
      file_out.write('ERROR: Execution timed out. timeout=%d\n'%self.timeout)

    if exceed_tol:
      file_out.write('ERROR: Difference between computed input and test_vals exceeded tolerance. TOL=%f\n'%self.tol)

    file_out.write('test_iter=%d, test_vals: '%self.test_iter)
    for j in self.test_vals:
      file_out.write('%f '%j)
    file_out.write('\n')

    file_out.write('delta_vals: ')
    for j in delta_vals:
      file_out.write('%f '%j)
    file_out.write('\n')

    file_out.write('SU2 output follows:\n')
    file_out.write('=========================================================\n\n')
    for j in output:
      file_out.write(j)
    file_out.close()

    # E-mail recipients if the test failed
    if not passed:
      self.email_failure() 

  def email_failure(self):

    # Create a MIME text file
    fp = open(self.tag, 'rb')
    msg = MIMEText(fp.read())
    fp.close()

    # E-mail header
    me             = 'su2_test@oscarthegrouch.stanford.edu'
    msg['Subject'] = 'Test failure for %s' % self.tag
    msg['From']    = me
    msg['To']      = ', '.join(self.recipients)

    # Send the e-mail
    s = smtplib.SMTP('localhost')
    s.sendmail(me, self.recipients, msg.as_string())
    s.quit()

  def do_adjust_iter(self):
  
    # Read the cfg file
    file_in = open(os.path.join(self.cfg_dir, self.cfg_file), 'r')
    lines   = file_in.readlines()
    file_in.close()
  
    # Rewrite the file with a .autotest extension
    self.cfg_file = "%s.autotest"%self.cfg_file
    file_out = open(os.path.join(self.cfg_dir, self.cfg_file),'w')
    file_out.write('%% This file automatically generated by autotest.py')
    file_out.write('%% Number of iterations changed to %d'%(self.test_iter+1))
    for line in lines:
      if line.find("EXT_ITER")==-1:
        file_out.write(line)
      else:
        file_out.write("EXT_ITER=%d"%(self.test_iter+1))
    file_out.close()

    

if __name__=="__main__":
  '''This program runs SU^2 and ensures that the output matches specified values. This will be used to do nightly checks to make sure nothing is broken. '''

  # Create a testcase object, giving it an identifying tag
  inv_NACA0012 = testcase('inv_NACA0012')

  # Specify the testcase configuration file and directory
  inv_NACA0012.cfg_dir  = "/home/su2_test/su2/trunk/TestCases/inv_NACA0012"
  inv_NACA0012.cfg_file = "inv_NACA0012.cfg"

  # Specify the iteration number and last 4 columns for the test
  inv_NACA0012.test_iter=230
  inv_NACA0012.test_vals=[-5.310748,0.158003,0.327863,0.021429]

  # How long (in seconds) is the program permitted to run for before timing out?
  inv_NACA0012.timeout  = 300

  # How much error can we tolerate for our check
  inv_NACA0012.tol      = 0.001

  # A 3D case
  turb_ONERAM6             = testcase('turb_ONERAM6')
  turb_ONERAM6.cfg_dir     = "/home/su2_test/su2/trunk/TestCases/turb_ONERAM6"
  turb_ONERAM6.cfg_file    = "turb_ONERAM6.cfg"
  turb_ONERAM6.test_iter   = 10
  turb_ONERAM6.test_vals   = [-4.331153,    -10.460434,       0.234600,       0.147112]
  turb_ONERAM6.timeout     = 300
  turb_ONERAM6.tol         = 0.001
  turb_ONERAM6.adjust_iter = True

  # Run the tests
  inv_NACA0012.run_test()
  turb_ONERAM6.run_test()
