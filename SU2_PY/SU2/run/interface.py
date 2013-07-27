## \file interface.py
#  \brief python package interfacing with the SU2 suite
#  \author Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 2.0.6
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import subprocess
from ..io import Config

# ------------------------------------------------------------
#  Setup
# ------------------------------------------------------------

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN )

# SU2 suite run command template
base_Command = os.path.join(SU2_RUN,'%s')

# check for slurm
slurm_job = os.environ.has_key('SLURM_JOBID')
    
# set mpi command
if slurm_job:
    mpi_Command = 'srun -n %i %s'
else:
    mpi_Command = 'mpirun -np %i %s'
    
# ------------------------------------------------------------
#  SU2 Suite Interface Functions
# ------------------------------------------------------------

def DDC(config):
    """ run SU2_DDC 
        partitions set by config.NUMBER_PART
        currently forced to run serially
    """
    # local copy
    konfig = copy.deepcopy(config)
    
    tempname = 'config_DDC.cfg'
    konfig.dump(tempname)
  
    processes = konfig['NUMBER_PART']
    
    the_Command = 'SU2_DDC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def CFD(config):
    """ run SU2_CFD
        partitions set by config.NUMBER_PART
    """
    
    konfig = copy.deepcopy(config)
    
    tempname = 'config_CFD.cfg'
    konfig.dump(tempname)
    
    processes = konfig['NUMBER_PART']
    
    the_Command = 'SU2_CFD ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def MAC(config):
    """ run SU2_MAC
        partitions set by config.NUMBER_PART
        currently forced to run serially
    """    
    konfig = copy.deepcopy(config)
    
    tempname = 'config_MAC.cfg'
    konfig.dump(tempname)
    
    # must run with rank 1
    processes = konfig['NUMBER_PART']
    processes = min([1,processes])    
    
    the_Command = 'SU2_MAC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def MDC(config):
    """ run SU2_MDC
        partitions set by config.NUMBER_PART
        forced to run in serial, expects merged mesh input
    """
    konfig = copy.deepcopy(config)
    
    tempname = 'config_MDC.cfg'
    konfig.dump(tempname) 
    
    # must run with rank 1
    processes = konfig['NUMBER_PART']
    #processes = min([1,processes])  # hack - TWL
    
    the_Command = 'SU2_MDC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def GPC(config):
    """ run SU2_GPC
        partitions set by config.NUMBER_PART
    """    
    konfig = copy.deepcopy(config)
    
    tempname = 'config_GPC.cfg'
    konfig.dump(tempname)   
    
    processes = konfig['NUMBER_PART']
    
    the_Command = 'SU2_GPC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def GDC(config):
    """ run SU2_GDC
        partitions set by config.NUMBER_PART
        forced to run in serial
    """    
    konfig = copy.deepcopy(config)
    
    tempname = 'config_GDC.cfg'
    konfig.dump(tempname)   
    
    # must run with rank 1
    processes = konfig['NUMBER_PART']
    processes = min([1,processes])       
        
    the_Command = 'SU2_GDC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def SMC(config):
    """ run SU2_SMC
        partitions set by config.NUMBER_PART
    """    
    konfig = copy.deepcopy(config)    
    
    tempname = 'config_SMC.cfg'
    konfig.dump(tempname)   
    
    # must run with rank 1
    processes = konfig['NUMBER_PART']
    processes = min([1,processes])       
    
    the_Command = 'SU2_SMC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def PBC(config):
    """ run SU2_PBC
        partitions set by config.NUMBER_PART
        currently forced to run serially
    """    
    konfig = copy.deepcopy(config)
    
    tempname = 'config_PBC.cfg'
    konfig.dump(tempname)
    
    # must run with rank 1
    processes = konfig['NUMBER_PART']
    processes = min([1,processes])      
    
    the_Command = 'SU2_PBC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return
        
def SOL(config):
    """ run SU2_SOL
      partitions set by config.NUMBER_PART
    """
  
    konfig = copy.deepcopy(config)
    
    tempname = 'config_SOL.cfg'
    konfig.dump(tempname)
  
    # must run with rank 1
    processes = konfig['NUMBER_PART']
    
    the_Command = 'SU2_SOL ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

# ------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------

def build_command( the_Command , processes=0 ):
    """ builds an mpi command for given number of processes """
    the_Command = base_Command % the_Command
    if processes > 0:
        the_Command = mpi_Command % (processes,the_Command)
    return the_Command

def run_command( Command ):
    """ runs os command with subprocess
        checks for errors from command
    """
    
    print "the command: %s" % Command
    print "the location: %s" % os.getcwd()
    sys.stdout.flush()
    
    proc = subprocess.Popen( Command, shell=True    ,
                             stdout=sys.stdout      , 
                             stderr=subprocess.PIPE  )
    return_code = proc.wait()
    message = proc.stderr.read()
    
    if return_code < 0:
        message = "SU2 process was terminated by signal '%s'\n%s" % (-return_code,message)
        raise SystemExit , message
    elif return_code > 0:
        message = "Path = %s\nCommand = %s\nSU2 process returned error '%s'\n%s" % (os.path.abspath(','),Command,return_code,message)
        raise Exception , message
    else:
        sys.stdout.write(message)
            
    return return_code

