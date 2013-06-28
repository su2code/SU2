
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
    mpi_Command = 'srun -n%i %s'
else:
    mpi_Command = 'mpirun -np %i %s'
    
# ------------------------------------------------------------
#  SU2 Suite Interface Functions
# ------------------------------------------------------------

def DDC(config):
    config = copy.deepcopy(config)
    
    tempname = 'config_DDC.cfg'
    config.dump(tempname)
    
    processes = config['NUMBER_PART']
    
    # must run with rank 1
    processes = min([1,processes])
    
    the_Command = 'SU2_DDC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def CFD(config):
    config = copy.deepcopy(config)
    
    tempname = 'config_CFD.cfg'
    config.dump(tempname)
    
    processes = config['NUMBER_PART']
    
    the_Command = 'SU2_CFD ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def MAC(config):
    config = copy.deepcopy(config)
    
    tempname = 'config_MAC.cfg'
    config.dump(tempname)
    
    # must run with rank 1
    processes = config['NUMBER_PART']
    processes = min([1,processes])    
    
    the_Command = 'SU2_MAC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def MDC(config):
    config = copy.deepcopy(config)
    
    tempname = 'config_MDC.cfg'
    config.dump(tempname) 
    
    processes = config['NUMBER_PART']
    
    the_Command = 'SU2_MDC ' + tempname
    the_Command = build_command( the_Command , processes )
    print 'the_Command: ', the_Command #-AA
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return
    
def GPC(config):
    config = copy.deepcopy(config)
    
    tempname = 'config_GPC.cfg'
    config.dump(tempname)   
    
    processes = config['NUMBER_PART']
    
    the_Command = 'SU2_GPC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def SMC(config):
    config = copy.deepcopy(config)    
    
    tempname = 'config_SMC.cfg'
    config.dump(tempname)   
    
    # must run with rank 1
    processes = config['NUMBER_PART']
    processes = min([1,processes])       
    
    the_Command = 'SU2_SMC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return

def PBC(config):
    config = copy.deepcopy(config)
    
    tempname = 'config_PBC.cfg'
    config.dump(tempname)
    
    # must run with rank 1
    processes = config['NUMBER_PART']
    processes = min([1,processes])      
    
    the_Command = 'SU2_PBC ' + tempname
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )
    
    #os.remove(tempname)
    
    return
        

# ------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------

def build_command( the_Command , processes=0 ):
    the_Command = base_Command % the_Command
    if processes > 0:
        the_Command = mpi_Command % (processes,the_Command)
    return the_Command

def run_command( Command ):
    
    proc = subprocess.Popen( Command, shell=True    ,
                             stdout=sys.stdout      , 
                             stderr=subprocess.PIPE  )
    return_code = proc.wait()
    message = proc.stderr.read()
    
    if return_code < 0:
        message = "SU2 process was terminated by signal '%s'\n%s" % (-return_code,message)
        raise SystemExit , message
    elif return_code > 0:
        message = "SU2 process returned error '%s'\n%s" % (return_code,message)
        raise Exception , message
    else:
        sys.stdout.write(message)
            
    return return_code

