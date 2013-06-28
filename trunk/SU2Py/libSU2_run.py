
import os, sys
import subprocess

SU2_RUN = os.environ['SU2_RUN'] 
sys.path.append( SU2_RUN )

base_Command = os.path.join(SU2_RUN,'%s')
mpi_Command = 'mpirun -np %i %s'

def SU2_DDC(config,partitions=0):
    
    # must run with rank 1
    partitions = min([1,partitions])
    
    the_Command = 'SU2_DDC ' + config
    the_Command = build_command( the_Command , partitions )
    run_command( the_Command )
    
    return

def SU2_CFD(config,partitions=0):
    
    the_Command = 'SU2_CFD ' + config
    the_Command = build_command( the_Command , partitions )
    run_command( the_Command )
    
    return

def SU2_MAC(config,partitions=0):
    
    # must run with rank 1
    partitions = min([1,partitions])    
    
    the_Command = 'SU2_MAC ' + config
    the_Command = build_command( the_Command , partitions )
    run_command( the_Command )
    
    return

def SU2_MDC(config,partitions=0):
    
    ## must run with rank 1
    #partitions = min([1,partitions])
    
    the_Command = 'SU2_MDC ' + config
    the_Command = build_command( the_Command , partitions )
    print 'the_Command: ', the_Command #-AA
    run_command( the_Command )
    
    return

def SU2_GPC(config,partitions=0):
    
    ## must run with rank 1
    #partitions = min([1,partitions])      
    
    the_Command = 'SU2_GPC ' + config
    the_Command = build_command( the_Command , partitions )
    run_command( the_Command )
    
    return

def SU2_SMC(config,partitions=0):
    
    # must run with rank 1
    partitions = min([1,partitions])      
    
    the_Command = 'SU2_SMC ' + config
    the_Command = build_command( the_Command , partitions )
    run_command( the_Command )
    
    return

def SU2_PBC(config,partitions=0):
    
    # must run with rank 1
    partitions = min([1,partitions])      
    
    the_Command = 'SU2_PBC ' + config
    the_Command = build_command( the_Command , partitions )
    run_command( the_Command )
    
    return
        

# ---------------------------------------------------------
#  Helper functions
# ---------------------------------------------------------

def build_command( the_Command , partitions=0 ):
    the_Command = base_Command % the_Command
    if partitions > 0:
        the_Command = mpi_Command % (partitions,the_Command)
    return the_Command

def run_command( Command ):
    
    try:
        return_code = subprocess.call( Command , shell=True ,
                                       stdout=sys.stdout    , 
                                       stderr=sys.stderr     )
        
        if return_code < 0:
            raise SystemExit , "SU2 process was terminated by signal '%s'" % -return_code
        elif return_code > 0:
            raise Exception , "SU2 process returned error '%s'" % return_code
            
    except OSError as the_error:
        raise OSError , "SU2 process failed: '%s'" % the_error
        
    return return_code

