#!/usr/bin/env python

## \file interface.py
#  \brief python package interfacing with the SU2 suite
#  \author T. Lukaczyk, F. Palacios
#  \version 7.5.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy
import subprocess
from ..io import Config
from ..util import which

# ------------------------------------------------------------
#  Setup
# ------------------------------------------------------------

SU2_RUN = os.environ['SU2_RUN']
sys.path.append( SU2_RUN )
quote = '"' if sys.platform == 'win32' else ''

# SU2 suite run command template
base_Command = os.path.join(SU2_RUN, '%s')

# check for slurm
slurm_job = 'SLURM_JOBID' in os.environ

# Check for custom mpi command
user_defined = 'SU2_MPI_COMMAND' in os.environ

# set mpi command
if user_defined:
    mpi_Command = os.environ['SU2_MPI_COMMAND']
elif slurm_job:
    mpi_Command = 'srun -n %i %s'
elif not which('mpirun') is None:
    mpi_Command = 'mpirun -n %i %s'
elif not which('mpiexec') is None:
    mpi_Command = 'mpiexec -n %i %s'
else:
    mpi_Command = ''

from .. import EvaluationFailure, DivergenceFailure
return_code_map = {
    1 : EvaluationFailure ,
    2 : DivergenceFailure ,
}

# ------------------------------------------------------------
#  SU2 Suite Interface Functions
# ------------------------------------------------------------

def CFD(config):
    """ run SU2_CFD
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)

    direct_diff = not konfig.get('DIRECT_DIFF',"") in ["NONE", ""]

    auto_diff = konfig.MATH_PROBLEM == 'DISCRETE_ADJOINT'

    if direct_diff:
        tempname = 'config_CFD_DIRECTDIFF.cfg'

        konfig.dump(tempname)

        processes = konfig['NUMBER_PART']

        the_Command = 'SU2_CFD_DIRECTDIFF%s %s' % (quote, tempname)

    elif auto_diff:
        tempname = 'config_CFD_AD.cfg'
        konfig.dump(tempname)

        processes = konfig['NUMBER_PART']

        the_Command = 'SU2_CFD_AD%s %s' % (quote, tempname)

    else:
        tempname = 'config_CFD.cfg'
        konfig.dump(tempname)

        processes = konfig['NUMBER_PART']

        the_Command = 'SU2_CFD%s %s' % (quote, tempname)

    the_Command = build_command( the_Command, processes )
    run_command( the_Command )

    #os.remove(tempname)

    return

def DEF(config):
    """ run SU2_DEF
        partitions set by config.NUMBER_PART
        forced to run in serial, expects merged mesh input
    """
    konfig = copy.deepcopy(config)

    tempname = 'config_DEF.cfg'
    konfig.dump(tempname)

    # must run with rank 1
    processes = konfig['NUMBER_PART']

    the_Command = 'SU2_DEF%s %s' % (quote, tempname)
    the_Command = build_command( the_Command, processes )
    run_command( the_Command )

    #os.remove(tempname)

    return

def DOT(config):
    """ run SU2_DOT
        partitions set by config.NUMBER_PART
    """
    konfig = copy.deepcopy(config)

    auto_diff = konfig.MATH_PROBLEM == 'DISCRETE_ADJOINT' or konfig.get('AUTO_DIFF','NO') == 'YES'

    if auto_diff:

        tempname = 'config_DOT_AD.cfg'
        konfig.dump(tempname)

        processes = konfig['NUMBER_PART']

        the_Command = 'SU2_DOT_AD%s %s' % (quote, tempname)
    else:

        tempname = 'config_DOT.cfg'
        konfig.dump(tempname)

        processes = konfig['NUMBER_PART']

        the_Command = 'SU2_DOT%s %s' % (quote, tempname)

    the_Command = build_command( the_Command, processes )
    run_command( the_Command )

    #os.remove(tempname)

    return

def GEO(config):
    """ run SU2_GEO
        partitions set by config.NUMBER_PART
        forced to run in serial
    """
    konfig = copy.deepcopy(config)

    tempname = 'config_GEO.cfg'
    konfig.dump(tempname)

    # must run with rank 1
    processes = konfig['NUMBER_PART']

    the_Command = 'SU2_GEO%s %s' % (quote, tempname)
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

    the_Command = 'SU2_SOL%s %s' % (quote, tempname)
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )

    #os.remove(tempname)

    return

def SOL_FSI(config):
    """ run SU2_SOL for FSI problems
      partitions set by config.NUMBER_PART
    """

    konfig = copy.deepcopy(config)

    tempname = 'config_SOL.cfg'
    konfig.dump(tempname)

    # must run with rank 1
    processes = konfig['NUMBER_PART']

    the_Command = 'SU2_SOL%s %s 2' % (quote, tempname)
    the_Command = build_command( the_Command , processes )
    run_command( the_Command )

    #os.remove(tempname)

    return


# ------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------

def build_command( the_Command , processes=0 ):
    """ builds an mpi command for given number of processes """
    the_Command = quote + (base_Command % the_Command)
    if processes > 1:
        if not mpi_Command:
            raise RuntimeError('could not find an mpi interface')
        the_Command = mpi_Command % (processes,the_Command)
    return the_Command

def run_command( Command ):
    """ runs os command with subprocess
        checks for errors from command
    """

    sys.stdout.flush()

    proc = subprocess.Popen( Command, shell=True    ,
                             stdout=sys.stdout      ,
                             stderr=subprocess.PIPE  )
    return_code = proc.wait()
    message = proc.stderr.read().decode()

    if return_code < 0:
        message = "SU2 process was terminated by signal '%s'\n%s" % (-return_code,message)
        raise SystemExit(message)
    elif return_code > 0:
        message = "Path = %s\nCommand = %s\nSU2 process returned error '%s'\n%s" % (os.path.abspath(','),Command,return_code,message)
        if return_code in return_code_map.keys():
            exception = return_code_map[return_code]
        else:
            exception = RuntimeError
        raise exception(message)
    else:
        sys.stdout.write(message)

    return return_code

