import sys
import shutil
from pathlib import Path

import torch

from mpi4py import MPI

import SU2
import pysu2
import pysu2ad

from su2torch.su2_function_mpi import RunCode, run_forward, run_adjoint


TEMP_CFG_BASENAME = 'tmp_cfg.cfg'


class SU2Function(torch.autograd.Function):
    """Function that uses the SU2 in-memory python wrapper."""

    def __init__(self, config_file, dims=2, num_zones=1, num_procs=1):
        assert num_zones == 1, 'Only supports 1 zone for now.'
        self.num_zones = num_zones
        self.dims = dims

        self.config_file = config_file
        base_config = SU2.io.Config(config_file)
        self.num_diff_inputs = len(base_config['DIFF_INPUTS'].split(','))
        self.num_diff_outputs = len(base_config['DIFF_OUTPUTS'].split(','))

        new_configs = {'MATH_PROBLEM': 'DIRECT'}
        self.forward_config = 'forward_{}'.format(TEMP_CFG_BASENAME)
        shutil.copy(config_file, self.forward_config)
        modify_config(base_config, new_configs, outfile=self.forward_config)
        # XXX create local single-rank forward_driver to get access to full mesh,
        #   could be slow for big meshes...
        self.forward_driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                                      self.dims, MPI.COMM_SELF)

        new_configs = {'MATH_PROBLEM': 'DISCRETE_ADJOINT'}
        self.adjoint_config = 'adjoint_{}'.format(TEMP_CFG_BASENAME)
        shutil.copy(config_file, self.adjoint_config)
        modify_config(base_config, new_configs, outfile=self.adjoint_config)

        self.num_procs = num_procs
        if num_procs > 1:
            intercomm = MPI.COMM_SELF.Spawn(sys.executable,
                                            args=str(Path(__file__).parent / 'su2_function_mpi.py'),
                                            maxprocs=num_procs)
            self.comm = intercomm
            assert self.comm.Get_rank() == 0, 'Rank is expected to be 0.'
            self.comm.bcast([self.num_zones, self.dims, self.num_diff_outputs], root=MPI.ROOT)
        else:
            self.comm = MPI.COMM_WORLD

    def forward(self, *inputs, **kwargs):
        # TODO Take batching into account
        if len(inputs) != self.num_diff_inputs:
            raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                            .format(len(inputs), self.config_file, self.num_diff_inputs))
        self.save_for_backward(*inputs)
        if self.num_procs > 1:
            self.comm.bcast(RunCode.RUN_FORWARD, root=MPI.ROOT)
            self.comm.bcast(self.forward_config, root=MPI.ROOT)
            self.comm.bcast(inputs, root=MPI.ROOT)  # TODO Any way to optimize? Potentially big
            outputs = self.comm.recv(source=0)
        else:
            driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones, self.dims, self.comm)
            outputs = run_forward(self.comm, driver, inputs, self.num_diff_outputs)
            driver.Postprocessing()
        return outputs

    def backward(self, *grad_outputs, **kwargs):
        if self.num_procs > 1:
            self.comm.bcast(RunCode.RUN_ADJOINT, root=MPI.ROOT)
            self.comm.bcast(self.adjoint_config, root=MPI.ROOT)
            self.comm.bcast(grad_outputs, root=MPI.ROOT)  # TODO Any way to optimize? Potentially big
            grads = self.comm.recv(source=0)
        else:
            driver = pysu2ad.CDiscAdjSinglezoneDriver(self.adjoint_config, self.num_zones, self.dims, self.comm)
            grads = run_adjoint(self.comm, driver, self.saved_tensors, grad_outputs)
            driver.Postprocessing()
        return grads

    def __del__(self):
        if self.num_procs > 1:
            self.comm.bcast(RunCode.STOP, root=MPI.ROOT)
            self.comm.Disconnect()
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()
        # super().__del__()


def modify_config(config, new_params, outfile=None):
    temp_config = config.copy()
    for k, v in new_params.items():
        temp_config[k] = v
    if outfile is not None:
        temp_config.write(outfile)
    return temp_config
