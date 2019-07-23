import os
import sys
import shutil
from pathlib import Path

import torch
import numpy as np

from mpi4py import MPI

import SU2
import pysu2
import pysu2ad

from su2_function_mpi import RunCode, run_forward, run_adjoint


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
        self.mpi_params = np.array([self.num_zones, self.dims, self.num_diff_outputs],
                                   dtype=np.int)

        new_configs = {'MATH_PROBLEM': 'DIRECT'}
        self.forward_config = 'forward_{}_{}'.format(os.getpid(), TEMP_CFG_BASENAME)
        shutil.copy(config_file, self.forward_config)
        modify_config(base_config, new_configs, outfile=self.forward_config)
        # XXX create local single-rank forward_driver to get access to full mesh,
        #   could be slow for big meshes...
        self.forward_driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                                      self.dims, MPI.COMM_SELF)

        new_configs = {'MATH_PROBLEM': 'DISCRETE_ADJOINT'}
        self.adjoint_config = 'adjoint_{}_{}'.format(os.getpid(), TEMP_CFG_BASENAME)
        shutil.copy(config_file, self.adjoint_config)
        modify_config(base_config, new_configs, outfile=self.adjoint_config)

        if num_procs > 1:
            self.intercomm = MPI.COMM_SELF.Spawn(sys.executable,
                                                 args=str(Path(__file__).parent / 'su2_function_mpi.py'),
                                                 maxprocs=num_procs-1)  # -1 because this process is also a worker
            self.comm = self.intercomm.Merge(high=False)

        else:
            self.comm = MPI.COMM_WORLD

        assert self.comm.Get_rank() == 0, 'Rank is expected to be 0.'
        self.comm.bcast([self.num_zones, self.dims, self.num_diff_outputs], root=0)

    def forward(self, *inputs, **kwargs):
        if len(inputs) != self.num_diff_inputs:
            raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                            .format(len(inputs), self.config_file, self.num_diff_inputs))
        self.save_for_backward(*inputs)

        self.comm.bcast(RunCode.RUN_FORWARD, root=0)
        self.comm.bcast(self.forward_config, root=0)
        self.comm.bcast(inputs, root=0)  # TODO Any way to optimize? Potentially big

        driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                         self.dims, self.comm)
        outputs = run_forward(self.comm, driver, inputs, self.num_diff_outputs)
        return outputs

    def backward(self, *grad_outputs, **kwargs):
        self.comm.bcast(RunCode.RUN_ADJOINT, root=0)
        self.comm.bcast(self.adjoint_config, root=0)
        self.comm.bcast(grad_outputs, root=0)  # TODO Any way to optimize? Potentially big

        adjoint_driver = pysu2ad.CDiscAdjSinglezoneDriver(self.adjoint_config, self.num_zones,
                                                          self.dims, self.comm)
        grads = run_adjoint(self.comm, adjoint_driver, self.saved_tensors, grad_outputs)
        return grads

    def __del__(self):
        os.remove(self.forward_config)
        os.remove(self.adjoint_config)
        if self.comm.Get_size() > 1:
            self.comm.bcast(RunCode.STOP, root=0)
            # TODO Disconnects hanging on cluster
            # self.comm.Disconnect()
            # self.intercomm.Disconnect()
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
