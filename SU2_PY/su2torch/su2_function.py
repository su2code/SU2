import sys
import shutil
from pathlib import Path

import torch

from mpi4py import MPI

import SU2
import pysu2

from su2torch.su2_function_mpi import RunCode


TEMP_CFG_BASENAME = 'tmp_cfg.cfg'


class SU2Function(torch.autograd.Function):
    """Function that uses the SU2 in-memory python wrapper."""

    def __init__(self, config_file, dims=2, num_zones=1, max_procs=-1):
        assert num_zones == 1, 'Only supports 1 zone for now.'
        self.num_zones = num_zones
        self.dims = dims
        self.max_procs = max_procs
        self.comms = []

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

    def forward(self, *inputs, **kwargs):
        if inputs[0].dim() < 2:
            raise TypeError('Input is expected to have first dimension for batch, '
                            'e.g. x[0, :] is first item in batch.')

        batch_size = inputs[0].shape[0]
        if 0 < self.max_procs < batch_size:
            raise TypeError('Batch size is larger than max_procs, not enough processes to run batch.')
        self.save_for_backward(*inputs)
        procs_per_example = self.max_procs // batch_size if self.max_procs > 0 else 1
        outputs = []
        for i in range(batch_size):
            x = [z[i] for z in inputs]
            if len(x) != self.num_diff_inputs:
                raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                                .format(len(x), self.config_file, self.num_diff_inputs))
            intercomm = MPI.COMM_SELF.Spawn(sys.executable,
                                            args=str(Path(__file__).parent / 'su2_function_mpi.py'),
                                            maxprocs=procs_per_example)
            self.comms.append(intercomm)
            assert intercomm.Get_rank() == 0, 'Rank is expected to be 0.'
            intercomm.bcast([self.num_zones, self.dims, self.num_diff_outputs], root=MPI.ROOT)
            intercomm.bcast(RunCode.RUN_FORWARD, root=MPI.ROOT)
            intercomm.bcast(self.forward_config, root=MPI.ROOT)
            intercomm.bcast(x, root=MPI.ROOT)  # TODO Any way to optimize? Potentially big
        for comm in self.comms:
            outputs.append(comm.recv(source=0))
        outputs = tuple(torch.cat([o[i].unsqueeze(0) for o in outputs]) for i in range(len(outputs[0])))
        return outputs

    def backward(self, *grad_outputs, **kwargs):
        batch_size = grad_outputs[0].shape[0]
        grads = []
        for i in range(batch_size):
            dl = [z[i] for z in grad_outputs]
            x = [z[i] for z in self.saved_tensors]
            self.comms[i].bcast(RunCode.RUN_ADJOINT, root=MPI.ROOT)
            self.comms[i].bcast(self.adjoint_config, root=MPI.ROOT)
            self.comms[i].bcast(dl, root=MPI.ROOT)  # TODO Any way to optimize? Potentially big
        for comm in self.comms:
            grads.append(comm.recv(source=0))
        grads = tuple(torch.cat([g[i].unsqueeze(0) for g in grads]) for i in range(len(grads[0])))
        return grads

    def __del__(self):
        for comm in self.comms:
            comm.bcast(RunCode.STOP, root=MPI.ROOT)
            comm.Disconnect()
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
