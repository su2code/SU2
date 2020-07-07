import math
import sys
from typing import Tuple

import os
import torch
from torch.nn.utils.rnn import pad_sequence

import SU2
from mpi4py import MPI  # Must be imported before pysu2 or else MPI error happens at some point
import pysu2

from su2torch.su2_function_mpi import RunCode, non_busy_post, non_busy_wait, modify_config

_global_max_ppe = -1


class SU2Module(torch.nn.Module):
    def __init__(self, config_file: str, mesh_file: str, dims: int = 2, num_zones: int = 1) -> None:
        """ Initialize the SU2 configurations for the provided config file.

        :param config_file: str - The SU2 configuration file name.
        :param mesh_file: str - Optional parameter, if not set defaults to the mesh filename set in the config file.
            Can be used to run a batch with different meshes for each sample.
            Passing in mesh_file with batch_index parameter in string format (e.g., 'b{batch_index}_mesh.su2')
            causes each element in batch to get assigned to the correct mesh file (0 indexed).
            If running multiple processes in parallel, take care to name each mesh file uniquely to avoid conflicts
            (e.g., unique = str(os.getpid()); mesh_file = 'b{batch_index}_' + unique + '_mesh.su2').
        :param dims: int - Number of dimensions for the problem (2D or 3D).
        :param num_zones: int - Number of zones in the simulation (only 1 supported currently).
        """
        super().__init__()
        assert num_zones == 1, 'Only supports 1 zone for now.'
        assert MPI.COMM_WORLD.Get_rank() == 0, 'Not rank 0 in comm'
        assert _global_max_ppe > 0, 'Before running SU2Function, a (single) call to activate_su2_mpi is needed.'
        self.num_zones = num_zones
        self.dims = dims
        self.mesh_file = mesh_file

        self.forward_config = config_file
        self.forward_driver = None

    def forward(self, *inputs: torch.Tensor) -> Tuple[torch.Tensor, ...]:
        return SU2Function.apply(*inputs, self.forward_config, self.mesh_file,
                                 self.num_zones, self.dims, self.set_forward_driver)

    def get_forward_driver(self):
        if self.forward_driver is None:
            raise AttributeError('Forward driver is only set after running forward()')
        return self.forward_driver

    def set_forward_driver(self, f):
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()
        self.forward_driver = f

    def __del__(self):
        """Close existing drivers and MPI communicators."""
        if hasattr(self, 'forward_driver') and self.forward_driver is not None:
            self.forward_driver.Postprocessing()


class SU2Function(torch.autograd.Function):
    num_params = 5

    @staticmethod
    def forward(ctx, *inputs):
        non_busy_post(MPI.COMM_WORLD)
        x = inputs[:-SU2Function.num_params]
        forward_config, mesh_file, num_zones, dims, set_forward_driver_hook = inputs[-SU2Function.num_params:]

        if x[0].dim() < 2:
            raise TypeError('Input is expected to have first dimension for batch, '
                            'e.g. x[0, :] is first item in batch.')
        batch_size = x[0].shape[0]
        max_ppe = _global_max_ppe
        workers = MPI.COMM_WORLD.Get_size() - 1
        if 0 <= workers < batch_size:
            raise TypeError('Batch size is larger than number of workers, not enough processes to run batch.')

        MPI.COMM_WORLD.bcast(RunCode.RUN_FORWARD, root=0)
        procs_per_example = min(max_ppe, math.ceil(workers / batch_size))
        MPI.COMM_WORLD.bcast([num_zones, dims, forward_config, mesh_file, procs_per_example, x], root=0)

        # instantiate forward_driver while workers work
        worker_forward_config = MPI.COMM_WORLD.recv(source=1)
        forward_driver = pysu2.CSinglezoneDriver(worker_forward_config, num_zones, dims, MPI.COMM_SELF)
        num_diff_inputs = forward_driver.GetnDiff_Inputs()
        num_diff_outputs = forward_driver.GetnDiff_Outputs()
        assert num_diff_inputs > 0 and num_diff_outputs > 0, \
            'Need to define at least one differentiable input and output. ' \
            'To run without differentiation, use the SU2Numpy class.'
        if len(x) != num_diff_inputs:
            raise TypeError(f'{len(x)} inputs were provided, but the config file '
                            f'({forward_config}) defines {num_diff_inputs} diff inputs.')
        set_forward_driver_hook(forward_driver)
        ctx.num_diff_inputs = num_diff_inputs

        outputs = []
        non_busy_wait(MPI.COMM_WORLD)
        for i in range(batch_size):
            output = MPI.COMM_WORLD.recv(source=1 + i * procs_per_example)
            outputs.append(output)
        outputs = tuple(pad_sequence([o[i] for o in outputs], batch_first=True)
                        for i in range(num_diff_outputs))
        return outputs

    @staticmethod
    def backward(ctx, *grad_outputs):
        non_busy_post(MPI.COMM_WORLD)
        max_ppe = _global_max_ppe
        workers = MPI.COMM_WORLD.Get_size() - 1
        MPI.COMM_WORLD.bcast(RunCode.RUN_ADJOINT, root=0)
        MPI.COMM_WORLD.bcast(grad_outputs, root=0)
        batch_size = grad_outputs[0].shape[0]
        procs_per_example = min(max_ppe, math.ceil(workers / batch_size))
        non_busy_wait(MPI.COMM_WORLD)
        grads = []
        for i in range(batch_size):
            grad = MPI.COMM_WORLD.recv(source=1 + i * procs_per_example)
            grads.append(grad)
        grads = tuple(pad_sequence([g[i] for g in grads], batch_first=True)
                      for i in range(ctx.num_diff_inputs))
        return grads + (None,) * SU2Function.num_params
