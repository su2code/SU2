import sys
import shutil
from pathlib import Path

import numpy as np
from mpi4py import MPI

import SU2
import pysu2

from su2torch.su2_function_mpi import RunCode
from su2torch.su2_function import modify_config


TEMP_CFG_BASENAME = 'tmp_cfg.cfg'


class SU2Numpy:
    """Class that uses the SU2 in-memory python wrapper
    to provide differentiable physics simulations.

    Usage example for scalar output case:

        # define differentiable inputs and outputs in the config
        # with DIFF_INPUTS and DIFF_OUTPUTS fields
        su2 = SU2Numpy('config.cfg')
        inputs = np.array([1.0])
        outputs = su2(inputs)
        # if output is a scalar, we can get the gradient of the output
        # with respect to the inputs by simply doing
        doutput_dinputs = loss.backward()
    """

    def __init__(self, config_file, dims=2, num_zones=1, max_procs=-1):
        """Initialize the SU2 configurations for the provided config file.

        :param config_file: str - The SU2 configuration file name.
        :param dims: int - Number of dimensions for the problem (2D or 3D).
        :param num_zones: int - Number of zones in the simulation (only 1 supported currently).
        :param max_procs: int - Maximum number of MPI processes to use for SU2. If set to -1 (default),
            number of processes will equal batch size. Otherwise, will use floor(max_procs / batch_size)
            processes per item in batch.
            In this case max_procs must be larger than the size of the batch passed in.
        """
        assert num_zones == 1, 'Only supports 1 zone for now.'
        self.num_zones = num_zones
        self.dims = dims
        self.max_procs = max_procs
        self.comms = []
        self.outputs_shape = None
        self.batch_size = -1

        self.config_file = config_file
        base_config = SU2.io.Config(config_file)

        new_configs = {'MATH_PROBLEM': 'DIRECT'}
        self.forward_config = 'forward_{}'.format(TEMP_CFG_BASENAME)
        shutil.copy(config_file, self.forward_config)
        modify_config(base_config, new_configs, outfile=self.forward_config)
        # XXX create local single-rank forward_driver to get access to full mesh,
        #   could be slow for big meshes...
        self.forward_driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                                      self.dims, MPI.COMM_SELF)
        self.num_diff_inputs = self.forward_driver.GetnDiff_Inputs()
        self.num_diff_outputs = self.forward_driver.GetnDiff_Outputs()

        new_configs = {'MATH_PROBLEM': 'DISCRETE_ADJOINT'}
        self.adjoint_config = 'adjoint_{}'.format(TEMP_CFG_BASENAME)
        shutil.copy(config_file, self.adjoint_config)
        modify_config(base_config, new_configs, outfile=self.adjoint_config)

    def __call__(self, *inputs):
        return self.forward(*inputs)

    def forward(self, *inputs):
        """ Runs a batch of SU2 simulations.

        :param inputs: The differentiable inputs for the batch of simulations.
            Number of inputs depends on the number of DIFF_INPUTS set in the configuration file.
            Each input is of shape BATCH_SIZE x SHAPE, where SHAPE is the shape of the given input.
            For example, a batch of 10 scalars would have input shape 10 x 1,
            a batch of 10 vectors of length N would have input shape 10 x N.
        :return: A tuple of tensors with the batch of differentiable outputs.
            Number of outputs depends on the number of DIFF_OUTPUTS set in the configuration file.
            As for the inputs, each output is of shape BATCH_SIZE x SHAPE,
            where SHAPE is the shape of the given output.
            Outputs are always either scalars or vectors.
        """
        if self.num_diff_inputs > 0 and inputs[0].ndim < 2:
            raise TypeError('Input is expected to have first dimension for batch, '
                            'e.g. x[0, :] is first item in batch.')

        self.batch_size = inputs[0].shape[0] if self.num_diff_inputs > 0 else 1
        if 0 <= self.max_procs < self.batch_size:
            raise TypeError('Batch size is larger than max_procs, not enough processes to run batch.')
        procs_per_example = self.max_procs // self.batch_size if self.max_procs > 0 else 1
        outputs = []
        for i in range(self.batch_size):
            x = [z[i] for z in inputs]
            if len(x) != self.num_diff_inputs:
                raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                                .format(len(x), self.config_file, self.num_diff_inputs))
            intercomm = MPI.COMM_SELF.Spawn(sys.executable,
                                            args=str(Path(__file__).parent / 'su2_function_mpi.py'),
                                            maxprocs=procs_per_example)
            self.comms.append(intercomm)
            assert intercomm.Get_rank() == 0, 'Rank is expected to be 0.'
            intercomm.bcast([self.num_zones, self.dims], root=MPI.ROOT)
            intercomm.bcast(RunCode.RUN_FORWARD, root=MPI.ROOT)
            intercomm.bcast(self.forward_config, root=MPI.ROOT)
            intercomm.bcast(x, root=MPI.ROOT)  # TODO Any way to optimize? Potentially big
        for comm in self.comms:
            outputs.append(comm.recv(source=0))
        outputs = tuple(np.concatenate([np.expand_dims(o[i], axis=0) for o in outputs])
                        for i in range(self.num_diff_outputs))
        self.outputs_shape = [o.shape for o in outputs]
        return outputs

    def backward(self, *grad_outputs):
        """Gives the gradient of some scalar loss with respect to the inputs of the previous
        forward call when provided the gradients of this loss with respect to the outputs of
        the forward call.

        :param grad_outputs: Gradients of a scalar loss with respect to the forward outputs.
            For example, if the loss is the sum of the outputs, the grad_outputs should be a all ones.
            This defaults to 1.0 when the output of the forward call is just a scalar (or batch of scalars).
        :return: The gradients of the loss with respect to the forward inputs.
        """
        if len(grad_outputs) == 0 and len(self.outputs_shape) == 1 and self.outputs_shape[0][1] == 1:
            # if no grad_outputs was provided and just one output scalar (or batch of scalars)
            # was used, then use a default grad outputs of 1.0
            grad_outputs = [np.ones(self.outputs_shape[0])]
        elif self.num_diff_outputs != len(grad_outputs):
            raise TypeError('To run backward() you need to provide the gradients of a scalar loss '
                            'with respect to the outputs of the forward pass')

        grads = []
        for i in range(self.batch_size):
            dl = [z[i] for z in grad_outputs]
            self.comms[i].bcast(RunCode.RUN_ADJOINT, root=MPI.ROOT)
            self.comms[i].bcast(self.adjoint_config, root=MPI.ROOT)
            self.comms[i].bcast(dl, root=MPI.ROOT)  # TODO Any way to optimize? Potentially big
        for comm in self.comms:
            grads.append(comm.recv(source=0))
        grads = tuple(np.concatenate([np.expand_dims(g[i], axis=0) for g in grads])
                      for i in range(self.num_diff_inputs))
        return grads

    def __del__(self):
        """Close existing drivers and MPI communicators."""
        for comm in self.comms:
            comm.bcast(RunCode.STOP, root=MPI.ROOT)
            comm.Disconnect()
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()
        # super().__del__()
