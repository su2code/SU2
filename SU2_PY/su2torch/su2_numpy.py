import math

import numpy as np
from mpi4py import MPI  # Must be imported before pysu2 or else MPI error happens at some point

import pysu2
from su2torch.su2_function_mpi import RunCode


# TODO Outdated, update to match the Torch version (su2_function.py)
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

    def __init__(self, config_file, dims=2, num_zones=1):
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
        assert MPI.COMM_WORLD.Get_rank() == 0, 'Not rank 0 in comm'
        self.comm = MPI.COMM_WORLD
        self.workers = self.comm.Get_size() - 1
        assert self.workers > 0, 'Need at least 1 master and 1 worker process.'
        self.num_zones = num_zones
        self.dims = dims
        self.outputs_shape = None
        self.batch_size = -1

        self.forward_config = config_file
        self.forward_driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                                      self.dims, MPI.COMM_SELF)
        self.num_diff_inputs = self.forward_driver.GetnDiff_Inputs()
        self.num_diff_outputs = self.forward_driver.GetnDiff_Outputs()

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
        if len(inputs) != self.num_diff_inputs:
            raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                            .format(len(inputs), self.forward_config, self.num_diff_inputs))
        if self.num_diff_inputs > 0 and inputs[0].ndim < 2:
            raise TypeError('Input is expected to have first dimension for batch, '
                            'e.g. x[0, :] is first item in batch.')
        self.batch_size = inputs[0].shape[0] if self.num_diff_inputs > 0 else 1
        if 0 <= self.workers < self.batch_size:
            raise TypeError('Batch size is larger than number of workers, not enough processes to run batch.')
        procs_per_example = math.ceil(self.workers / self.batch_size)

        self.comm.bcast(RunCode.RUN_FORWARD, root=0)
        self.comm.bcast([self.num_zones, self.dims, self.forward_config, inputs], root=0)
        outputs = []
        for i in range(self.batch_size):
            output = self.comm.recv(source=1 + i * procs_per_example)
            outputs.append(output)
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

        procs_per_example = math.ceil(self.workers / self.batch_size)
        self.comm.bcast(RunCode.RUN_ADJOINT, root=0)
        self.comm.bcast(grad_outputs, root=0)
        grads = []
        for i in range(self.batch_size):
            grad = self.comm.recv(source=1 + i * procs_per_example)
            grads.append(grad)
        grads = tuple(np.concatenate([np.expand_dims(g[i], axis=0) for g in grads])
                      for i in range(self.num_diff_inputs))
        return grads

    def __del__(self):
        """Close existing drivers and MPI communicators."""
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()

