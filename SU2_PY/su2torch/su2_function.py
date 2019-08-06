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
    """PyTorch function that uses the SU2 in-memory python wrapper
    to provide differentiable physics simulations.

    Usage example:

        # define differentiable inputs and outputs in the config
        # with DIFF_INPUTS and DIFF_OUTPUTS fields
        su2 = SU2Function('config.cfg')
        inputs = torch.tensor([1.0], requires_grad=True)
        outputs = su2(inputs)
        loss = some_torch_operations(outputs)
        loss.backward()
        # now we have gradient of loss with respect to inputs here:
        inputs.grad
    """

    def __init__(self, config_file, dims=2, num_zones=1, max_procs=-1):
        """ Initialize the SU2 configurations for the provided config file.

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
        if inputs[0].dim() < 2:
            raise TypeError('Input is expected to have first dimension for batch, '
                            'e.g. x[0, :] is first item in batch.')

        batch_size = inputs[0].shape[0]
        if 0 <= self.max_procs < batch_size:
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

    def backward(self, *grad_outputs):
        """Called implicitly by PyTorch on the backward pass.

        Takes in the gradients of some scalar loss with respect to the outputs of the previous forward
        call and returns the gradients of this loss with respect to the inputs of the forward call.

        :param grad_outputs: Gradients of a scalar loss with respect to the forward outputs.
        :return: The gradients of the loss with respect to the forward inputs.
        """
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
        """Finish existing drivers and MPI communicators."""
        for comm in self.comms:
            comm.bcast(RunCode.STOP, root=MPI.ROOT)
            comm.Disconnect()
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()
        # super().__del__()


def modify_config(config, new_params, outfile=None):
    """Modify a config, saving the modifications to outfile if provided."""
    temp_config = config.copy()
    for k, v in new_params.items():
        temp_config[k] = v
    if outfile is not None:
        temp_config.write(outfile)
    return temp_config
