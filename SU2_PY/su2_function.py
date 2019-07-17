import sys
import shutil
from pathlib import Path
from optparse import OptionParser

import torch
import numpy as np
from mpi4py import MPI

import SU2
import pysu2ad as pysu2
from su2_function_mpi import RunCode, run_forward, run_adjoint


TEMP_CFG_BASENAME = 'tmp_cfg.cfg'


class SU2MPIFunction(torch.autograd.Function):
    """Function that uses the SU2 in-memory python wrapper."""

    def __init__(self, config_file, dims=2, num_zones=1, num_procs=4):
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
        self.forward_config = 'forward_' + TEMP_CFG_BASENAME
        shutil.copy(config_file, self.forward_config)
        modify_config(base_config, new_configs, outfile=self.forward_config)
        self.forward_driver = None

        new_configs = {'MATH_PROBLEM': 'DISCRETE_ADJOINT'}
        self.adjoint_config = 'adjoint_' + TEMP_CFG_BASENAME
        shutil.copy(config_file, self.adjoint_config)
        modify_config(base_config, new_configs, outfile=self.adjoint_config)
        self.adjoint_driver = None

        intercomm = MPI.COMM_SELF.Spawn(sys.executable,
                                        args=str(Path(__file__).parent / 'su2_function_mpi.py'),
                                        maxprocs=num_procs-1)  # XXX -1 because this process is also a worker
        self.comm = intercomm.Merge(high=False)
        assert self.comm.Get_rank() == 0, 'Rank is expected to be 0.'
        self.comm.bcast([self.num_zones, self.dims, self.num_diff_outputs], root=0)

    def forward(self, *inputs, **kwargs):
        if len(inputs) != self.num_diff_inputs:
            raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                            .format(len(inputs), self.config_file, self.num_diff_inputs))

        self.save_for_backward(*inputs)
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()

        self.comm.bcast(self.forward_config, root=0)
        self.comm.bcast(RunCode.RUN_FORWARD, root=0)
        self.comm.bcast(inputs, root=0)  # TODO Any way to optimize? Potentially big

        self.forward_driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                                      self.dims, self.comm)
        run_forward(self.comm, self.forward_driver, inputs)
        shutil.move("./restart_flow.dat", "./solution_flow.dat")

        outputs = tuple(torch.tensor(self.forward_driver.GetDiff_Outputs_Vars(i))
                        for i in range(self.num_diff_outputs))

        # TODO Way to get results in-memory, without writing to file?

        return outputs

    def backward(self, *grad_outputs, **kwargs):
        if self.adjoint_driver is not None:
            self.adjoint_driver.Postprocessing()

        self.comm.bcast(self.adjoint_config, root=0)
        self.comm.bcast(RunCode.RUN_ADJOINT, root=0)
        self.comm.bcast(grad_outputs, root=0)  # TODO Any way to optimize? Potentially big

        self.adjoint_driver = pysu2.CDiscAdjSinglezoneDriver(self.adjoint_config, self.num_zones,
                                                             self.dims, self.comm)
        run_adjoint(self.comm, self.adjoint_driver, self.saved_tensors, grad_outputs)

        grads = tuple(torch.tensor(self.adjoint_driver.GetTotal_Sens_Diff_Inputs(i))
                      for i in range(self.adjoint_driver.GetnDiff_Inputs()))
        return grads


class SU2Function(torch.autograd.Function):
    """Function that uses the SU2 in-memory python wrapper."""

    def __init__(self, config_file, dims=2, num_zones=1):
        assert num_zones == 1, 'Only supports 1 zone for now.'
        self.num_zones = num_zones
        self.dims = dims

        self.config_file = config_file
        base_config = SU2.io.Config(config_file)
        self.num_diff_inputs = len(base_config['DIFF_INPUTS'].split(','))
        self.num_diff_outputs = len(base_config['DIFF_OUTPUTS'].split(','))

        rand_vec = torch.ones(2)
        rand_vec[1] = 2 * torch.rand(1) - 1
        rand_vec = rand_vec / rand_vec.norm()
        self.rand_vec = rand_vec

        new_configs = {'MATH_PROBLEM': 'DIRECT'}
        self.forward_config = 'forward_' + TEMP_CFG_BASENAME
        shutil.copy(config_file, self.forward_config)
        modify_config(base_config, new_configs, outfile=self.forward_config)
        self.forward_driver = None

        new_configs = {'MATH_PROBLEM': 'DISCRETE_ADJOINT'}
        self.adjoint_config = 'adjoint_' + TEMP_CFG_BASENAME
        shutil.copy(config_file, self.adjoint_config)
        modify_config(base_config, new_configs, outfile=self.adjoint_config)
        self.adjoint_driver = None

    def forward(self, *inputs, **kwargs):
        if len(inputs) != self.num_diff_inputs:
            raise TypeError('{} inputs were provided, but the config file ({}) defines {} diff inputs.'
                            .format(len(inputs), self.config_file, self.num_diff_inputs))

        self.save_for_backward(*inputs)
        # TODO How to reset driver from previous run? Or is it necessary to reinstatiate?
        if self.forward_driver is not None:
            self.forward_driver.Postprocessing()
        self.forward_driver = pysu2.CSinglezoneDriver(self.forward_config, self.num_zones,
                                                      self.dims, MPI.COMM_WORLD)
        # TODO Add checks to make sure these are run before Preprocess?
        for i, x in enumerate(inputs):
            self.forward_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
        self.forward_driver.ApplyDiff_Inputs_Vars()

        self.forward_driver.Preprocess(0)
        self.forward_driver.Run()
        outputs = tuple(torch.tensor(self.forward_driver.GetDiff_Outputs_Vars(i))
                        for i in range(self.num_diff_outputs))

        # TODO Are all these necessary?
        self.forward_driver.Update()
        self.forward_driver.Monitor(0)
        self.forward_driver.Output(0)
        # XXX Postprocessing cleans out objects, don't do this for now so theyre accessible after run
        # self.forward_driver.Postprocessing()
        # TODO Way to get results in-memory, without writing to file?
        shutil.move("./restart_flow.dat", "./solution_flow.dat")

        return outputs

    def backward(self, *grad_outputs, **kwargs):
        # TODO How to reset driver from previous run? Or is it necessary to reinstatiate?
        if self.adjoint_driver is not None:
            self.adjoint_driver.Postprocessing()
        self.adjoint_driver = pysu2.CDiscAdjSinglezoneDriver(self.adjoint_config, self.num_zones,
                                                             self.dims, MPI.COMM_WORLD)
        # TODO Add checks to make sure these are run before Preprocess
        inputs = self.saved_tensors
        for i, x in enumerate(inputs):
            self.adjoint_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
        self.adjoint_driver.ApplyDiff_Inputs_Vars()
        for i, g in enumerate(grad_outputs):
            self.adjoint_driver.SetBackprop_Derivs(g.flatten().tolist(), i)

        self.adjoint_driver.Preprocess(0)
        self.adjoint_driver.Run()
        grads = tuple(torch.tensor(self.adjoint_driver.GetTotal_Sens_Diff_Inputs(i))
                      for i in range(self.adjoint_driver.GetnDiff_Inputs()))
        return grads


class SU2IOFunction(torch.autograd.Function):
    """Function that uses the SU2 IO wrapper."""
    def __init__(self, config_file, partitions=0, step=1e-4, nzones=1, grad_method='CONTINUOUS_ADJOINT'):
        self.step = step

        if config_file is SU2.io.Config:
            config = config_file
        else:
            config = SU2.io.Config(config_file)
        config.NUMBER_PART = partitions
        config.NZONES = int(nzones)
        config['GRADIENT_METHOD'] = grad_method
        # Force CSV output in order to compute gradients
        config.WRT_CSV_SOL = 'YES'

        assert config.OPT_COMBINE_OBJECTIVE == 'NO', 'Only works without combined objective.'

        state = SU2.io.State()
        state.FILES.MESH = config.MESH_FILENAME

        project = SU2.opt.Project(config, state, warn=False)

        # TODO Is this the appropriate way to get the list of values?
        self.output_vars = config.OBJECTIVE_FUNCTION.replace(' ', '').split(',')
        self.scales = torch.tensor([v['SCALE'] for v in config['OPT_OBJECTIVE'].values()])

        self.project = project
        # self.config = config
        # self.state = state

    def forward(self, dvs, **kwargs):
        self.save_for_backward(dvs)
        dvs = dvs.tolist()
        self.project.obj_f(dvs)
        state = self.project.latest_design.state
        funcs = [v for k, v in state['FUNCTIONS'].items() if k in self.output_vars]
        return torch.tensor(funcs)

        # XXX Other version, doesnt apply dvs
        # info = SU2.run.direct(self.config)
        # self.state.update(info)
        # SU2.io.restart2solution(self.config, self.state)
        #
        # # XXX return only values for which we will compute adjoint
        # # (i.e., specified in obj function in config)
        # values = [v for k, v in self.state['FUNCTIONS'].items() if k in self.output_vars]
        # return torch.tensor(values)
        # # XXX return all values
        # # return torch.tensor(self.state['FUNCTIONS'].values())

    def backward(self, *grad_outputs):
        dvs = self.saved_tensors[0].tolist()
        grads = self.project.obj_df(dvs)
        grads = torch.tensor(grads).t() / self.scales  # XXX De-scale here or keep scale from config applied?
        # Adjust grad sign because some are flipped in SU2
        for i in range(len(self.output_vars)):
            grads[i] *= self.get_grad_sign(self.output_vars[i])

        return grads @ grad_outputs[0]

        # XXX Other version, returns combo grad
        # info = SU2.run.adjoint(self.config)
        # self.state.update(info)
        # SU2.io.restart2solution(self.config, self.state)
        #
        # # Gradient Projection
        # info = SU2.run.projection(self.config, self.step)
        # self.state.update(info)
        # #
        # # # XXX Concatenate values for all output variables from the forward pass
        # values = [v for k, v in self.state['GRADIENTS'].items() if k in self.output_vars]
        # return torch.tensor(values)

    @staticmethod
    def get_grad_sign(var):
        sign = 1
        if (
            var == "LIFT" or
            var == "EFFICIENCY" or
            var == "THRUST" or
            var == "FIGURE_OF_MERIT" or
            var == "SURFACE_TOTAL_PRESSURE" or
            var == "SURFACE_STATIC_PRESSURE" or
            var == "SURFACE_MASSFLOW" or
            var == "SURFACE_MACH" or
            var == "TOTAL_STATIC_EFFICIENCY"
        ):
            sign = -1
        return sign


def modify_config(config, new_params, outfile=None):
    temp_config = config.copy()
    for k, v in new_params.items():
        temp_config[k] = v
    if outfile is not None:
        temp_config.write(outfile)
    return temp_config


if __name__ == '__main__':
    # Command Line Options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename",
                      help="read config from FILE", metavar="FILE")
    parser.add_option("-n", "--partitions", dest="partitions", default=1,
                      help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_option("-s", "--step", dest="step", default=1E-4,
                      help="DOT finite difference STEP", metavar="STEP")
    parser.add_option("-v", "--validate", dest="validate", default="False",
                      help="Validate the gradient using direct diff. mode", metavar="VALIDATION")
    parser.add_option("-z", "--zones", dest="nzones", default="1",
                      help="Number of Zones", metavar="ZONES")

    (options, args) = parser.parse_args()
    options.partitions = int(options.partitions)
    options.step = float(options.step)
    options.validate = options.validate.upper() == 'TRUE'
    options.nzones = int(options.nzones)

    f = SU2Function(options.filename,
                    options.partitions,
                    options.step,
                    options.nzones)

    n_dvs = sum(f.project.config['DEFINITION_DV']['SIZE'])
    direct = f(torch.zeros(n_dvs, requires_grad=True))
    dummy_loss = direct.sum()
    dummy_loss.backward()
    print('Done')
