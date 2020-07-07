import os
import time
import shutil
from enum import IntEnum

import torch
import numpy as np
from mpi4py import MPI

import SU2
import pysu2
import pysu2ad

from typing import Sequence, Union, Tuple, Dict, TypeVar
GenTensor = TypeVar('GenTensor', torch.Tensor, np.ndarray)


_non_busy_wait_max_time = 0.1


class RunCode(IntEnum):
    """Run codes for communication with worker processes."""
    STOP = -1
    RUN_FORWARD = 0
    RUN_ADJOINT = 1


def run_forward(comm: MPI.Intracomm, forward_driver: pysu2.CSinglezoneDriver,
                inputs: Sequence[GenTensor]) -> Tuple[GenTensor, ...]:
    """Runs a simulation with the provided driver, using the inputs to set the values
    defined in DIFF_INPUTS in the config file.

    :param comm: The communicator for the processes running the simulation.
    :param forward_driver: The driver for the simulation, created using the same comm as passed into this function.
    :param inputs: The inputs used to set the DIFF_INPUTS as defined in the configuration file.
    :return: The outputs of the simulation, as defined in DIFF_OUTPUTS in the config file.
    """
    rank = comm.Get_rank()
    for i, x in enumerate(inputs):
        forward_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
    forward_driver.ApplyDiff_Inputs_Vars()

    forward_driver.StartSolver()
    comm.Barrier()

    # are we using numpy or torch
    is_numpy = len(inputs) == 0 or type(inputs[0]) is np.ndarray
    if is_numpy:
        array_func = np.array
        cat_func = np.concatenate
    else:
        import torch
        array_func = inputs[0].new_tensor
        cat_func = torch.cat

    num_diff_outputs = forward_driver.GetnDiff_Outputs()
    outputs = [array_func(forward_driver.GetDiff_Outputs_Vars(i))
               for i in range(num_diff_outputs)]

    for i in range(num_diff_outputs):
        if outputs[i].shape[0] > 1:
            # if dealing with full-grid, reorder according to GlobalIndex
            if comm.Get_size() > 1:
                # gather outputs in rank 0 if more than one rank
                outputs[i] = comm.gather(outputs[i], root=0)
                global_inds = comm.gather(forward_driver.GetAllGlobalIndices(), root=0)
                if rank == 0:
                    outputs[i] = cat_func(outputs[i])
                    global_inds = list(sum(global_inds, tuple()))  # join tuples
            else:
                global_inds = list(forward_driver.GetAllGlobalIndices())

            if rank == 0:
                # TODO Make the list integers on the C side
                global_inds = np.array(global_inds, dtype=np.long)
                assert outputs[i].shape[0] == len(global_inds), \
                    'Only full grid outputs supported by now (besides scalars).'
                # order by global_inds
                outputs[i][global_inds] = outputs[i].copy() if is_numpy else outputs[i].clone()
            else:
                outputs[i] = None
    return tuple(outputs)


def run_adjoint(comm: MPI.Intracomm, adjoint_driver: pysu2ad.CDiscAdjSinglezoneDriver,
                inputs: Sequence[GenTensor], grad_outputs: Sequence[GenTensor]) -> Tuple[GenTensor, ...]:
    """Runs a simulation with the provided driver, using the inputs to set the values
    defined in DIFF_INPUTS in the config file.

    :param comm: The communicator for the processes running the simulation.
    :param adjoint_driver: The driver for the adjoint computation, created using the same comm as passed into this function.
    :param inputs: The same inputs used to set the DIFF_INPUTS in the forward pass.
    :param grad_outputs: Gradients of a scalar loss with respect to the forward outputs, see SU2Function's backward() method.
    :return: The gradients of the loss with respect to the inputs.
    """
    rank = comm.Get_rank()
    for i, x in enumerate(inputs):
        adjoint_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
    adjoint_driver.ApplyDiff_Inputs_Vars()
    for i, g in enumerate(grad_outputs):
        adjoint_driver.SetBackprop_Derivs(g.flatten().tolist(), i)

    adjoint_driver.StartSolver()

    # are we using numpy or torch
    is_numpy = len(inputs) == 0 or type(inputs[0]) is np.ndarray
    if is_numpy:
        array_func = np.array
        cat_func = np.concatenate
    else:
        import torch
        array_func = inputs[0].new_tensor
        cat_func = torch.cat

    num_diff_inputs = adjoint_driver.GetnDiff_Inputs()
    grads = [array_func(adjoint_driver.GetTotal_Sens_Diff_Inputs(i))
             for i in range(num_diff_inputs)]
    for i in range(num_diff_inputs):
        if grads[i].shape[0] > 1:
            # if dealing with full-grid, reorder according to GlobalIndex
            if comm.Get_size() > 1:
                # gather outputs in rank 0 if more than one rank
                grads[i] = comm.gather(grads[i], root=0)
                global_inds = comm.gather(adjoint_driver.GetAllGlobalIndices(), root=0)
                if rank == 0:
                    grads[i] = cat_func(grads[i])
                    global_inds = list(sum(global_inds, tuple()))  # join tuples
            else:
                global_inds = list(adjoint_driver.GetAllGlobalIndices())

            if rank == 0:
                global_inds = np.array(global_inds, dtype=np.long)
                assert grads[i].shape[0] == len(global_inds), \
                    'Only full grid outputs supported by now (besides scalars).'
                # order by global_inds
                grads[i][global_inds] = grads[i].copy() if is_numpy else grads[i].clone()
            else:
                grads[i] = None
    return tuple(grads)


def modify_config(config: SU2.io.Config, new_params: Dict[str, str],
                  outfile: Union[str, os.PathLike, None] = None) -> SU2.io.Config:
    """Modify a config, saving the modifications to outfile if provided."""
    temp_config = config.copy()
    for k, v in new_params.items():
        temp_config[k] = v
    if outfile is not None:
        temp_config.write(outfile)
    return temp_config


def activate_su2_mpi(remove_temp_files: bool = True, max_procs_per_example: int = 1,
                     non_busy_wait_max_time: float = 0.1) -> None:
    assert MPI.COMM_WORLD.Get_size() > 1, 'Need at least 1 master and 1 worker process, run with "mpirun -np ...'

    if MPI.COMM_WORLD.Get_rank() != 0:
        global _non_busy_wait_max_time
        _non_busy_wait_max_time = non_busy_wait_max_time
        main(remove_temp_files=remove_temp_files)
        exit(0)

    # Only rank 0 from here on
    def stop():
        non_busy_post(MPI.COMM_WORLD)
        MPI.COMM_WORLD.bcast(RunCode.STOP, root=0)
    import atexit
    atexit.register(stop)

    import su2torch.su2_function
    su2torch.su2_function._global_max_ppe = max_procs_per_example


def non_busy_wait(comm: MPI.Intracomm) -> None:
    b = comm.Ibarrier()
    start = time.time()
    while not b.Get_status():
        time.sleep(min((time.time() - start) / 2, _non_busy_wait_max_time))


def non_busy_post(comm: MPI.Intracomm) -> None:
    comm.Ibarrier()


def main(remove_temp_files: bool = True) -> None:
    """Runs a loop for the worker processes.
    Can be signaled to run either a forward simulation or an adjoint computation
    using RunCodes.
    """
    local_comm = MPI.COMM_WORLD.Create_group(MPI.Group.Excl(MPI.COMM_WORLD.Get_group(), [0]))
    local_rank = local_comm.Get_rank()
    local_size = local_comm.Get_size()
    ppid = str(os.getppid())

    x = inputs = batch_comm = batch_index = batch_rank = forward_config = None
    num_zones = dims = batch_solution_filename = batch_restart_filename = None
    batch_size = procs_per_example = 1
    while True:
        non_busy_wait(MPI.COMM_WORLD)
        run_type = MPI.COMM_WORLD.bcast(None, root=0)
        if run_type == RunCode.STOP:
            # remove temporary files
            if local_rank == 0 and remove_temp_files:
                import sys
                os.system(f'rm b*_{ppid}_* 2> /dev/null')
            break

        if run_type == RunCode.RUN_FORWARD:
            if procs_per_example != 1 and procs_per_example != local_size:
                # disconnect batch_comm from previous run, if it was created
                batch_comm.Disconnect()
            num_zones, dims, forward_config, mesh_file, procs_per_example, inputs = MPI.COMM_WORLD.bcast(None, root=0)
            batch_size = inputs[0].shape[0]
            batch_index = local_rank // procs_per_example
            if procs_per_example == 1:
                batch_comm = MPI.COMM_SELF
            elif procs_per_example == local_size:
                batch_comm = local_comm
            else:
                batch_comm = local_comm.Split(batch_index, local_rank)
            if local_rank >= batch_size * procs_per_example:
                # these procs wont be used
                non_busy_post(MPI.COMM_WORLD)
                continue
            batch_rank = batch_comm.Get_rank()
            x = [z[batch_index] for z in inputs]

            batch_forward_config = f'b{batch_index}_{ppid}_{forward_config}'
            if batch_rank == 0:
                old_config = SU2.io.Config(forward_config)
                restart_filename = old_config['RESTART_FLOW_FILENAME']
                batch_restart_filename = f'b{batch_index}_{ppid}_{restart_filename}'
                mesh_file = mesh_file.format(batch_index=batch_index) if mesh_file else old_config['MESH_FILENAME']
                new_config = {'RESTART_FLOW_FILENAME': batch_restart_filename,
                              'MESH_FILENAME': mesh_file}
                shutil.copy(forward_config, batch_forward_config)
                modify_config(old_config, new_config, outfile=batch_forward_config)
            if local_rank == 0:
                MPI.COMM_WORLD.send(batch_forward_config, dest=0)
            batch_comm.Barrier()

            forward_driver = pysu2.CSinglezoneDriver(batch_forward_config, num_zones, dims, batch_comm)
            # TODO SetRestart_FlowFileName is not necessary anymore, remove from C++
            # forward_driver.SetRestart_FlowFileName(batch_restart_filename)
            outputs = run_forward(batch_comm, forward_driver, x)
            output_lengths = [o.shape[0] for o in outputs]
            non_busy_post(MPI.COMM_WORLD)
            if batch_rank == 0:
                MPI.COMM_WORLD.send(outputs, dest=0)
                # TODO Way to get results in-memory, without writing to file?
                batch_solution_filename = batch_restart_filename.replace('restart', 'solution')
                shutil.move(batch_restart_filename, batch_solution_filename)
            forward_driver.Postprocessing()

        elif run_type == RunCode.RUN_ADJOINT:
            assert inputs is not None, 'Run forward simulation before running the adjoint.'
            inputs = None
            grad_outputs = MPI.COMM_WORLD.bcast(None, root=0)
            if local_rank >= batch_size * procs_per_example:
                # these procs wont be used
                non_busy_post(MPI.COMM_WORLD)
                continue
            dl = [z[batch_index, :output_lengths[i]] for i, z in enumerate(grad_outputs)]

            batch_adjoint_config = 'b{}_{}_adjoint_{}'.format(batch_index, str(os.getppid()), forward_config)
            if batch_rank == 0:
                old_config = SU2.io.Config(forward_config)
                mesh_file = mesh_file.format(batch_index=batch_index) if mesh_file else old_config['MESH_FILENAME']
                new_config = {'MATH_PROBLEM': 'DISCRETE_ADJOINT',
                              'SOLUTION_FLOW_FILENAME': batch_solution_filename,
                              'RESTART_ADJ_FILENAME': batch_restart_filename.replace('flow', 'adj'),
                              'MESH_FILENAME': mesh_file}
                shutil.copy(forward_config, batch_adjoint_config)
                modify_config(old_config, new_config, outfile=batch_adjoint_config)
            batch_comm.Barrier()
            adjoint_driver = pysu2ad.CDiscAdjSinglezoneDriver(batch_adjoint_config, num_zones, dims, batch_comm)
            grads = run_adjoint(batch_comm, adjoint_driver, x, dl)
            non_busy_post(MPI.COMM_WORLD)
            if batch_rank == 0:
                MPI.COMM_WORLD.send(grads, dest=0)
            adjoint_driver.Postprocessing()
