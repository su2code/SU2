from enum import IntEnum

import pysu2ad as pysu2
from mpi4py import MPI


class RunCode(IntEnum):
    # Run codes for worker processes
    STOP = -1
    RUN_FORWARD = 0
    RUN_ADJOINT = 1


def run_forward(comm, forward_driver, inputs):
    # TODO Add checks to make sure these are run before Preprocess?
    for i, x in enumerate(inputs):
        forward_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
    forward_driver.ApplyDiff_Inputs_Vars()

    forward_driver.Preprocess(0)
    forward_driver.Run()
    forward_driver.Output(0)
    comm.Barrier()


def run_adjoint(comm, adjoint_driver, inputs, grad_outputs):
    # TODO Add checks to make sure these are run before Preprocess
    for i, x in enumerate(inputs):
        adjoint_driver.SetDiff_Inputs_Vars(x.flatten().tolist(), i)
    adjoint_driver.ApplyDiff_Inputs_Vars()
    for i, g in enumerate(grad_outputs):
        adjoint_driver.SetBackprop_Derivs(g.flatten().tolist(), i)

    adjoint_driver.Preprocess(0)
    adjoint_driver.Run()


def main():
    intercomm = MPI.Comm.Get_parent()
    comm = intercomm.Merge(high=True)

    num_zones, dims, num_diff_outputs = comm.bcast(None, root=0)

    inputs = None
    while True:
        config = comm.bcast(None, root=0)
        run_type = comm.bcast(None, root=0)
        if run_type == RunCode.STOP:
            break

        if run_type == RunCode.RUN_FORWARD:
            inputs = comm.bcast(None, root=0)
            forward_driver = pysu2.CSinglezoneDriver(config, num_zones, dims, comm)
            run_forward(comm, forward_driver, inputs)
            forward_driver.Postprocessing()

        elif run_type == RunCode.RUN_ADJOINT:
            assert inputs is not None, 'Run forward simulation before running the adjoint.'
            grad_outputs = comm.bcast(None, root=0)
            adjoint_driver = pysu2.CDiscAdjSinglezoneDriver(config, num_zones, dims, comm)
            run_adjoint(comm, adjoint_driver, inputs, grad_outputs)
            adjoint_driver.Postprocessing()
            inputs = None

    comm.Disconnect()
    intercomm.Disconnect()


if __name__ == '__main__':
    # will run this when run with mpirun
    main()
