import sys
import numpy as np
from mpi4py import MPI
import SU2, pysu2


def setup_mpi():
    if 'mpi4py.MPI' in sys.modules:
        have_mpi = True
        comm = MPI.COMM_WORLD
        status = MPI.Status()
        myid = comm.Get_rank()
    else:
        have_mpi = False
        comm = None
        status = None
        myid = 0
    return have_mpi, comm, status, myid


class SU2Interface:
    def __init__(self, config_filename):
        self.have_mpi, self.comm, self.status, self.myid = setup_mpi()
        self.FluidSolver = pysu2.CSinglezoneDriver(config_filename, 1, self.comm)

    def run_solver(self):
        self.comm.barrier()
        # run solver
        self.FluidSolver.Preprocess(0)
        self.FluidSolver.Run()
        self.FluidSolver.Postprocess()
        # write outputs
        self.FluidSolver.Monitor(0)
        self.FluidSolver.Output(0)
        self.comm.barrier()

    def save_forces(self):
        solver_all_moving_markers = np.array(self.FluidSolver.GetDeformableMarkerTags())
        solver_markers = self.FluidSolver.GetMarkerTags()
        solver_marker_ids = self.FluidSolver.GetMarkerIndices()
        # The surface marker and the partitioning of the solver usually don't agree.
        # Thus, it is necessary to figure out if the partition of the current mpi process has
        # a node that belongs to a moving surface marker.
        has_moving_marker = [marker in solver_markers for marker in solver_all_moving_markers]

        f = open('forces_'+str(self.myid)+'.csv','w')

        for marker in solver_all_moving_markers[has_moving_marker]:
            solver_marker_id = solver_marker_ids[marker]
            n_vertices = self.FluidSolver.GetNumberMarkerNodes(solver_marker_id)
            for i_vertex in range(n_vertices):
                fxyz = self.FluidSolver.GetMarkerFlowLoad(solver_marker_id, i_vertex)
                iPoint = self.FluidSolver.GetMarkerNode(solver_marker_id, i_vertex)
                GlobalIndex = self.FluidSolver.GetNodeGlobalIndex(iPoint)
                f.write('{}, {:.2f}, {:.2f}, {:.2f}\n'.format(GlobalIndex, fxyz[0], fxyz[1], 0.0))
        f.close()

cfd_interface = SU2Interface('config.cfg')
cfd_interface.run_solver()
cfd_interface.save_forces()

