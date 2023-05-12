import sys
import numpy as np
from mpi4py import MPI
import pysu2

"""
This test case updates the translation and rotation rates given in the config file via the python wrapper
using the python functions SetTranslationRate() and SetRotationRate(). 
Expected results:
            |  Inner_Iter|    rms[Rho]|   rms[RhoU]|   rms[RhoV]|   rms[RhoE]|         CFx|         CFy|         CFz|         CMx|         CMy|         CMz|
A0A=0.0 deg: no lift and moments for the symmetrical airfoil (for reference)
             |          51|   -9.205443|   -9.406638|   -9.542742|   -8.868505|    0.008431|   -0.000001|    0.000000|    0.000000|    0.000000|   -0.000037|
AoA=0.5 deg: increase in lift (CFy), small pitching moment (CMz) because origin at the leading edge
             |          65|   -9.016705|   -9.182525|   -9.424299|   -8.656815|    0.009455|    0.133624|    0.000000|    0.000000|    0.000000|    0.045316|
AoA=0.5 deg and rot_z=-30 deg/s: small increase in pitching moment (CMz) and with the origin at the leading edge, the pitch rate adds a little lift (CFy)
             |          50|   -9.031666|   -9.303253|   -9.288828|   -8.652367|    0.009733|    0.153685|    0.000000|    0.000000|    0.000000|    0.053223|
Thus, a check of the final surface forces is an indicator that the update of the moving frame worked properly.
"""

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
        
    def update_moving_frame(self):
        # Reduce AoA from 0.0 to 0.5 deg
        self.FluidSolver.SetTranslationRate(xDot=-265.0469838528219, yDot=-2.3130299864229564, zDot=0.0)
        # Add a pitch rate (in the 2D case, nose up pitch is a negative rotation about the z-axis)
        self.FluidSolver.SetRotationRate(rot_x=0.0, rot_y=0.0, rot_z=-30.0/180.0*np.pi)

    def save_forces(self):
        solver_all_moving_markers = np.array(self.FluidSolver.GetDeformableMarkerTags())
        solver_markers = self.FluidSolver.GetMarkerTags()
        solver_marker_ids = self.FluidSolver.GetMarkerIndices()
        # The surface marker and the partitioning of the solver usually don't agree.
        # Thus, it is necessary to figure out if the partition of the current mpi process has
        # a node that belongs to a moving surface marker.
        has_moving_marker = [marker in solver_markers for marker in solver_all_moving_markers]

        with open('forces_'+str(self.myid)+'.csv','w') as f:
            for marker in solver_all_moving_markers[has_moving_marker]:
                solver_marker_id = solver_marker_ids[marker]
                n_vertices = self.FluidSolver.GetNumberMarkerNodes(solver_marker_id)
                for i_vertex in range(n_vertices):
                    fxyz = self.FluidSolver.GetMarkerFlowLoad(solver_marker_id, i_vertex)
                    iPoint = self.FluidSolver.GetMarkerNode(solver_marker_id, i_vertex)
                    GlobalIndex = self.FluidSolver.GetNodeGlobalIndex(iPoint)
                    f.write('{}, {:.2f}, {:.2f}, {:.2f}\n'.format(GlobalIndex, fxyz[0], fxyz[1], 0.0))

cfd_interface = SU2Interface('config.cfg')
cfd_interface.update_moving_frame()
cfd_interface.run_solver()
cfd_interface.save_forces()

