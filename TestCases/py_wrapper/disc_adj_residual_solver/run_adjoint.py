import sys
import pysu2ad
from optparse import OptionParser
import numpy as np

def matrix_to_numpy(mat):
    """Convert a CPyWrapperMatrixView to a numpy array."""
    num_rows, num_cols = mat.Shape()
    arr = np.zeros((num_rows, num_cols))
    for i in range(num_rows):
        row = mat.Get(i)
        for j in range(num_cols):
            arr[i, j] = row[j]

    return arr


def main():
    # command line options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
    parser.add_option("--parallel", action="store_true", help="Specify if we need to initialize MPI", dest="with_MPI",
                      default=False)

    (options, args) = parser.parse_args()
    options.nDim  = 2
    options.nZone = 1

    # initialize MPI
    if options.with_MPI:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.Get_rank()
    else:
        comm = 0
        rank = 0

    # initialize and preprocess solver
    adj_driver = pysu2ad.CDiscAdjSinglezoneDriver(options.filename, options.nZone, comm)

    sys.stdout.flush()
    if options.with_MPI:
        comm.Barrier()

    adj_driver.Preprocess(0)
    adj_driver.Run()
    adj_driver.Postprocess()
    adj_driver.Update()
    adj_driver.Monitor(0)
    adj_driver.Output(0)

    # extract sensitivities
    if rank == 0:
        print("\n------------------------------ Sensitivities -----------------------------\n")

    # farfield variable sensitivities
    obj_farfield = adj_driver.GetObjectiveFarfieldVariablesSensitivities()

    if rank == 0:
        print("Farfield variable sensitivities (dObjective/dVariables):")
        print(f"    dI/dMach = {obj_farfield[0]:.6e}")
        print(f"    dI/dAoA  = {obj_farfield[1]:.6e}")

    # dI/dq (nPoint, nVar)
    obj_states = adj_driver.ObjectiveStatesSensitivities()
    obj_states_np = matrix_to_numpy(obj_states)
    if rank == 0:
        print(f"\ndObjective/dStates: shape = {obj_states_np.shape}")
        print(obj_states_np)

    # dA/dq^T * psi (nPoint, nVar)
    res_states = adj_driver.ResidualsStatesSensitivities()
    res_states_np = matrix_to_numpy(res_states)
    if rank == 0:
        print(f"\ndResiduals/dStates^T * psi: shape = {res_states_np.shape}")
        print(res_states_np)

    # dI/dxv (nPoint, nDim)
    obj_coords = adj_driver.ObjectiveCoordinatesSensitivities()
    obj_coords_np = matrix_to_numpy(obj_coords)
    if rank == 0:
        print(f"\ndObjective/dCoordinates: shape = {obj_coords_np.shape}")
        print(obj_coords_np)

    # dA/dxv^T * psi (nPoint, nDim)
    res_coords = adj_driver.ResidualsCoordinatesSensitivities()
    res_coords_np = matrix_to_numpy(res_coords)
    if rank == 0:
        print(f"\ndResiduals/dCoordinates^T * psi: shape = {res_coords_np.shape}")
        print(res_coords_np)

    adj_driver.Finalize()


if __name__ == '__main__':
    main()
