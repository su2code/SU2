import os
import sys
import numpy as np
from optparse import OptionParser
import shutil
import pysu2ad

try:
    from mpi4py import MPI
    COMM = MPI.COMM_WORLD
    RANK = COMM.Get_rank()
except ImportError:
    COMM = 0
    RANK = 0

NDIM  = 2
NZONE = 1


def matrix_to_numpy(mat):
    """Convert a CPyWrapperMatrixView to a numpy array."""
    num_rows, num_cols = mat.Shape()
    arr = np.zeros((num_rows, num_cols))
    for i in range(num_rows):
        row = mat.Get(i)
        for j in range(num_cols):
            arr[i, j] = row[j]

    return arr


def parse_surface_file(filename):
    """Parse surface file for sensitivities."""
    coords    = []
    sens_x    = []
    sens_y    = []
    sens_surf = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(('TITLE', 'VARIABLES', 'ZONE')):
                continue
            vals = line.split('\t')

            if len(vals) > 5:
                coords.append((float(vals[0]), float(vals[1])))
                sens_x.append(float(vals[-3]))
                sens_y.append(float(vals[-2]))
                sens_surf.append(float(vals[-1]))

    return np.array(coords), np.array(sens_x), np.array(sens_y), np.array(sens_surf)


def run_primal(config_file):
    """Run primal forward flow solution."""
    driver = pysu2ad.CSinglezoneDriver(config_file, NZONE, COMM)

    sys.stdout.flush()
    if COMM != 0:
        COMM.Barrier()

    driver.Preprocess(0)
    driver.Run()
    driver.Postprocess()
    driver.Update()
    driver.Monitor(0)
    driver.Output(0)
    driver.Finalize()


def run_residuals_adjoint(config_file):
    """Run residuals discrete adjoint solver."""
    adj_driver = pysu2ad.CDiscAdjSinglezoneDriver(config_file, NZONE, COMM)

    sys.stdout.flush()
    if COMM != 0:
        COMM.Barrier()

    adj_driver.Preprocess(0)
    adj_driver.Run()
    adj_driver.Postprocess()
    adj_driver.Update()
    adj_driver.Monitor(0)
    adj_driver.Output(0)

    # farfield variable sensitivities
    obj_farfield = adj_driver.GetObjectiveFarfieldVariablesSensitivities()

    if RANK == 0:
        print("Farfield variable sensitivities (dObjective/dVariables):")
        print(f"    dI/dMach = {obj_farfield[0]:.6e}")
        print(f"    dI/dAoA  = {obj_farfield[1]:.6e}")

    # dI/dq (nPoint, nVar)
    obj_states = adj_driver.ObjectiveSolutionSensitivities()
    obj_states_np = matrix_to_numpy(obj_states)
    if RANK == 0:
        print(f"\ndObjective/dStates: shape = {obj_states_np.shape}")
        print(obj_states_np)

    # dA/dq^T * psi (nPoint, nVar)
    res_states = adj_driver.ResidualsSolutionSensitivities()
    res_states_np = matrix_to_numpy(res_states)
    if RANK == 0:
        print(f"\ndResiduals/dStates^T * psi: shape = {res_states_np.shape}")
        print(res_states_np)

    # dI/dxv (nPoint, nDim)
    obj_coords = adj_driver.ObjectiveCoordinatesSensitivities()
    obj_coords_np = matrix_to_numpy(obj_coords)
    if RANK == 0:
        print(f"\ndObjective/dCoordinates: shape = {obj_coords_np.shape}")
        print(obj_coords_np)

    # dA/dxv^T * psi (nPoint, nDim)
    res_coords = adj_driver.ResidualsCoordinatesSensitivities()
    res_coords_np = matrix_to_numpy(res_coords)
    if RANK == 0:
        print(f"\ndResiduals/dCoordinates^T * psi: shape = {res_coords_np.shape}")
        print(res_coords_np)

    adj_driver.Finalize()


def run_fixed_point_adjoint(config_file):
    """Run fixed-point discrete adjoint solver."""
    adj_driver = pysu2ad.CDiscAdjSinglezoneDriver(config_file, NZONE, COMM)

    sys.stdout.flush()
    if COMM != 0:
        COMM.Barrier()

    adj_driver.Preprocess(0)
    adj_driver.Run()
    adj_driver.Postprocess()
    adj_driver.Update()
    adj_driver.Monitor(0)
    adj_driver.Output(0)
    adj_driver.Finalize()


def compare_surface_sensitivities(res_file, fp_file):
    """Compare surface sensitivity outputs between residual and fixed-point adjoint solvers."""
    if RANK != 0:
        return

    print('\n---------------------------------------------------------------------------')
    print('SURFACE SENSITIVITY COMPARISON: RESIDUAL-BASED VS FIXED-POINT')
    print(f'    Residual adjoint file   : {res_file}')
    print(f'    Fixed-point adjoint file: {fp_file}')
    print('---------------------------------------------------------------------------')

    coords_r, sx_r, sy_r, ss_r = parse_surface_file(res_file)
    coords_f, sx_f, sy_f, ss_f = parse_surface_file(fp_file)

    n_points = len(coords_r)

    if n_points != len(coords_f):
        raise RuntimeError('Mismatch in number of surface points between adjoint solver outputs.')

    header = f"{'(x, y)':>13s} {'Sens (res)':>14s} {'Sens (fp)':>14s} {'ratio':>10s}"
    print(header)
    print('-' * len(header))

    max_rel = 0.0

    for i in range(n_points):
        x, y = coords_r[i]
        r = ss_r[i] / ss_f[i] if abs(ss_f[i]) > 1e-20 else float('inf')

        if abs(ss_f[i]) > 1e-20:
            max_rel = max(max_rel, abs(r - 1.0))
        print(f'({x:5.1f},{y:5.1f}) {ss_r[i]:>14.6e} {ss_f[i]:>14.6e} {r:>10.4f}')

    print(f'\nMax relative error for surface sensitivity: {max_rel:.4e}')

    # print results for regression test
    d_sx = np.max(np.abs(sx_r - sx_f))
    d_sy = np.max(np.abs(sy_r - sy_f))
    d_ss = np.max(np.abs(ss_r - ss_f))
    l2   = np.linalg.norm(np.concatenate([sx_r - sx_f, sy_r - sy_f, ss_r - ss_f]))

    print(f'\nPrintout for regression test...')
    print('\n------------------------------ Begin Solver -----------------------------')
    print(f'{77777:>8d} {d_sx:.6f} {d_sy:.6f} {d_ss:.6f} {l2:.6f}')

def main():
    # command line options
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="filename", help="Read config from FILE", metavar="FILE")
    parser.add_option("--parallel", action="store_true", help="Specify if we need to initialize MPI", dest="with_MPI",
                      default=False)

    (options, args) = parser.parse_args()

    if options.filename is None:
        parser.error("Config file must be specified with -f or --file.")

    # run primal solver
    if RANK == 0:
        print("\n-------------------- FORWARD FLOW SOLVER --------------------")
    run_primal(options.filename)

    shutil.copy('restart.dat', 'solution.dat')

    # run residual adjoint solver
    res_config = 'adjoint_res.cfg'
    with open(options.filename, 'r') as f:
        cfg_text = f.read()
    cfg_text = cfg_text.replace('MATH_PROBLEM = DIRECT', 'MATH_PROBLEM = DISCRETE_ADJOINT')
    with open(res_config, 'w') as f:
        f.write(cfg_text)

    if RANK == 0:
        print("\n-------------------- RESIDUAL ADJOINT SOLVER --------------------")
    run_residuals_adjoint(res_config)

    res_surface = 'surface_adjoint_residual.dat'
    if RANK == 0:
        os.rename('surface_adjoint.dat', res_surface)

    # fixed-point adjoint solver
    fp_config = 'adjoint_fp.cfg'
    with open(res_config, 'r') as f:
        cfg_text = f.read()
    cfg_text = cfg_text.replace('KIND_DISC_ADJ = RESIDUALS', 'KIND_DISC_ADJ = FIXED_POINT')
    with open(fp_config, 'w') as f:
        f.write(cfg_text)

    if options.with_MPI:
        COMM.Barrier()

    print("\n-------------------- FIXED-POINT ADJOINT SOLVER --------------------")
    run_fixed_point_adjoint(fp_config)

    fp_surface = 'surface_adjoint_fixed_point.dat'
    if RANK == 0:
        os.rename('surface_adjoint.dat', fp_surface)

    # compare surface sensitivities
    compare_surface_sensitivities(res_surface, fp_surface)

    # clean up temp configs
    if RANK == 0 and os.path.exists(res_config):
        os.remove(res_config)
    if RANK == 0 and os.path.exists(fp_config):
        os.remove(fp_config)


if __name__ == '__main__':
    main()
