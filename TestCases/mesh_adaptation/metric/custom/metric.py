from mpi4py import MPI

import pysu2
from SU2.metric import CustomSensorRegistry

comm = MPI.COMM_WORLD


def total_pressure(driver) -> list[float]:
    """Compute local total pressure at all nodes for compressible flow."""
    prim_idx = driver.GetPrimitiveIndices()
    nDim = driver.GetNumberDimensions()
    primVars = driver.Primitives()
    nNodes = driver.GetNumberNodes() - driver.GetNumberHaloNodes()

    gamma = 1.4
    coeff = 0.5 * (gamma - 1.0)
    exp_fac = gamma / (gamma - 1.0)
    press_col = prim_idx["PRESSURE"]
    vel_cols = [prim_idx["VELOCITY_X"], prim_idx["VELOCITY_Y"]]
    a_col = prim_idx["SOUND_SPEED"]
    if nDim == 3:
        vel_cols.append(prim_idx["VELOCITY_Z"])

    p_tot = [0.0] * nNodes
    for iNode in range(nNodes):
        row = primVars(iNode)
        p = row[press_col]
        a = row[a_col]
        vel2 = sum(row[k] ** 2 for k in vel_cols)
        mach2 = vel2 / (a * a)

        p_tot[iNode] = p * pow(1.0 + coeff * mach2, exp_fac)
    return p_tot


def main():
    # Intialize driver
    driver = pysu2.CSinglezoneDriver("rans_naca0012.cfg", 1, comm)

    # Initialize custom metric sensors
    custom_sensors = CustomSensorRegistry({"TOTAL_PRESSURE": total_pressure})

    # Initialize map of Sensor: idx
    custom_sensors.initialize(driver)

    driver.Preprocess(0)
    driver.Run()

    # Store custom sensor
    custom_sensors.populate(driver)

    # Postprocess solution/metric and output
    driver.Postprocess()
    driver.Update()
    driver.Monitor(0)
    driver.Output(0)

    driver.Finalize()


if __name__ == "__main__":
    main()
