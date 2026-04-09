import math

from mpi4py import MPI

import pysu2
from SU2.metric import CustomSensorRegistry

comm = MPI.COMM_WORLD


def mach_number(driver, iNode: int) -> float:
    """Compute local Mach number at a single node for compressible flow."""
    prim_idx = driver.GetPrimitiveIndices()
    nDim = driver.GetNumberDimensions()
    primVars = driver.Primitives()

    vel_cols = [prim_idx["VELOCITY_X"], prim_idx["VELOCITY_Y"]]
    if nDim == 3:
        vel_cols.append(prim_idx["VELOCITY_Z"])
    a_col = prim_idx["SOUND_SPEED"]

    vel2 = sum(primVars.Get(iNode, k) ** 2 for k in vel_cols)
    a = primVars.Get(iNode, a_col)
    return math.sqrt(vel2) / max(a, 1e-20)


def main():
    # Intialize driver
    driver = pysu2.CSinglezoneDriver("rans_naca0012.cfg", 1, comm)

    # Initialize custom metric sensors
    custom_sensors = CustomSensorRegistry({"MACH": mach_number})

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
