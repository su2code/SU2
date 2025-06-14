import pysu2ad
import numpy as np

try:
    from mpi4py import MPI
except ImportError:
    comm = None
    rank = 0
else:
    comm = MPI.COMM_WORLD
    rank = MPI.COMM_WORLD.Get_rank()

filename = 'flow_adjoint.cfg'
nzone = 1

adj_driver = pysu2ad.CDiscAdjSinglezoneDriver(filename, nzone, comm)
