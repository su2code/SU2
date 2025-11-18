import pysu2			            # imports the SU2 wrapped module
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
rotation_vector = np.linspace(0,20,10)
SU2Driver = pysu2.CSinglezoneDriver("spinning_cylinder.cfg",1, comm)

for i, rate in enumerate(rotation_vector):
    SU2Driver.SetMarkerRotationRate(0,0,0,rate)
    SU2Driver.Preprocess(i)
    SU2Driver.Run()
    SU2Driver.Postprocess()
    SU2Driver.Output(i)
    SU2Driver.Update()
SU2Driver.Finalize()

   