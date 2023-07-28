import sys
from optparse import OptionParser	# use a parser for configuration
import pysu2			            # imports the SU2 wrapped module
from math import *
import os
import pysu2ad
from mpi4py import MPI
import numpy as np
import shutil
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ad_com = MPI.COMM_WORLD
rank_ad = ad_com.Get_rank() 
n_of_steps = 10
rotation_vector = np.zeros(n_of_steps)+20
rotation_vector[5:9] = -rotation_vector[5:9]

SU2Driver = pysu2.CSinglezoneDriver("spinning_cylinder.cfg",1, comm)

for i in range(n_of_steps):
    
    SU2Driver.SetMarkerRotationRate(0,0,0,rotation_vector[n_of_steps-i-1])
    SU2Driver.Preprocess(i)
    SU2Driver.Run()
    SU2Driver.Postprocess()
    SU2Driver.Output(i)
    SU2Driver.Update()
SU2Driver.Finalize()

   