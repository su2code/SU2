import sys
import pysu2ad
from optparse import OptionParser
import numpy as np


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

    if rank == 0:
        print("\n------------------------------ Begin Solver -----------------------------\n")
    sys.stdout.flush()
    if options.with_MPI:
        comm.Barrier()

    adj_driver.Preprocess(0)
    adj_driver.Run()
    adj_driver.Postprocess()
    adj_driver.Update()
    adj_driver.Monitor(0)
    adj_driver.Output(0)
    adj_driver.Finalize()


if __name__ == '__main__':
    main()
