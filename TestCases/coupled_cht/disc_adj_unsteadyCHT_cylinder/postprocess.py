# Compute and print absolute difference between Discrete Adjoint
# and Finite Difference gradient. Prints also percentage difference.

import pandas as pd

# load files
DAgrad = pd.read_csv("of_grad.csv")
FDvals = pd.read_csv("doe.his")

# additional values
FDstep = 1e-4
FDstring = '  tavgT'
DAstring = 'AVG_TEMPERATURE gradient '

# create FD gradient
FDgrad = (FDvals[FDstring].iloc[1:] - FDvals[FDstring].iloc[0]) / FDstep
# The above operation creates a pd.series that starts at 1 because the first FDval is the baseline.
# This is reset to match the DAgrad
FDgrad = FDgrad.reset_index(inplace=False)

# absolute difference
absoluteDiff = DAgrad[DAstring] - FDgrad[FDstring]
print("DAgrad - FDgrad\n", absoluteDiff)

# relative difference in percent
relDiffPercent = abs(DAgrad[DAstring] - FDgrad[FDstring])/DAgrad[DAstring] * 100
print("(DAgrad - FDgrad) / DAgrad * 100\n", relDiffPercent)

