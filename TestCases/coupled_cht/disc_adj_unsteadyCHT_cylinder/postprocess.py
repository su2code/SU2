# Compute and print absolute difference between Discrete Adjoint
# and Finite Difference gradient. Prints also percentage difference.
#
# Run this script after `python gradient_validation.py` successfully finished

import pandas as pd

# load files
DAgrad = pd.read_csv("DOE/DOT/of_grad.csv")
FDvals = pd.read_csv("doe.his")

# additional values
FDstep = 1e-4
FDstring = '  tavgT'
DAstring = 'AVG_TEMPERATURE gradient '

# create FD gradient
FDgrad = (FDvals[FDstring].iloc[:18] - FDvals[FDstring].iloc[18]) / FDstep

# absolute difference
absoluteDiff = DAgrad[DAstring] - FDgrad
print("DAgrad - FDgrad\n", absoluteDiff)

# relative difference in percent
relDiffPercent = abs(DAgrad[DAstring] - FDgrad)/DAgrad[DAstring] * 100
print("(DAgrad - FDgrad) / DAgrad * 100\n", relDiffPercent)
