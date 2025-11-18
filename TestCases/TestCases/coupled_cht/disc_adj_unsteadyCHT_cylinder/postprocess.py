# Compute and print absolute difference between Discrete Adjoint
# and Finite Difference gradient. Prints also percentage difference.
#
# Run this script after `python gradient_validation.py` successfully finished

import pandas as pd

# load files
DAgrad = pd.read_csv("DOE/DOT/of_grad.csv")
FDvals = pd.read_csv("doe.his")

# additional values
nDV = 18
FDstep = 1e-4
FDstring = '  tavgT'
DAstring = 'AVG_TEMPERATURE gradient '

# create FD gradient
FDgrad = (FDvals[FDstring].iloc[:nDV] - FDvals[FDstring].iloc[nDV]) / FDstep

# absolute difference
absoluteDiff = DAgrad[DAstring] - FDgrad
print("absolute diff = DAgrad - FDgrad")

# relative difference in percent
relDiffPercent = (DAgrad[DAstring] - FDgrad)/abs(DAgrad[DAstring]) * 100
print("relative diff = (DAgrad - FDgrad) / abs(DAgrad) * 100")


print('')
print('+-----------+-------------------+-------------------+-------------------+-------------------+')
print('| DV number |       DA gradient |       FD gradient |     absolute diff | relative diff [%] |')
print('+-----------+-------------------+-------------------+-------------------+-------------------+')

for i_dv in range(0, nDV,1):
    print('|{0:10d} |{1:18.10f} |{2:18.10f} |{3:18.10f} |{4:18.10f} |'.format(i_dv, DAgrad[DAstring].iloc[i_dv], FDgrad[i_dv], absoluteDiff[i_dv], relDiffPercent[i_dv]))

print('+-----------+-------------------+-------------------+-------------------+-------------------+')
print('')
