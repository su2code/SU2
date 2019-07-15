# --------------------------------------------------------------------------- #
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from matplotlib.colors import LightSource

# --------------------------------------------------------------------------- #
# implort .dat surface file solution
data = pd.read_csv("surface_flow.dat", nrows=4264, skiprows=3, sep='\t', header=None)
x = data[0][:]
y = data[1][:]
vel_z = data[6][:]

# create surface triangulation
points2D = np.vstack([x,y]).T
tri = Delaunay(points2D)

# --------------------------------------------------------------------------- #
analytic_sol = -1/(4*1.8e-5) * (-0.001/5e-4) * (5e-3**2 - ((x**2 + y**2)**(0.5))**2 )
# plot the percentage of deviation '(analytic - sim)/sim*100' for each point
perc_devi_from_anal = abs(analytic_sol - vel_z) / max(analytic_sol) * 100

# get absolute maximum of dataset
maxvel = max(abs(perc_devi_from_anal))
# --------------------------------------------------------------------------- #

fig, ax = plt.subplots(2,2)
# --------------------------------------------------------------------------- #
# 1. analytical solution
ax[0,0].set_title("Analytical solution")
ax[0,0].set_aspect('equal')
tcf1 = ax[0,0].tricontourf(x, y, abs(analytic_sol))
ax[0,0].scatter(x,y, s=0.1, color='black', marker='.')

print(min(analytic_sol))

fig.colorbar(tcf1, ax=ax[0,0])
# --------------------------------------------------------------------------- #
# 2. simulated solution
ax[0,1].set_title("Simulated solution")
ax[0,1].set_aspect('equal')
tcf = ax[0,1].tricontourf(x, y, vel_z)
ax[0,1].scatter(x,y, s=0.1, color='black', marker='.')

fig.colorbar(tcf, ax=ax[0,1])
# --------------------------------------------------------------------------- #
# 3. absolute value deviation between analytic and simulated
ax[1,0].set_title("abs(analytic-simulated)")
ax[1,0].set_aspect('equal')
tcf = ax[1,0].tricontourf(x, y, abs(analytic_sol-vel_z), cmap=plt.cm.Greys)
ax[1,0].scatter(x,y, s=0.1, color='black', marker='.')

fig.colorbar(tcf, ax=ax[1,0])
# --------------------------------------------------------------------------- #a
# 4. percentual deviation scaled by the maximal value
ax[1,1].set_title("abs(analytic-simulated) / max(analytic) * 100")
#tcf = ax.tricontourf(x, y, perc_devi_from_anal, cmap=plt.cm.seismic, vmin=-maxvel, vmax=maxvel)
ax[1,1].set_aspect('equal')
tcf = ax[1,1].tricontourf(x, y, perc_devi_from_anal, cmap=plt.cm.Greys, vmin=0.0, vmax=maxvel)
ax[1,1].scatter(x,y, s=0.1, color='black', marker='.')

fig.colorbar(tcf, ax=ax[1,1])
# --------------------------------------------------------------------------- #a

#plt.savefig('foo.png', dpi=500)
plt.show()

