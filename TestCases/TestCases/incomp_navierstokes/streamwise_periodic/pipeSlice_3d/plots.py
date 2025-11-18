#! /usr/bin/python3
# --------------------------------------------------------------------------- #
# Kattmann, 16.07.2019
# This python script provides some plots to test the match between analytical
# and simulated solution for a 3D circular laminar pipe flow, either from 
# streamwise periodic simulation or the outlet of a suitable long pipe.
#
# requires: surface_flow.dat (SURFACE_TECPLOT_ASCII) in current directory
#
# output: plots (opened in separate window, not saved)
#
# optional: which plots to show
showLineplot       = True
show2Dsurfaceplots = False
show3Dplots        = False
# --------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

# --------------------------------------------------------------------------- #
# Import data from surface_flow.dat into pandas dataframe
data = pd.read_csv("surface_flow.dat", nrows=4264, skiprows=3, sep='\t', header=None)
x = data[0][:]
y = data[1][:]
vel_z = data[6][:]

# Create Delaunay surface triangulation from scatterd dataset
points2D = np.vstack([x,y]).T
tri = Delaunay(points2D)

# --------------------------------------------------------------------------- #
# Create analytic solution vector on the same points as the imported data
dynanmic_vsicosity = 1.8e-5
pressure_drop = 1e-3
domain_length = 5e-4
radius = 5e-3

analytic_sol = -1/(4*dynanmic_vsicosity) * (-pressure_drop/domain_length) * \
    (radius**2 - ((x**2 + y**2)**(0.5))**2 )

perc_devi_from_anal = abs(analytic_sol - vel_z) / max(analytic_sol) * 100
maxvel = max(abs(perc_devi_from_anal)) # get absolute maximum of dataset

# --------------------------------------------------------------------------- #
# Plot velocity on line from domain midpoint to wall
if showLineplot:
    plt.close()

    # interpolator (ip) for simulated and analytical dataset
    ip_sim = LinearNDInterpolator(tri, vel_z)
    ip_ana = LinearNDInterpolator(tri, analytic_sol)
    # line (which lies on the x-axis) where values will be interpolated
    n_sample_points = 30
    x_line = np.linspace(0, radius-5e-6, n_sample_points)
    y_line = np.zeros(n_sample_points)
    ip_pos = np.vstack((x_line,y_line)).T

    ax = plt.axes()
    plt.plot(ip_sim(ip_pos), x_line, color='b', marker='', linestyle='--', linewidth=3, label='simulated')
    plt.plot(ip_ana(ip_pos), x_line, color='r', marker='', linestyle=':' , linewidth=3, label='analytical')
    plt.legend()
    plt.title('Velocity profile: analytic vs simulated (interpolated values)')
    plt.xlabel('velocity [m/s]')
    plt.ylabel('radius [m]')
    ax.set_aspect(aspect=max(ip_sim(ip_pos)) / max(x_line)) # make plot square
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.grid(True, linestyle='--')
    plt.show()

# --------------------------------------------------------------------------- #
# Plot various 2D surface plots of sim. and analy. data
if show2Dsurfaceplots:
    plt.close()

    fig, ax = plt.subplots(2,2)

    # 1. analytical solution
    ax_tmp = ax[0,0]

    tcf = ax_tmp.tricontourf(x, y, abs(analytic_sol))
    ax_tmp.scatter(x,y, s=0.1, color='black', marker='.')

    ax_tmp.set_title("Analytical solution")
    ax_tmp.set_aspect('equal')
    fig.colorbar(tcf, ax=ax_tmp)

    # 2. simulated solution
    ax_tmp = ax[1,0]

    tcf = ax_tmp.tricontourf(x, y, vel_z)
    ax_tmp.scatter(x,y, s=0.1, color='black', marker='.')

    ax_tmp.set_title("Simulated solution")
    ax_tmp.set_aspect('equal')
    fig.colorbar(tcf, ax=ax_tmp)

    # 3. absolute value deviation between analytic and simulated
    ax_tmp = ax[0,1]

    tcf = ax_tmp.tricontourf(x, y, abs(analytic_sol-vel_z), cmap=plt.cm.Greys)
    ax_tmp.scatter(x,y, s=0.1, color='black', marker='.')

    ax_tmp.set_title("abs(analytic-simulated)")
    ax_tmp.set_aspect('equal')
    fig.colorbar(tcf, ax=ax_tmp)

    # 4. percentual deviation scaled by the maximal value
    ax_tmp = ax[1,1]

    tcf = ax_tmp.tricontourf(x, y, perc_devi_from_anal, cmap=plt.cm.Greys, vmin=0.0, vmax=maxvel)
    ax_tmp.scatter(x,y, s=0.1, color='black', marker='.')

    ax_tmp.set_title("abs(analytic-simulated) / max(analytic) * 100")
    ax_tmp.set_aspect('equal')
    fig.colorbar(tcf, ax=ax_tmp)

    plt.show()

# --------------------------------------------------------------------------- #
if show3Dplots:
    # Plot 3D surfaces of sim. and analy. data
    plt.close()

    # Scatter plot deviation
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.scatter(x, y, perc_devi_from_anal)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z-Velocity deviation [%]')

    plt.show()

    # Surface plot deviation
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_trisurf(x, y, perc_devi_from_anal, triangles=tri.simplices, cmap='jet', linewidth=0)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z-Velocity deviation [%]')
    fig.colorbar(surf)

    plt.show()

    # Surface plot of velocity
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    surf = ax.plot_trisurf(x, y, vel_z, triangles=tri.simplices, cmap='jet', linewidth=0)
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z-Velocity [m/s]')
    fig.colorbar(surf)

    plt.show()
