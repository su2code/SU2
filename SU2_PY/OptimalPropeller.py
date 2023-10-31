## \file OptimalPropeller.py
#  \brief Python script for generating the ActuatorDisk.dat file.
#  \author E. Saetta, L. Russo, R. Tognaccini
#  \version 8.0.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.
# ==============================================================================================
# Name        : OptimalPropeller
# Author      : Ettore Saetta, Lorenzo Russo, Renato Tognaccini
#               Theoretical and Applied Aerodynamic Research Group (TAARG),
#               University of Naples Federico II.
# Version     : 1.0.0 - Python
# Date        : 01/09/2020
# Copyright   :
# Description : Compute the optimal load distribution along the propeller radius using
#               the inviscid theory of the optimal propeller.
# Reference   : Glauert H., Airplane Propellers, in Aerodynamic Theory, Ed. Durand W. F.,
#               Vol. IV, pp. 169 - 360, Springer, 1935.
# Input       : Interactive.
# Output      : ActuatorDisk.cfg, containing part of SU2 .cfg file.
#               ActuatorDisk.dat, containing propeller load distribution to be read by SU2_CFD.
# Note        : Python 3 or higher needed.
# ==============================================================================================

import math
import numpy as np
import pylab as pl

##########################
###     Functions      ###
##########################


def a_distribution(w0, Chi):
    """Function used to compute the value of the axial interference factor using the inviscid theory of the optimal propeller."""

    return (w0 * pow(Chi, 2)) / (pow(Chi, 2) + pow((1 + (w0)), 2))


def write_su2_config_file():
    """Write the actuator disk configuration file"""

    with open("ActuatorDisk.cfg", "w") as f:
        f.write("% Automatic generated actuator disk configuration file.\n")
        f.write("%\n")
        f.write("% The first two elements of MARKER_ACTDISK must be filled.\n")
        f.write("% An example of this file can be found in the TestCases directory.\n")
        f.write("%\n")
        f.write("% Author: Ettore Saetta, Lorenzo Russo, Renato Tognaccini.\n")
        f.write("% Theoretical and Applied Aerodynamic Research Group (TAARG),\n")
        f.write("% University of Naples Federico II\n")
        f.write("\n")
        f.write("ACTDISK_TYPE = VARIABLE_LOAD\n")
        f.write("ACTDISK_FILENAME = ActuatorDisk.dat\n")
        f.write("MARKER_ACTDISK = ( , , 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)\n")

    print("SU2 file generated!")


def write_external_file(CTrs, CPrs):
    """Function to write the actuator disk input data file"""

    with open("ActuatorDisk.dat", "w") as f:
        f.write(
            "# Automatic generated actuator disk input data file using the Optimal Propeller code.\n"
        )
        f.write("# Data file needed for the actuator disk VARIABLE_LOAD type.\n")
        f.write(
            "# The load distribution is obtained using the inviscid theory of the optimal propeller\n"
        )
        f.write("# using global data.\n")
        f.write("#\n")
        f.write("# The first three lines must be filled.\n")
        f.write("# An example of this file can be found in the TestCases directory.\n")
        f.write("#\n")
        f.write("# Author: Ettore Saetta, Lorenzo Russo, Renato Tognaccini.\n")
        f.write("# Theoretical and Applied Aerodynamic Research Group (TAARG),\n")
        f.write("# University of Naples Federico II\n")
        f.write(
            "# -------------------------------------------------------------------------------------\n"
        )
        f.write("#\n")
        f.write("MARKER_ACTDISK= \n")
        f.write("CENTER= \n")
        f.write("AXIS= \n")
        f.write("RADIUS= " + str(R) + "\n")
        f.write("ADV_RATIO= " + str(J) + "\n")
        f.write("NROW= " + str(stations) + "\n")
        f.write("# rs=r/R        dCT/drs       dCP/drs       dCR/drs\n")

        for i in range(0, stations):
            f.write(f"  {r[i]:.7f}     {CTrs[i]:.7f}     {CPrs[i]:.7f}     0.0\n")


##########################
###        Main        ###
##########################

print("------------------ Optimal Propeller vsn 7.0.6 ------------------")
print("| Computation of the optimal dCT/dr and dCP/dr distributions.   |")
print("| Based on the inviscid theory of the optimal propeller.        |")
print("|                                                               |")
print("| This code is used to generate the actuator disk input data    |")
print("| file needed for the VARIABLE_LOAD actuator disk type          |")
print("| implemented in SU2 7.0.6.                                     |")
print("|                                                               |")
print("| Author: Ettore Saetta, Lorenzo Russo, Renato Tognaccini.      |")
print("| Theoretical and Applied Aerodynamic Research Group (TAARG),   |")
print("| University of Naples Federico II.                             |")
print("-----------------------------------------------------------------")
print("")
print("Warning: present version requires input in SI units.")
print("")

# Number of radial stations in input.
stations = int(input("Number of radial stations: "))

# Resize the vectors using the number of radial stations.
r = np.empty(stations)
dCp = np.empty(stations)
a_new = np.empty(stations)
a_old = np.empty(stations)
a_0 = np.empty(stations)
a_optimal = np.empty(stations)
ap_optimal = np.empty(stations)

# Thrust coefficient in input.
Ct = float(input("\nCT (Renard definition): "))

# Propeller radius in input.
R = float(input("\nR (propeller radius [m]): "))

# Hub radius in input.
rhub = float(input("\nr_hub (hub radius [m]): "))

# Advance ratio in input.
J = float(input("\nJ (advance ratio): "))

# Freestream velocity in input.
Vinf = float(input("\nVinf (m/s): "))

# Asking if the tip loss Prandtl correction function needs to be used.
prandtl_input = input("\nUsing tip loss Prandtl correction? (<y>/n): ")

if prandtl_input.lower() in ["yes", "y", ""]:
    # Number of propeller blades in input.
    N = int(input("\nN (number of propeller blades): "))
    prandtl_correction = True
else:
    prandtl_correction = False

# Computation of the non-dimensional hub radius.
rs_hub = rhub / R

# Computation of the non-dimensional radial stations.
for i in range(1, stations + 1):
    r[i - 1] = i / float(stations)
    if r[i - 1] <= rs_hub:
        i_hub = i - 1

# Computation of the propeller diameter.
D = 2 * R
# Computation of the propeller angular velocity (Rounds/s).
n = Vinf / (D * J)
# Computation of the propeller angular velocity (Rad/s).
Omega = n * 2 * math.pi

# Computation of the tip loss Prandtl correction function F.
if prandtl_correction:
    F = (2 / math.pi) * np.arccos(
        np.exp(-0.5 * N * (1 - r) * np.sqrt(1 + pow(Omega * R / Vinf, 2)))
    )
else:
    F = np.ones((stations))

# Computation of the non-dimensional radius chi=Omega*r/Vinf.
chi = Omega * r * R / Vinf


eps = 5e-20
# Computation of the propeller radial stations spacing.
h = 1.0 / stations

# Computation of the first try induced velocity distribution.
w = (2 / np.power(Vinf, 2)) * (
    (-1 / Vinf)
    + np.sqrt(
        1 + ((np.power(D, 4) * (Ct) * np.power(n, 2)) / (np.power(Vinf, 2) * np.pi * r))
    )
)

# Computation of the first try Lagrange moltiplicator.
w_0 = sum(w) / (Vinf * stations)

# Computation of the first try axial interference factor distribution.
for i in range(0, stations):
    a_0[i] = a_distribution(w_0 * F[i], chi[i])

# Computation of the thrust coefficient distribution
dCt_0 = math.pi * J**2 * r * (1 + a_0) * a_0

# Computation of the total thrust coefficient.
Ct_0 = sum(h * dCt_0[i_hub:])

# Compute the error with respect to the thrust coefficient given in input.
err_0 = Ct_0 - Ct
print("\n\nCONVERGENCE HISTORY:")
print(err_0)

# Computation of the second try Lagrange moltiplicator.
w_old = w_0 + 0.1

# Computation of the second try axial interference factor distribution.
for i in range(0, stations):
    a_old[i] = a_distribution(w_old * F[i], chi[i])

# Computation of the thrust coefficient distribution
dCt_old = math.pi * J**2 * r * (1 + a_old) * a_old

# Computation of the total thrust coefficient.
Ct_old = sum(h * dCt_old[i_hub:])

# Compute the error with respect to the thrust coefficient given in input.
err_old = Ct_old - Ct
print(err_old)

##########################
###     Iterations     ###
##########################
# Iterate using the false position methods.
# Based on the error from the thrust coefficient given in input.
iteration = 2
err_new = err_old
while math.fabs(err_new) >= eps and err_0 != err_old:
    iteration += 1

    # Computation of the new Lagrange moltiplicator value based on the false position method.
    w_new = (w_old * err_0 - w_0 * err_old) / (err_0 - err_old)

    # Computation of the new axial interference factor distribution.
    for i in range(0, stations):
        a_new[i] = a_distribution(w_new * F[i], chi[i])

    # Computation of the new thrust coefficient distribution.
    dCt_new = math.pi * J**2 * r * (1 + a_new) * a_new

    # Computation of the new total thrust coefficient.
    Ct_new = sum(h * dCt_new[i_hub:])

    # Computation of the total thrust coefficient error with respect to the input value.
    err_new = Ct_new - Ct
    print(err_new)

    # Updating the stored values for the next iteration.
    err_0 = err_old
    err_old = err_new

    w_0 = w_old
    w_old = w_new

# Computation of the correct axial and rotational interference factors (a and ap).
for i in range(0, stations):
    a_optimal[i] = a_distribution(w_new * F[i], chi[i])
    ap_optimal[i] = (w_new * F[i]) * (
        (1 + w_new * F[i]) / (chi[i] * chi[i] + math.pow(1 + w_new * F[i], 2))
    )

# Computation of the correct thrust coefficient distribution.
dCt_optimal = math.pi * J**2 * r * (1 + a_optimal) * a_optimal

# Computation of the correct power coefficient distribution.
for i in range(0, stations):
    dCp[i] = (R * 4 * math.pi / (math.pow(n, 3) * math.pow(D, 5))) * (
        math.pow(Vinf, 3) * math.pow(1 + a_optimal[i], 2) * a_optimal[i] * r[i] * R
        + math.pow(Omega, 2)
        * Vinf
        * (1 + a_optimal[i])
        * math.pow(ap_optimal[i], 2)
        * math.pow(r[i] * R, 3)
    )

##########################
###   Check Results    ###
##########################
# Computation of the total power coefficient.
Cp = sum(h * dCp[i_hub:])

# Computation of the total thrust coefficient.
Ct_optimal = sum(h * dCt_optimal[i_hub:])

# Computation of the static pressure jump distribution.
DeltaP = dCt_optimal * (2 * Vinf**2) / (J**2 * math.pi * r)

# Computation of the thrust over density (T) using the static pressure jump distribution.
T = sum(2 * math.pi * r[i_hub:] * math.pow(R, 2) * h * DeltaP[i_hub:])

# Computation of the thrust coefficient using T.
Ct_Renard = T / (math.pow(n, 2) * math.pow(D, 4))

# Computation of the efficiency.
eta = J * (Ct_optimal / Cp)

# Screen output used to check that everything worked correcty.
print("%%%%%%%%%%%%%%%%%%%%%%%%% CHECK OUTPUT VALUES %%%%%%%%%%%%%%%%%%%%%%%%%")
print(f"       dCT distribution integral: {Ct_optimal:.4f}")
print(f"       dCT computed using the static pressure jump: {Ct_Renard:.4f}")
print(f"       dCP distribution integral: {Cp:.4f}")
print(f"       Thrust over Density (T/rho): {T:.4f} [N*m^3/kg]")
print(f"       Efficiency eta: {eta:.4f}")
print(f"       w0/Vinf: {w_new:.4f}")
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

##########################
###    File Writing    ###
##########################

# Write the corresponding SU2 configuration file
write_su2_config_file()

# Write the actuator disk data file. This is the actuator disk input data file.
write_external_file(dCt_optimal, dCp)

##########################
###        Plots       ###
##########################
# Automatically plot the computed propeller performance.

pl.figure(1)
pl.plot(r, dCt_optimal, "r", markersize=4, label="$\\frac{dCT}{d\overline{r}}$")
pl.plot(r, dCp, "k", markersize=4, label="$\\frac{dCP}{d\overline{r}}$")
pl.grid(True)
pl.legend(numpoints=3)
pl.xlabel("$\overline{r}$")
pl.ylabel("")
pl.title("Load Distribution")

pl.figure(2)
pl.plot(chi, a_optimal, "r", markersize=4, label="$a$")
pl.plot(chi, ap_optimal, "k", markersize=4, label="$a^1$")
pl.grid(True)
pl.legend(numpoints=3)
pl.xlabel("$\chi$")
pl.ylabel("")
pl.title("Interference Factors")

if prandtl_correction:
    pl.figure(3)
    pl.plot(r, F, "k", markersize=4)
    pl.grid(True)
    pl.xlabel("$\overline{r}$")
    pl.ylabel("$F(\overline{r})$")
    pl.title("Tip Loss Prandtl Correction Function")

pl.show()
