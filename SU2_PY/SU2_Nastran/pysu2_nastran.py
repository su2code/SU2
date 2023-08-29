#!/usr/bin/env python

## \file pysu2_nastran.py
#  \brief Structural solver using Nastran models
#  \authors Nicola Fonzi, Vittorio Cavalieri, based on the work of David Thomas
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import numpy as np
import scipy.linalg as linalg
from math import *

# ----------------------------------------------------------------------
#  Config class
# ----------------------------------------------------------------------


class ImposedMotionClass:
    def __init__(self, time0, typeOfMotion, parameters, mode):

        self.time0 = time0
        self.typeOfMotion = typeOfMotion
        self.mode = mode

        self.amplitude = parameters["AMPLITUDE"]
        self.timeStart = parameters["TIME_START"]
        if "TIME_STOP" in parameters.keys():
            self.timeStop = parameters["TIME_STOP"]
        else:
            self.timeStop = inf

        if self.typeOfMotion == "SINUSOIDAL":
            self.bias = parameters["BIAS"]
            self.frequency = parameters["FREQUENCY"]

        elif self.typeOfMotion == "BLENDED_STEP":
            self.kmax = parameters["K_MAX"]
            self.vinf = parameters["V_INF"]
            self.lref = parameters["L_REF"]
            self.tmax = 2 * pi / self.kmax * self.lref / self.vinf
            self.omega0 = 1 / 2 * self.kmax

        else:
            raise Exception(
                "Imposed function {} not found, please implement it in pysu2_nastran.py".format(
                    self.tipo
                )
            )

    def GetDispl(self, time):
        time = time - self.time0 - self.timeStart
        if self.typeOfMotion == "SINUSOIDAL":
            if (time < 0.0) or (time > self.timeStop):
                return 0.0
            return self.bias + self.amplitude * sin(2 * pi * self.frequency * time)

        if self.typeOfMotion == "BLENDED_STEP":
            if (time < 0.0) or (time > self.timeStop):
                return 0.0
            if time < self.tmax:
                return (
                    self.amplitude
                    / 2.0
                    * (1.0 - cos(self.omega0 * time * self.vinf / self.lref))
                )
            return self.amplitude

    def GetVel(self, time):
        time = time - self.time0 - self.timeStart

        if self.typeOfMotion == "SINUSOIDAL":
            if (time < 0.0) or (time > self.timeStop):
                return 0.0
            return (
                self.amplitude
                * cos(2 * pi * self.frequency * time)
                * 2
                * pi
                * self.frequency
            )

        if self.typeOfMotion == "BLENDED_STEP":
            if (time < 0.0) or (time > self.timeStop):
                return 0.0
            if time < self.tmax:
                return (
                    self.amplitude
                    / 2.0
                    * sin(self.omega0 * time * self.vinf / self.lref)
                    * (self.omega0 * self.vinf / self.lref)
                )
            return 0.0

    def GetAcc(self, time):
        time = time - self.time0 - self.timeStart

        if self.typeOfMotion == "SINUSOIDAL":
            if (time < 0.0) or (time > self.timeStop):
                return 0.0
            return (
                -self.amplitude
                * sin(2 * pi * self.frequency * time)
                * (2 * pi * self.frequency) ** 2
            )

        if self.typeOfMotion == "BLENDED_STEP":
            if (time < 0.0) or (time > self.timeStop):
                return 0.0
            if time < self.tmax:
                return (
                    self.amplitude
                    / 2.0
                    * cos(self.omega0 * time * self.vinf / self.lref)
                    * (self.omega0 * self.vinf / self.lref) ** 2
                )
            return 0.0


class RefSystem:
    def __init__(self):
        self.CID = 0
        self.RID = 0
        self.Origin = np.array([[0.0], [0.0], [0.0]])
        self.Rot = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

    def SetOrigin(self, A):
        AX, AY, AZ = A
        self.Origin[0] = AX
        self.Origin[1] = AY
        self.Origin[2] = AZ

    def SetRotMatrix(self, x, y, z):
        self.Rot = np.array(
            [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]]
        )

    def SetCID(self, CID):
        self.CID = CID

    def SetRID(self, RID):
        self.RID = RID

    def GetOrigin(self):
        return self.Origin

    def GetRotMatrix(self):
        return self.Rot

    def GetRID(self):
        return self.RID

    def GetCID(self):
        return self.CID


class Point:
    """
    Class containing data regarding all the structural nodes.
    Coord0: Coordinates at the initial time iteration.
    Coord: Coordinates at the current time iteration.
    Coord_n: Coordinates at the previous time iteration.
    Vel: Velocity at the current time iteration.
    Vel_n: Velocity at the previous time iteration.
    Force: Nodal force provided by the aerodynamics.
    ID: ID of the node.
    CP: Coordinate system definition of the position.
    CD: Coordinate system definition of the output coming from Nastran.
    """

    def __init__(self):
        self.Coord0 = np.zeros((3, 1))
        self.Coord = np.zeros((3, 1))
        self.Coord_n = np.zeros((3, 1))
        self.Vel = np.zeros((3, 1))
        self.Vel_n = np.zeros((3, 1))
        self.Force = np.zeros((3, 1))
        self.ID = 0
        self.CP = 0
        self.CD = 0

    def GetCoord0(self):
        return self.Coord0

    def GetCoord(self):
        return self.Coord

    def GetCoord_n(self):
        return self.Coord_n

    def GetVel(self):
        return self.Vel

    def GetVel_n(self):
        return self.Vel_n

    def GetForce(self):
        return self.Force

    def GetID(self):
        return self.ID

    def GetCP(self):
        return self.CP

    def GetCD(self):
        return self.CD

    def SetCoord0(self, val_Coord):
        x, y, z = val_Coord
        self.Coord0[0] = x
        self.Coord0[1] = y
        self.Coord0[2] = z

    def SetCoord(self, val_Coord):
        x, y, z = val_Coord
        self.Coord[0] = x
        self.Coord[1] = y
        self.Coord[2] = z

    def SetCoord_n(self, val_Coord):
        x, y, z = val_Coord
        self.Coord_n[0] = x
        self.Coord_n[1] = y
        self.Coord_n[2] = z

    def SetVel(self, val_Vel):
        vx, vy, vz = val_Vel
        self.Vel[0] = vx
        self.Vel[1] = vy
        self.Vel[2] = vz

    def SetVel_n(self, val_Vel):
        vx, vy, vz = val_Vel
        self.Vel_n[0] = vx
        self.Vel_n[1] = vy
        self.Vel_n[2] = vz

    def SetForce(self, val_Force):
        fx, fy, fz = val_Force
        self.Force[0] = fx
        self.Force[1] = fy
        self.Force[2] = fz

    def SetID(self, ID):
        self.ID = ID

    def SetCP(self, CP):
        self.CP = CP

    def SetCD(self, CD):
        self.CD = CD

    def updateCoordVel(self):
        self.Coord_n = np.copy(self.Coord)
        self.Vel_n = np.copy(self.Vel)


class Solver:
    """
    Structural solver main class.
    It contains all the required methods for the coupling with SU2.
    """

    def __init__(self, config_fileName, ImposedMotion):
        """
        Constructor of the structural solver class.
        """

        self.Config_file = config_fileName
        self.Config = {}

        print("\n")
        print(
            " Configuring the structural tester solver for FSI simulation ".center(
                80, "-"
            )
        )
        self.__readConfig()

        self.Mesh_file = self.Config["MESH_FILE"]
        self.Punch_file = self.Config["PUNCH_FILE"]
        self.FSI_marker = self.Config["MOVING_MARKER"]
        self.Unsteady = self.Config["TIME_MARCHING"] == "YES"
        self.ImposedMotion = ImposedMotion
        if self.Unsteady:
            print("Dynamic computation.")
        self.nDof = self.Config["NMODES"]
        print("Reading number of modes from file")

        # Structural properties
        print("Reading the modal and stiffnes matrix from file")
        self.ModalDamping = self.Config["MODAL_DAMPING"]
        if self.ModalDamping == 0:
            print("The structural model is undamped")
        else:
            print("Assuming {}% of modal damping".format(self.ModalDamping * 100))

        self.deltaT = self.Config["DELTA_T"]
        self.rhoAlphaGen = self.Config["RHO"]

        self.nPoint = int()
        self.nMarker = int()
        self.nRefSys = int()
        self.node = []
        self.markers = {}
        self.refsystems = []
        self.ImposedMotionToSet = True
        self.ImposedMotionFunction = []

        print("\n")
        print(" Reading the mesh ".center(80, "-"))
        self.__readNastranMesh()

        print("\n")
        print(" Creating the structural model ".center(80, "-"))
        self.__setStructuralMatrices()

        print("\n")
        print(" Setting the integration parameters ".center(80, "-"))
        self.__setIntegrationParameters()
        self.__setInitialConditions()

        # Prepare the output file
        if self.Config["RESTART_SOL"] == "NO":
            histFile = open("StructHistoryModal.dat", "w")
            header = "Time\t" + "Time Iteration\t" + "FSI Iteration\t"
            for imode in range(self.nDof):
                header = (
                    header
                    + "q"
                    + str(imode + 1)
                    + "\t"
                    + "qdot"
                    + str(imode + 1)
                    + "\t"
                    + "qddot"
                    + str(imode + 1)
                    + "\t"
                )
            header = header + "\n"
            histFile.write(header)
            histFile.close()
        else:
            self.__setRestart()

    def __readConfig(self):
        """
        This methods obtains the configuration options from the structural solver input
        file.
        """

        with open(self.Config_file) as configfile:
            while 1:
                line = configfile.readline()
                if not line:
                    break

                # remove line returns
                line = line.strip("\r\n")
                # make sure it has useful data
                if (not "=" in line) or (line[0] == "%"):
                    continue
                # split across equal sign
                line = line.split("=", 1)
                this_param = line[0].strip()
                this_value = line[1].strip()

                # integer values
                if (this_param == "NMODES") or (this_param == "RESTART_ITER"):
                    self.Config[this_param] = int(this_value)

                # float values
                elif (
                    (this_param == "DELTA_T")
                    or (this_param == "MODAL_DAMPING")
                    or (this_param == "RHO")
                ):
                    self.Config[this_param] = float(this_value)

                # string values
                elif (
                    (this_param == "TIME_MARCHING")
                    or (this_param == "MESH_FILE")
                    or (this_param == "PUNCH_FILE")
                    or (this_param == "RESTART_SOL")
                    or (this_param == "MOVING_MARKER")
                ):
                    self.Config[this_param] = this_value

                # lists values
                elif (
                    (this_param == "INITIAL_MODES")
                    or (this_param == "IMPOSED_MODES")
                    or (this_param == "IMPOSED_PARAMETERS")
                ):
                    self.Config[this_param] = eval(this_value)

                else:
                    raise Exception("{} is an invalid option !".format(this_param))

    def __readNastranMesh(self):
        """
        This method reads the nastran 3D mesh.
        """

        def nastran_float(s):
            if s.find("E") == -1:
                s = s.replace("-", "e-")
                s = s.replace("+", "e+")
                if s[0] == "e":
                    s = s[1:]
            return float(s)

        self.nMarker = 0
        self.nPoint = 0
        self.nRefSys = 0

        with open(self.Mesh_file, "r") as meshfile:
            print("Opened mesh file " + self.Mesh_file + ".")
            while 1:
                line = meshfile.readline()
                if not line:
                    break

                pos = line.find("GRID")
                if pos == 30:
                    line = line.strip("\r\n")
                    self.node.append(Point())
                    line = line[30:]
                    ID = int(line[8:16])
                    CP = self.__checkBlankField(line[16:24])
                    x = nastran_float(line[24:32])
                    y = nastran_float(line[32:40])
                    z = nastran_float(line[40:48])
                    if CP != 0:
                        for iRefSys in range(self.nRefSys):
                            if self.refsystems[iRefSys].GetCID() == CP:
                                break
                        if self.refsystems[iRefSys].GetCID() != CP:
                            raise Exception(
                                "Definition reference {} system not found".format(CP)
                            )
                        DeltaPos = self.refsystems[iRefSys].GetOrigin()
                        RotatedPos = (
                            self.refsystems[iRefSys]
                            .GetRotMatrix()
                            .dot(np.array([[x], [y], [z]]))
                        )
                        x = RotatedPos[0] + DeltaPos[0]
                        y = RotatedPos[1] + DeltaPos[1]
                        z = RotatedPos[2] + DeltaPos[2]
                    CD = self.__checkBlankField(line[48:56])
                    self.node[self.nPoint].SetCoord((x, y, z))
                    self.node[self.nPoint].SetID(ID)
                    self.node[self.nPoint].SetCP(CP)
                    self.node[self.nPoint].SetCD(CD)
                    self.node[self.nPoint].SetCoord0((x, y, z))
                    self.node[self.nPoint].SetCoord_n((x, y, z))
                    self.nPoint += 1
                    continue

                pos = line.find("CORD2R")
                if pos == 30:
                    line = line.strip("\r\n")
                    self.refsystems.append(RefSystem())
                    line = line[30:]
                    CID = int(line[8:16])
                    self.refsystems[self.nRefSys].SetCID(CID)
                    RID = int(line[16:24])
                    if RID != 0:
                        raise Exception(
                            "ERROR: Reference system {} must be defined with respect to global reference system".format(
                                CID
                            )
                        )
                    self.refsystems[self.nRefSys].SetRID(RID)
                    AX = nastran_float(line[24:32])
                    AY = nastran_float(line[32:40])
                    AZ = nastran_float(line[40:48])
                    BX = nastran_float(line[48:56])
                    BY = nastran_float(line[56:64])
                    BZ = nastran_float(line[64:72])
                    z_direction = np.array([BX - AX, BY - AY, BZ - AZ])
                    z_direction = z_direction / linalg.norm(z_direction)
                    line = meshfile.readline()
                    line = line.strip("\r\n")
                    line = line[30:]
                    CX = nastran_float(line[8:16])
                    CY = nastran_float(line[16:24])
                    CZ = nastran_float(line[24:32])
                    y_direction = np.cross(z_direction, [CX - AX, CY - AY, CZ - AZ])
                    y_direction = y_direction / linalg.norm(y_direction)
                    x_direction = np.cross(y_direction, z_direction)
                    x_direction = x_direction / linalg.norm(x_direction)
                    self.refsystems[self.nRefSys].SetRotMatrix(
                        x_direction, y_direction, z_direction
                    )
                    self.refsystems[self.nRefSys].SetOrigin((AX, AY, AZ))
                    self.nRefSys += 1
                    continue

                pos = line.find("SET1")
                if pos == 30:
                    line = line.strip("\r\n")
                    line = line[37:]
                    line = line.split()
                    existValue = True
                    markerTag = line.pop(0)
                    self.markers[markerTag] = []
                    while existValue:
                        if line[0] == "+":
                            line = meshfile.readline()
                            line = line.strip("\r\n")
                            line = line[37:]
                            line = line.split()
                        ID = int(line.pop(0))
                        for iPoint in range(self.nPoint):
                            if self.node[iPoint].GetID() == ID:
                                break
                        if (iPoint == (self.nPoint - 1)) and (
                            self.node[iPoint].GetID() != ID
                        ):
                            raise Exception(
                                "Point {} in the set {} was not found in the mesh".format(
                                    ID, markerTag
                                )
                            )
                        self.markers[markerTag].append(iPoint)
                        existValue = len(line) >= 1
                    self.nMarker += 1
                    continue

        if not any(self.FSI_marker in key for key in self.markers.keys()):
            raise Exception("The FSI marker was not found in the available sets")

        self.markers[self.FSI_marker].sort()

        print("Number of points: {}".format(self.nPoint))
        print("Number of markers: {}".format(self.nMarker))
        print("Number of reference systems: {}".format(self.nRefSys))
        print("Moving marker: {}".format(self.FSI_marker))
        print(
            "Number of points in the moving marker".format(
                len(self.markers[self.FSI_marker])
            )
        )

    def __checkBlankField(self, string):
        """
        This method considers that Nastran apply 0 when the reference system is not specified
        """

        if string == " " * 8:
            return int(0)
        return int(string)

    def __setStructuralMatrices(self):
        """
        This method reads the punch file and obtains the modal shapes and modal stiffnesses.
        """

        self.M = np.zeros((self.nDof, self.nDof))
        self.K = np.zeros((self.nDof, self.nDof))
        self.C = np.zeros((self.nDof, self.nDof))

        self.q = np.zeros((self.nDof, 1))
        self.qdot = np.zeros((self.nDof, 1))
        self.qddot = np.zeros((self.nDof, 1))
        self.a = np.zeros((self.nDof, 1))

        self.q_n = np.zeros((self.nDof, 1))
        self.qdot_n = np.zeros((self.nDof, 1))
        self.qddot_n = np.zeros((self.nDof, 1))
        self.a_n = np.zeros((self.nDof, 1))

        self.F = np.zeros((self.nDof, 1))

        self.Ux = np.zeros((self.nPoint, self.nDof))
        self.Uy = np.zeros((self.nPoint, self.nDof))
        self.Uz = np.zeros((self.nPoint, self.nDof))

        with open(self.Punch_file, "r") as punchfile:
            print("Opened punch file " + self.Punch_file + ".")
            while 1:
                line = punchfile.readline()
                if not line:
                    break

                pos = line.find("MODE ")
                if pos != -1:
                    line = line.strip("\r\n").split()
                    n = int(line[5])
                    imode = n - 1
                    k_i = float(line[2])
                    self.M[imode][imode] = 1
                    self.K[imode][imode] = k_i
                    w_i = sqrt(k_i)
                    self.C[imode][imode] = 2 * self.ModalDamping * w_i
                    iPoint = 0
                    for indexIter in range(self.nPoint):
                        line = punchfile.readline()
                        line = line.strip("\r\n").split()
                        if line[1] == "G":
                            ux = float(line[2])
                            uy = float(line[3])
                            uz = float(line[4])
                            if self.node[iPoint].GetCD() != 0:
                                for iRefSys in range(self.nRefSys):
                                    if (
                                        self.refsystems[iRefSys].GetCID()
                                        == self.node[iPoint].GetCD()
                                    ):
                                        break
                                if (
                                    self.refsystems[iRefSys].GetCID()
                                    != self.node[iPoint].GetCD()
                                ):
                                    raise Exception(
                                        "Output reference {} system not found".format(
                                            self.node[iPoint].GetCD()
                                        )
                                    )
                                RotatedOutput = (
                                    self.refsystems[iRefSys]
                                    .GetRotMatrix()
                                    .dot(np.array([[ux], [uy], [uz]]))
                                )
                                ux = RotatedOutput[0]
                                uy = RotatedOutput[1]
                                uz = RotatedOutput[2]
                            self.Ux[iPoint][imode] = ux
                            self.Uy[iPoint][imode] = uy
                            self.Uz[iPoint][imode] = uz
                            iPoint = iPoint + 1
                            line = punchfile.readline()
                        if line[1] == "S":
                            line = punchfile.readline()

                    if n == self.nDof:
                        break

        self.__setNonDiagonalStructuralMatrices()

        self.UxT = self.Ux.transpose()
        self.UyT = self.Uy.transpose()
        self.UzT = self.Uz.transpose()

        if n < self.nDof:
            raise Exception(
                "ERROR: available {} degrees of freedom instead of {} as requested".format(
                    n, self.nDof
                )
            )
        else:
            print("Using {} degrees of freedom".format(n))

    def __setNonDiagonalStructuralMatrices(self):
        """
        This method is part of an advanced feature of this solver that allows to set
        nondiagonal matrices for the structural modes.
        """

        K_updated = self.__readNonDiagonalMatrix("NDK")
        M_updated = self.__readNonDiagonalMatrix("NDM")
        C_updated = self.__readNonDiagonalMatrix("NDC")
        if K_updated and M_updated and (not C_updated):
            print("Setting modal damping")
            self.__setNonDiagonalDamping()
        elif (not K_updated) and (not M_updated):
            print("Modal stiffness and mass matrices are diagonal")
        elif (not K_updated) and M_updated:
            raise Exception("Non-Diagonal stiffness matrix is missing")
        elif (not M_updated) and K_updated:
            raise Exception("Non-Diagonal mass matrix is missing")

    def __readNonDiagonalMatrix(self, keyword):
        """
        This method reads from the punch file the definition of nondiagonal structural
        matrices.
        """

        matrixUpdated = False

        with open(self.Punch_file, "r") as punchfile:

            while 1:
                line = punchfile.readline()
                if not line:
                    break

                pos = line.find(keyword)
                if pos != -1:
                    while 1:
                        line = punchfile.readline()
                        line = line.strip("\r\n").split()
                        if line[0] != "-CONT-":
                            i = int(line[0]) - 1
                            j = 0
                            el = line[1:]
                            ne = len(el)
                        elif line[0] == "-CONT-":
                            el = line[1:]
                            ne = len(el)
                        if keyword == "NDK":
                            self.K[i][j : j + ne] = np.array(el)
                        elif keyword == "NDM":
                            self.M[i][j : j + ne] = np.array(el)
                        elif keyword == "NDC":
                            self.C[i][j : j + ne] = np.array(el)
                        j = j + ne
                        if i + 1 == self.nDof and j == self.nDof:
                            matrixUpdated = True
                            break

        return matrixUpdated

    def __setNonDiagonalDamping(self):

        D, V = linalg.eig(self.K, self.M)
        D = D.real
        D = np.sqrt(D)
        Mmodal = ((V.transpose()).dot(self.M)).dot(V)
        Mmodal = np.diag(Mmodal)
        C = 2 * self.ModalDamping * np.multiply(D, Mmodal)
        C = np.diag(C)
        Vinv = linalg.inv(V)
        C = C.dot(Vinv)
        VinvT = Vinv.transpose()
        self.C = VinvT.dot(C)

    def __setIntegrationParameters(self):
        """
        This method uses the time step size to define the integration parameters.
        """

        self.alpha_m = (2.0 * self.rhoAlphaGen - 1.0) / (self.rhoAlphaGen + 1.0)
        self.alpha_f = (self.rhoAlphaGen) / (self.rhoAlphaGen + 1.0)
        self.gamma = 0.5 + self.alpha_f - self.alpha_m
        self.beta = 0.25 * (self.gamma + 0.5) ** 2

        self.gammaPrime = self.gamma / (self.deltaT * self.beta)
        self.betaPrime = (1.0 - self.alpha_m) / (
            (self.deltaT**2) * self.beta * (1.0 - self.alpha_f)
        )

        print("Time integration with the alpha-generalized algorithm.")
        print("rho : {}".format(self.rhoAlphaGen))
        print("alpha_m : {}".format(self.alpha_m))
        print("alpha_f : {}".format(self.alpha_f))
        print("gamma : {}".format(self.gamma))
        print("beta : {}".format(self.beta))
        print("gammaPrime : {}".format(self.gammaPrime))
        print("betaPrime : {}".format(self.betaPrime))

    def __setInitialConditions(self):
        """
        This method uses the list of initial modal amplitudes to set the initial conditions
        """

        print("Setting initial conditions.")

        print("Using modal amplitudes from config file")
        for imode in range(self.nDof):
            if imode in self.Config["INITIAL_MODES"].keys():
                self.q[imode] = float(self.Config["INITIAL_MODES"][imode])
                self.q_n[imode] = float(self.Config["INITIAL_MODES"][imode])

        RHS = np.zeros((self.nDof, 1))
        RHS += self.F
        RHS -= self.C.dot(self.qdot)
        RHS -= self.K.dot(self.q)
        self.qddot = linalg.solve(self.M, RHS)
        self.qddot_n = np.copy(self.qddot)
        self.a = np.copy(self.qddot)
        self.a_n = np.copy(self.qddot)

    def __reset(self, vector):
        """
        This method set to zero any vector.
        """

        for ii in range(vector.shape[0]):
            vector[ii] = 0.0

    def __computeInterfacePosVel(self, initialize):
        """
        This method uses the mode shapes to compute, based on the modal velocities, the
        nodal velocities at the interface.
        """

        # Multiply the modal matrices with modal amplitudes
        X_vel = self.Ux.dot(self.qdot)
        Y_vel = self.Uy.dot(self.qdot)
        Z_vel = self.Uz.dot(self.qdot)

        X_disp = self.Ux.dot(self.q)
        Y_disp = self.Uy.dot(self.q)
        Z_disp = self.Uz.dot(self.q)

        for iPoint in range(self.nPoint):
            coord0 = self.node[iPoint].GetCoord0()
            self.node[iPoint].SetCoord(
                (
                    X_disp[iPoint] + coord0[0],
                    Y_disp[iPoint] + coord0[1],
                    Z_disp[iPoint] + coord0[2],
                )
            )
            self.node[iPoint].SetVel((X_vel[iPoint], Y_vel[iPoint], Z_vel[iPoint]))

            if initialize:
                self.node[iPoint].SetCoord_n(
                    (
                        X_disp[iPoint] + coord0[0],
                        Y_disp[iPoint] + coord0[1],
                        Z_disp[iPoint] + coord0[2],
                    )
                )
                self.node[iPoint].SetVel_n(
                    (X_vel[iPoint], Y_vel[iPoint], Z_vel[iPoint])
                )

    def __setRestart(self):
        """
        This method sets all the variables needed for the correct restart.
        """

        # read the Structhistory to obtain the mode amplitudes
        nM1Set = False
        nSet = False
        firstLineRead = False
        couplingLineRead = False

        with open("StructHistoryModal.dat", "r") as file:
            print("Opened history file StructHistoryModal.dat.")
            line = file.readline()
            while 1:
                line = file.readline()
                if not line:
                    break
                line = line.strip("\r\n").split()

                # The old time_0 for imposed motion can either be the first line of the StructHistoryModal, if TimeIterTreshold was -1 (immediate coupling), or the second line. In the former case, time_0 is 0.0, so it is easy to recognize it
                if not firstLineRead:
                    firstLineRead = True
                    if float(line[0]) == 0.0:
                        couplingLineRead = True
                        self.timeStartCoupling = 0.0
                else:
                    if couplingLineRead:
                        pass
                    else:
                        self.timeStartCoupling = float(line[0])
                        couplingLineRead = True

                if int(line[1]) == (self.Config["RESTART_ITER"] - 2):
                    index = 0
                    for index_mode in range(self.nDof):
                        self.q[index_mode] = float(line[index + 3])
                        self.qdot[index_mode] = float(line[index + 4])
                        self.qddot[index_mode] = float(line[index + 5])
                        index += 3
                    del index
                    # push back the mode amplitudes velocities and accelerations
                    self.__computeInterfacePosVel(True)
                    self.q_n = np.copy(self.q)
                    self.qdot_n = np.copy(self.qdot)
                    self.qddot_n = np.copy(self.qddot)
                    self.a_n = np.copy(self.a)
                    nM1Set = True
                if int(line[1]) == (self.Config["RESTART_ITER"] - 1):
                    index = 0
                    for index_mode in range(self.nDof):
                        self.q[index_mode] = float(line[index + 3])
                        self.qdot[index_mode] = float(line[index + 4])
                        self.qddot[index_mode] = float(line[index + 5])
                        index += 3
                    del index
                    self.__computeInterfacePosVel(False)
                    nSet = True
                    break

        if (not nM1Set) or (not nSet):
            raise Exception(
                "The restart iteration was not found in the structural history"
            )

    def __temporalIteration(self, time):
        """
        This method integrates in time the solution.
        """

        self.__reset(self.q)
        self.__reset(self.qdot)
        self.__reset(self.qddot)
        self.__reset(self.a)

        if not self.ImposedMotion:
            eps = 1e-6

            self.__SetLoads()

            # Prediction step

            self.a += (self.alpha_f) / (1 - self.alpha_m) * self.qddot_n
            self.a -= (self.alpha_m) / (1 - self.alpha_m) * self.a_n

            self.q = np.copy(self.q_n)
            self.q += self.deltaT * self.qdot_n
            self.q += (0.5 - self.beta) * self.deltaT * self.deltaT * self.a_n
            self.q += self.deltaT * self.deltaT * self.beta * self.a

            self.qdot = np.copy(self.qdot_n)
            self.qdot += (1 - self.gamma) * self.deltaT * self.a_n
            self.qdot += self.deltaT * self.gamma * self.a

            # Correction step
            res = self.__ComputeResidual()

            while linalg.norm(res) >= eps:
                St = self.__TangentOperator()
                Deltaq = -1 * (linalg.solve(St, res))
                self.q += Deltaq
                self.qdot += self.gammaPrime * Deltaq
                self.qddot += self.betaPrime * Deltaq
                res = self.__ComputeResidual()

            self.a += (1 - self.alpha_f) / (1 - self.alpha_m) * self.qddot
        else:
            if self.ImposedMotionToSet:
                if self.Config["RESTART_SOL"] == "NO":
                    # If yes we already set it in the __setRestart function
                    self.timeStartCoupling = time
                iImposedFunc = 0
                for imode in self.Config["IMPOSED_MODES"].keys():
                    for isuperposed in range(len(self.Config["IMPOSED_MODES"][imode])):
                        typeOfMotion = self.Config["IMPOSED_MODES"][imode][isuperposed]
                        parameters = self.Config["IMPOSED_PARAMETERS"][imode][
                            isuperposed
                        ]
                        self.ImposedMotionFunction.append(
                            ImposedMotionClass(
                                self.timeStartCoupling, typeOfMotion, parameters, imode
                            )
                        )
                        iImposedFunc += 1
                self.ImposedMotionToSet = False
            for iImposedFunc in range(len(self.ImposedMotionFunction)):
                imode = self.ImposedMotionFunction[iImposedFunc].mode
                self.q[imode] += self.ImposedMotionFunction[iImposedFunc].GetDispl(time)
                self.qdot[imode] += self.ImposedMotionFunction[iImposedFunc].GetVel(
                    time
                )
                self.qddot[imode] += self.ImposedMotionFunction[iImposedFunc].GetAcc(
                    time
                )
            self.a = np.copy(self.qddot)

    def __SetLoads(self):
        """
        This method uses the nodal forces and the mode shapes to obtain the modal forces.
        """
        nodeList = self.markers[self.FSI_marker]
        FX = np.zeros((self.nPoint, 1))
        FY = np.zeros((self.nPoint, 1))
        FZ = np.zeros((self.nPoint, 1))
        for iPoint in nodeList:
            Force = self.node[iPoint].GetForce()
            FX[iPoint] = float(Force[0])
            FY[iPoint] = float(Force[1])
            FZ[iPoint] = float(Force[2])
        self.F = self.UxT.dot(FX) + self.UyT.dot(FY) + self.UzT.dot(FZ)

    def __ComputeResidual(self):
        """
        This method computes the residual for integration.
        """

        res = (
            self.M.dot(self.qddot) + self.C.dot(self.qdot) + self.K.dot(self.q) - self.F
        )

        return res

    def __TangentOperator(self):
        """
        This method computes the tangent operator for solution.
        """

        # The problem is linear, so the tangent operator is straightforward.
        St = self.betaPrime * self.M + self.gammaPrime * self.C + self.K

        return St

    def exit(self):
        """
        This method cleanly exits the structural solver.
        """

        print(
            "\n**************** Exiting the structural tester solver ****************"
        )

    def run(self, time):
        """
        This method is the main function for advancing the solution of one time step.
        """
        self.__temporalIteration(time)
        header = "Time\t"
        for imode in range(min([self.nDof, 5])):
            header = (
                header
                + "q"
                + str(imode + 1)
                + "\t"
                + "qdot"
                + str(imode + 1)
                + "\t"
                + "qddot"
                + str(imode + 1)
                + "\t"
            )
        header = header + "\n"
        print(header)
        line = "{:6.4f}".format(time) + "\t"
        for imode in range(min([self.nDof, 5])):
            line = (
                line
                + "{:6.4f}".format(float(self.q[imode]))
                + "\t"
                + "{:6.4f}".format(float(self.qdot[imode]))
                + "\t"
                + "{:6.4f}".format(float(self.qddot[imode]))
                + "\t"
            )
        line = line + "\n"
        print(line)
        self.__computeInterfacePosVel(False)

    def activateMode(self, iMode):
        """
        This method is used to artificially set only one mode activated, thus
        with non zero amplitude.
        """
        self.__reset(self.q)
        self.q[iMode] = 1.0
        self.__computeInterfacePosVel(True)

    def setInitialDisplacements(self):
        """
        This method provides public access to the method __computeInterfacePosVel and
        sets velocities for previous time steps.
        """

        self.__computeInterfacePosVel(True)

    def writeSolution(self, time, timeIter, FSIIter):
        """
        This method is the main function for output. It writes the file StructHistoryModal.dat
        """

        # Modal History
        histFile = open("StructHistoryModal.dat", "a")
        line = str(time) + "\t" + str(timeIter) + "\t" + str(FSIIter) + "\t"
        for imode in range(self.nDof):
            line = (
                line
                + str(float(self.q[imode]))
                + "\t"
                + str(float(self.qdot[imode]))
                + "\t"
                + str(float(self.qddot[imode]))
                + "\t"
            )
        line = line + "\n"
        histFile.write(line)
        histFile.close()

    def updateSolution(self):
        """
        This method updates the solution.
        """

        self.q_n = np.copy(self.q)
        self.qdot_n = np.copy(self.qdot)
        self.qddot_n = np.copy(self.qddot)
        self.a_n = np.copy(self.a)
        self.__reset(self.q)
        self.__reset(self.qdot)
        self.__reset(self.qddot)
        self.__reset(self.a)

        for iPoint in range(self.nPoint):
            self.node[iPoint].updateCoordVel()

    def applyload(self, iVertex, fx, fy, fz):
        """
        This method can be accessed from outside to set the nodal forces.
        """
        iPoint = self.getVertexGlobalIndex(self.FSI_marker, iVertex)
        self.node[iPoint].SetForce((fx, fy, fz))

    def getNumberOfModes(self):
        """
        This method provides the number of degrees of freedom used in
        the structural solver.
        """
        return self.nDof

    def getFSIMarkerID(self):
        """
        This method provides the ID of the interface marker
        """
        return self.FSI_marker

    def getNumberOfSolidInterfaceNodes(self, markerID):

        return len(self.markers[markerID])

    def getVertexGlobalIndex(self, markerID, iVertex):

        # This solver is serial, thus global=local
        return self.markers[markerID][iVertex]

    def getInterfaceNodePosInit(self, markerID, iVertex):

        iPoint = self.markers[markerID][iVertex]
        Coord0 = self.node[iPoint].GetCoord0()
        return Coord0

    def getInterfaceNodeDisp(self, markerID, iVertex):

        iPoint = self.markers[markerID][iVertex]
        Coord = self.node[iPoint].GetCoord()
        Coord0 = self.node[iPoint].GetCoord0()
        return Coord - Coord0

    def getInterfaceNodeVel(self, markerID, iVertex):

        iPoint = self.markers[markerID][iVertex]
        Vel = self.node[iPoint].GetVel()
        return Vel

    def getInterfaceNodeVelNm1(self, markerID, iVertex):

        iPoint = self.markers[markerID][iVertex]
        Vel = self.node[iPoint].GetVel_n()
        return Vel

    def IsAHaloNode(self, markerID, iVertex):

        # There are no halo nodes in this solver as it is serial
        iPoint = self.markers[markerID][iVertex]
        halo = False
        return halo
