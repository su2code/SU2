#!/usr/bin/env python

## \file SolidSolverTester.py
#  \brief Structural solver tester (one or two degree of freedom) used for testing the Py wrapper for external FSI coupling.
#  \author David Thomas
#  \version 7.0.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

import os, sys, shutil, copy
import numpy as np
import scipy as sp
import scipy.linalg as linalg
from math import *
from util import switch

# ----------------------------------------------------------------------
#  Config class
# ----------------------------------------------------------------------

class Point:
  """ Description. """

  def __init__(self):
    self.Coord0 = np.zeros((3,1))
    self.Coord = np.zeros((3,1))
    self.Coord_n = np.zeros((3,1))
    self.Vel = np.zeros((3,1))
    self.Vel_n = np.zeros((3,1))
    self.Force = np.zeros((3,1))

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

  def updateCoordVel(self):
    self.Coord_n = np.copy(self.Coord)
    self.Vel_n = np.copy(self.Vel)

class Solver:
  """Description"""
  
  def __init__(self, config_fileName):
    """ Description. """

    self.Config_file = config_fileName
    self.Config = {}    

    print("\n------------------------------ Configuring the structural tester solver for FSI simulation ------------------------------")
    self.__readConfig()

    self.Mesh_file = self.Config['MESH_FILE']
    self.FSI_marker = self.Config['MOVING_MARKER']
    self.Unsteady = (self.Config['TIME_MARCHING']=="YES")
    if self.Unsteady:
      print('Dynamic computation.')
    if self.Config['STRUCT_TYPE'] == "AIRFOIL":
      self.nDof = 2
      print("Structural model : pitching-plunging airfoil.")
    else:
      self.nDof = 0
      

    # Structural properties
    self.m = self.Config['SPRING_MASS'] # airfoil mass [kg]
    self.Kh = self.Config['SPRING_STIFFNESS'] # plunging stiffness [N/m]
    self.Ka = self.Config['TORSIONAL_STIFFNESS'] # pitching stiffness [N]
    self.If = self.Config['INERTIA_FLEXURAL'] # inertia around flexural axis [kg m^2]
    self.Ch = self.Config['SPRING_DAMPING'] # plunging damping [Ns/m]
    self.Ca = self.Config['TORSIONAL_DAMPING'] # pitching damping [Ns]
    self.c = self.Config['CORD'] # airfoil cord [m]
    self.b = self.c/2.0 # airfoil semi cord [m]
    self.xf = self.Config['FLEXURAL_AXIS'] # position of the flexural axis [m]
    self.xCG = self.Config['GRAVITY_CENTER'] # position of the center of gravity [m]
    self.S = self.m*(self.xCG - self.xf)
    self.h0 = self.Config['INITIAL_DISP']
    self.a0 = self.Config['INITIAL_ANGLE']
    self.startTime = self.Config['START_TIME']
    self.stopTime = self.Config['STOP_TIME']
    self.deltaT = self.Config['DELTA_T']
    self.rhoAlphaGen = self.Config['RHO']
    
    self.nDim= int()
    self.nElem = int()
    self.nPoint = int()
    self.nMarker = int()
    self.node = []
    self.markers = {}

    print("\n------------------------------ Reading the SU2 mesh ------------------------------")
    self.__readSU2Mesh()

    print("\n------------------------------ Creating the structural model ------------------------------")
    self.__setStructuralMatrices()
  
    print("\n------------------------------ Setting the integration parameters ------------------------------")
    self.__setIntegrationParameters()
    self.__setInitialConditions()

  def __readConfig(self):
    """ Description. """

    with open(self.Config_file) as configfile:
      while 1:
	line = configfile.readline()
	if not line:
	  break

        # remove line returns
        line = line.strip('\r\n')
        # make sure it has useful data
        if (not "=" in line) or (line[0] == '%'):
          continue
        # split across equal sign
        line = line.split("=",1)
        this_param = line[0].strip()
        this_value = line[1].strip()

        for case in switch(this_param):
	  #integer values
	  #if case("NB_FSI_ITER")		:
	    #self.Config[this_param] = int(this_value)
	    #break

	  #float values
	  if case("DELTA_T")			: pass
	  if case("START_TIME")		      	: pass
	  if case("STOP_TIME")		      	: pass
	  if case("SPRING_MASS")		: pass
	  if case("INERTIA_FLEXURAL")		: pass
	  if case("SPRING_STIFFNESS")		: pass
	  if case("SPRING_DAMPING")		: pass
	  if case("TORSIONAL_STIFFNESS")	: pass
	  if case("TORSIONAL_DAMPING")		: pass
	  if case("CORD")		      	: pass
	  if case("FLEXURAL_AXIS")	      	: pass
	  if case("GRAVITY_CENTER")	      	: pass
	  if case("INITIAL_DISP")	      	: pass
	  if case("INITIAL_ANGLE")	      	: pass
	  if case("RHO")	      		: 
	    self.Config[this_param] = float(this_value)
	    break

	  #string values
	  if case("TIME_MARCHING")	: pass
	  if case("MESH_FILE")			: pass
	  if case("CSD_SOLVER")		      	: pass
	  if case("MOVING_MARKER")		: pass
	  if case("STRUCT_TYPE")		:
	    self.Config[this_param] = this_value
	    break

 	  if case():
	    print(this_param + " is an invalid option !")
            break

  def __readSU2Mesh(self):
    """ Description. """
    
    with open(self.Mesh_file, 'r') as meshfile:
      print('Opened mesh file ' + self.Mesh_file + '.')
      while 1:
        line = meshfile.readline()
	if not line:
	  break
        
	pos = line.find('NDIM')
	if pos != -1:
	  line = line.strip('\r\n')
          line = line.split("=",1)
	  self.nDim = int(line[1])
	  continue
	
	pos = line.find('NELEM')
	if pos != -1:
	  line = line.strip('\r\n')
          line = line.split("=",1)
	  self.nElem = int(line[1])
	  continue

	pos = line.find('NPOIN')
	if pos != -1:
	  line = line.strip('\r\n')
          line = line.split("=",1)
	  self.nPoint = int(line[1])
          for iPoint in range(self.nPoint):
	    self.node.append(Point())
	    line = meshfile.readline()
	    line = line.strip('\r\n')
	    line = line.split(' ',self.nDim)
	    x = float(line[0])
	    y = float(line[1])
            z = 0.0
	    if self.nDim == 3:
	      z = float(line[2])
	    self.node[iPoint].SetCoord((x,y,z))
            self.node[iPoint].SetCoord0((x,y,z))
	    self.node[iPoint].SetCoord_n((x,y,z))
	  continue

	pos = line.find('NMARK')
	if pos != -1:
	  line = line.strip('\r\n')
          line = line.split("=",1)
	  self.nMarker = int(line[1])
	  continue

	pos = line.find('MARKER_TAG')
	if pos != -1:
	  line = line.strip('\r\n')
	  line = line.replace(" ", "")
          line = line.split("=",1)
	  markerTag = line[1]
	  if markerTag == self.FSI_marker:
	    self.markers[markerTag] = []
	    line = meshfile.readline()
	    line = line.strip('\r\n')
	    line = line.split("=",1)
	    nElem = int(line[1])
	    for iElem in range(nElem):
	      line = meshfile.readline()
	      line = line.strip('\r\n')
	      line = line.split(' ',1)
	      elemType = int(line[0])
	      if elemType == 3:
	        nodes = line[1].split(' ', 1)
		if not int(nodes[0]) in self.markers[markerTag]:
		   self.markers[markerTag].append(int(nodes[0]))
		if not int(nodes[1]) in self.markers[markerTag]:
		   self.markers[markerTag].append(int(nodes[1]))
	      else:
		print("Element type {} is not recognized !!".format(elemType))
	    continue
	  else:
	    continue

    print("Number of dimensions: {}".format(self.nDim))
    print("Number of elements: {}".format(self.nElem))
    print("Number of point: {}".format(self.nPoint))
    print("Number of markers: {}".format(self.nMarker))
    if len(self.markers) > 0:
      print("Moving marker(s):")
      for mark in self.markers.keys():
        print(mark)

  def __setStructuralMatrices(self):
    """ Descriptions. """
    
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

    if self.Config['STRUCT_TYPE'] == "AIRFOIL":
      print('Setting pitching-plunging airfoil system')
      print('Number of DOF : ')
      
      self.centerOfRotation = np.zeros((3,1))
      self.centerOfRotation[0] = self.xf
      self.centerOfRotation_n = np.zeros((3,1))
      self.centerOfRotation_n[0] = self.xf
      self.M[0][0] = self.m
      self.M[0][1] = self.S
      self.M[1][0] = self.S
      self.M[1][1] = self.If
      self.C[0][0] = self.Ch
      self.C[1][1] = self.Ca
      self.K[0][0] = self.Kh
      self.K[1][1] = self.Ka
      
      print('Airfoil mass : {} [kg]'.format(self.m))
      print('Airfoil cord : {} [m]'.format(self.c))
      print('Position of the flexural axis (from the leading edge) : {} [m]'.format(self.xf))
      print('Position of the center of gravity (from the leading edge)  : {} [m]'.format(self.xCG))
      print('Inertia around the flexural axis : {} [kg m^2]'.format(self.If))
      print('Static unbalance : {} [kg m]'.format(self.S))
      print('Plunging stiffness : {} [N/m]'.format(self.Kh))
      print('Plunging dampimg : {} [Ns/m]'.format(self.Ch))
      print('Pitching stiffness : {} [N]'.format(self.Ka))
      print('Pitching dampimg : {} [Ns]'.format(self.Ca))

  def __setIntegrationParameters(self):
    """ Description. """
    
    self.alpha_m = (2.0*self.rhoAlphaGen-1.0)/(self.rhoAlphaGen+1.0)
    self.alpha_f = (self.rhoAlphaGen)/(self.rhoAlphaGen+1.0)
    self.gamma = 0.5+self.alpha_f-self.alpha_m
    self.beta = 0.25*(self.gamma+0.5)**2

    self.gammaPrime = self.gamma/(self.deltaT*self.beta)
    self.betaPrime = (1.0-self.alpha_m)/((self.deltaT**2)*self.beta*(1.0-self.alpha_f))

    print('Time integration with the alpha-generalized algorithm.')
    print('rho : {}'.format(self.rhoAlphaGen))
    print('alpha_m : {}'.format(self.alpha_m))
    print('alpha_f : {}'.format(self.alpha_f))
    print('gamma : {}'.format(self.gamma))
    print('beta : {}'.format(self.beta))
    print('gammaPrime : {}'.format(self.gammaPrime))
    print('betaPrime : {}'.format(self.betaPrime))

  def __setInitialConditions(self):
    """ Description. """

    print('Setting initial conditions.')
    self.__reset(self.F)
    self.__reset(self.q_n)

    print('Initial plunge displacement : {} [m]'.format(self.h0))
    self.q[0] = self.h0
    if self.nDof == 2:
      print('Initial pitch angle : {} [rad]'.format(self.a0))
      self.q[1] = self.a0

    RHS = np.zeros((self.nDof,1))
    RHS += self.F
    RHS -= self.C.dot(self.qdot) 
    RHS -= self.K.dot(self.q)
    self.qddot = linalg.solve(self.M, RHS) 
    self.a = np.copy(self.qddot)

    self.centerOfRotation[1] = self.q[0]

  def __reset(self, vector):
    """ Description. """

    for ii in range(vector.shape[0]):
      vector[ii] = 0.0

  def __computeInterfacePosVel(self, initialize):
    """ Description. """

    dTheta = 0.0
    dPhi = 0.0
    newCenter = np.zeros((3,1))
    Centerdot = np.zeros((3,1))
    newVel = np.zeros((3,1))
    
    if self.nDof == 2:
      dPsi = -(self.q[1] - self.q_n[1])
      newCenter[0] = self.centerOfRotation[0]
      newCenter[1] = -self.q[0]
      newCenter[2] = self.centerOfRotation[2]
      Centerdot[0] = 0.0
      Centerdot[1] = -self.qdot[0]
      Centerdot[2] = 0.0
      psidot = self.qdot[1]

    cosPsi = cos(dPsi)
    sinPsi = sin(dPsi)

    rotMatrix = np.array([[cosPsi, -sinPsi, 0.0],[sinPsi, cosPsi, 0.0],[0.0, 0.0, 1.0]])

    for iMarker in self.markers.keys():
      vertexList = self.markers[iMarker]
      for iPoint in vertexList:
        Coord = self.node[iPoint].GetCoord()
        Coord_n = self.node[iPoint].GetCoord_n()

        if self.Unsteady:
	  r = Coord_n - self.centerOfRotation_n
	else:
	  r = Coord - self.centerOfRotation

	rotCoord = rotMatrix.dot(r)

        newCoord = newCenter + rotCoord
        newVel[0] = Centerdot[0]+psidot*(newCoord[1]-newCenter[1])
	newVel[1] = Centerdot[1]-psidot*(newCoord[0]-newCenter[0])
	newVel[2] = Centerdot[2]+0.0

        self.node[iPoint].SetCoord((newCoord[0], newCoord[1], newCoord[2]))
        self.node[iPoint].SetVel((newVel[0], newVel[1], newVel[2]))

	if initialize:
	  self.node[iPoint].SetCoord_n((newCoord[0], newCoord[1], newCoord[2]))
	  self.node[iPoint].SetVel_n((newVel[0], newVel[1], newVel[2]))

    self.centerOfRotation = np.copy(newCenter)

  def __temporalIteration(self):
    """ Description. """
    
    eps = 1e-6 

    self.__SetLoads()

    # Prediction step
    self.__reset(self.qddot)
    self.__reset(self.a)

    self.a += (self.alpha_f)/(1-self.alpha_m)*self.qddot_n
    self.a -= (self.alpha_m)/(1-self.alpha_m)*self.a_n

    self.q = np.copy(self.q_n)
    self.q += self.deltaT*self.qdot_n
    self.q += (0.5-self.beta)*self.deltaT*self.deltaT*self.a_n
    self.q += self.deltaT*self.deltaT*self.beta*self.a

    self.qdot = np.copy(self.qdot_n)
    self.qdot += (1-self.gamma)*self.deltaT*self.a_n
    self.qdot += self.deltaT*self.gamma*self.a

    # Correction step
    res = self.__ComputeResidual()

    while linalg.norm(res) >= eps:
      St = self.__TangentOperator()
      Deltaq = -1*(linalg.solve(St,res))
      self.q += Deltaq
      self.qdot += self.gammaPrime*Deltaq
      self.qddot += self.betaPrime*Deltaq
      res = self.__ComputeResidual()

    self.a += (1-self.alpha_f)/(1-self.alpha_m)*self.qddot
      

  def __SetLoads(self):
    """ Description """
    
    makerID = self.markers.keys()[0]
    nodeList = self.markers[makerID]

    FX = 0.0
    FY = 0.0
    FZ = 0.0
    MZ = 0.0

    for iPoint in nodeList:
      Force = self.node[iPoint].GetForce()
      Coord = self.node[iPoint].GetCoord()
      FX += float(Force[0])
      FY += float(Force[1])
      FZ += float(Force[2])
      MZ += float(Force[1]*(Coord[0]-self.centerOfRotation[0])-Force[0]*(Coord[1]-self.centerOfRotation[1]))

    self.F[0] = -FY
    self.F[1] = -MZ

  def __ComputeResidual(self):
    """ Description. """

    res = self.M.dot(self.qddot) + self.C.dot(self.qdot) + self.K.dot(self.q) - self.F

    return res

  def __TangentOperator(self):
    """ Description. """

    # The problem is linear, so the tangent operator is straightforward.
    St = self.betaPrime*self.M + self.gammaPrime*self.C + self.K

    return St

  def exit(self):
    """ Description. """

    print("\n**************** Exiting the structural tester solver ****************")

  def run(self,t0,t1):
    """ Description. """

    self.__temporalIteration()

    print("Time\tDisp 1\tDisp2\tVel 1\tVel2\tAcc 1\tAcc 2")
    print(str(t1) + '\t' + str(float(self.q[0])) + '\t' + str(float(self.q[1])) + '\t' + str(float(self.qdot[0])) + '\t' + str(float(self.qdot[1])) + '\t' + str(float(self.qddot[0])) + '\t' + str(float(self.qddot[1])))

    self.__computeInterfacePosVel(False)

  def setInitialDisplacements(self):
    """ Description. """

    self.__computeInterfacePosVel(True)

  def writeSolution(self, time, FSIIter, TimeIter, NbTimeIter):
    """ Description. """

    if time == 0:
      histFile = open('StructHistory.dat', "w")
      histFile.write("Time\tDisp 1\tDisp2\tVel 1\tVel2\tAcc 1\tAcc 2\tAccVar 1\tAccVar 2\n")
    else:
      histFile = open('StructHistory.dat', "a")
    histFile.write(str(time) + '\t' + str(float(self.q[0])) + '\t' + str(float(self.q[1])) + '\t' + str(float(self.qdot[0])) + '\t' + str(float(self.qdot[1])) + '\t' + str(float(self.qddot[0])) + '\t' + str(float(self.qddot[1])) + '\t' + str(float(self.a[0])) + '\t' + str(float(self.a[1])) + '\n')
    histFile.close()

  def updateSolution(self):
    """ Description. """

    self.q_n = np.copy(self.q)
    self.qdot_n = np.copy(self.qdot)
    self.qddot_n = np.copy(self.qddot)
    self.a_n = np.copy(self.a)
    self.__reset(self.q)
    self.__reset(self.qdot)
    self.__reset(self.qddot)
    self.__reset(self.a)

    makerID = self.markers.keys()[0]
    nodeList = self.markers[makerID]

    for iPoint in nodeList:
      self.node[iPoint].updateCoordVel()

    self.centerOfRotation_n = np.copy(self.centerOfRotation)

  def applyload(self, iVertex, fx, fy, fz, time):
    """ Description """

    makerID = self.markers.keys()[0]
    iPoint = self.getInterfaceNodeGlobalIndex(makerID, iVertex)
    self.node[iPoint].SetForce((fx,fy,fz))

  def getFSIMarkerID(self):
    """ Description. """
    
    list = self.markers.keys()
    return list[0]

  def getNumberOfSolidInterfaceNodes(self, markerID):
    """ Description. """

    return len(self.markers[markerID])

  def getInterfaceNodeGlobalIndex(self, markerID, iVertex):
    """ Description. """

    return self.markers[markerID][iVertex]

  def getInterfaceNodePosX(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    return float(Coord[0])

  def getInterfaceNodePosY(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    return float(Coord[1])

  def getInterfaceNodePosZ(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    return float(Coord[2])

  def getInterfaceNodeDispX(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    Coord0 = self.node[iPoint].GetCoord0()
    return float(Coord[0]-Coord0[0])

  def getInterfaceNodeDispY(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    Coord0 = self.node[iPoint].GetCoord0()
    return float(Coord[1]-Coord0[1])

  def getInterfaceNodeDispZ(self, markerID, iVertex):
    """ Desciption. """

    iPoint = self.markers[markerID][iVertex]
    Coord = self.node[iPoint].GetCoord()
    Coord0 = self.node[iPoint].GetCoord0()
    return float(Coord[2]-Coord0[2])

  def getInterfaceNodeVelX(self, markerID, iVertex):
    """ Description """

    iPoint = self.markers[markerID][iVertex]
    Vel = self.node[iPoint].GetVel()
    return float(Vel[0])

  def getInterfaceNodeVelY(self, markerID, iVertex):
    """ Description """

    iPoint = self.markers[markerID][iVertex]
    Vel = self.node[iPoint].GetVel()
    return float(Vel[1])

  def getInterfaceNodeVelZ(self, markerID, iVertex):
    """ Description """

    iPoint = self.markers[markerID][iVertex]
    Vel = self.node[iPoint].GetVel()
    return float(Vel[2])

  def getInterfaceNodeVelXNm1(self, markerID, iVertex):
    """ Description """

    iPoint = self.markers[markerID][iVertex]
    Vel = self.node[iPoint].GetVel_n()
    return float(Vel[0])

  def getInterfaceNodeVelYNm1(self, markerID, iVertex):
    """ Description """

    iPoint = self.markers[markerID][iVertex]
    Vel = self.node[iPoint].GetVel_n()
    return float(Vel[1])

  def getInterfaceNodeVelZNm1(self, markerID, iVertex):
    """ Description """

    iPoint = self.markers[markerID][iVertex]
    Vel = self.node[iPoint].GetVel_n()
    return float(Vel[2])

  def getRotationCenterPosX(self):
    """ Description. """

    return float(self.centerOfRotation[0])

  def getRotationCenterPosY(self):
    """ Description. """

    return float(self.centerOfRotation[1])

  def getRotationCenterPosZ(self):
    """ Description. """

    return float(self.centerOfRotation[2])
