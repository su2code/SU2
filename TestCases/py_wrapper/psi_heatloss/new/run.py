#!/usr/bin/env python

## \file run.py
#  \brief turbulent premixed dump combustor simulation (PSI flame) 
# phi=0.5, methane-air, U=40 m/s
#  \version 8.1.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2024, SU2 Contributors (cf. AUTHORS.md)
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

import sys
import pysu2
import numpy as np
from mpi4py import MPI
import statistics as st

# unburnt temperature of the propane-air mixture
# flame temperature of the methane-air mixture (phi=0.5, P=5)
Tf = 1777

Tu = 673.0
T_ref= 298
Pu = 5.0
phi = 0.5
# unburnt density at P=5
rho_u = 2.52
# unburnt thermal conductivity of methane-air at phi=0.5 (P=5) 
k_u = 0.0523
# unburnt heat capacity of methane-air at phi=0.5 (P=1)
cp_u = 1144.0

# P = rho*R*T
# 5 = 2.55 * R * 673
# R = 0.0029

# ################################################################## #
# create a function for the initial progress variable                # 
# ################################################################## #
def initC(coord):
    x = coord[0]
    y = coord[1]
    #z = coord[2]
    #print("x,y = ",x," ",y) 
    C = [0.0,276000]
    # location where the flame should be
    flame_x = 0.012
    if (x < flame_x):
      C[0] = 0.0
      #C[1] = 276000
    else:
      C[0] = 1.0
      #C[1] = 0.0

    return C   

# ################################################################## #
# loop over all vertices and set the species progress variable       # 
# ################################################################## #
def SetInitialSpecies(SU2Driver):
    allCoords = SU2Driver.Coordinates()
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    print("index of species solver = ",iSPECIESSOLVER)
    nVarsSpecies = SU2Driver.GetNumberSolverVars(iSPECIESSOLVER)
    print("number of species solver variables:",nVarsSpecies)
    for iPoint in range(SU2Driver.GetNumberNodes() - SU2Driver.GetNumberHaloNodes()):
      coord = allCoords.Get(iPoint) 
      C = initC(coord)
      # now update the initial condition 
      SU2Driver.SetSolutionVector(iSPECIESSOLVER, iPoint, C)

def SetInitialVelocity(SU2Driver):
    allCoords = SU2Driver.Coordinates()
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    print("index of FLOW solver = ",iFLOWSOLVER)
    nVarsFlow = SU2Driver.GetNumberSolverVars(iFLOWSOLVER)
    print("number of flow solver variables:",nVarsFlow)
    for iPoint in range(SU2Driver.GetNumberNodes() - SU2Driver.GetNumberHaloNodes()):
      coord = allCoords.Get(iPoint) 
      C = initC(coord)
      # now update the initial condition 
      SU2Driver.SetSolutionVector(iFLOWSOLVER, iPoint,[C])

def update_temperature(SU2Driver, iPoint):
    # first, get the progress variable
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    # returns a list
    C = SU2Driver.GetSolutionVector(iSPECIESSOLVER, iPoint)
    T = Tu*(1-C[0]) + Tf*C[0]
    P0 = 101325.0
    # kg/kmol.K
    R = 8.314
    # kg/kmol
    M = 28.5 
    RHO = P0/(R*T/M)
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    solvar = list(SU2Driver.GetSolutionVector(iFLOWSOLVER, iPoint))
    primvar = list(SU2Driver.GetPrimitiveVector(iFLOWSOLVER, iPoint))
    # the list with names
    solindex = getsolvar(SU2Driver)
    primindex = SU2Driver.GetPrimitiveIndices()
    #print("primindex = ",primindex)
    # the actual values
    #print("solindex = ",solindex)
    #print("primvar = ",primvar)
    #print("solvar=",solvar)
    iTEMP = solindex.get("TEMPERATURE")
    iDENSITY = primindex.get("DENSITY")
    #iDIFFUSIVITY = primindex.get("DENSITY")
    #print("itemp=",iTEMP)
    #print("irho=",iDENSITY)
    #print("T=",T)
    solvar[iTEMP] = T
    #
    #
    SU2Driver.SetSolutionVector(iFLOWSOLVER, iPoint, solvar)

    # also set the primitive variable list 
    #primvar[iDENSITY] = 1.225 #RHO 
    #primvar[iTEMP] = T
    #SU2Driver.SetPrimitiveVector(iFLOWSOLVER, iPoint, primvar)
    
    # how do we get for the scalar solver the diffusion coefficient?



def zimont(SU2Driver, iPoint):

    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    iSSTSOLVER = SU2Driver.GetSolverIndices()['SST']
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    primindex = SU2Driver.GetPrimitiveIndices()
    primvar = list(SU2Driver.GetPrimitiveVector(iFLOWSOLVER, iPoint))

    solindex = getsolvar(SU2Driver)
    T = solindex.get("TEMPERATURE")
    iDENSITY = primindex.get("DENSITY")
    iMU = primindex.get("LAMINAR_VISCOSITY")
    C = SU2Driver.GetSolutionVector(iSPECIESSOLVER, iPoint)
  
   # Xf=(C[1]/1000000)
    Slu = 0.232
   # if Xf>-0.3:
   # Slu = -0.3168*(Xf**5)+0.8312*(Xf**4)+1.3077*(Xf**3)+0.8763*(Xf**2)+0.3087*(Xf)+0.0445
   # else:
      # Slu=0  
    

     # Slu= -0.3168*(h_tot**5)+0.8312*(h_tot**4)+1.3077*(h_tot**3)+0.8763*(h_tot**2)+0.3087*(h_tot)+0.0445
   # elif (-0.6 < h_tot) and (h_tot < -0.5):
     # Slu= 0.0791*h_tot+0.0247
   # else :
     # Slu=0

    rho = primvar[iDENSITY]  
    mu = primvar[iMU]  
    nu=mu/rho
    tke, dissipation = SU2Driver.GetSolutionVector(iSSTSOLVER,iPoint)
    gradc = SU2Driver.GetGradient(iSPECIESSOLVER,iPoint,0)
    norm_gradc = np.sqrt(gradc[0]*gradc[0] + gradc[1]*gradc[1])

    up = np.sqrt((2.0/3.0) * tke )
    lt = (0.09**0.75) * (tke**1.5) / dissipation
    Re = up*lt/nu
    Le = 1.0
    Ut = Slu * (1.0 + (0.46/Le)*np.power(Re,0.25)*np.power(up/Slu,0.3)*np.power(Pu,0.2))

    norm_gradc = np.sqrt(gradc[0]*gradc[0] + gradc[1]*gradc[1])

    Sc = rho_u * Ut * norm_gradc
    
    # print("Sc = ",Sc)
    # print("Slu = ",Slu)
    return Sc

def getsolvar(SU2Driver):
    primindex = SU2Driver.GetPrimitiveIndices()
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    nVars = SU2Driver.GetNumberSolverVars(iFLOWSOLVER)
    varindex = primindex.copy()
    for prim in varindex.copy():
      if varindex[prim] >=nVars:
        del varindex[prim]
    varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))
    return varindex

def getcp(T,coeff):
    #Low
    a1=coeff[0]
    a2=coeff[1]
    a3=coeff[2]
    a4=coeff[3]
    a5=coeff[4]
    #High
    a8=coeff[7]
    a9=coeff[8]
    a10=coeff[9]
    a11=coeff[10]
    a12=coeff[11]

    if T<=1000:
       cp=8.314*(a1+T*a2+a3*pow(T,2)+a4*pow(T,3)+a5*pow(T,4))
    elif T>1000:   
       cp=8.314*(a8+T*a9+a10*pow(T,2)+a11*pow(T,3)+a12*pow(T,4))

    return cp   


def getenthalpy(SU2Driver, iPoint):
    
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    C = SU2Driver.GetSolutionVector(iSPECIESSOLVER, iPoint)
    x_CH4= 0.049900+C[0]*(0-0.049900)
    x_O2=0.199601+C[0]*(0.105042-0.199601)
    x_N2=0.750499+C[0]*(0.789916-0.750499)
    x_CO2=0+C[0]*(0.052521-0)
    x_H2O=0+C[0]*(0.052521-0)
    
    h=7856.742599
    Twall=300

    h_ref_CH4=-x_CH4*74800
    h_ref_CO2=-x_CO2*393500
    h_ref_H2O=-x_H2O*241800
    h_chem=h_ref_CH4+h_ref_CO2+h_ref_H2O

    T_guess=1500
    tol=1e-3
    alpha=0.3

    CH4_coeff=[5.14911468,-0.0136622009,0.0000491453921,-0.0000000484246767,0.0000000000166603441,-10246.5983,-4.63848842,1.65326226,0.01002631,-3.31661E-06,5.36483E-10,-3.14697E-14,-10009.5936,9.90506283]
    O2_coeff=[3.78245636,-0.002996734,9.8473E-06,-9.6813E-09,3.24373E-12,-1063.94356,3.65767573,3.66096083,0.000656366,-1.41149E-07,2.05798E-11,-1.29913E-15,-1215.97725,3.41536184]
    N2_coeff=[3.53100528,-0.000123661,-5.02999E-07,2.43531E-09,-1.40881E-12,-1046.97628,2.96747038,2.95257637,0.0013969,-4.92632E-07,7.8601E-11,-4.60755E-15,-923.948688,5.87188762]
    CO2_coeff=[2.356813,0.00898413,-7.12206E-06,2.4573E-09,-1.42885E-13,-48371.971,9.9009035,4.6365111,0.002741457,-9.95898E-07,1.60387E-10,-9.16199E-15,-49024.904,-1.9348955]
    H2O_coeff=[4.1986352,-0.002036402,6.52034E-06,-5.48793E-09,1.77197E-12,-30293.726,-0.84900901,2.6770389,0.002973182,-7.73769E-07,9.44335E-11,-4.269E-15,-29885.894,6.88255]

    def f(T_guess,h_chem):
      T_ref=298
      T_inlet=673
      h_sens=x_CH4*st.mean([getcp(T_guess,CH4_coeff),getcp(T_inlet,CH4_coeff)])*(T_guess-T_inlet)+x_O2*st.mean([getcp(T_guess,O2_coeff),getcp(T_inlet,O2_coeff)])*(T_guess-T_inlet)+x_N2*st.mean([getcp(T_guess,N2_coeff),getcp(T_inlet,N2_coeff)])*(T_guess-T_inlet)+x_CO2*st.mean([getcp(T_guess,CO2_coeff),getcp(T_ref,CO2_coeff)])*(T_guess-T_ref)+x_H2O*st.mean([getcp(T_guess,H2O_coeff),getcp(T_ref,H2O_coeff)])*(T_guess-T_ref)
      ht=h_sens+h_chem
      return (ht-h)

    def fprime(T_guess,h_chem):
      return (f(T_guess+1e-6,h_chem)-f(T_guess,h_chem))/1e-6
    
    if (C[0]!=0):
      while(abs(f(T_guess,h_chem))>tol):
          T_guess=T_guess-alpha*((f(T_guess,h_chem)/fprime(T_guess,h_chem)))
    elif (C[0]==0):
          T_guess=673
    
    h_loss=x_CH4*st.mean([getcp(Twall,CH4_coeff),getcp(T_guess,CH4_coeff)])*(Twall-T_guess)+x_O2*st.mean([getcp(Twall,O2_coeff),getcp(T_guess,O2_coeff)])*(Twall-T_guess)+x_N2*st.mean([getcp(Twall,N2_coeff),getcp(T_guess,N2_coeff)])*(Twall-T_guess)+x_CO2*st.mean([getcp(Twall,CO2_coeff),getcp(T_guess,CO2_coeff)])*(Twall-T_guess)+x_H2O*st.mean([getcp(Twall,H2O_coeff),getcp(T_guess,H2O_coeff)])*(Twall-T_guess)
    h_bc=((h+h_loss)*1000)/28.2

    return h_bc


# def ApplyScalar(driver,marker_ids):
#     for marker_id in marker_ids:
#       if marker_id < 0:
#         continue

#       #h1=[]
#       #h2=[]
#       #print('Marker id',marker_id)
#       for i_vertex in range(driver.GetNumberMarkerNodes(marker_id)):
        
#         h1=0
#         h2=-130213
#         hf=[h1,h2]
#         driver.SetMarkerCustomScalar(marker_id, i_vertex, hf)
        

def main():
  """
  Run the flow solver with a custom inlet (function of time and space).
  """
  comm = MPI.COMM_WORLD
  #comm = 0

  # Initialize the primal driver of SU2, this includes solver preprocessing.
  try:
    driver = pysu2.CSinglezoneDriver('psi.cfg', 1, comm)
  except TypeError as exception:
    print('A TypeError occured in pysu2.CSinglezoneDriver : ', exception)
    raise

  print("\n------------------------------ Begin Solver -----------------------------")
  sys.stdout.flush()

  nDim = driver.GetNumberDimensions()

  # index to the flow solver
  # C.FLOW
  # INC.FLOW
  # HEAT
  # FLAMELET
  # SPECIES
  # SA
  # SST
  iFLOWSOLVER = driver.GetSolverIndices()['INC.FLOW']
  print("index of flow solver = ",iFLOWSOLVER)
  iSPECIESSOLVER = driver.GetSolverIndices()['SPECIES']
  print("index of species solver = ",iSPECIESSOLVER)
  iSSTSOLVER = driver.GetSolverIndices()['SST']
  print("index of turbulence solver = ",iSSTSOLVER)


  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()
  print("indices of primitives=",primindex)
  print("number of primitives:",len(primindex))

  nElem = driver.GetNumberElements()
  print("number of elements:",nElem)

  nVars = driver.GetNumberSolverVars(iFLOWSOLVER)
  print("number of flow solver variables:",nVars)

  nVarsSpecies = driver.GetNumberSolverVars(iSPECIESSOLVER)
  print("number of species solver variables:",nVarsSpecies)
  nVarsTurb = driver.GetNumberSolverVars(iSSTSOLVER)
  print("number of turbulence solver variables:",nVarsTurb)

  varindex = primindex.copy()
  for prim in varindex.copy():
    if varindex[prim] >=nVars:
      del varindex[prim]
  varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))


# it is possible to get the solver type by doing
#  iFLOWSOLVER = driver.GetSolverIndices()['INC.FLOW']
# with the solver type we then get the solver variables using:

#  nVars = driver.GetNumberSolverVars(iFLOWSOLVER)
# and the solver variable names:
#  print("solver variable names:",varindex)

# we can overwrite the solution using:
# driver.SetSolutionVector(iSolver, iPoint, solutionVector)
#

  print("solver variable names:",varindex)
  iDENSITY = primindex.get("DENSITY")
  print("index of density = ",iDENSITY)

  index_Vel = varindex.get("VELOCITY_X")
  print("index of velocity = ",index_Vel)
  index_Temp = varindex.get("TEMPERATURE")
  custom_source_vector = [0.0 for i in range(nVars)]
  S = [0.0 for i in range(nVarsSpecies)]
  print("custom source vector = ", custom_source_vector)

  #print("max. number of inner iterations: ",driver.GetNumberInnerIter());
  #print("max nr of outer iterations: ",driver.GetNumberOuterIter());

  # set initial condition
  print("Start calling SetInitialSpecies")
  SetInitialSpecies(driver)
  print("End calling SetInitialSpecies")

  # super important to actually push the commands.
  sys.stdout.flush()
  
  all_marker_ids = driver.GetMarkerIndices()
  marker_names = ['wall_top','wall_side']
  #print('Marker id trial',all_marker_ids)
  marker_ids = []
  for name in marker_names:
    marker_ids.append(all_marker_ids[name] if name in all_marker_ids else -1)
 
  # run 500 iterations
  for inner_iter in range(701):

   # ApplyScalar(driver,marker_ids)
    driver.Preprocess(inner_iter)
    driver.Run()

    # set the source term, per point
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      #print("inode=",i_node)
      # add source term:
      # default TFC of Zimont: rho*Sc = rho_u * U_t * grad(c)
      S = [zimont(driver,i_node),0.0]
      driver.SetPointCustomSource(iSPECIESSOLVER, i_node,S)
      
      custom_source_vector[index_Temp] = -0.92*0.0284*50.5*1e6*(zimont(driver,i_node)) #should be 50.5 instead of 32.9
      #S_c = [x * constant for x in S]
      driver.SetPointCustomSource(iFLOWSOLVER, i_node,custom_source_vector)
      # at this point we also need to update the temperature based on the progress variable:
      # sete the temperature to T = c*Tf + (1-c)*Tu
      # update_temperature(driver, i_node)
      

    driver.Postprocess()
    driver.Update()



    #driver.Monitor(inner_iter)

    driver.Output(inner_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
