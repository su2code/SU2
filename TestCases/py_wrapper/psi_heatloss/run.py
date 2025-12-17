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

# unburnt temperature of the propane-air mixture
# flame temperature of the methane-air mixture (phi=0.5, P=5)
Tf = 1803
Twall=1200
Tu = 673.0
Pu = 5.0
phi = 0.5
# unburnt density at P=5
rho_u = 2.50
# unburnt thermal conductivity of methane-air at phi=0.5 (P=5) 
k_u = 0.0523
# unburnt heat capacity of methane-air at phi=0.5 (P=1)
cp_u = 1311.0

NASA_COEFFICIENTS = {
    'CH4': {
        'low_coeffs': [5.14911468, -0.013662201, 4.91454E-05, -4.84247E-08, 1.66603E-11, -10246.5983, -4.63848842],
        'high_coeffs': [1.65326226, 0.01002631, -3.31661E-06, 5.36483E-10, -3.14697E-14, -10009.5936, 9.90506283]
    },
    'O2': {
        'low_coeffs': [3.78245636, -0.002996734, 9.8473E-06, -9.6813E-09, 3.24373E-12, -1063.94356, 3.65767573],
        'high_coeffs': [3.66096083, 0.000656366, -1.41149E-07, 2.05798E-11, -1.29913E-15, -1215.97725, 3.41536184]
    },
    'N2': {
        'low_coeffs': [3.53100528, -0.000123661, -5.02999E-07, 2.43531E-09, -1.40881E-12, -1046.97628, 2.96747038],
        'high_coeffs': [2.95257637, 0.0013969, -4.92632E-07, 7.8601E-11, -4.60755E-15, -923.948688, 5.87188762]
    },
    'CO2': {
        'low_coeffs': [2.356813, 0.00898413, -7.12206E-06, 2.4573E-09, -1.42885E-13, -48371.971, 9.9009035],
        'high_coeffs': [4.6365111, 0.002741457, -9.95898E-07, 1.60387E-10, -9.16199E-15, -49024.904, -1.9348955]
    },
    'H2O': {
        'low_coeffs': [4.1986352, -0.002036402, 6.52034E-06, -5.48793E-09, 1.77197E-12, -30293.726, -0.84900901],
        'high_coeffs': [2.6770389, 0.002973182, -7.73769E-07, 9.44335E-11, -4.269E-15, -29885.894, 6.88255]
    },
    'H2': {
        'low_coeffs': [2.34433112, 0.007980521, -1.94782E-05, 2.01572E-08, -7.37612E-12, -917.935173, 0.683010238],
        'high_coeffs': [2.93286575, 0.000826608, -1.46402E-07, 1.541E-11, -6.88805E-16, -813.065581, -1.02432865]
    }
}

# Temperature range boundaries
T_LOW = 300.0  # K
T_HIGH = 1000.0  # K

def get_nasa_coefficients(species, T):
    """
    Return the appropriate coefficients for the given temperature
    """
    if species not in NASA_COEFFICIENTS:
        raise ValueError(f"Species {species} not found in database")
    
    if T <= T_HIGH:
        return NASA_COEFFICIENTS[species]['low_coeffs']
    else:
        return NASA_COEFFICIENTS[species]['high_coeffs']
    
        
def calculate_cp_nasa(T, coeffs):
    """
    Calculate Cp in kJ/kmol-K using NASA polynomial coefficients
    Cp/R = a1 + a2*T + a3*T² + a4*T³ + a5*T⁴
    """
    R = 8.314462618  # kJ/kmol-K (universal gas constant)
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    cp_over_R = a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4
    return cp_over_R * R

def calculate_enthalpy_nasa(T, coeffs):
    """
    Calculate enthalpy H(T) in kJ/kmol using NASA polynomial coefficients
    H/RT = a1 + a2/2*T + a3/3*T² + a4/4*T³ + a5/5*T⁴ + a6/T
    """
    R = 8.314462618  # kJ/kmol-K
    a1, a2, a3, a4, a5, a6, a7 = coeffs
    H_over_RT = a1 + a2/2*T + a3/3*T**2 + a4/4*T**3 + a5/5*T**4 + a6/T
    return H_over_RT * R * T

def calculate_adiabatic_flame_temperature(phi, T_in, c, tolerance=1e-3, max_iterations=100):
    """
    Calculate adiabatic flame temperature for methane-air mixture with H2 addition
    
    Parameters:
    phi: equivalence ratio
    T_in: inlet temperature (K)
    c: progress variable (0 to 1)
    h2_fraction: fraction of H2 in the fuel mixture (0 to 1)
    tolerance: convergence tolerance
    max_iterations: maximum number of iterations
    """
    # Calculate composition based on progress variable c
    # For fuel mixture: (1 - h2_fraction) CH4 + h2_fraction H2
    n_CH4 = c*(0-0.6)+0.6  # Linear interpolation
    n_H2 = c*(0-0.4)+0.4             # Linear interpolation
    
    # Oxygen requirement: CH4 needs 2O2, H2 needs 0.5O2
    
    n_O2_air =c*(1.433-2.8)+2.8
    n_N2_air = 10.533  # Nitrogen scales with oxygen
    
    # Products composition
    # CO2 from CH4 combustion only
    n_CO2 = c*(0.605-0)+0
    # H2O from both CH4 and H2 combustion
    n_H2O = c*(1.612-0)+0
    n_O2 = max(0, n_O2_air - (2 * n_CH4 + 0.5 * n_H2))  # Excess oxygen
    n_N2 = n_N2_air

    # Enthalpy of formation (kJ/kmol) at 298.15K
    hf_CH4 = -74850
    hf_H2 = 0  # H2 enthalpy of formation is 0 (element)
    hf_O2 = 0
    hf_N2 = 0
    hf_CO2 = -393520
    hf_H2O = -241820

    # Reference temperature
    T_ref = 298.15  # K

    # Initial guess for adiabatic flame temperature
    T_guess = 1800  # K

    for iteration in range(max_iterations):
        # Get coefficients for current temperature guess (products)
        coeffs_CO2 = get_nasa_coefficients('CO2', T_guess)
        coeffs_H2O = get_nasa_coefficients('H2O', T_guess)
        coeffs_O2_prod = get_nasa_coefficients('O2', T_guess)
        coeffs_N2_prod = get_nasa_coefficients('N2', T_guess)
        
        # Get coefficients for inlet temperature (reactants)
        coeffs_CH4_in = get_nasa_coefficients('CH4', T_in)
        coeffs_H2_in = get_nasa_coefficients('H2', T_in)
        coeffs_O2_in = get_nasa_coefficients('O2', T_in)
        coeffs_N2_in = get_nasa_coefficients('N2', T_in)
        
        # Get coefficients for reference temperature
        coeffs_CH4_ref = get_nasa_coefficients('CH4', T_ref)
        coeffs_H2_ref = get_nasa_coefficients('H2', T_ref)
        coeffs_O2_ref = get_nasa_coefficients('O2', T_ref)
        coeffs_N2_ref = get_nasa_coefficients('N2', T_ref)
        coeffs_CO2_ref = get_nasa_coefficients('CO2', T_ref)
        coeffs_H2O_ref = get_nasa_coefficients('H2O', T_ref)

        # 1. Calculate Enthalpy of Reactants at inlet temperature
        H_reactants = 0
        H_reactants += n_CH4 * (hf_CH4 + (calculate_enthalpy_nasa(T_in, coeffs_CH4_in) - calculate_enthalpy_nasa(T_ref, coeffs_CH4_ref)))
        H_reactants += n_H2 * (hf_H2 + (calculate_enthalpy_nasa(T_in, coeffs_H2_in) - calculate_enthalpy_nasa(T_ref, coeffs_H2_ref)))
        H_reactants += n_O2_air * (hf_O2 + (calculate_enthalpy_nasa(T_in, coeffs_O2_in) - calculate_enthalpy_nasa(T_ref, coeffs_O2_ref)))
        H_reactants += n_N2_air * (hf_N2 + (calculate_enthalpy_nasa(T_in, coeffs_N2_in) - calculate_enthalpy_nasa(T_ref, coeffs_N2_ref)))

        # 2. Calculate Enthalpy of Products at guessed temperature
        H_products = 0
        H_products += n_CO2 * (hf_CO2 + (calculate_enthalpy_nasa(T_guess, coeffs_CO2) - calculate_enthalpy_nasa(T_ref, coeffs_CO2_ref)))
        H_products += n_H2O * (hf_H2O + (calculate_enthalpy_nasa(T_guess, coeffs_H2O) - calculate_enthalpy_nasa(T_ref, coeffs_H2O_ref)))
        H_products += n_O2 * (hf_O2 + (calculate_enthalpy_nasa(T_guess, coeffs_O2_prod) - calculate_enthalpy_nasa(T_ref, coeffs_O2_ref)))
        H_products += n_N2 * (hf_N2 + (calculate_enthalpy_nasa(T_guess, coeffs_N2_prod) - calculate_enthalpy_nasa(T_ref, coeffs_N2_ref)))

        # 3. Energy Balance Residual: H_reactants - H_products = 0
        residual = H_reactants - H_products

        # 4. Calculate total heat capacity of the product mixture at T_guess
        cp_mix = (n_CO2 * calculate_cp_nasa(T_guess, coeffs_CO2) +
                 n_H2O * calculate_cp_nasa(T_guess, coeffs_H2O) +
                 n_O2 * calculate_cp_nasa(T_guess, coeffs_O2_prod) +
                 n_N2 * calculate_cp_nasa(T_guess, coeffs_N2_prod))

        # 5. Newton-Raphson update: T_new = T_guess + residual / cp_mix
        T_new = T_guess + residual / cp_mix

        # 6. Check for convergence
        error = abs(T_new - T_guess)

        if error < tolerance:
            return T_new

        T_guess = T_new

    return T_guess      

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
    C = 0.0 
    # location where the flame should be
    flame_x = 0.012
    if (x < flame_x):
      C = 0.0 
    else:
      C = 1.0

    return C   

# ################################################################## #
# loop over all vertices and set the species progress variable       # 
# ################################################################## #
def SetInitialSpecies(SU2Driver):
    allCoords = SU2Driver.Coordinates()
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    for iPoint in range(SU2Driver.GetNumberNodes() - SU2Driver.GetNumberHaloNodes()):
      coord = allCoords.Get(iPoint)
      C = initC(coord)
      # now update the initial condition for the species
      SU2Driver.Solution(iSPECIESSOLVER).Set(iPoint,0,C)

def update_temperature(SU2Driver, iPoint):
    # first, get the progress variable
    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    # Note: returns a list
    C = SU2Driver.Solution(iSPECIESSOLVER)(iPoint,0)
    h= SU2Driver.Solution(iSPECIESSOLVER)(iPoint,1)
    #print("enthalpy",h)
    
    Tad=calculate_adiabatic_flame_temperature(phi,Tu,C)
    #print("Adiabatic temp=",Tad)
    y=-86.48*(C**4)+321.51*(C**3)-706.92*(C**2)+548.63*C-38.021
    Tad_f=Tad - y
    T=(h-304721)/(1161) + Tad_f

    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    # the list with names
    solindex = getsolvar(SU2Driver)
    iTEMP = solindex.get("TEMPERATURE")
    #print("Adiabatic temp=",T)
    SU2Driver.Solution(iFLOWSOLVER).Set(iPoint,iTEMP,T)



def zimont(SU2Driver, iPoint):

    iSSTSOLVER = SU2Driver.GetSolverIndices()['SST']
    tke, dissipation = SU2Driver.Solution(iSSTSOLVER)(iPoint)

    iSPECIESSOLVER = SU2Driver.GetSolverIndices()['SPECIES']
    # get the gradient of species_0
    gradc = SU2Driver.Gradient(iSPECIESSOLVER)(iPoint,0)
    primindex = SU2Driver.GetPrimitiveIndices()
    iDENSITY = primindex.get("DENSITY")
    iMU = primindex.get("LAMINAR_VISCOSITY")
    
    h= SU2Driver.Solution(iSPECIESSOLVER)(iPoint,1)
    h_scaled=h//1000000
    # laminar burning velocity of methane-air at phi=0.5, P=5

    #Slu = 0.33745

    rho = SU2Driver.Primitives()(iPoint,iDENSITY)
    mu = SU2Driver.Primitives()(iPoint,iMU)
    nu=mu/rho
    # Turbulent Flamespeed Closure with Dinkelacker correction
    up = np.sqrt((2.0/3.0) * tke )
    lt = (0.09**0.75) * (tke**1.5) / dissipation
    Re = up*lt/nu
    Le = 0.498
    if h_scaled>=-0.3:
      Slu= -1.6508*(h_scaled**5)+1.3050*(h_scaled**4)+1.7592*(h_scaled**3)+1.1562*(h_scaled**2)+0.3947*(h_scaled)+0.0545
      Ut = Slu * (1.0 + (0.46/Le)*np.power(Re,0.25)*np.power(up/Slu,0.3)*np.power(Pu,0.2))
    else :
      Ut = 0
    
    norm_gradc = np.sqrt(gradc[0]*gradc[0] + gradc[1]*gradc[1])
    Sc = rho_u * Ut * norm_gradc

    return Sc

def getsolvar(SU2Driver):
    primindex = SU2Driver.GetPrimitiveIndices()
    iFLOWSOLVER = SU2Driver.GetSolverIndices()['INC.FLOW']
    nVars = SU2Driver.Solution(iFLOWSOLVER).Shape()[1]
    varindex = primindex.copy()
    for prim in varindex.copy():
      if varindex[prim] >=nVars:
        del varindex[prim]
    varindex = dict(sorted(varindex.items(), key=lambda item: item[1]))
    return varindex

def ApplyScalar(driver,marker_ids):
    for marker_id in marker_ids:
      if marker_id < 0:
        continue

      #h1=[]
      #h2=[]
      #print('Marker id',marker_id)
      for i_vertex in range(driver.GetNumberMarkerNodes(marker_id)):
        iSPECIESSOLVER = driver.GetSolverIndices()['SPECIES']
        C=driver.Solution(iSPECIESSOLVER)(i_vertex,0)
        Tad=calculate_adiabatic_flame_temperature(phi,Tu,C)
        y=-86.48*(C**4)+321.51*(C**3)-706.92*(C**2)+548.63*C-38.021   #Correction factor
        Tad_f=Tad - y
        h1=1
        h2=304721+1161*(Twall-Tad_f)
        hf=[h1,h2]
        driver.SetMarkerCustomScalar(marker_id, i_vertex, hf)

def main():
  """
  Run the flow solver with a custom inlet (function of time and space).
  """
  comm = MPI.COMM_WORLD
  #comm = 0
  rank = comm.Get_rank()
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
  iSPECIESSOLVER = driver.GetSolverIndices()['SPECIES']
  iSSTSOLVER = driver.GetSolverIndices()['SST']
  # all the indices and the map to the names of the primitives
  primindex = driver.GetPrimitiveIndices()
  nElem = driver.GetNumberElements()
  nVars = driver.Solution(iFLOWSOLVER).Shape()[1]
  nVarsSpecies = driver.Solution(iSPECIESSOLVER).Shape()[1]
  nVarsTurb = driver.Solution(iSSTSOLVER).Shape()[1]

  solindex = getsolvar(driver)
  iTEMP = solindex.get("TEMPERATURE")
 
  if rank == 0:
    print("Dimensions of the problem = ",nDim)
    print("index of flow solver = ",iFLOWSOLVER)
    print("index of turbulence solver = ",iSSTSOLVER)
    print("indices of primitives=",primindex)
    print("number of primitives:",len(primindex))
    print("number of elements:",nElem)
    print("number of flow solver variables:",nVars)
    print("number of species solver variables:",nVarsSpecies)
    print("number of turbulence solver variables:",nVarsTurb)
    sys.stdout.flush()


  with open('psi.cfg') as f:
    if 'RESTART_SOL= YES' in f.read():
      if rank == 0:
        print("restarting from file")
    else:
        # We can set an initial condition by calling this function:
        if rank == 0:
          print("Using user defined initial condition.")
        SetInitialSpecies(driver) 

  sys.stdout.flush()
  

  all_marker_ids = driver.GetMarkerIndices()
  marker_names = ['wall_top']
  #print('Marker id trial',all_marker_ids)
  marker_ids = []
  for name in marker_names:
    marker_ids.append(all_marker_ids[name] if name in all_marker_ids else -1)

  # run 500 iterations
  for inner_iter in range(2000):

    if (rank==0):
      print("python iteration ", inner_iter)
    #ApplyScalar(driver,marker_ids)
    driver.Preprocess(inner_iter)
    driver.Run()

    Source = driver.UserDefinedSource(iSPECIESSOLVER)

    # set the source term, per point
    for i_node in range(driver.GetNumberNodes() - driver.GetNumberHaloNodes()):
      # add source term:
      # default TFC of Zimont: rho*Sc = rho_u * U_t * grad(c)
      S = zimont(driver, i_node)
      Source.Set(i_node, 0, S)

    # for the update of temperature, we need to update also the halo nodes
    #for i_node in range(driver.GetNumberNodes()):
    #  # set the temperature to T = c*Tf + (1-c)*Tu
      update_temperature(driver, i_node)

    driver.Postprocess()
    driver.Update()
    # Monitor the solver and output solution to file if required.
    #driver.Monitor(inner_iter)
    # Output the solution to file
    driver.Output(inner_iter)

  # Finalize the solver and exit cleanly.
  driver.Finalize()

if __name__ == '__main__':
  main()
