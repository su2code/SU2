#!/usr/bin/env python 
## \file downstream_function.py


import os, sys, shutil, copy
import numpy as np
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

def downstream_function(config, state ):
  nvar = 5
  for iDV in range(len(config.DV_KIND)):
    if config.DV_KIND[iDV] == 'OTHER':
      nvar = nvar+1
  d_in = [0.0]*nvar
  obj = objective(config, state,d_in)
  print " objective ", obj
  return obj

def objective(config, state, d_in ):
  # Values from config object
  P0=float(config.FREESTREAM_PRESSURE)
  T0=float(config.FREESTREAM_TEMPERATURE)

  # Values from history file
  rho3 = state['HISTORY']['DIRECT']['FLUXAVG_OUTLET_DENSITY'][-1]+d_in[0];
  M3 = float(state['HISTORY']['DIRECT']['AVG_OUTLET_MACH'][-1])+d_in[1]
  P3 = state['HISTORY']['DIRECT']['AVG_OUTLET_PRESSURE'][-1]+d_in[4];
  
  # CUSTOM DV: Always set to their initial values within this file as well as within the config file
  CustomDV = 3.0
  CustomDV2 = 3.0
  FFD_occured = 0

  for iDV in range(len(config.DV_KIND)):
    # This section ensures that the custom design variable is read in
    #if (config.MATH_PROBLEM == 'CONTINUOUS_ADJOINT'):
    #  if (config.DV_KIND[iDV] == 'FFD_CONTROL_POINT_2D') or (config.DV_KIND[iDV] == 'FFD_CONTROL_POINT'):
    #    FFD_occured = 1;
    if config.DV_KIND[iDV] == 'CUSTOM':
      CustomDV = config.DV_PARAM['PARAM'][iDV][0]
      if( len(config.DV_PARAM['FFDTAG'][iDV])>0):
        CustomDV = float(config.DV_PARAM['FFDTAG'][iDV][0])
      if len(config.DV_VALUE)>=len(config.DV_KIND):
        CustomDV = CustomDV+config.DV_VALUE[iDV]
  

  if (len(d_in)>5):
    CustomDV = CustomDV + d_in[5]
    CustomDV2 = CustomDV2 + d_in[6]

  a0 = np.sqrt(T0*1.4*287.87)
  # Here is a random function to act a placeholder.
  obj_val = (P3/P0)+(M3-rho3)+CustomDV-rho3/T0+rho3*P3-CustomDV2/2

  return obj_val

def downstream_gradient(config,state,step=1e-8):

  nvar = 7

  if type(step)==list:
    step = step[0]
  d_in = [0.0]*nvar

  J0 = objective(config,state,d_in)
  d_in[4]=step
  dJdP = (objective(config,state,d_in)-J0)/step
  d_in[4]=0.0
  d_in[0]=step
  dJdrho = (objective(config,state,d_in)-J0)/step
  d_in[0]=0.0
  d_in[1]=step
  dJdu = (objective(config,state,d_in)-J0)/step
  d_in[1]=0.0
  d_in[5]=step
  dJd1 = (objective(config,state,d_in)-J0)/step
  d_in[5]=0.0
  d_in[6]=step
  dJd2 = (objective(config,state,d_in)-J0)/step
  gradient = (dJdrho,dJdu,0.0,0.0,dJdP,dJd1, dJd2)

  return gradient
