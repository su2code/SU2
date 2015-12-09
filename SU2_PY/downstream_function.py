#!/usr/bin/env python 
## \file downstream_function.py


import os, sys, shutil, copy
import numpy as np
from optparse import OptionParser
sys.path.append(os.environ['SU2_RUN'])
import SU2

def downstream_function(config, state ):

  return objective(config, state,(0.0,0.0,0.0,0.0,0.0))

def objective(config, state, d_in ):
  # Values from config object
  gamma=1.4
  R = 287.15

  P0=float(config.FREESTREAM_PRESSURE)
  T0=float(config.FREESTREAM_TEMPERATURE)
  M0=float(config.MACH_NUMBER)
  a0=np.sqrt(gamma*R*T0)
  rho0=P0/T0/R
  V0 = M0*a0;

  # Values from history file at the monitored outlet
  T1 = float(state['HISTORY']['DIRECT']['AVG_OUTLET_TEMPERATURE'][-1])
  P1 = state['HISTORY']['DIRECT']['AVG_OUTLET_PRESSURE'][-1]+d_in[4]
  M1 = float(state['HISTORY']['DIRECT']['AVG_OUTLET_MACH'][-1])
  rho1=P1/T1/R+d_in[0]
  V1 = M1/np.sqrt(gamma*R*T1)+d_in[1]
  
  obj_val = P1+0.5*V1*V1*rho1;
  
  return obj_val

def downstream_gradient(config,state,step=1e-4):
  print config.keys()
  J0 = objective(config,state,(0.0,0.0,0.0,0.0,0.0))
  dJdP = (objective(config,state,(0.0,0.0,0.0,0.0,step))-J0)/step
  dJdrho = (objective(config,state,(step,0.0,0.0,0.0,0.0))-J0)/step
  dJdu = (objective(config,state,(0.0,step,0.0,0.0,0.0))-J0)/step
  print dJdrho, dJdu, dJdP
  gradient = (dJdrho,dJdu,0.0,0.0,dJdP)

  return gradient
