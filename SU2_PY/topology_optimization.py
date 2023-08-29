#!/usr/bin/env python

## \file topology_optimization.py
#  \brief Python script to drive SU2 in topology optimization.
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
#
#########################################################################################
#                                                                                       #
# This script is provided to show how the feature can be used.                          #
# It is not meant to be generic, topology optimization calls for some parameters to be  #
# ramped and the strategy to do so is hard coded.                                       #
# A bit of hacking will be required if you want to deviate from what is done here.      #
# The hard coded bits are explained as they appear.                                     #
#                                                                                       #
#########################################################################################

import os
import sys
import math
import time
import shutil
import subprocess as sp
import numpy as np
import scipy.optimize

####### SETUP #######

obj_scale = 1 / 1.25e-3  # scale the objective so that it starts at 1-4
con_scale = 1 / 0.5  # 1 over upper bound (e.g. max volume)
var_scale = 1.0  # variable scale

# maximum number of iterations
maxJev_t = 1000
# max iters for gray initialization, i.e. soft filter settings
maxJev_i = 200
# num iters between updates of the filter settings and constraint penalty factor
nJev_u = 40

# tolerances
ftol_u = 1e-5  # during updates
ftol_f = 1e-7  # final iteration
# the exterior penalty method is used to impose the constraint,
# this is the maximum constraint violation, below it the penalty factor is not increased
htol = 5e-3

# general options for L-BFGS-B
options = {"disp": True, "maxcor": 10, "ftol": ftol_u, "gtol": 1e-18}

# these are the commands for the direct and adjoint runs, modify to run parallel
commands = ["SU2_CFD ", "SU2_CFD_AD "]

# file through which SU2 gets the design densities
inputFile = "element_properties.dat"

# names of the output files [objective value, objective gradient, constraint value, ...]
outputFiles = ["grad_compliance.dat", "grad_vol_frac.dat"]

# settings for direct run and adjoint of the objective and constraint
fnames = ["settings.cfg", "settings_compliance.cfg", "settings_volfrac.cfg"]

# use the DILATE, ERODE, (DILATE,ERODE), or (ERODE,DILATE) filters, the first value is
# for gray initialization, then it is ramped until a solid-void topology is obtained
filterParam = [0.01, 1, 4, 16, 64, 200]


####### SU2 Driver #######


class Driver:
    def __init__(self, commands, inputFile, configFiles, outputFiles):
        self._inputFile = inputFile
        self._objValFile = "history.csv"
        self._objDerFile = outputFiles[0]
        self._conValFile = "history.csv"
        self._conDerFile = outputFiles[1]
        self._objValCommand = commands[0] + configFiles[0] + " > objval.stdout"
        self._objDerCommand = commands[1] + configFiles[1] + " > objder.stdout"
        self._conDerCommand = commands[1] + configFiles[2] + " > conval.stdout"

    # end

    def _assert_isfinite(self, val):
        if math.isinf(val) or math.isnan(val):
            raise ValueError

    # end

    def _write_input(self, x):
        fid = open(self._inputFile, "w")
        lines = ["\n"]
        for val in x:
            lines.append("0  0  0  0  0  " + str(val / var_scale) + "\n")
        # end
        fid.writelines(lines)
        fid.close()

    # end

    def obj_val(self, x):
        # write inputs
        self._write_input(x)

        # clear previous output and run direct solver
        try:
            os.remove(self._objValFile)
        except:
            pass

        try:
            sp.call(self._objValCommand, shell=True)
            with open(self._objValFile, "r") as fid:
                lines = fid.readlines()
            for col, name in enumerate(lines[0].split(",")):
                if "TopComp" in name:
                    val = float(lines[1].split(",")[col])
                    break
            # the return code of mpirun is useless, we test the value of the function
            self._assert_isfinite(val)
        except:
            raise RuntimeError("Objective function evaluation failed")
        # end

        return val * obj_scale

    # end

    def obj_der(self, x):
        # inputs written in obj_val_driver

        # clear previous output and run direct solver
        try:
            os.remove(self._objDerFile)
        except:
            pass
        N = x.shape[0]
        y = np.ndarray((N,))

        try:
            # main command
            sp.call(self._objDerCommand, shell=True)

            fid = open(self._objDerFile, "r")
            lines = fid.readlines()
            fid.close()
            for i in range(N):
                val = float(lines[i][0:-1])
                self._assert_isfinite(val)
                y[i] = val * obj_scale / var_scale
            # end
        except:
            raise RuntimeError("Objective gradient evaluation failed")
        # end

        return y

    # end

    def con_val(self, x):
        # inputs written in obj_val_driver

        try:
            with open(self._conValFile, "r") as fid:
                lines = fid.readlines()
            for col, name in enumerate(lines[0].split(",")):
                if "VolFrac" in name:
                    val = float(lines[1].split(",")[col])
                    break
            self._assert_isfinite(val)
        except:
            raise RuntimeError("Constraint function evaluation failed")
        # end

        return val * con_scale - 1

    # end

    def con_der(self, x):
        # inputs written in obj_val_driver

        # clear previous output and run solver
        try:
            os.remove(self._conDerFile)
        except:
            pass
        N = x.shape[0]
        y = np.ndarray((N,))

        # read result
        try:
            sp.call(self._conDerCommand, shell=True)

            fid = open(self._conDerFile, "r")
            lines = fid.readlines()
            fid.close()
            for i in range(N):
                val = float(lines[i][0:-1])
                self._assert_isfinite(val)
                y[i] = val * con_scale / var_scale
            # end
        except:
            raise RuntimeError("Constraint function evaluation failed")
        # end

        return y

    # end


# end


####### Helpers #######

# updates the parameters in the config files
def update_settings(fnames, params):
    for fname in fnames:
        fid = open(fname, "r")
        lines = fid.readlines()
        fid.close()

        for param in params:
            for i in range(len(lines)):
                if lines[i].startswith(param.name()):
                    lines[i] = param.name() + "= " + repr(param.value()) + "\n"
                    break
                # end
            # end
        # end

        fid = open(fname, "w")
        fid.writelines(lines)
        fid.close()
    # end


# end


# use a list as a function
class ValueList:
    def __init__(self, values):
        self._values = values
        self._ub = len(values) - 1

    def val(self, idx):
        return self._values[min(idx, self._ub)]


# end


# helper class to hold parameters that are ramped
class IncrParam:
    def __init__(self, name, init, incr, maxi, func=None):
        self._name = name
        self._init = init
        self._incr = incr
        self._maxi = maxi
        self._func = func
        self._value = 0
        self.reset()

    # end

    def name(self):
        return self._name

    def reset(self):
        self._value = self._init

    def update(self):
        self._value = min(self._value + self._incr, self._maxi)

    def finished(self):
        return self._value == self._maxi

    def value(self):
        if self._func == None:
            return self._value
        else:
            return self._func(self._value)
        # end

    # end


# end


# Exterior penalty method wrapper
class ExteriorPenaltyMethod:
    def __init__(self, driver, r0=8, rmax=1024, c=2):
        self._driver = driver
        self._r = r0
        self._c = c
        self._rmax = rmax
        self._fval = 0
        self._hval = 0
        # timers
        self._funTime = 0
        self._jacTime = 0

    # end

    def fun(self, x):
        self._funTime -= time.time()
        f = self._driver.obj_val(x)
        h = self._driver.con_val(x)
        self._funTime += time.time()
        self._fval = f
        self._hval = h
        return f + self._r * max(0.0, h) * h

    # end

    def jac(self, x):
        self._jacTime -= time.time()
        df = self._driver.obj_der(x)
        dh = self._driver.con_der(x)
        self._jacTime += time.time()

        # log current values of f and h
        hisfile.write(repr(self._fval) + "  " + repr(self._hval) + "\n")
        hisfile.flush()

        return df + 2 * self._r * max(0.0, self._hval) * dh

    # end

    def update(self):
        self._r = min(self._r * self._c, self._rmax)

    # end


# end


####### RUN OPTIMIZATION #######

paramValues = ValueList(filterParam)
params = [
    IncrParam("TOPOL_OPTIM_KERNEL_PARAM", 0, 1, len(filterParam) - 1, paramValues.val)
]

obj = ExteriorPenaltyMethod(Driver(commands, inputFile, fnames, outputFiles))

logfile = open("optimization.log", "w")
hisfile = open("optimization.his", "w")

line = "### Optimization Started ###\n"
print(line)
logfile.write(line + "\n")
logfile.flush()

totTime = -time.time()
nJacEval = 0
nFunEval = 0
itCount = 0

# initial values and bounds
fid = open(inputFile, "r")
N = len(fid.readlines()) - 1
fid.close()
x = np.ones((N,)) * var_scale / con_scale
lb = np.zeros((N,))
ub = np.ones((N,)) * var_scale
bounds = np.array((lb, ub), float).transpose()

## 1st Phase: Run with "gray" filter settings ##
# get the constraint and function within some tolerance
line = "1: Gray filter (initialization)"
print(line)
logfile.write(line + "\n")
logfile.flush()

update_settings(fnames, params)
success = False

while nJacEval < maxJev_i:
    options["maxiter"] = min(nJev_u, maxJev_i - nJacEval)

    optimum = scipy.optimize.minimize(
        obj.fun, x, method="L-BFGS-B", jac=obj.jac, bounds=bounds, options=options
    )
    itCount += 1
    x = optimum.x
    nJacEval += optimum.nit
    nFunEval += optimum.nfev

    line = " Iter {:d}: f= {:f}  h= {:e}  r= {:f}  nfev= {:d}  njev= {:d}".format(
        itCount, obj._fval, obj._hval, obj._r, optimum.nfev, optimum.nit
    )
    print(line)
    logfile.write(line + "\n")
    logfile.flush()

    if obj._hval > htol:  # increase penalty
        obj.update()
    elif optimum.success:  # check convergence
        success = True
        break
    else:  # continue until convergence or maxJev_i
        pass
    # end
# end
tmp = inputFile.split(".")
shutil.copy(inputFile, tmp[0] + "_gray." + tmp[1])

if not (success):
    line = " Initialization did not converge to desired tolerances"
    print(line)
    logfile.write(line + "\n")
    logfile.flush()
# end

## 2nd Phase: Make filter more "black-white" ##
line = "\n2: Black-White filter"
print(line)
logfile.write(line + "\n")
logfile.flush()

options["maxiter"] = nJev_u
finalIter = False

while nJacEval < maxJev_t and not (finalIter):
    finalIter = True
    for i in range(len(params)):
        params[i].update()
        finalIter &= params[i].finished()
    # end
    update_settings(fnames, params)
    if obj._hval > htol:
        obj.update()

    options["ftol"] = (ftol_u, ftol_f)[int(finalIter)]
    options["maxiter"] = max(nJev_u, (maxJev_t - nJacEval) * int(finalIter))

    optimum = scipy.optimize.minimize(
        obj.fun, x, method="L-BFGS-B", jac=obj.jac, bounds=bounds, options=options
    )
    itCount += 1
    x = optimum.x
    nJacEval += optimum.nit
    nFunEval += optimum.nfev

    line = " Iter {:d}: f= {:f}  h= {:e}  r= {:f}  nfev= {:d}  njev= {:d}".format(
        itCount, obj._fval, obj._hval, obj._r, optimum.nfev, optimum.nit
    )
    print(line)
    logfile.write(line + "\n")
    logfile.flush()
# end
tmp = inputFile.split(".")
shutil.copy(inputFile, tmp[0] + "_bw." + tmp[1])

success = finalIter and obj._hval < htol and optimum.success
totTime += time.time()

line = (
    "\n### Optimization Finished ###\n"
    + "Summary: "
    + ("Failure\n", "Success\n")[int(success)]
    + "  fval: {:f}  hval: {:e}\n".format(obj._fval, obj._hval)
    + "Details:\n"
    + "  iter: {:d}  ttot: {:f}s\n".format(itCount, totTime)
    + "  nfev: {:d}  tfev: {:f}s\n".format(nFunEval, obj._funTime)
    + "  njev: {:d}  tjev: {:f}s\n".format(nJacEval, obj._jacTime)
)
print(line)
logfile.write(line + "\n")
logfile.close()
hisfile.close()
