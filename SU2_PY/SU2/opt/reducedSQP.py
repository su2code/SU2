#!/usr/bin/env python
# This Python file uses the following encoding: utf-8

## \file ReducedSQP.py
#  \brief Python script for performing the reducedSQP optimization w.
#  \author T. Dick

import sys
import numpy as np
from scipy import optimize
from .project import Project

def reduced_sqp(x0, func, f_eqcons, f_ieqcons, fprime, fprime_eqcons, fprime_ieqcons, fdotdot, project, acc ):
    "This is the implementation of the reduced SQP optimizer for smoothed derivatives"

    # preprocessing before the first optimization run

    # set the inout parameters
    p = x0
    nu = 1.0
    err = 1.0
    step = 1

    # main loop
    while ( err > acc ):

        sys.stdout.write('Optimizer iteration: ' + str(step) + ' current err: ' + str(err) + '\n')

        # evaluate the function
        F = func(p, project)
        E = f_eqcons(p, project)
        D_F = fprime(p, project)
        D_E = fprime_eqcons(p, project)
        H_F = fdotdot(p, project)

        sys.stdout.write('   objective function: ' + str(F) + ' , constrain: ' + str(E) + '\n')

        # assemble NLES
        Jac = sqp_jacobian(H_F, D_E)
        Rhs = sqp_rhs(D_F, E)

        # solve the Newton step
        sol = np.linalg.solve(Jac, Rhs)

        #update the design
        p += sol[0:project.n_dv]
        nu = sol[len(sol)-1]
        err = np.linalg.norm(D_F,2)
        step += 1

        sys.stdout.write('   current design: ' + str(p) + ' , Lagrange multiplier: ' + str(nu) + '\n')

    return 0


def sqp_jacobian(H_F, D_E):
    """ This function assembles the Jacobian for the Newton step in the reduced SQP method.
        input: objective function Hessian approximation, equality constrain gradient.
        output: the LHS matrix for the Newton step
    """
    Jac = np.block([ [H_F, np.transpose(D_E)], [D_E, 0] ])
    return Jac


def sqp_rhs(D_F, E):
    """ This function assembles the right hand side for the Newton step in the reduced SQP method.
        input: objective function gradient, equality constrain.
        output: the RHS vector for the Newton step
    """
    Rhs = np.block([ -D_F, -E ])
    return Rhs


"""
left over code snippets
    F = func(x0, project)
    E = f_eqcons(x0, project)
    C = f_ieqcons(x0, project)
    D_F = fprime(x0, project)
    D_E = fprime_eqcons(x0, project)
    D_C = fprime_ieqcons(x0, project)
    H_F = fdotdot(x0, project)
    Jac = sqp_jacobian(H_F,D_E, project)
    stop = 1

"""
