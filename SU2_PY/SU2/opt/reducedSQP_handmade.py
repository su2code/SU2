#!/usr/bin/env python
# This Python file uses the following encoding: utf-8

## \file reducedSQP_handmade.py
#  \brief Python script for performing the SQP optimization with a semi handwritten optimizer.
#  \author T. Dick

import sys
import numpy as np
import cvxopt


def SQPhandimplementation(x0, func, f_eqcons, f_ieqcons, fprime, fprime_eqcons, fprime_ieqcons, fdotdot, parameter, iter, acc, xb=None):
    """ This is a implementation of a SQP optimizer
        It is written for smoothed derivatives
        Accepts:
            - equality and inequality constraints
            - approximated hessians
            - additional parameters
        can be applied to analytical functions for test purposes"""

    sys.stdout.write('Using the hand implemented version. Setting up the environment...' + '\n')

    # start by evaluating the the functions
    p = x0
    F = func(p, parameter)
    E = f_eqcons(p, parameter)
    C = f_ieqcons(p, parameter)
    D_F = fprime(p, parameter)
    D_E = fprime_eqcons(p, parameter)
    D_C = fprime_ieqcons(p, parameter)
    H_F = fdotdot(p, parameter)

    # set the inout parameters
    # nu1 = [1.0]*len(E)
    # nu2 = [1.0]*len(C)
    err = 2*acc+1
    step = 1

    # main optimizer loop
    while (err > acc and step <= iter):

        sys.stdout.write('Optimizer iteration: ' + str(step) +
                         ' current err: ' + str(err) + '\n')

        if step > 1:
            # reevaluate the functions
            F = func(p, parameter)
            E = f_eqcons(p, parameter)
            C = f_ieqcons(p, parameter)
            D_F = fprime(p, parameter)
            D_E = fprime_eqcons(p, parameter)
            D_C = fprime_ieqcons(p, parameter)
            H_F = fdotdot(p, parameter)

        sys.stdout.write('   objective function: ' + str(F) +
                         ' , equality constrain: ' + str(E) +
                         ' , inequality constrain: ' + str(C) + '\n')

        # assemble equality constraints
        if np.size(E) > 0:
            A = cvxopt.matrix(D_E)
            b = cvxopt.matrix(E)

        # expand inequality constraints by bounds
        if xb is None:
            xb = [-1e-1, 1e-1]
        Id = np.identity(len(p))
        if np.size(C) > 0:
            G = cvxopt.matrix(np.block([[D_C], [-Id], [Id]]))
            h = cvxopt.matrix(np.append(-C, np.append([-xb[0]]*len(p), [xb[1]]*len(p))))
        else:
            G = cvxopt.matrix(np.block([[-Id], [Id]]))
            h = cvxopt.matrix(np.append([-xb[0]]*len(p), [xb[1]]*len(p)))

        # pack objective function
        P = cvxopt.matrix(H_F)
        q = cvxopt.matrix(D_F)

        # solve the interior quadratic problem
        if np.size(E) > 0:
            sol = cvxopt.solvers.qp(P, q, G, h, A, b)
        else:
            sol = cvxopt.solvers.qp(P, q, G, h)

        # line search
        criteria = True
        delta_p = np.transpose(np.array(sol['x']))
        while (criteria):
            p_temp = p + delta_p
            p_temp = p_temp[0]
            sys.stdout.write("descend step: " + str(delta_p) + "\n")
            F_new = func(p_temp, parameter)
            # test for optimization progress
            if (F_new > F):
                sys.stdout.write("not a reduction in objective function. \n")
                delta_p = 0.5*delta_p
            elif ((np.size(C) > 0) & np.any(C >= 0)):
                sys.stdout.write("constrained not fullfilled. \n")
                delta_p = 0.5*delta_p
            else:
                sys.stdout.write("descend step accepted. \n")
                criteria = False

            if (np.linalg.norm(delta_p, 2) < 10*acc):
                sys.stdout.write("can't find a good step. \n")
                criteria = False

        # update the design
        p = p_temp
        err = np.linalg.norm(delta_p, 2)
        step += 1

        sys.stdout.write('   current design: ' + str(p) + '\n')

    return 0

# end SQPhandimplementation
