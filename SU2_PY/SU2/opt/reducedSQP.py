#!/usr/bin/env python
# This Python file uses the following encoding: utf-8

## \file ReducedSQP.py
#  \brief Python script for performing the reducedSQP optimization w.
#  \author T. Dick

import sys
import numpy as np
from scipy import optimize
from .. import eval as su2eval
from .project import Project
import cvxopt

global glob_project

def reduced_sqp(x0, func, f_eqcons, f_ieqcons, fprime, fprime_eqcons, fprime_ieqcons, fdotdot, project, iter, acc, xb ):
    "This is the implementation of the reduced SQP optimizer for smoothed derivatives"

    SCIPY=False
    HANDMADE=True

    # use the SciPy optimizer
    if SCIPY:

        # preprocessing before the optimization run
        opt = { 'maxiter' : iter, 'verbose' : 3 }

        #pack the constaints in form for trust-constr
        equal   = optimize.NonlinearConstraint( feq, 0.0 , 0.0)
        inequal = optimize.NonlinearConstraint( fieq, 0.0 , np.inf)

        # some global function voodoo
        global glob_project
        glob_project = project

        # call the optimizer from Scipy minimize
        outputs = optimize.minimize( fun          = f                 ,
                                     x0           = x0                ,
                                     method       = 'trust-constr'    ,
                                     jac          = df            ,
                                     hess         = ddf           ,
                                     constraints  = [equal, inequal]  ,
                                     bounds       = xb                ,
                                     tol          = acc               ,
                                     options      = opt               )


    # use the self implemented SQP optimizer
    elif HANDMADE:

        sys.stdout.write('Using the hand implemented version. Setting up the environment...' + '\n')

        # evaluate the functions
        p = x0
        F = func(p, project)
        E = f_eqcons(p, project)
        C = f_ieqcons(p,project)
        D_F = fprime(p, project)
        D_E = fprime_eqcons(p, project)
        D_C = fprime_ieqcons(p,project)
        H_F = fdotdot(p, project)

        # set the inout parameters
        nu1 = [1.0]*len(E)
        nu2 = [1.0]*len(C)
        err = np.linalg.norm(D_F,2)
        step = 1

        # main loop
        while ( err > acc and step <= iter ):

            sys.stdout.write('Optimizer iteration: ' + str(step) + ' current err: ' + str(err) + '\n')

            if step>1:
                # reevaluate the functions
                F = func(p, project)
                E = f_eqcons(p, project)
                C = f_ieqcons(p,project)
                D_F = fprime(p, project)
                D_E = fprime_eqcons(p, project)
                D_C = fprime_ieqcons(p,project)
                H_F = fdotdot(p, project)

            sys.stdout.write('   objective function: ' + str(F) + ' , equality constrain: ' + str(E) + ' , inequality constrain: ' + str(C) + '\n')

            # assemble NLES
            #Jac = np.block([ [H_F, np.transpose(D_E)], [D_E, 0] ])
            #Rhs = np.block([ -D_F, -E ])

            #expand inequality constraints by bounds
            ub       = float ( project.config.OPT_BOUND_UPPER )
            lb       = float ( project.config.OPT_BOUND_LOWER )
            Id = np.identity(len(p))
            if len(C)>0:
                G = cvxopt.matrix(np.block([ [-D_C], [Id], [-Id]]))
                h = cvxopt.matrix(np.append( C, np.append([ub]*len(p), [lb]*len(p))))
            else:
                G = np.block([ [Id], [-Id]])
                h = np.append([ub]*len(p), [lb]*len(p))

            # for debugging
            P = cvxopt.matrix(H_F)
            q = cvxopt.matrix(D_F)


            # solve the interior quadratic problem
            sol = cvxopt.solvers.qp(P,
                                    q,
                                    G,
                                    h)

            #update the design
            delta_p = np.transpose(np.array(sol['x']))
            p += delta_p
            p = p[0]
            err = np.linalg.norm(delta_p,2)
            step += 1

            sys.stdout.write('   current design: ' + str(p) + ' , Lagrange multiplier: ' + str(0) + '\n')

        return 0

    elif DESCEND:

        sys.stdout.write('Using gradient descend \n')
        # evaluate the functions
        p = x0
        F = func(p, project)
        E = f_eqcons(p, project)
        C = f_ieqcons(p,project)
        D_F = fprime(p, project)
        D_E = fprime_eqcons(p, project)
        D_C = fprime_ieqcons(p,project)
        H_F = fdotdot(p, project)




    #return the results
    return outputs

#
# end of treduced_sqp
#


def lagrangian_interface(x, project):
    """ This function evaluates the Lagrangian of the problem
        This takes the form of func + multiplier*f_eqcons
        We assume the flow solver to be converged
    """
    p = x[0:len(x)-1]
    nu = x[len(x)-1]
    F = func(p, project)
    E = f_eqcons(p,project)
    return F + nu*E


def hessian_interface(x, project):
    """ This function assembles the Jacobian for the Newton step in the reduced SQP method.
        input: objective function Hessian approximation, equality constrain gradient.
        output: the LHS matrix for the Newton step
    """
    p = x[0:len(x)-1]
    H_F = fdotdot(p, project)
    D_E = fprime_eqcons(p. project)
    Jac = np.block([ [H_F, np.transpose(D_E)], [D_E, 0] ])
    return Jac


def jacobian_interface(x, project):
    """ This function assembles the right hand side for the Newton step in the reduced SQP method.
        input: objective function gradient, equality constrain.
        output: the RHS vector for the Newton step
    """
    p = x[0:len(x)-1]
    D_F = fprime(p, project)
    E = f_eqcons(p, project)
    Rhs = np.block([ -D_F, -E ])
    return Rhs




# this is a stupid idea
def f(x):
    global glob_project
    obj_list = glob_project.obj_f(x)
    obj = 0
    for this_obj in obj_list:
        obj = obj+this_obj

    return obj

def df(x):
    global glob_project
    dobj_list = glob_project.obj_df(x)
    dobj=[0.0]*len(dobj_list[0])

    for this_dobj in dobj_list:
        idv=0
        for this_dv_dobj in this_dobj:
            dobj[idv] = dobj[idv]+this_dv_dobj;
            idv+=1
    dobj = np.array( dobj )

    return dobj

def ddf(x):
    global glob_project
    dobj_list = glob_project.obj_ddf(x)
    dobj=[0.0]*len(dobj_list[0])

    for this_dobj in dobj_list:
        idv=0
        for this_dv_dobj in this_dobj:
            dobj[idv] = dobj[idv]+this_dv_dobj;
            idv+=1
    dobj = np.array( dobj )

    return dobj


def feq(x):
    global glob_project
    cons = glob_project.con_ceq(x)

    if cons: cons = np.array(cons)
    else:    cons = np.zeros([0])

    return cons

def fieq(x):

    global glob_project
    cons = glob_project.con_cieq(x)

    if cons: cons = np.array(cons)
    else:    cons = np.zeros([0])

    return -cons

"""
left over code snippets

"""
