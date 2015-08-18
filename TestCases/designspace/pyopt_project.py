#!/usr/bin/env python 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   PROJECT: NACA 0012 Optimization with pyOpt
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, numpy
import pyOpt
import libSU2, tasks_su2
from tasks_project import Project
import scipy.io

# -------------------------------------------------------------------
#  Setup Project
# -------------------------------------------------------------------    

# filenames
config_filename  = 'config_NACA0012.cfg'  # SU2 config file
design_filename  = 'design_NACA0012.pkl'  # design data (objectives,gradients,variables)
project_filename = 'project_NACA0012.pkl' # project data (a class structure)
design_matname   = os.path.splitext( design_filename )[0] + '.mat'

# load old project
if os.path.exists(project_filename):
    The_Project = libSU2.load_data(project_filename)
    The_Project.folder_self = os.getcwd()        

# start new project
else:
    # new design data
    design_init = { 'VARIABLES'  : [] ,
                    'OBJECTIVES' : {} ,
                    'GRADIENTS'  : {}  }
    libSU2.save_data(design_filename,design_init)
    # start project
    The_Project = Project( config_name = config_filename ,
                           design_name = design_filename  )        
#: if load/start project


# -------------------------------------------------------------------
#  Setup Optimizer
# -------------------------------------------------------------------     

n_DV = len( The_Project.config_current['DEFINITION_DV']['KIND'] )

# Initial guess
x0    = numpy.zeros(n_DV)    

# Bounds
xb    = numpy.array([-0.10,0.10]) * numpy.ones([n_DV,2])

# Max Iterations
its   = 20

# -------------------------------------------------------------------     
# Objective Function 

def objective_function(X,Prj):
        
    nDV = len( Prj.config_current['DEFINITION_DV']['KIND'] )
    X = X[0:nDV]
    print 'X = ' + str(X)
    
    f    = tasks_su2.eval_f(X,Prj)
    ceq  = tasks_su2.eval_ceq(X,Prj)
    cieq = tasks_su2.eval_cieq(X,Prj)
    
    g = numpy.vstack([ ceq , -cieq ])
    
    fail = 0
    
    return f,g,fail

#: def objective_function()

# -------------------------------------------------------------------     
# Gradient Function 

def gradient_function(X,Prj):
    
    nDV = len( Prj.config_current['DEFINITION_DV']['KIND'] )
    X = X[0:nDV]    
    
    df    = tasks_su2.eval_df(X,Prj)
    dceq  = tasks_su2.eval_dceq(X,Prj)
    dcieq = tasks_su2.eval_dcieq(X,Prj)
    
    df = numpy.array([ df ])
    dg = numpy.vstack([ dceq , -dcieq ])
    
    fail = 0
    
    return df,dg,fail

#: def gradient_function()

# function handles
the_Objective = lambda X: objective_function(X,The_Project)
the_Gradient  = lambda X,f,g: gradient_function(X,The_Project)

# -------------------------------------------------------------------     
# Instantiate Optimization Problem 

The_Problem = pyOpt.Optimization('NACA0012 Constrained Drag Minimization',the_Objective)

# variables
Var_Names = The_Project.config_current['DEFINITION_DV']['KIND']
for i,name in enumerate(Var_Names):
    var_name = 'X%i - %s' % (i,name)
    The_Problem.addVar(var_name,'c',lower=xb[i,0],upper=xb[i,1],value=x0[i])

# objective
Obj_Name = The_Project.config_current['OPT_OBJFUNC']['OBJECTIVE']
The_Problem.addObj(Obj_Name)

# equality constraint
Con_Names = The_Project.config_current['OPT_CONSTR']['EQUALITY'].keys()
for i,name in enumerate(Con_Names):
    The_Problem.addCon(name,'e')

# inequality constraint
Con_Names = The_Project.config_current['OPT_CONSTR']['INEQUALITY'].keys()
for i,name in enumerate(Con_Names):
    The_Problem.addCon(name,'i')

print The_Problem


# -------------------------------------------------------------------
#  Run Optimizer
# -------------------------------------------------------------------     

# SLSQP
The_Optimizer = pyOpt.SLSQP()
The_Optimizer.setOption('IPRINT',-1)
The_Optimizer.setOption('MAXIT',its)
The_Optimizer.setOption('ACC',1e-8)

# CONMIN
#The_Optimizer = pyOpt.CONMIN()
#The_Optimizer.setOption('IPRINT',1)
#The_Optimizer.setOption('ITMAX',its)
#The_Optimizer.setOption('DABFUN',1e-6)
#The_Optimizer.setOption('DELFUN',1e-6)

# Call the Optimizer
The_Optimizer( The_Problem , sens_type=the_Gradient )

print The_Problem.solution(0)


# -------------------------------------------------------------------
#  Save Data
# -------------------------------------------------------------------     

# pull all design data
design_data = The_Project.design_current

# save project data
libSU2.save_data( project_filename,The_Project)

# save python data
libSU2.save_data( design_filename,design_data)

# save matlab data
libSU2.save_data( design_matname,design_data)




