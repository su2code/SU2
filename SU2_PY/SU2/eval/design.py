#!/usr/bin/env python

## \file design.py
#  \brief python package for designs
#  \author T. Lukaczyk, F. Palacios
#  \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

import os, sys, shutil, copy, glob, re
from .. import io   as su2io
from .  import func as su2func
from .  import grad as su2grad
from ..io import redirect_folder, save_data

# todo:
# shouldnt be needed, but self.append_state() (ie after initialization)


# ----------------------------------------------------------------------
#  Design Class
# ----------------------------------------------------------------------

class Design(object):
    """ SU2.eval.Design(config,state=None,folder='DESIGNS/DSN_*')
    
        Starts a design class, which manages a config and state.
        Will run design in folder, and with self indexing name if '*' is
        included in the folder name.
        Methods are wrappers for SU2.eval.func() and SU2.eval.grad()
       
        Attributes:
            state  - design state
            config - design config
            files  - design files
            folder - design folder
            funcs  - design function value bunch
            grads  - design gradient values bunch
        
        Methods:
            Optimizer Interface
            The following methods take a design vector for input
            as a list (shape n) or numpy array (shape n or nx1 or 1xn).
            Values are returned as floats or lists or lists of lists.
            See SU2.eval.obj_f, etc for more detail.
            
            obj_f(dvs)     - objective function              : float
            obj_df(dvs)    - objective function derivatives  : list
            con_ceq(dvs)   - equality constraints            : list
            con_dceq(dvs)  - equality constraint derivatives : list[list]
            con_cieq(dvs)  - inequality constraints          : list
            con_dcieq(dvs) - inequality constraint gradients : list[list]
            
            Fucntional Interface
            The following methods take an objective function name for input.
            func(func_name)                  - function of specified name
            grad(func_name,method='CONTINUOUS_ADJOINT') - gradient of specified name
    """
    
    def __init__(self, problem, state=None, folder='DESIGNS/DSN_*'):
        """ Initializes an SU2 Design """
        
        if '*' in folder: folder = su2io.next_folder(folder)
        
#        print "New Design: %s" % folder
        
        problem= copy.deepcopy(problem)
        config = problem.config
        state  = copy.deepcopy(state)
        state  = su2io.State(state)
        # At the design level, we need to add the files that contain the design variables
        state.find_files(problem.physics, problem.DV_KIND)

        self.problem= problem
        self.config = config
        self.state  = state
        self.files  = state.FILES
        self.funcs  = state.FUNCTIONS
        self.grads  = state.GRADIENTS
        self.folder = folder
        
        self.filename = 'design.pkl'
            
        # initialize folder with files
        pull,link = state.pullnlink(problem)
        with redirect_folder(folder, pull, link, force=True):
            # save design, config
            save_data(self.filename,self)
            config.dump('config_DSN.cfg')
            problem.dump_dv_new('./')   # dumps the values of the design variables
        
    def _eval(self,eval_func,*args):
        """ Evaluates an SU2 Design 
            always adds config and state to the inputs list
        """

        problem= self.problem
        config = self.config
        state  = self.state
        files  = self.files
        folder = self.folder
        
        filename = self.filename

        # Sends files from base folder (or from previous design) into ./DESIGNS/DSN_003

        # check folder
        assert os.path.exists(folder) , 'cannot find design folder %s' % folder
        
        # list files to pull and link
        pull, link = state.pullnlink(problem)
        
        # output redirection, don't re-pull files
        with redirect_folder(folder,pull,link,force=False) as push:
            
            # get timestamp
            timestamp = state.tic()
            
            # run 
            inputs = args + (problem, state)
            vals = eval_func(*inputs)
            
            # save design
            if state.toc(timestamp):
                save_data(filename,self)
            
        #: with redirect folder
        
        # update files
        files.update(state['FILES'])
        
        return vals
    
    def obj_f(self,dvs):
        """ Evaluates SU2 Design Objectives """
        return self._eval(obj_f,dvs)
    
    def obj_df(self,dvs):
        """ Evaluates SU2 Design Objective Gradients """
        return self._eval(obj_df,dvs)

    def con_ceq(self,dvs):
        """ Evaluates SU2 Design Equality Constraints """
        return self._eval(con_ceq,dvs)
    
    def con_dceq(self,dvs):
        """ Evaluates SU2 Design Equality Constraint Gradients """
        return self._eval(con_dceq,dvs)
    
    def con_cieq(self,dvs):
        """ Evaluates SU2 Design Inequality Constraints """
        return self._eval(con_cieq,dvs)
    
    def con_dcieq(self,dvs):
        """ Evaluates SU2 Design Inequality Constraint Gradients """
        return self._eval(con_dcieq,dvs) 

    def func(self,func_name):
        """ Evaluates SU2 Design Functions by Name """
        return self._eval(su2func,func_name)
    
    def grad(self,func_name,method='CONTINUOUS_ADJOINT'):
        """ Evaluates SU2 Design Gradients by Name """
        return self._eval(su2grad,func_name,method)
    
    def touch(self):
        return self._eval(touch)
    
    def skip(self,*args,**kwarg):
        return self._eval(skip)
        
    
    def __repr__(self):
        return '<Design> %s' % self.folder
    def __str__(self):
        output = self.__repr__()
        output += '\n%s' % self.state
        return output
    
#: class Design()


# ----------------------------------------------------------------------
#  Optimization Interface Functions
# ----------------------------------------------------------------------
        
def obj_f(dvs,problem,state=None):
    """ val = SU2.eval.obj_f(dvs,config,state=None)
    
        Evaluates SU2 Objectives 
        Wraps SU2.eval.func()
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Outputs a float.
    """

    # unpack config and state 
    problem.unpack_dvs(dvs)
    config = problem.config
    state = su2io.State(state)

    def_objs = problem['OBJECTIVE_FUNCTION']
    objectives = def_objs.keys()

#    if objectives: print('Evaluate Objectives')
    # evaluate each objective
    vals_out = []
    func = 0.0
    for i_obj,this_obj in enumerate(objectives):
        scale = def_objs[this_obj]['SCALE']
        sign  = su2io.get_objectiveSign(this_obj)
        
        # Evaluate Objective Function
        # scaling and sign
        func += su2func(this_obj, problem, state) * sign * scale
        
    vals_out.append(func)
    #: for each objective
    
    return vals_out

#: def obj_f()

def obj_df(dvs,problem,state=None):
    """ vals = SU2.eval.obj_df(dvs,config,state=None)
    
        Evaluates SU2 Objective Gradients
        Wraps SU2.eval.grad()
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Outputs a list of gradients.
    """    
    
    # unpack config and state
    problem.unpack_dvs(dvs)
    config      = problem.config_adj
    state       = su2io.State(state)
    grad_method = problem['GRADIENT_METHOD']

    def_objs = problem['OBJECTIVE_FUNCTION']
    objectives = def_objs.keys()
    n_obj = len( objectives )
    multi_objective = (config['OPT_COMBINE_OBJECTIVE']=="YES")
     
    dv_scales = config['DEFINITION_DV']['SCALE']
    dv_size   = config['DEFINITION_DV']['SIZE']

    vals_out = []
    for i_obj, this_obj in enumerate(objectives):
        config.GRADIENT_METHOD = grad_method
        grad = su2grad(this_obj, grad_method, problem, state)
        # scaling : obj scale and sign are accounted for in combo gradient, dv scale now applied
        k = 0
        for i_dv, dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] / dv_scl
                k = k + 1
        vals_out.append(grad)

    '''
    #  if objectives: print('Evaluate Objective Gradients')
    # evaluate each objective
    vals_out = []
    if (multi_objective and n_obj>1):
        scale = [1.0]*n_obj
        for i_obj,this_obj in enumerate(objectives):
            sign = su2io.get_objectiveSign(this_obj)
            scale[i_obj] = def_objs[this_obj]['SCALE']*sign
            
        config['OBJECTIVE_WEIGHT']=','.join(map(str,scale))
        
        grad= su2grad(objectives,grad_method,config,state)
        # scaling : obj scale  adn sign are accounted for in combo gradient, dv scale now applied
        k = 0
        for i_dv,dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] / dv_scl
                k = k + 1

        vals_out.append(grad)
    elif (objectives[0] == 'REFERENCE_GEOMETRY' or objectives[0] == 'REFERENCE_NODE'):
        for i_obj,this_obj in enumerate(objectives):
          print "GRAD METHOD", grad_method
          config.GRADIENT_METHOD = 'DISCRETE_ADJOINT'
          grad = su2grad(this_obj, grad_method, config, state)
          # scaling : obj scale and sign are accounted for in combo gradient, dv scale now applied
          k = 0                  
          for i_dv,dv_scl in enumerate(dv_scales):
              for i_grd in range(dv_size[i_dv]):
                  grad[k] = grad[k] / dv_scl
                  k = k + 1
          vals_out.append(grad)
    else:
        marker_monitored = config['MARKER_MONITORING']
        for i_obj,this_obj in enumerate(objectives):
            scale = def_objs[this_obj]['SCALE']
            sign  = su2io.get_objectiveSign(this_obj)
            # Correct marker monitoring for case where multiple objectives are evaluated separately
            if n_obj>1 and len(marker_monitored)>1:
                config['MARKER_MONITORING'] = marker_monitored[i_obj]

            
            # Evaluate Objective Gradient
    #        sys.stdout.write('  %s... ' % this_obj.title())
            grad = su2grad(this_obj,grad_method,config,state)
    #        sys.stdout.write('done\n')
            
            # scaling and sign
            k = 0
            for i_dv,dv_scl in enumerate(dv_scales):
                for i_grd in range(dv_size[i_dv]):
                    grad[k] = grad[k] * sign * scale / dv_scl
                    k = k + 1
            
            vals_out.append(grad)

    '''
    #: for each objective
    
    return vals_out

#: def obj_df()

def con_ceq(dvs,problem,state=None):
    """ vals = SU2.eval.con_ceq(dvs,config,state=None)
    
        Evaluates SU2 Equality Constraints
        Wraps SU2.eval.func()
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Returns a list of constraint values, ordered 
        by the OPT_CONSTRAINT config parameter.
    """
    
    # unpack state and config
    problem.unpack_dvs(dvs)
    state = su2io.State(state)
    
    def_cons = problem['CONSTRAINTS']['EQUALITY']
    constraints = def_cons.keys()
    
 #   if constraints: sys.stdout.write('Evalaute Equality Constraints')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        
        # Evaluate Constraint Function
#        sys.stdout.write('  %s... ' % this_con.title())
        func = su2func(this_con, problem, state)
#        sys.stdout.write('done: %.6f\n' % func)
        
        # scaling and centering
        func = (func - value) * scale
        
        vals_out.append(func)
        
    #: for each constraint
    
    return vals_out

#: def obj_ceq()

def con_dceq(dvs, problem, state=None):
    """ vals = SU2.eval.con_dceq(dvs,config,state=None)
    
        Evaluates SU2 Equality Constraint Gradients
        Wraps SU2.eval.grad()
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Returns a list of lists of constraint gradients,
        ordered by the OPT_CONSTRAINT config parameter.
    """
    
    # unpack state and config
    problem.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = problem['GRADIENT_METHOD']
    
    def_cons = problem['CONSTRAINTS']['EQUALITY']
    constraints = def_cons.keys()

    # config = problem.config
    # dv_scales = config['DEFINITION_DV']['SCALE']
    # dv_size   = config['DEFINITION_DV']['SIZE']

#    if constraints: sys.stdout.write('Evaluate Equality Constraint Gradients ...')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        
        # Evaluate Constraint Gradient
#        sys.stdout.write('  %s... ' % this_con.title())
        grad = su2grad(this_con,grad_method,config,state)
#        sys.stdout.write('done\n')
        
        # scaling
        k = 0
        for i_dv,dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] * scale / dv_scl
                k = k + 1

        vals_out.append(grad)
        
    #: for each constraint
    
    return vals_out

#: def obj_dceq()

def con_cieq(dvs, problem, state=None):
    """ vals = SU2.eval.con_cieq(dvs,config,state=None)
    
        Evaluates SU2 Inequality Constraints
        Wraps SU2.eval.func()
        Convention is con(x)<=0
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Returns a list of constraint gradients, ordered 
        by the OPT_CONSTRAINT config parameter.
    """    
    
    # unpack state and config    
    problem.unpack_dvs(dvs)
    state = su2io.State(state)
    
    def_cons = problem['CONSTRAINTS']['INEQUALITY']
    constraints = def_cons.keys()
    
#    if constraints: sys.stdout.write('Evaluate Inequality Constraints')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        sign  = def_cons[this_con]['SIGN']
        sign  = su2io.get_constraintSign(sign)
        
        # Evaluate Constraint Function
#        sys.stdout.write('  %s... ' % this_con.title())
        func = su2func(this_con,config,state)
#        sys.stdout.write('done: %s\n' % func)
        
        # scaling and centering
        func = (func - value) * scale * sign
        
        vals_out.append(func)
    
    #: for each constraint
    
    return vals_out

#: def obj_cieq()

def con_dcieq(dvs, problem, state=None):
    """ vals = SU2.eval.con_dceq(dvs,config,state=None)
    
        Evaluates SU2 Inequality Constraint Gradients
        Wraps SU2.eval.grad()
        Convention is con(x)<=0
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Returns a list of lists of constraint gradients,
        ordered by the OPT_CONSTRAINT config parameter.
    """    
    
    # unpack state and config
    problem.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = problem['GRADIENT_METHOD']
    
    def_cons = problem['CONSTRAINTS']['INEQUALITY']
    constraints = def_cons.keys()

    config = problem.config

    # dv_scales = config['DEFINITION_DV']['SCALE']
    # dv_size   = config['DEFINITION_DV']['SIZE']

#    if constraints: sys.stdout.write('Evaluate Inequality Constraint Gradients')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        sign  = def_cons[this_con]['SIGN']
        sign  = su2io.get_constraintSign(sign)        
        
        # Evaluate Constraint Gradient
#        sys.stdout.write('  %s... ' % this_con.title())
        grad = su2grad(this_con,grad_method,config,state)
#        sys.stdout.write('done\n')
        
        # scaling and sign
        k = 0
        for i_dv,dv_scl in enumerate(dv_scales):
            for i_grd in range(dv_size[i_dv]):
                grad[k] = grad[k] * sign * scale / dv_scl
                k = k + 1

        vals_out.append(grad)
        
    #: for each constraint
    
    return vals_out
             
    
#: def obj_dcieq()

def touch(config,state):
    """ SU2.eval.touch(config,state)
        resets state timestamp 
    """
    state.set_timestamp()
    
def skip(config,state):
    """ SU2.eval.skip(config,state)
        does nothing
    """
    pass
