#!/usr/bin/env python

## \file design.py
#  \brief python package for designs
#  \author T. Lukaczyk, F. Palacios
#  \version 4.1.2 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
    
    def __init__(self, config, state=None, folder='DESIGNS/DSN_*'):
        """ Initializes an SU2 Design """
        
        ## ???: Move to Project, no next folder here
        
        if '*' in folder: folder = su2io.next_folder(folder)
        
#        print "New Design: %s" % folder
        
        config = copy.deepcopy(config)
        state  = copy.deepcopy(state)
        state  = su2io.State(state)
        state.find_files(config)
        
        self.config = config
        self.state  = state
        self.files  = state.FILES
        self.funcs  = state.FUNCTIONS
        self.grads  = state.GRADIENTS
        self.folder = folder
        
        self.filename = 'design.pkl'
            
        # initialize folder with files
        pull,link = state.pullnlink(config)
        with redirect_folder(folder,pull,link,force=True):
            # save design, config
            save_data(self.filename,self)
            config.dump('config_DSN.cfg')
        
    def _eval(self,eval_func,*args):
        """ Evaluates an SU2 Design 
            always adds config and state to the inputs list
        """
        
        config = self.config
        state  = self.state
        files  = self.files
        folder = self.folder
        
        filename = self.filename

        # check folder
        assert os.path.exists(folder) , 'cannot find design folder %s' % folder
        
        # list files to pull and link
        pull,link = state.pullnlink(config)
        
        # output redirection, don't re-pull files
        with redirect_folder(folder,pull,link,force=False) as push:
            
            # get timestamp
            timestamp = state.tic()
            
            # run 
            inputs = args + (config,state)
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
        
def obj_f(dvs,config,state=None):
    """ val = SU2.eval.obj_f(dvs,config,state=None)
    
        Evaluates SU2 Objectives 
        Wraps SU2.eval.func()
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Outputs a float.
    """
    
    # unpack config and state 
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    
    def_objs = config['OPT_OBJECTIVE']
    objectives = def_objs.keys()
    n_obj = len( objectives )
    assert n_obj == 1 , 'SU2 currently only supports one objective'
    
#    if objectives: print('Evaluate Objectives')
    
    # evaluate each objective
    vals_out = []
    for i_obj,this_obj in enumerate(objectives):
        scale = def_objs[this_obj]['SCALE']
        sign  = su2io.get_objectiveSign(this_obj)
        
        # Evaluate Objective Function
#        sys.stdout.write('  %s... ' % this_obj.title())
        func = su2func(this_obj,config,state)
#        sys.stdout.write('done: %.6f\n' % func)
        
        # scaling and sign
        func = func * sign * scale
        
        vals_out.append(func)
    
    #: for each objective
    
    return vals_out

#: def obj_f()

def obj_df(dvs,config,state=None):
    """ vals = SU2.eval.obj_df(dvs,config,state=None)
    
        Evaluates SU2 Objective Gradients
        Wraps SU2.eval.grad()
        
        Takes a design vector for input as a list (shape n) 
        or numpy array (shape n or nx1 or 1xn), a config
        and optionally a state.
        
        Outputs a list of gradients.
    """    
    
    # unpack config and state
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = config.get('GRADIENT_METHOD','CONTINUOUS_ADJOINT')
    
    def_objs = config['OPT_OBJECTIVE']
    objectives = def_objs.keys()
    n_obj = len( objectives )
    assert n_obj == 1 , 'SU2 currently only supports one objective'
    
    dv_scales = config['DEFINITION_DV']['SCALE']
    dv_size   = config['DEFINITION_DV']['SIZE']
    
#    if objectives: print('Evaluate Objective Gradients')
    
    # evaluate each objective
    vals_out = []
    for i_obj,this_obj in enumerate(objectives):
        scale = def_objs[this_obj]['SCALE']
        sign  = su2io.get_objectiveSign(this_obj)
        
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
    
    #: for each objective
    
    return vals_out

#: def obj_df()

def con_ceq(dvs,config,state=None):
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
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    
    def_cons = config['OPT_CONSTRAINT']['EQUALITY']
    constraints = def_cons.keys()
    
 #   if constraints: sys.stdout.write('Evalaute Equality Constraints')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        
        # Evaluate Constraint Function
#        sys.stdout.write('  %s... ' % this_con.title())
        func = su2func(this_con,config,state)
#        sys.stdout.write('done: %.6f\n' % func)
        
        # scaling and centering
        func = (func - value) * scale
        
        vals_out.append(func)
        
    #: for each constraint
    
    return vals_out

#: def obj_ceq()

def con_dceq(dvs,config,state=None):
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
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = config.get('GRADIENT_METHOD','CONTINUOUS_ADJOINT')
    
    def_cons = config['OPT_CONSTRAINT']['EQUALITY']
    constraints = def_cons.keys()
    
    dv_scales = config['DEFINITION_DV']['SCALE']
    dv_size   = config['DEFINITION_DV']['SIZE']

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

def con_cieq(dvs,config,state=None):
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
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    
    def_cons = config['OPT_CONSTRAINT']['INEQUALITY']
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

def con_dcieq(dvs,config,state=None):
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
    config.unpack_dvs(dvs)
    state = su2io.State(state)
    grad_method = config.get('GRADIENT_METHOD','CONTINUOUS_ADJOINT')
    
    def_cons = config['OPT_CONSTRAINT']['INEQUALITY']
    constraints = def_cons.keys()
    
    dv_scales = config['DEFINITION_DV']['SCALE']
    dv_size   = config['DEFINITION_DV']['SIZE']

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
