
import os, sys, shutil, copy, glob, re
from .. import io   as su2io
from .  import func as su2func
from .  import grad as su2grad
from ..io import redirect_folder, save_data

# todo: 
# save self in design folder
# setup working folder on initialization

# nice:
# shouldnt be needed, but self.append_state() (ie after initialization)

class Design(object):
    
    def __init__( self, config, state=None, 
                  folder='DESIGNS/DSN_*'   ):
        
        if '*' in folder: folder = su2io.next_folder(folder)
        
        print "New Design: %s" % folder
        
        config = copy.deepcopy(config)
        state  = copy.deepcopy(state)
        state  = su2io.state.default_state(state)
        state.find_files(config)
        
        self.config = config
        self.state  = state
        self.files  = state.FILES
        self.funcs  = state.FUNCTIONS
        self.grads  = state.GRADIENTS
        self.folder = folder
            
        # initialize folder with files
        pull,link = state.pullnlink(config)
        with redirect_folder(folder,pull,link,force=True):
            pass
        
    def _eval(self,eval_func,*args):
        
        config = self.config
        state  = self.state
        files  = self.files
        folder = self.folder
        
        # list files to pull and link
        pull,link = state.pullnlink(config)
        
        # output redirection, don't re-pull files
        with redirect_folder(folder,pull,link,force=False) as push:
            
            # run 
            inputs = args + (config,state)
            vals = eval_func(*inputs)
            
            # save project
            save_data('design.pkl',self)
            
        #: with redirect folder
        
        # update files
        files.update(state['FILES'])
        
        return vals
    
    def obj_f(self,dvs):
        return self._eval(obj_f,dvs)
    
    def obj_df(self,dvs):
        return self._eval(obj_df,dvs)

    def con_ceq(self,dvs):
        return self._eval(con_ceq,dvs)
    
    def con_dceq(self,dvs):
        return self._eval(con_dceq,dvs)
    
    def con_cieq(self,dvs):
        return self._eval(con_cieq,dvs)
    
    def con_dcieq(self,dvs):
        return self._eval(con_dcieq,dvs) 

    def func(self,func_name):
        return self._eval(su2func,func_name)
    
    def grad(self,func_name,method='ADJOINT'):
        return self._eval(su2grad,func_name,method)
    
    def __repr__(self):
        return '<Design> %s' % self.folder
    def __str__(self):
        output = self.__repr__()
        output += '\n%s' % self.state
        return output
            
            
            
            
            
	
def obj_f(dvs,config,state=None):
    
    # unpack config and state 
    config.unpack_dvs(dvs)
    state = su2io.state.default_state(state)
    
    def_objs = config['OPT_OBJECTIVE']
    objectives = def_objs.keys()
    n_obj = len( objectives )
    assert n_obj == 1 , 'SU2 currently only supports one objective'
    
    if objectives: print('Evaluate Objectives')
    
    # evaluate each objective
    vals_out = []
    for i_obj,this_obj in enumerate(objectives):
        scale = def_objs[this_obj]['SCALE']
        sign  = su2io.get_objectiveSign(this_obj)
        
        # Evaluate Objective Function
        sys.stdout.write('  %s... ' % this_obj.title())
        func = su2func(this_obj,config,state)
        sys.stdout.write('done: %.6f\n' % func)
        
        # scaling and sign
        func = func * sign * scale
        
        vals_out.append(func)
    
    #: for each objective
    
    return vals_out

def obj_df(dvs,config,state=None):
    
    # unpack config and state
    config.unpack_dvs(dvs)
    state = su2io.state.default_state(state)
    
    # assumption
    grad_method = 'ADJOINT'
    
    def_objs = config['OPT_OBJECTIVE']
    objectives = def_objs.keys()
    n_obj = len( objectives )
    assert n_obj == 1 , 'SU2 currently only supports one objective'
    
    dv_scales = config['DEFINITION_DV']['SCALE']
    
    if objectives: print('Evaluate Objective Gradients')
    
    # evaluate each objective
    vals_out = []
    for i_obj,this_obj in enumerate(objectives):
        scale = def_objs[this_obj]['SCALE']
        sign  = su2io.get_objectiveSign(this_obj)
        
        # Evaluate Objective Gradient
        sys.stdout.write('  %s... ' % this_obj.title())
        grad = su2grad(this_obj,grad_method,config,state)
        sys.stdout.write('done\n')
        
        # scaling and sign
        for i_grd,dv_scl in enumerate(dv_scales):
            grad[i_grd] = grad[i_grd] * sign * scale / dv_scl
        
        vals_out.append(grad)
    
    #: for each objective
    
    return vals_out

def con_ceq(dvs,config,state=None):
    
    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.state.default_state(state)
    
    def_cons = config['OPT_CONSTRAINT']['EQUALITY']
    constraints = def_cons.keys()
    
    if constraints: print('Evalaute Equality Constraints')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        
        # Evaluate Constraint Function
        sys.stdout.write('  %s... ' % this_con.title())
        func = su2func(this_con,config,state)
        sys.stdout.write('done: %.6f\n' % func)
        
        # scaling and centering
        func = (func - value) * scale
        
        vals_out.append(func)
        
    #: for each constraint
    
    return vals_out

def con_dceq(dvs,config,state=None):
    
    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.state.default_state(state)
    
    # assumption
    grad_method = 'ADJOINT'
    
    def_cons = config['OPT_CONSTRAINT']['EQUALITY']
    constraints = def_cons.keys()
    
    dv_scales = config['DEFINITION_DV']['SCALE']
    
    if constraints: print('Evaluate Equality Constraint Gradients ...')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        
        # Evaluate Constraint Gradient
        sys.stdout.write('  %s... ' % this_con.title())
        grad = su2grad(this_con,grad_method,config,state)
        sys.stdout.write('done\n')
        
        # scaling
        for i_grd,dv_scl in enumerate(dv_scales):
            grad[i_grd] = grad[i_grd] * scale / dv_scl     
        
        vals_out.append(grad)
        
    #: for each constraint
    
    return vals_out

def con_cieq(dvs,config,state=None):
    
    # unpack state and config    
    config.unpack_dvs(dvs)
    state = su2io.state.default_state(state)
    
    def_cons = config['OPT_CONSTRAINT']['INEQUALITY']
    constraints = def_cons.keys()
    
    if constraints: print('Evaluate Inequality Constraints')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        sign  = def_cons[this_con]['SIGN']
        sign  = su2io.get_constraintSign(sign)
        
        # Evaluate Constraint Function
        sys.stdout.write('  %s... ' % this_con.title())
        func = su2func(this_con,config,state)
        sys.stdout.write('done: %s\n' % func)
        
        # scaling and centering
        func = (func - value) * scale * sign
        
        vals_out.append(func)
    
    #: for each constraint
    
    return vals_out

def con_dcieq(dvs,config,state=None):
    
    # unpack state and config
    config.unpack_dvs(dvs)
    state = su2io.state.default_state(state)
    
    # assumption
    grad_method = 'ADJOINT'
    
    def_cons = config['OPT_CONSTRAINT']['INEQUALITY']
    constraints = def_cons.keys()
    
    dv_scales = config['DEFINITION_DV']['SCALE']
    
    if constraints: print('Evaluate Inequality Constraint Gradients')
    
    # evaluate each constraint
    vals_out = []
    for i_obj,this_con in enumerate(constraints):
        scale = def_cons[this_con]['SCALE']
        value = def_cons[this_con]['VALUE']
        sign  = def_cons[this_con]['SIGN']
        sign  = su2io.get_constraintSign(sign)        
        
        # Evaluate Constraint Gradient
        sys.stdout.write('  %s... ' % this_con.title())
        grad = su2grad(this_con,grad_method,config,state)
        sys.stdout.write('done\n')
        
        # scaling and sign
        for i_grd,dv_scl in enumerate(dv_scales):
            grad[i_grd] = grad[i_grd] * sign * scale / dv_scl          

        vals_out.append(grad)
        
    #: for each constraint
    
    return vals_out
    
