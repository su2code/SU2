
import os, sys, shutil, copy
import numpy as np
from .. import io   as su2io
from .. import eval as su2eval
from .. import util as su2util
from ..io import redirect_folder


inf = 1.0e20

# todo: 
# dump design data

# nice:
# shouldnt be needed, but self.append_state() (ie after initialization)
# append_design()

class Project(object):
    
    def __init__( self, config, state=None , 
                  designs=[], folder='.'  ):
        
        folder = folder.rstrip('/')+'/'
        if '*' in folder: folder = su2io.next_folder(folder)        
        
        print 'New Project: %s' % (folder)
        
        config = copy.deepcopy(config)
        state  = copy.deepcopy(state)
        state  = su2io.state.default_state(state)
        state.find_files(config)
        
        self.config  = config      # base config
        self.state   = state       # base state
        self.files   = state.FILES # base files
        self.designs = designs     # design list
        self.folder  = folder      # project folder
        self.data    = su2util.ordered_bunch() # project design data
        
        # initialize folder with files
        pull,link = state.pullnlink(config)
        with redirect_folder(folder,pull,link,force=True):
            pass
            
        return
    
    def _eval(self,config,func,*args):
        
        konfig = copy.deepcopy(config) # design config
        config = self.config           # project config
        state  = self.state            # project state
        folder = self.folder           # project folder
        
        # list project files to pull and link
        pull,link = state.pullnlink(config)
        
        # project folder redirection, don't overwrite files
        with redirect_folder(folder,pull,link,force=False) as push:        
        
            # find closest design
            closest,delta = self.find_design(konfig)
            # found existing design
            if delta == 0.0 and closest:
                design = closest
            # start new design
            else:
                design = self.new_design(konfig,closest)
            #: if new design
            
            # run design
            vals = design._eval(func,*args)
            
            # recompile design data
            self.compile_data()
            
            # plot data
            self.plot_data()
            
            # save project
            su2io.save_data('project.pkl',self)
            
        #: with redirect folder
        
        # done, return output
        return vals
    
    def _unpack_dvs(self,dvs):
        dvs = copy.deepcopy(dvs)
        konfig = copy.deepcopy( self.config )
        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        konfig.unpack_dvs(dvs)
        return konfig, dvs
    
    def obj_f(self,dvs):
        func = su2eval.obj_f
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
        
    def obj_df(self,dvs):
        func = su2eval.obj_df
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_ceq(self,dvs):
        func = su2eval.con_ceq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_dceq(self,dvs):
        func = su2eval.con_dceq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_cieq(self,dvs):
        func = su2eval.con_cieq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def con_dcieq(self,dvs):
        func = su2eval.con_dcieq
        konfig,dvs = self._unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)
    
    def func(self,func_name,config):
        func = su2eval.func
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name)
    
    def grad(self,func_name,method,config):
        func = su2eval.grad
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func, func_name,method)
    
    def user(self,user_func,config,*args):
        raise NotImplementedError
        #return self._eval(config, user_func,*args) 
        
        
    def find_design(self,config):
                
        designs = self.designs
        
        keys_check = ['DV_VALUE_NEW']
        
        if not designs: 
            return [] , inf
        
        diffs = []
        for this_design in designs:
            this_config = this_design.config
            distance = config.dist(this_config,keys_check)
            diffs.append(distance) 
                        
        #: for each design 
        
        # pick closest design
        i_min = np.argmin(diffs)
        delta  = diffs[i_min]
        closest = designs[i_min]
        
        return closest, delta 
    
    def new_design(self,config,closest=[]):
        
        konfig = copy.deepcopy(config)
        ztate  = copy.deepcopy(self.state)
        
        # use closest design as seed
        if closest:
            # copy useful state info
            seed_folder = closest.folder
            seed_files  = closest.files
            for key in seed_files.keys():
                # ignore mesh
                if key == 'MESH': continue 
                # build file path
                name = seed_files[key]
                name = os.path.join(seed_folder,name)
                # update pull files
                ztate.FILES[key] = name

        # start new design (pulls files to folder)
        design = su2eval.Design(konfig,ztate)
        
        # update local state filenames
        for key in design.files:
            name = design.files[key]
            name = os.path.split(name)[-1]
            design.files[key] = name
        
        # add design to project 
        self.designs.append(design)        
        
        return design
    
    def compile_data(self,default=np.nan):
        """
        """
        
        data = su2io.State()
        data.VARIABLES = []
        del data.FILES
        
        n_dv = 0
        
        # populate fields
        for i,design in enumerate(self.designs):
            for key in design.state.FUNCTIONS.keys():
                data.FUNCTIONS[key] = []
            for key in design.state.GRADIENTS.keys():
                data.GRADIENTS[key] = []
            for key in design.state.HISTORY.keys():
                data.HISTORY[key] = []
            this_ndv = len( design.state.design_vector() )
            if i == 0:
                n_dv = this_ndv
            else:
                if n_dv != this_ndv:
                    print 'warning - different dv vector length during reduce_data()'
                
        # populate data
        for design in self.designs:
            this_designvector = design.state.design_vector()
            data.VARIABLES.append( this_designvector )
            for key in data.FUNCTIONS.keys():
                if design.state.FUNCTIONS.has_key(key):
                    new_func = design.state.FUNCTIONS[key]
                else:
                    new_func = default
                data.FUNCTIONS[key].append(new_func)
            for key in data.GRADIENTS.keys():
                if design.state.GRADIENTS.has_key(key):
                    new_grad = design.state.GRADIENTS[key]
                else:
                    new_grad = [default] * len( this_designvector )
                data.GRADIENTS[key].append(new_grad)
            for key in data.HISTORY.keys():
                if design.state.HISTORY.has_key(key):
                    new_func = design.state.HISTORY[key].ITERATION[-1]
                else:
                    new_func = 0.0
                data.HISTORY[key].append(new_func)                
        
        self.data = data
            
        return self.data
    
    def plot_data(self):
        output_format = self.config.OUTPUT_FORMAT
        functions     = self.data.FUNCTIONS
        history       = self.data.HISTORY
        
        data_plot     = su2util.ordered_bunch()
        data_plot.EVALUATION = range(1,len(self.designs)+1)
        data_plot.update(functions)
        data_plot.update(history)
        
        su2util.write_plot('project.plt',output_format,data_plot)
        
    def __repr__(self):
        return '<Project> with %i <Design>' % len(self.designs)
    def __str__(self):
        output = self.__repr__()
        return output    