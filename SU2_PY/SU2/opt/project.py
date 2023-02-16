#!/usr/bin/env python

## \file project.py
#  \brief package for optimization projects
#  \author T. Lukaczyk, F. Palacios
#  \version 7.5.1 "Blackbird"
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

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

# -------------------------------------------------------------------
#  Imports
# -------------------------------------------------------------------

import os, sys, shutil, copy, glob, time
import numpy as np
from .. import io   as su2io
from .. import eval as su2eval
from .. import util as su2util
from ..io import redirect_folder
from ..io import historyOutFields
from warnings import warn, simplefilter
#simplefilter(Warning,'ignore')

inf = 1.0e20


# -------------------------------------------------------------------
#  Project Class
# -------------------------------------------------------------------

class Project(object):
    """ project = SU2.opt.Project(self,config,state=None,
                                  designs=[],folder='.')

        Starts a project class to manage multiple designs

        Runs multiple design classes, avoiding redundancy
        Looks for closest design on restart
        Currently only based on DV_VALUE_NEW
        Exposes all methods of SU2.eval.design

        Attributes:
             config  - base config
             state   - base state
             files   - base files
             designs - list of designs
             folder  - project working folder
             results - project design results

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

            Functional Interface
            The following methods take an objective function name for input.
            func(func_name,config)        - function of specified name
            grad(func_name,method,config) - gradient of specified name,
                                            where method is 'CONTINUOUS_ADJOINT' or 'FINDIFF'
            setup config for given dvs with
            config = project.unpack_dvs(dvs)
    """

    _design_folder = 'DESIGNS/DSN_*'
    _design_number = '%03d'


    def __init__( self, config, state=None ,
                  designs=None, folder='.' ,
                  warn = True                ):

        folder = folder.rstrip('/')+'/'
        if '*' in folder: folder = su2io.next_folder(folder)
        if designs is None: designs = []

        print('New Project: %s' % (folder))

        # setup config
        config = copy.deepcopy(config)

        # data_dict creation does not preserve the ordering of the config file.
        # This section ensures that the order of markers and objectives match
        # It is only needed when more than one objective is used.
        def_objs = config['OPT_OBJECTIVE']
        if len(def_objs)>1:
            objectives = def_objs.keys()
            marker_monitoring = []
            weights = []
            for i_obj, this_obj in enumerate(objectives):
                marker_monitoring+=[def_objs[this_obj]['MARKER']]
                weights+=[str(def_objs[this_obj]['SCALE'])]
            config['MARKER_MONITORING'] = marker_monitoring
            config['OBJECTIVE_WEIGHT'] = ",".join(weights)
            config['OBJECTIVE_FUNCTION'] = ",".join(objectives)

        for this_obj in def_objs:
            if this_obj in su2io.optnames_multi:
                this_obj = this_obj.split('_')[1]
            group = historyOutFields[this_obj]['GROUP']
            if not group in config.HISTORY_OUTPUT:
                config.HISTORY_OUTPUT.append(group)

        # setup state
        if state is None:
            state = su2io.State()
        else:
            state  = copy.deepcopy(state)
            state  = su2io.State(state)
        state.find_files(config)

        if 'MESH' not in state.FILES:
            raise Exception('Could not find mesh file: %s' % config.MESH_FILENAME)

        self.config  = config      # base config
        self.state   = state       # base state
        self.files   = state.FILES # base files
        self.designs = designs     # design list
        self.folder  = folder      # project folder
        self.results = su2util.ordered_bunch() # project design results

        # output filenames
        self.filename = 'project.pkl'
        self.results_filename = 'results.pkl'

        # initialize folder with files
        pull,link = state.pullnlink(config)

        with redirect_folder(folder,pull,link,force=True):

            # look for existing designs
            folders = glob.glob(self._design_folder)
            if len(folders)>0:
                sys.stdout.write('Removing old designs in 10s.')
                sys.stdout.flush()
                if warn: time.sleep(10)
                sys.stdout.write(' Done!\n\n')
                for f in folders: shutil.rmtree(f)
            #: if existing designs

            # save project
            su2io.save_data(self.filename,self)

        return

    def _eval(self,config,func,*args):
        """ evalautes a config, checking for existing designs
        """

        konfig = copy.deepcopy(config) # design config
        config = self.config           # project config
        state  = self.state            # project state
        folder = self.folder           # project folder
        filename = self.filename

        # check folder
        assert os.path.exists(folder) , 'cannot find project folder %s' % folder

        # list project files to pull and link
        pull,link = state.pullnlink(config)

        # project folder redirection, don't overwrite files
        with redirect_folder(folder,pull,link,force=False) as push:

            # start design
            design = self.new_design(konfig)

            if config.get('CONSOLE','VERBOSE') == 'VERBOSE':
                print(os.path.join(self.folder,design.folder))
            timestamp = design.state.tic()

            # set right option in design config.
            if konfig.get('TIME_DOMAIN', 'NO') == 'YES' and konfig.get('RESTART_SOL', 'NO') == 'YES':
                design.config['RESTART_SOL'] = 'YES'

            # run design+
            vals = design._eval(func,*args)

            # check for update
            if design.state.toc(timestamp):

                # recompile design results
                self.compile_results()

                # plot results
                self.plot_results()

                # save data
                su2io.save_data(filename,self)

            #: if updated

        #: with redirect folder

        # done, return output
        return vals

    def unpack_dvs(self,dvs):
        dvs = copy.deepcopy(dvs)
        konfig = copy.deepcopy( self.config )
        if isinstance(dvs, np.ndarray): dvs = dvs.tolist()
        konfig.unpack_dvs(dvs)
        return konfig, dvs

    def obj_f(self,dvs):
        func = su2eval.obj_f
        konfig,dvs = self.unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)

    def obj_df(self,dvs):
        func = su2eval.obj_df
        konfig,dvs = self.unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)

    def con_ceq(self,dvs):
        func = su2eval.con_ceq
        konfig,dvs = self.unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)

    def con_dceq(self,dvs):
        func = su2eval.con_dceq
        konfig,dvs = self.unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)

    def con_cieq(self,dvs):
        func = su2eval.con_cieq
        konfig,dvs = self.unpack_dvs(dvs)
        return self._eval(konfig, func,dvs)

    def con_dcieq(self,dvs):
        func = su2eval.con_dcieq
        konfig,dvs = self.unpack_dvs(dvs)
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

    def add_design(self,config):
        #func = su2eval.touch # hack - TWL
        func = su2eval.skip
        konfig = copy.deepcopy(config)
        return self._eval(konfig, func)

    def new_design(self,config):
        """ finds an existing design for given config
            or starts a new design with a closest design
            used for restart data
        """
         # local konfig
        konfig = copy.deepcopy(config)
        # find closest design
        closest,delta = self.closest_design(konfig)
        # found existing design
        if delta == 0.0 and closest:
            design = closest
        # start new design
        else:
            design = self.init_design(konfig,closest)
        #: if new design
        return design

    def get_design(self,config):
        konfig = copy.deepcopy(config)
        closest,delta = self.closest_design(konfig)
        if delta == 0.0 and closest:
            design = closest
        else:
            raise Exception('design not found for this config')
        return design

    def closest_design(self,config):
        """ looks for an existing or closest design
            given a config
        """

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

    def init_design(self,config,closest=None):
        """ starts a new design
            works in project folder
        """

        konfig = copy.deepcopy(config)
        ztate  = copy.deepcopy(self.state)
        if closest is None: closest = []

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
                if isinstance(name,list):
                    built_name = []
                    for elem in name:
                        built_name.append(os.path.join(seed_folder,elem))
                    ztate.FILES[key] = built_name
                else:
                    name = os.path.join(seed_folder,name)
                    # update pull files
                    ztate.FILES[key] = name

        # name new folder
        folder = self._design_folder.replace('*',self._design_number)
        folder = folder % (len(self.designs) + 1)

        # start new design (pulls files to folder)
        design = su2eval.Design(konfig,ztate,folder)

        # update local state filenames ( ??? why not in Design() )
        for key in design.files:
            name = design.files[key]
            if isinstance(name,list):
                built_name = []
                for elem in name:
                    built_name.append(os.path.split(elem)[-1])
                design.files[key] = built_name
            else:
                name = os.path.split(name)[-1]
                design.files[key] = name

        # add design to project
        self.designs.append(design)

        return design

    def compile_results(self,default=np.nan):
        """ results = SU2.opt.Project.compile_results(default=np.nan)
            builds a Bunch() of design results

            Inputs:
                default - value for missing values

            Outputs:
                results - state with items filled with list of
                values ordered by each design iteration.

                results.VARIABLES
                results.FUNCTIONS
                results.GRADIENTS
                results.HISTORY.DIRECT
                results.HISTORY.ADJOINT_*

        """

        results = su2io.State()
        results.VARIABLES = []
        del results.FILES
        filename = self.results_filename

        n_dv = 0

        # populate fields
        for i,design in enumerate(self.designs):
            for key in design.state.FUNCTIONS.keys():
                results.FUNCTIONS[key] = []
            for key in design.state.GRADIENTS.keys():
                results.GRADIENTS[key] = []
            for TYPE in design.state.HISTORY.keys():
                if not TYPE in results.HISTORY:
                    results.HISTORY[TYPE] = su2util.ordered_bunch()
                for key in design.state.HISTORY[TYPE].keys():
                    results.HISTORY[TYPE][key] = []
            this_ndv = len( design.state.design_vector() )

            # check design vectors are of same length
            if i == 0:
                n_dv = this_ndv
            else:
                if n_dv != this_ndv:
                    warn('different dv vector length during compile_results()')
        #: for each design

        # populate results
        for design in self.designs:
            this_designvector = design.state.design_vector()
            results.VARIABLES.append( this_designvector )
            for key in results.FUNCTIONS.keys():
                if key in design.state.FUNCTIONS:
                    new_func = design.state.FUNCTIONS[key]
                else:
                    new_func = default
                results.FUNCTIONS[key].append(new_func)
            for key in results.GRADIENTS.keys():
                if key in design.state.GRADIENTS:
                    new_grad = design.state.GRADIENTS[key]
                else:
                    new_grad = [default] * len( this_designvector )
                results.GRADIENTS[key].append(new_grad)
            for TYPE in results.HISTORY.keys():
                for key in results.HISTORY[TYPE].keys():
                    if key in results.FUNCTIONS.keys():
                        new_func = results.FUNCTIONS[key][-1]
                    elif ( TYPE in design.state.HISTORY.keys() and
                            key in design.state.HISTORY[TYPE].keys() ):
                        new_func = design.state.HISTORY[TYPE][key][-1]
                    else:
                        new_func = default
                    results.HISTORY[TYPE][key].append(new_func)
        #: for each design

        # save
        self.results = results
        su2io.save_data(filename,results)

        return self.results

    def deep_compile(self):
        """ Project.deep_compile()
            recompiles project using design files saved in each design folder
            useful if designs were run outside of project class
        """

        project_folder = self.folder
        designs = self.designs

        with su2io.redirect_folder(project_folder):
            for i_dsn,design in enumerate(designs):
                design_filename = os.path.join(design.folder,design.filename)
                self.designs[i_dsn] = su2io.load_data(design_filename)

            self.compile_results()
            su2io.save_data(self.filename,self)

        return

    def plot_results(self):
        """ writes a tecplot file for plotting design results
        """
        output_format = self.config.TABULAR_FORMAT
        functions     = self.results.FUNCTIONS
        history       = self.results.HISTORY

        results_plot = su2util.ordered_bunch()
        results_plot.EVALUATION = range(1,len(self.designs)+1)
        results_plot.update(functions)
        results_plot.update(history.get('DIRECT',{}))

        if (output_format == 'CSV'):
          su2util.write_plot('history_project.csv',output_format,results_plot)
        else:
          su2util.write_plot('history_project.dat',output_format,results_plot)

    def save(self):
        with su2io.redirect_folder(self.folder):
            su2io.save_data(self.filename,self)

    def __repr__(self):
        return '<Project> with %i <Design>' % len(self.designs)
    def __str__(self):
        output = self.__repr__()
        return output
