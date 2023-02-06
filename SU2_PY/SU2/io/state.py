#!/usr/bin/env python

## \file state.py
#  \brief python package for state
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, time
from ..io   import expand_part, expand_zones, expand_time, get_adjointSuffix, add_suffix, \
                   get_specialCases, Config, expand_multipoint, optnames_multi
from ..util import bunch
from ..util import ordered_bunch


# ----------------------------------------------------------------------
#  State Factory
# ----------------------------------------------------------------------

def State_Factory(state=None,config=None):
    """ state = SU2.io.State()

        Starts a state class, an extension of ordered_bunch().
        Stores data generated while traversing SU2 tool chain

        Fields:
            FUNCTIONS - ordered bunch of objective function values
            GRADIENTS - ordered bunch of gradient value lists
            VARIABLES - ordered bunch of variables
            FILES     - ordered bunch of file types
            HISTORY   - ordered bunch of history information

        Fields can be accessed by item or attribute
        ie: state['FUNCTIONS'] or state.FUNCTIONS

        Methods:
            update()        - updates self with another state
            pullnlink()     - returns files to pull and link
            design_vector() - vectorizes design variables
            find_files()    - finds existing mesh and solutions

        Example of a filled state:
        FUNCTIONS:
            LIFT: 0.2353065809
            DRAG: 0.042149736
            SIDEFORCE: 0.0
            MOMENT_X: 0.0
            MOMENT_Y: 0.0
            MOMENT_Z: 0.046370243
            FORCE_X: 0.0370065195
            FORCE_Y: 0.2361700759
            FORCE_Z: 0.0
            EFFICIENCY: 5.5826347517
        GRADIENTS:
            DRAG: [0.133697, 0.41473, 0.698497, (...)
        VARIABLES:
            DV_VALUE_NEW: [0.002, 0.002, 0.002, (...)
        FILES:
            MESH: mesh.su2
            DIRECT: solution_flow.dat
            ADJOINT_DRAG: solution_adj_cd.dat
            FLOW_META: flow.meta
            MULTIPOINT_DIRECT: [solution_flow_point0.dat solution_flow_point1.dat, ...]
            MULTIPOINT_ADJOINT_DRAG: [solution_adj_point0_cd.dat solution_adj_point1_cd.dat, ...]
            MULTIPOINT_MESH_FILENAME: [mesh_0.su2, mesh_1.su2, ... ]
            MULTIPOINT_FLOW_META: [flow_point0.meta, flow_point1.meta, ...]
        HISTORY:
            DIRECT: {ITERATION=[1.0, 2.0, 3.0, (...)
            ADJOINT_DRAG: {ITERATION=[1.0, 2.0, 3.0, (...)
        WND_CAUCHY_DATA:
            TIME_ITER
            UNST_ADJOINT_ITER
            ITER_AVERAGE_OBJ

    """

    if isinstance(state,Config) and not config:
        config = state
        state = None

    if not state is None:
        assert isinstance(state,State) , 'input is must be a state instance'
        return state

    NewClass = State()

    for key in ['FUNCTIONS','GRADIENTS','VARIABLES','FILES','HISTORY','WND_CAUCHY_DATA']:
        NewClass[key] = ordered_bunch()

    if config:
        NewClass.find_files(config)
        # WND_Convergence Data
        NewClass['WND_CAUCHY_DATA'] = {'TIME_ITER': config['TIME_ITER'],
                                       'UNST_ADJOINT_ITER': config['UNST_ADJOINT_ITER'],
                                       'ITER_AVERAGE_OBJ': config['ITER_AVERAGE_OBJ']}

    return NewClass


# ----------------------------------------------------------------------
#  State Class
# ----------------------------------------------------------------------

class State(ordered_bunch):
    """ state = SU2.io.state.State()

        This is the State class that should be generated with the
        Factory Function SU2.io.state.State_Factory()

        Parameters:
            none, should be loaded with State_Factory()

        Methods:
            update()        - updates self with another state
            pullnlink()     - returns files to pull and link
            design_vector() - vectorizes design variables
            find_files()    - finds existing mesh and solutions

    """

    _timestamp = 0

    def update(self,ztate):
        """ Updates self given another state
        """

        if not ztate: return
        assert isinstance(ztate,State) , 'must update with another State-type'
        for key in self.keys():
            if isinstance(ztate[key],dict):
                self[key].update( ztate[key] )
            elif ztate[key]:
                self[key] = ztate[key]

        self.set_timestamp()


    def __repr__(self):
        return self.__str__()

    def __str__(self):
        output = 'STATE:'
        for k1, v1 in self.items():
            output += '\n    %s:' % k1
            if isinstance(v1,dict):
                for k2, v2 in v1.items():
                    output += '\n        %s: %s' % (k2,v2)
            else:
                output += '\n        %s' % v1
        return output

    def pullnlink(self,config):
        """ pull,link = SU2.io.State.pullnlink(config)
            returns lists pull and link of files for folder
            redirection, based on a given config
        """

        pull = []; link = []

        # choose files to pull and link
        for key, value in self.FILES.items():
            # link big files
            if key == 'MESH':
                # mesh (merged or partitioned)
                value = expand_part(value,config)
                link.extend(value)
            elif key == 'DIRECT':
                # direct solution
                value = expand_zones(value,config)
                value = expand_time(value,config)
                link.extend(value)
            elif 'ADJOINT_' in key and  (not 'MULTIPOINT' in key):
                # adjoint solution
                value = expand_zones(value,config)
                value = expand_time(value,config)
                link.extend(value)
            elif 'MULTIPOINT' in key:
                # multipoint files
                if key != 'MULTIPOINT_MESH_FILENAME':
                    # DIRECT and ADJOINT files
                    value = expand_zones(value,config)
                    value = expand_time(value,config)
                for elem in value:
                    if elem:
                        link.append(elem)
            #elif key == 'STABILITY':
                #pass
            # copy all other files
            else:
                pull.append(value)

        #: for each filename

        return pull,link

    def design_vector(self):
        """ vectorizes State.VARIABLES
        """
        vector = []
        for value in self.VARIABLES.values():
            if isinstance(value,dict):
                for v in value.values():
                    vector.append(v)
            elif not isinstance(value,list):
                value = [value]
            vector.extend(value)
        return vector

    def find_files(self,config):
        """ SU2.io.State.find_files(config)
            finds mesh and solution files for a given config.
            updates state.FILES with filenames.
            files already logged in state are not overridden.
            will ignore solutions if config.RESTART_SOL == 'NO'.
        """

        files = self.FILES

        mesh_name     = config.MESH_FILENAME
        if config.get('READ_BINARY_RESTART', 'YES') == 'NO':
            if not 'RESTART_ASCII' in config.get('OUTPUT_FILES',['RESTART']):
                print ('RESTART_ASCII must be in OUTPUT_FILES if READ_BINARY_RESTART is set to NO')
                sys.exit()

        direct_name   = config.SOLUTION_FILENAME
        adjoint_name  = config.SOLUTION_ADJ_FILENAME

        if 'RESTART_ASCII' in config.get('OUTPUT_FILES', ['RESTART']):
            direct_name = direct_name.split('.')[0] + '.csv'
            adjoint_name = adjoint_name.split('.')[0] + '.csv'
        else:
            direct_name = direct_name.split('.')[0] + '.dat'
            adjoint_name = adjoint_name.split('.')[0] + '.dat'

        targetea_name = 'TargetEA.dat'
        targetcp_name = 'TargetCp.dat'
        targetheatflux_name = 'TargetHeatFlux.dat'

        adj_map = get_adjointSuffix()
        restart = config.RESTART_SOL == 'YES'
        special_cases = get_specialCases(config)

        if config.get('OPT_OBJECTIVE'):
            def_objs = config['OPT_OBJECTIVE']
            objectives = def_objs.keys()
            multipoint = any(elem in optnames_multi for elem in objectives)
        else:
            multipoint = False

        def register_file(label,filename):
            if not label in files:
                if label.split('_')[0] in ['DIRECT', 'ADJOINT']:
                    names = expand_zones(filename, config)
                    found = False
                    for name in names:
                        if os.path.exists(name):
                            found = True
                        else:
                            found = False
                            break

                    if found:
                        files[label] = filename
                        print('Found: %s' % filename)

                elif label.split('_')[0] in ['MULTIPOINT']:
                    # if multipoint, list of files needs to be added
                    file_list= [];
                    for name in filename:
                        if os.path.exists(name):
                            file_list.append(name)
                            print('Found: %s' % name)
                        else:
                            # if file doesn't exist, enter empty string as placeholder
                            file_list.append('')
		            # If even one of the multipoint files is found, add the list
                    if any(file for file in file_list):
                        files[label] = file_list
                else:
                    if os.path.exists(filename):
                        files[label] = filename
                        print('Found: %s' % filename)
            else:
                if label.split("_")[0] in ['DIRECT', 'ADJOINT']:
                    for name in expand_zones(files[label], config):
                        assert os.path.exists(name), 'state expected file: %s' % filename
                elif label.split('_')[0] in ['MULTIPOINT']:
                    for name in expand_zones(files[label], config):
                        if name:
                            if not os.path.exists(name):
                                raise AssertionError('state expected file: %s' % name)
                else:
                    assert os.path.exists(files[label]) , 'state expected file: %s' % filename
        #: register_file()

        # mesh
        if multipoint:
            mesh_list = [elem.strip() for elem in config['MULTIPOINT_MESH_FILENAME'].replace("(", "").replace(")", "").split(',')]
            if len(set(mesh_list)) > 1:
                # Only register MULTIPOINT_MESH_FILENAME if multiple meshes are specified
                register_file('MULTIPOINT_MESH_FILENAME', mesh_list)
            mesh_name = mesh_list[0]

        register_file('MESH',mesh_name)

        # old style restart
        if not 'RESTART_FILE_1' in files.keys():
            # direct solutions
            if restart:
                register_file('DIRECT',direct_name)
                if multipoint:
                    name_list = expand_multipoint(direct_name,config)
                    name_list = expand_zones(name_list,config)
                    register_file('MULTIPOINT_DIRECT',name_list)

            # flow meta data file
            if restart:
                register_file('FLOW_META','flow.meta')
                if multipoint:
                    name_list = expand_multipoint('flow.meta',config)
                    register_file('MULTIPOINT_FLOW_META',name_list)

            # adjoint solutions
            if restart:
                for obj, suff in adj_map.items():
                    ADJ_LABEL = 'ADJOINT_' + obj
                    adjoint_name_suffixed = add_suffix(adjoint_name,suff)
                    register_file(ADJ_LABEL,adjoint_name_suffixed)
                    if multipoint:
                        name_list = expand_zones(add_suffix(expand_multipoint(adjoint_name,config), suff), config)
                        multipoint_adj_name = 'MULTIPOINT_' + ADJ_LABEL
                        register_file(multipoint_adj_name, name_list)

        # equivalent area
        if 'EQUIV_AREA' in special_cases:
            register_file('TARGET_EA',targetea_name)

        # pressure inverse design
        if 'INV_DESIGN_CP' in special_cases:
            register_file('TARGET_CP',targetcp_name)

        # heat flux inverse design
        if 'INV_DESIGN_HEATFLUX' in special_cases:
            register_file('TARGET_HEATFLUX',targetheatflux_name)

        return

    def __setitem__(self,k,v):
        if self._initialized:
            self.set_timestamp()
        super(State,self).__setitem__(k,v)

    def set_timestamp(self):
        self._timestamp = time.time()

    def tic(self):
        """ timestamp = State.tic()
            returns the time that this state was last modified
        """
        return self._timestamp

    def toc(self,timestamp):
        """ updated = State.toc(timestamp)
            returns True if state was modified since last timestamp
        """
        return self._timestamp > timestamp


#: def State
