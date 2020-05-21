#!/usr/bin/env python3

## \file FSI_config.py
#  \brief Python class for handling configuration file for FSI computation.
#  \author David Thomas, Rocco Bombardieri
#  \version 7.0.2 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2020, SU2 Contributors (cf. AUTHORS.md)
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

from .switch import switch

# ----------------------------------------------------------------------
#  FSI Configuration Class
# ----------------------------------------------------------------------

class FSIConfig:
    """
    Class that contains all the parameters coming from the FSI configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self, FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.readConfig()

    def __str__(self):
        tempString = str()
        for key, value in self._ConfigContent.items():
            tempString += "{} = {}\n".format(key, value)
        return tempString

    def __getitem__(self, key):
        return self._ConfigContent[key]

    def __setitem__(self, key, value):
        self._ConfigContent[key] = value

    def readConfig(self):
        input_file = open(self.ConfigFileName)
        while 1:
            line = input_file.readline()
            if not line:
                break
            # remove line returns
            line = line.strip('\r\n')
            # make sure it has useful data
            if (not "=" in line) or (line[0] == '%'):
                continue
            # split across equal sign
            line = line.split("=", 1)
            this_param = line[0].strip()
            this_value = line[1].strip()

            for case in switch(this_param):
                # integer values
                if case("RESTART_ITER"): pass
                if case("NDIM"): pass
                if case("NB_EXT_ITER"): pass
                if case("NB_FSI_ITER"):
                    self._ConfigContent[this_param] = int(this_value)
                    break

                # float values
                if case("FSI_TOLERANCE"):
                    self._ConfigContent[this_param] = float(this_value)
                    break
                if case("RELAX_PARAM"):
                    self._ConfigContent[this_param] = float(this_value)
                    break

                # string values
                if case("SU2_CONFIG"): pass
                if case("PYBEAM_CONFIG"): pass
                if case("MLS_CONFIG_FILE_NAME"): pass
                if case("INTERNAL_FLOW"):
                    self._ConfigContent[this_param] = this_value
                    break

                if case():
                    print(this_param + " is an invalid option !")
            # end for
           
    # def dump()

# ----------------------------------------------------------------------
#  FSI optimization Configuration Class
# ----------------------------------------------------------------------

class OptConfig:
    """
    Class that contains all the parameters coming from the FSI optimization configuration file.
    Read the file and store all the options into a dictionary.
    """

    def __init__(self, FileName):
        self.ConfigFileName = FileName
        self._ConfigContent = {}
        self.readConfig()

    def __str__(self):
        tempString = str()
        for key, value in self._ConfigContent.items():
            tempString += "{} = {}\n".format(key, value)
        return tempString

    def __getitem__(self, key):
        return self._ConfigContent[key]

    def __setitem__(self, key, value):
        self._ConfigContent[key] = value

    def readConfig(self):
        input_file = open(self.ConfigFileName)
        while 1:
            line = input_file.readline()
            if not line:
                break
            # remove line returns
            line = line.strip('\r\n')
            # make sure it has useful data
            if (not "=" in line) or (line[0] == '%'):
                continue
            # split across equal sign
            line = line.split("=", 1)
            this_param = line[0].strip()
            this_value = line[1].strip()

            for case in switch(this_param):
                # integer values
                if case("OPT_ITERATIONS"): 
                    self._ConfigContent[this_param] = int(this_value)
                    break

                # float values
                if case("OPT_CONSTRAINT_SCALE"): pass
                if case("OPT_GRADIENT_FACTOR"): pass
                if case("OPT_RELAX_FACTOR"): pass
                if case("OPT_ACCURACY"): pass
                if case("OPT_BOUND_UPPER"): pass
                if case("OPT_BOUND_LOWER"): 
                    self._ConfigContent[this_param] = float(this_value)
                    break

                # string values
                if case("CONFIG_PRIMAL"): pass
                if case("CONFIG_ADJOINT"): pass
                if case("CONFIG_GEO"): pass
                if case("CONFIG_DEF"):
                    self._ConfigContent[this_param] = this_value
                    break
                    
                # OPT_CONSTRAINT values
                if case("OPT_CONSTRAINT"):
                   constraints_str = this_value.split(";") 
                   x = len(constraints_str)
                   constraint = []
                   for i in range(x):
                      constraint_value = [None] * 4
                      if constraints_str[i].find('*') != -1:
                         string = constraints_str[i].split('*')
                         scale = float(string[1].strip())
                         value = string[0].strip(') (')
                      else:
                         scale = None   
                         value = constraints_str[i].strip(') (')
                      if value.find('>') != -1:
                         value_s = value.split('>')
                         constraint_value[0] = value_s[0].strip()
                         constraint_value[1] = '>'
                         constraint_value[2] = float(value_s[1].strip()) 
                      elif value.find('<') != -1:
                         value_s = value.split('<')
                         constraint_value[0] = value_s[0].strip()         
                         constraint_value[1] = '<'
                         constraint_value[2] = float(value_s[1].strip())  
                      elif value.find('=') != -1:
                         value_s = value.split('=')
                         constraint_value[0] = value_s[0].strip()
                         constraint_value[1] = '='
                         constraint_value[2] = float(value_s[1].strip()) 
                      else:
                         raise SystemExit(constraint_value[1] + '  is an constraint option') 
                      constraint_value[3] = scale                  
                      constraint.append(constraint_value)   
                      
                   self._ConfigContent[this_param] = constraint
                   break
                    
                    
                if case():
                    print(this_param + " is an invalid option !")
            # end for

    # def dump()