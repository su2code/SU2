#!/usr/bin/env python 

## \file updateHistoryMap.py
#  \brief Python script for updating the historyMap.py file.
#  \author T. Albring
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
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
import os, pprint

su2_home = os.environ['SU2_HOME'] 

fileList = ['CFlowOutput.cpp', 
'CFlowIncOutput.cpp', 
'CFlowCompOutput.cpp',
'CHeatOutput.cpp',
'CFlowCompFEMOutput.cpp',
'CElasticityOutput.cpp',
'CAdjHeatOutput.cpp',
'CAdjFlowIncOutput.cpp',
'CAdjFlowCompOutput.cpp',
'CAdjElasticityOutput.cpp']

fileList = [os.path.join(su2_home, 'SU2_CFD/src/output/' + i) for i in fileList]

def parse_output(files):
    outputFields = dict()

    for file in files:
        print('Parsing ' + file)
        f = open(file,'r')
        while(1):
            s = f.readline().strip(' ')
            if not s:
                break
            if s.startswith('AddHistoryOutput('):
                s = s.replace('AddHistoryOutput', '').strip('()').split(',')
                curOutputField = dict()
                name = s[0].strip(' ()"\n;')
                curOutputField['HEADER']  = s[1].strip(' ()"\n;')
                curOutputField['GROUP']   = s[3].strip(' ()"\n;')
                curOutputField['DESCRIPTION'] = s[4].strip(' ()"\n;')
                if len(s) == 6:
                    curOutputField['TYPE'] = s[5].strip(' ()"\n;').split('::')[1]
                else:
                    curOutputField['TYPE'] = 'DEFAULT'
                outputFields[name] = curOutputField
        f.close()

    addedOutputFields = dict()

    for field in outputFields:
        if outputFields[field]['TYPE'] == 'COEFFICIENT':
            curOutputField = dict()
            name = 'D_' + field
            curOutputField['HEADER'] = 'd[' + outputFields[field]['HEADER'] + ']'
            curOutputField['GROUP'] = 'D_' + outputFields[field]['GROUP'] 
            curOutputField['TYPE'] = 'D_COEFFICIENT'
            curOutputField['DESCRIPTION'] = 'Derivative value'
            addedOutputFields[name] = curOutputField

    outputFields.update(addedOutputFields)
    f = open(os.path.join(su2_home) + 'SU2_PY/SU2/io/historyMap.py', 'w')
    f.write('history_header_map = ')
    pprint.pprint(outputFields, f)
    f.close()

parse_output(fileList)
