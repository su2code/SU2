#!/usr/bin/env python 

## \file change_version_number.py
#  \brief Python script for updating the version number of the SU2 suite.
#  \author A. Aranake
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# Run the script from the base directory (ie $SU2HOME). Grep will search directories recursively for matches in version number
import os,sys

oldvers = '4.3.0 "Cardinal"'
newvers = '5.0.0 "Raven"'

os.system('rm -rf version.txt')

# Grep flag cheatsheet:
# -I : Ignore binary files
# -F : Match exact pattern (instead of regular expressions)
# -w : Match whole word
# -r : search directory recursively
# -v : Omit search string (.svn omitted, line containing ISC is CGNS related)
os.system("grep -IFwr '%s' *|grep -vF '.svn' |grep -v ISC > version.txt"%oldvers)

# Create a list of files to adjust
filelist = []
f = open('version.txt','r')
for line in f.readlines():
  candidate = line.split(':')[0]
  if not candidate in filelist and candidate.find(sys.argv[0])<0:
    filelist.append(candidate)
f.close()
print filelist

# Prompt user before continuing 
yorn = ''
while(not yorn.lower()=='y'):
  yorn = raw_input('Replace %s with %s in the listed files? [Y/N]: '%(oldvers,newvers))
  if yorn.lower()=='n':
    print 'The file version.txt contains matches of oldvers'
    sys.exit()

# Loop through and correct all files
for fname in filelist:
  s = open(fname,'r').read()
  s_new = s.replace(oldvers,newvers)

  f = open(fname,'w')
  f.write(s_new)
  f.close()

os.system('rm -rf version.txt')
