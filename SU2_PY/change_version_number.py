#!/usr/bin/env python 

## \file change_version_number.py
#  \brief Python script for updating the version number of the SU2 suite.
#  \author A. Aranake
#  \version 7.0.0 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
# 
# The SU2 Project is maintained by the SU2 Foundation 
# (http://su2foundation.org)
#
# Copyright 2012-2019, SU2 Contributors (cf. AUTHORS.md)
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

# Run the script from the base directory (ie $SU2HOME). Grep will search directories recursively for matches in version number
import os,sys

#oldvers = '2012-2018'
#newvers = '2012-2019'
oldvers = '6.2.0 "Falcon"'
newvers = '7.0.0 "Blackbird"'

if sys.version_info[0] > 2:
  # In PY3, raw_input is replaced with input.
  # For original input behaviour, just write eval(input())
  raw_input = input


if os.path.exists('version.txt'):
  os.remove('version.txt')

# Grep flag cheatsheet:
# -I : Ignore binary files
# -F : Match exact pattern (instead of regular expressions)
# -w : Match whole word
# -r : search directory recursively
# -v : Omit search string (.svn omitted, line containing ISC is CGNS related)

#TODO: replace with portable instructions. This works only on unix systems
os.system("grep -IFwr '%s' *|grep -vF '.svn' |grep -v ISC > version.txt"%oldvers)

# Create a list of files to adjust
filelist = []
f = open('version.txt','r')
for line in f.readlines():
  candidate = line.split(':')[0]
  if not candidate in filelist and candidate.find(sys.argv[0])<0:
    filelist.append(candidate)
f.close()
print(filelist)

# Prompt user before continuing 
yorn = ''
while(not yorn.lower()=='y'):
  yorn = raw_input('Replace %s with %s in the listed files? [Y/N]: '%(oldvers,newvers))
  if yorn.lower()=='n':
    print('The file version.txt contains matches of oldvers')
    sys.exit()

# Loop through and correct all files
for fname in filelist:
  s = open(fname,'r').read()
  s_new = s.replace(oldvers,newvers)

  f = open(fname,'w')
  f.write(s_new)
  f.close()

os.system('rm -rf version.txt')
