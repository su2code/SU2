#!/usr/bin/env python 

## \file change_version_number.py
#  \brief Python script for doing the continuous adjoint computation using the SU2 suite.
#  \author Aniket Aranake
#
# Stanford University Unstructured (SU2) Code
# Copyright (C) 2012 Aerospace Design Laboratory
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Run the script from the base directory (ie $SU2HOME). Grep will search directories recursively for matches in version number
import os,sys

oldvers = '2.0.5'
newvers = '2.0.5'

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
