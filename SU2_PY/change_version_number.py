#!/usr/bin/env python

## \file change_version_number.py
#  \brief Python script for updating the version number of the SU2 suite.
#  \author A. Aranake
#  \version 8.0.0 "Harrier"
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
from optparse import OptionParser

# Run the script from the base directory (ie $SU2HOME). Grep will search directories recursively for matches in version number
import os, sys

parser = OptionParser()
parser.add_option(
    "-v", "--version", dest="version", help="the new version number", metavar="VERSION"
)
parser.add_option(
    "-r",
    "--releasename",
    dest="releasename",
    help="Name of the new release",
    metavar="RELEASENAME",
)
parser.add_option(
    "-y",
    action="store_true",
    dest="yes",
    help="Answer yes to all questions",
    metavar="YES",
)
(options, args) = parser.parse_args()

if not options.version:
    parser.error("new version number must be provided with -v option")

oldvers = '8.0.0 "Harrier"'
oldvers_q = r"8.0.0 \"Harrier\""
newvers = str(options.version) + ' "' + str(options.releasename) + '"'
newvers_q = str(options.version) + ' \\"' + str(options.releasename) + '\\"'
# oldvers = 'Copyright 2012-2023, SU2'
# oldvers_q = oldvers
# newvers = 'Copyright 2012-2023, SU2'
# newvers_q = newvers

if sys.version_info[0] > 2:
    # In PY3, raw_input is replaced with input.
    # For original input behaviour, just write eval(input())
    raw_input = input


if os.path.exists("version.txt"):
    os.remove("version.txt")

# Grep flag cheatsheet:
# -I : Ignore binary files
# -F : Match exact pattern (instead of regular expressions)
# -w : Match whole word
# -r : search directory recursively
# -v : Omit search string (.svn omitted, line containing ISC is CGNS related)

# TODO: replace with portable instructions. This works only on unix systems
os.system("grep -IFwr '%s' *|grep -vF '.svn' |grep -v ISC > version.txt" % oldvers)
os.system(
    "grep -IFwr '%s' --exclude='version.txt' *|grep -vF '.svn' |grep -v ISC >> version.txt"
    % oldvers_q
)

# Create a list of files to adjust
filelist = []
f = open("version.txt", "r")
for line in f.readlines():
    candidate = line.split(":")[0]
    if not candidate in filelist:
        filelist.append(candidate)
f.close()
print(filelist)

# Prompt user before continuing
yorn = ""
while not yorn.lower() == "y" and not options.yes:
    yorn = raw_input(
        "Replace %s with %s and %s with %s in the listed files? [Y/N]: "
        % (oldvers, newvers, oldvers_q, newvers_q)
    )
    if yorn.lower() == "n":
        print("The file version.txt contains matches of oldvers")
        sys.exit()

# Loop through and correct all files
for fname in filelist:
    with open(fname, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        for old, new in zip((oldvers, oldvers_q), (newvers, newvers_q)):
            # Avoid breaking the formating of some headers
            if old + "  " in line:
                n = len(new) - len(old)
                lines[i] = line.replace(old + " " * max(0, n), new + " " * max(0, -n))
            elif old in line:
                lines[i] = line.replace(old, new)
    with open(fname, "w") as f:
        f.writelines(lines)

os.system("rm -rf version.txt")
