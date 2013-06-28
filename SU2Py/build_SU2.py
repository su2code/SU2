#!/usr/bin/python

## \file build_SU2.py
#  \brief Python script for automated building of the SU2 suite.
#  \author Thomas D. Economon
#  \version 1.0.2
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


# Import packages.
import os, time, sys, string
from optparse import OptionParser, OptionGroup

# Parse command line options.
parser=OptionParser()
parser.add_option("-o", "--os", dest="operatingsys",
                  help="specify build OS (macosx, redhat)", metavar="OS")
parser.add_option("--compiler", dest="compiler", default=False,
                  help="specify C++ compiler other than default", metavar="CXX")
parser.add_option("-p", "--parallel",
                  action="store_true", dest="parallel", default=False,
                  help="build parallel version")
parser.add_option("-c", "--clean",
                  action="store_true", dest="clean", default=False,
                  help="clean all directories")
group = OptionGroup(parser, "CGNS Options",
                    "Build with CGNS using using these options.")
group.add_option("--with-cgns",
                 action="store_true", dest="cgns", default=False,
                 help="build with CGNS support")
group.add_option("--cgns-inc-path", dest="inc_cgns", default="-I/usr/local/include",
                 help="specify path to CGNS header", metavar="-I/path/to/cgns/include")
group.add_option("--cgns-lib-path", dest="lib_cgns", default="-L/usr/local/lib -lcgns",
                 help="specify path to CGNS library", metavar="/path/to/libcgns.a")
parser.add_option_group(group)
group = OptionGroup(parser, "SU2 Suite Components",
                    "Build or clean individual components rather"
                    " than the entire suite by specifying one or"
                    " more of the following build flags.")
group.add_option("--SU2_CFD",
                 action="store_true", dest="SU2_CFD", default=False,
                 help="build SU2_CFD")
group.add_option("--SU2_DDC",
                 action="store_true", dest="SU2_DDC", default=False,
                 help="build SU2_DDC")
group.add_option("--SU2_GDC",
                 action="store_true", dest="SU2_GDC", default=False,
                 help="build SU2_GDC")
group.add_option("--SU2_GPC",
                 action="store_true", dest="SU2_GPC", default=False,
                 help="build SU2_GPC")
group.add_option("--SU2_MAC",
                 action="store_true", dest="SU2_MAC", default=False,
                 help="build SU2_MAC")
group.add_option("--SU2_MDC",
                 action="store_true", dest="SU2_MDC", default=False,
                 help="build SU2_MDC")
group.add_option("--SU2_PBC",
                 action="store_true", dest="SU2_PBC", default=False,
                 help="build SU2_PBC")
parser.add_option_group(group)
(options, args)=parser.parse_args()

# Perform some logic for the options.
if options.SU2_CFD or options.SU2_DDC or options.SU2_GDC or options.SU2_GPC or options.SU2_MAC or options.SU2_MDC or options.SU2_PBC:
  if options.SU2_CFD:
    build_cfd = True
  else:
    build_cfd = False
  if options.SU2_DDC:
    build_ddc = True
  else:
    build_ddc = False
  if options.SU2_GDC:
    build_gdc = True
  else:
    build_gdc = False
  if options.SU2_GPC:
    build_gpc = True
  else:
    build_gpc = False
  if options.SU2_MAC:
    build_mac = True
  else:
    build_mac = False
  if options.SU2_MDC:
    build_mdc = True
  else:
    build_mdc = False
  if options.SU2_PBC:
    build_pbc = True
  else:
    build_pbc = False
else:
  build_cfd = True
  build_ddc = True
  build_gdc = True
  build_gpc = True
  build_mac = True
  build_mdc = True
  build_pbc = True

# Set the SU2_HOME environment variable, if not set already.
try:
  su2_home = os.environ['SU2_HOME']
  os.chdir('../')
  if su2_home != os.getcwd():
    su2_home = os.getcwd()
    os.environ['SU2_HOME'] = su2_home
except:
  os.chdir('../')
  su2_home = os.getcwd()
  os.environ['SU2_HOME'] = su2_home

# Function for setting the compiler in the makefiles.
def set_compiler(compiler, makefile_path):
  try:
    makefile = open(makefile_path,'r')
    newfile  = open('makefile.new','w')
    for line in makefile:
      if line.find("=") > 0:
        iEquals = line.index("=")
        key     = (line[:iEquals].lower()).split()
        value   = (line[iEquals+1:]).split()
        # Strings are automatically converted to lower case.
        if key[0] == 'cxx':
          iStart = line.rindex("=")
          newfile.write(line[:iStart] + '= ' + ( "%s" % compiler ) + '\n')
        else:
          newfile.write(line)
      else:
        newfile.write(line)
    os.system('mv makefile.new ' + makefile_path)
  except:
    print '\n!!! Error: Could not set specified compiler !!!\n\n'
    exit(1)

# Function for handling CGNS changes to the makefile.
def set_cgns(makefile_path, cgns_inc_path, cgns_lib_path):
  try:
    makefile = open(makefile_path,'r')
    newfile  = open('makefile.new','w')
    for line in makefile:
      if line.find("-DNO_CGNS") > 0 and line.find("$") < 0:
        iStart   = line.index("-DNO_CGNS")
        iEnd     = iStart + len("-DNO_CGNS")
        newfile.write(line[:iStart] + line[iEnd:])
      elif line.find("INC_CGNS") > 0:
        iEquals = line.index("=")
        key     = (line[:iEquals].lower()).split()
        value   = (line[iEquals+1:]).split()
        # Strings are automatically converted to lower case.
        if key[0] == 'inc_cgns':
          iStart = line.rindex("=")
          newfile.write(line[:iStart] + '= ' + ( "%s" % cgns_inc_path ) + '\n')
        else:
          newfile.write(line)
      elif line.find("LIB_CGNS") > 0:
        iEquals = line.index("=")
        key     = (line[:iEquals].lower()).split()
        value   = (line[iEquals+1:]).split()
        # Strings are automatically converted to lower case.
        if key[0] == 'lib_cgns':
          iStart = line.rindex("=")
          newfile.write(line[:iStart] + '= ' + ( "%s" % cgns_lib_path ) + '\n')
        else:
          newfile.write(line)
      else:
        newfile.write(line)
    os.system('mv makefile.new ' + makefile_path)
  except:
    print '\n!!! Error: Could not set the CGNS requirements !!!\n\n'
    exit(1)

# Print message to the console.
print "\n"
print "-------------------------------------------------------------------------"
print "|                   SU2 Suite (Automated Build Code)                    |"
print "-------------------------------------------------------------------------"

# Build instructions based on operating system.
if options.operatingsys == "macosx":
  
  if options.clean:
    
    # Perform cleaning for Mac OS X.
    print '\n  Cleaning SU2 on Mac OS X.'
    print '  Setting environment variable SU2_HOME=' + su2_home
    
    # Clean SU2_CFD
    if build_cfd:
      print '\n  Cleaning SU2_CFD...\n'
      os.chdir(su2_home + '/SU2_CFD')
      os.system('cp ./config/makefile.darwin.gcc4.2.1.in ./makefile.in')
      os.system('make clean')
      os.chdir(su2_home + '/SU2Py')
      os.system('rm -f SU2_CFD')
    
    # Clean SU2_DDC
    if build_ddc:
      print '\n  Cleaning METIS & SU2_DDC...\n'
      os.chdir(su2_home + '/metis-4.0.3')
      os.system('make clean')
      os.chdir(su2_home + '/SU2_DDC')
      os.system('cp ./config/makefile.darwin.gcc4.2.1.in ./makefile.in')
      os.system('make clean')
      os.chdir(su2_home + '/SU2Py')
      os.system('rm -f SU2_DDC')
  
  elif options.parallel:
    
    # Build parallel version for Mac OS X.
    print '\n  Building a parallel version of SU2 on Mac OS X.'
    print '  Setting environment variable SU2_HOME=' + su2_home
    if options.cgns:
      print '  Building against the CGNS library.'
    
    # Build SU2_CFD
    if build_cfd:
      print '\n  Compiling SU2_CFD...\n'
      os.chdir(su2_home + '/SU2_CFD')
      os.system('cp ./config/makefile.darwin.gcc4.2.1.mpi.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        print '\n!!! Warning: SU2_CFD does not yet support parallel CGNS. Option ignored. !!!\n'
      os.system('make all')
    
    # Build SU2_DDC
    if build_ddc:
      print '\n  Compiling METIS & SU2_DDC...\n'
      os.chdir(su2_home + '/metis-4.0.3')
      os.system('make')
      os.chdir(su2_home + '/SU2_DDC')
      os.system('cp ./config/makefile.darwin.gcc4.2.1.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        print '\n!!! Warning: SU2_DDC does not yet support CGNS. Option ignored. !!!\n'
      os.system('make all')    
  
  else:
    
    # Build serial version for Mac OS X.
    print '\n  Building a serial version of SU2 on Mac OS X.'
    print '  Setting environment variable SU2_HOME=' + su2_home
    if options.cgns:
      print '  Building against the CGNS library.'
    
    # Build SU2_CFD
    if build_cfd:
      print '\n  Compiling SU2_CFD...\n'
      os.chdir(su2_home + '/SU2_CFD')
      os.system('cp ./config/makefile.darwin.gcc4.2.1.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        set_cgns('makefile.in',options.inc_cgns, options.lib_cgns)
      os.system('make all')
    
    # Build SU2_DDC
    if build_ddc:
      print '\n  Compiling METIS & SU2_DDC...\n'
      os.chdir(su2_home + '/metis-4.0.3')
      os.system('make')
      os.chdir(su2_home + '/SU2_DDC')
      os.system('cp ./config/makefile.darwin.gcc4.2.1.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        print '\n!!! Warning: SU2_DDC does not yet support CGNS. Option ignored. !!!\n'
      os.system('make all')

elif options.operatingsys == "redhat":
  
  if options.clean:
    
    # Perform cleaning for Linux Redhat.
    print '\n  Cleaning SU2 on Linux Redhat.'
    print '  Setting environment variable SU2_HOME=' + su2_home
    
    # Clean SU2_CFD
    if build_cfd:
      print '\n  Cleaning SU2_CFD...\n'
      os.chdir(su2_home + '/SU2_CFD')
      os.system('cp ./config/makefile.redhat.gcc4.1.2.in ./makefile.in')
      os.system('make clean')
      os.chdir(su2_home + '/SU2Py')
      os.system('rm -f SU2_CFD')
    
    # Clean SU2_DDC
    if build_ddc:
      print '\n  Cleaning METIS & SU2_DDC...\n'
      os.chdir(su2_home + '/metis-4.0.3')
      os.system('make clean')
      os.chdir(su2_home + '/SU2_DDC')
      os.system('cp ./config/makefile.redhat.gcc4.1.2.in ./makefile.in')
      os.system('make clean')
      os.chdir(su2_home + '/SU2Py')
      os.system('rm -f SU2_DDC')
  
  elif options.parallel:
    
    # Build parallel version for Linux Redhat.
    print '\n  Building a parallel version of SU2 on Linux Redhat.'
    print '  Setting environment variable SU2_HOME=' + su2_home
    if options.cgns:
      print '  Building against the CGNS library.'
    
    # Build SU2_CFD
    if build_cfd:
      print '\n  Compiling SU2_CFD...\n'
      os.chdir(su2_home + '/SU2_CFD')
      os.system('cp ./config/makefile.redhat.gcc4.1.2.mpi.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        print '\n!!! Warning: SU2_CFD does not yet support parallel CGNS. Option ignored. !!!\n'
      os.system('make all')
    
    # Build SU2_DDC
    if build_ddc:
      print '\n  Compiling METIS & SU2_DDC...\n'
      os.chdir(su2_home + '/metis-4.0.3')
      os.system('make')
      os.chdir(su2_home + '/SU2_DDC')
      os.system('cp ./config/makefile.redhat.gcc4.1.2.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        print '\n!!! Warning: SU2_DDC does not yet support CGNS. Option ignored. !!!\n'
      os.system('make all')
  
  else:
    
    # Build serial version for Linux Redhat.
    print '\n  Building a serial version of SU2 on Linux Redhat.'
    print '  Setting environment variable SU2_HOME=' + su2_home
    
    # Build SU2_CFD
    if build_cfd:
      print '\n  Compiling SU2_CFD...\n'
      os.chdir(su2_home + '/SU2_CFD')
      os.system('cp ./config/makefile.redhat.gcc4.1.2.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        set_cgns('makefile.in',options.inc_cgns, options.lib_cgns)
      os.system('make all')
    
    # Build SU2_DDC
    if build_ddc:
      print '\n  Compiling METIS & SU2_DDC...\n'
      os.chdir(su2_home + '/metis-4.0.3')
      os.system('make')
      os.chdir(su2_home + '/SU2_DDC')
      os.system('cp ./config/makefile.redhat.gcc4.1.2.in ./makefile.in')
      if options.compiler != False:
        set_compiler(options.compiler, 'makefile.in')
      if options.cgns:
        print '\n!!! Warning: SU2_DDC does not yet support CGNS. Option ignored. !!!\n'
      os.system('make all')

else:
  
  print '\n!!! Error: Unrecognized operating system !!!'
  print '      Options are: macosx, redhat           \n'
  exit(1)


# Build/Clean all other codes which don't depend on OS.

if options.clean:
  
  # Clean SU2_GDC
  if build_gdc:
    print '\n  Cleaning SU2_GDC...\n'
    os.chdir(su2_home + '/SU2_GDC')
    os.system('make clean')
    os.chdir(su2_home + '/SU2Py')
    os.system('rm -f SU2_GDC')
  
  # Clean SU2_GPC
  if build_gpc:
    print '\n  Cleaning SU2_GPC...\n'
    os.chdir(su2_home + '/SU2_GPC')
    os.system('make clean')
    os.chdir(su2_home + '/SU2Py')
    os.system('rm -f SU2_GPC')
  
  # Clean SU2_MAC
  if build_mac:
    print '\n  Cleaning SU2_MAC...\n'
    os.chdir(su2_home + '/SU2_MAC')
    os.system('make clean')
    os.chdir(su2_home + '/SU2Py')
    os.system('rm -f SU2_MAC')
  
  # Clean SU2_MDC
  if build_mdc:
    print '\n  Cleaning SU2_MDC...\n'
    os.chdir(su2_home + '/SU2_MDC')
    os.system('make clean')
    os.chdir(su2_home + '/SU2Py')
    os.system('rm -f SU2_MDC')
  
  # Clean SU2_PBC
  if build_pbc:
    print '\n  Cleaning SU2_PBC...\n'
    os.chdir(su2_home + '/SU2_PBC')
    os.system('make clean')
    os.chdir(su2_home + '/SU2Py')
    os.system('rm -f SU2_PBC')
  
  print '\n  Clean complete.\n\n'

else:
  
  # Build SU2_GDC
  if build_gdc:
    print '\n  Compiling SU2_GDC...\n'
    os.chdir(su2_home + '/SU2_GDC')
    os.system('cp ./config/makefile.in ./makefile.in')
    if options.compiler != False:
      set_compiler(options.compiler, 'makefile.in')
    if options.cgns:
      print '\n!!! Warning: SU2_GDC does not yet support CGNS. Option ignored. !!!\n'
    os.system('make all')
  
  # Build SU2_GPC
  if build_gpc:
    print '\n  Compiling SU2_GPC...\n'
    os.chdir(su2_home + '/SU2_GPC')
    os.system('cp ./config/makefile.in ./makefile.in')
    if options.compiler != False:
      set_compiler(options.compiler, 'makefile.in')
    if options.cgns:
      print '\n!!! Warning: SU2_GPC does not yet support CGNS. Option ignored. !!!\n'
    os.system('make all')
  
  # Build SU2_MAC
  if build_mac:
    print '\n  Compiling SU2_MAC...\n'
    os.chdir(su2_home + '/SU2_MAC')
    os.system('cp ./config/makefile.in ./makefile.in')
    if options.compiler != False:
      set_compiler(options.compiler, 'makefile.in')
    if options.cgns:
      print '\n!!! Warning: SU2_MAC does not yet support CGNS. Option ignored. !!!\n'
    os.system('make all')
  
  # Build SU2_MDC
  if build_mdc:
    print '\n  Compiling SU2_MDC...\n'
    os.chdir(su2_home + '/SU2_MDC')
    os.system('cp ./config/makefile.in ./makefile.in')
    if options.compiler != False:
      set_compiler(options.compiler, 'makefile.in')
    if options.cgns:
      print '\n!!! Warning: SU2_MDC does not yet support CGNS. Option ignored. !!!\n'
    os.system('make all')
  
  # Build SU2_PBC
  if build_pbc:
    print '\n  Compiling SU2_PBC...\n'
    os.chdir(su2_home + '/SU2_PBC')
    os.system('cp ./config/makefile.in ./makefile.in')
    if options.compiler != False:
      set_compiler(options.compiler, 'makefile.in')
    if options.cgns:
      print '\n!!! Warning: SU2_PBC does not yet support CGNS. Option ignored. !!!\n'
    os.system('make all')
  
  print '\n  Build complete. Products can be found in ' + su2_home + '/SU2Py.\n\n'

