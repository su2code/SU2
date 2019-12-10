#!/usr/bin/env python

## \file configure.py
#  \brief An extended configuration script.
#  \author T. Albring
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

from __future__ import print_function, division, absolute_import
from optparse import OptionParser, BadOptionError
import sys,time, os, subprocess, os.path, glob, re, shutil, fileinput
from subprocess import call

# "Pass-through" option parsing -- an OptionParser that ignores
# unknown options and lets them pile up in the leftover argument
# list.  Useful to pass unknown arguments to the automake configure.
class PassThroughOptionParser(OptionParser):

    def _process_long_opt(self, rargs, values):
        try:
            OptionParser._process_long_opt(self, rargs, values)
        except BadOptionError as err:
            self.largs.append(err.opt_str)

    def _process_short_opts(self, rargs, values):
        try:
            OptionParser._process_short_opts(self, rargs, values)
        except BadOptionError as err:
            self.largs.append(err.opt_str)

def main():

    # Command Line Options

    usage = './preconfigure.py [options]' \
            '\nNote: Options not listed below are passed to the automake configure.' \
            '\n      Compiler flags must be set with \'export CFLAGS=...\' or \'export CXXFLAGS=...\' ' \
            '\n      before calling this script.'

    parser = PassThroughOptionParser(usage = usage)

    parser.add_option("--enable-direct-diff", action="store_true",
                      help="Enable direct differentiation mode support", dest="directdiff", default=False)
    parser.add_option("--enable-autodiff", action="store_true",
                      help="Enable Automatic Differentiation support", dest="ad_support", default=False)
    parser.add_option("--with-ad", action="store",  type = "string",  help="AD Tool, CODI/ADOLC", default="CODI", dest="adtool")
    parser.add_option("--enable-mpi", action="store_true",
                      help="Enable mpi support", dest="mpi_enabled", default=False)
    parser.add_option("--enable-PY_WRAPPER", action="store_true",
                      help="Enable Python wrapper compilation", dest="py_wrapper_enabled", default=False)
    parser.add_option("--disable-tecio", action="store_true",
                      help="Disable Tecplot binary support", dest="tecio_disabled", default=False)
    parser.add_option("--disable-normal", action="store_true",
                      help="Disable normal mode support", dest="normal_mode", default=False)
    parser.add_option("-c" , "--check", action="store_true",
                      help="Check the source code for potential problems", dest="check", default=False)
    parser.add_option("-r" , "--replace", action="store_true",
                      help="Do a search and replace of necessary symbols. Creates back up of source files.", dest="replace", default=False)
    parser.add_option("-d" , "--delete", action="store_true",
                      help="Removes the back up files.", dest="remove", default=False)
    parser.add_option("-v" , "--revert", action="store_true",
                      help="Revert files to original state.", dest="revert", default=False)
    parser.add_option("-u", "--update", action="store_true",
                      help="Update and recompile submodules.", dest="update", default=False)
    (options, args)=parser.parse_args()

    options.adtool = options.adtool.upper()

    if options.directdiff == False:
        adtool_dd = ""
    else:
        adtool_dd = options.adtool
    if options.ad_support == False:
        adtool_da = ""
    else:
        adtool_da = options.adtool

    conf_environ = os.environ

    made_adolc = False
    made_codi = False

    header()

    modes =  {'SU2_BASE'     : not options.normal_mode == True,
              'SU2_DIRECTDIFF' : adtool_dd ,
              'SU2_AD'    : adtool_da }

    # Create a dictionary from the arguments
    argument_dict = dict(zip(args[::2],args[1::2]))

    # Set the default installation path (if not set with --prefix)
    argument_dict['--prefix'] = argument_dict.get('--prefix', os.getcwd().rstrip())


    if not options.check:
        if any([modes["SU2_AD"] == 'CODI',  modes["SU2_DIRECTDIFF"] == 'CODI']):
            conf_environ, made_codi  = init_codi(argument_dict,modes,options.mpi_enabled, options.update)

        configure(argument_dict,
                  conf_environ,
                  options.mpi_enabled,
                  options.py_wrapper_enabled,
                  options.tecio_disabled,
                  modes,
                  made_adolc,
                  made_codi)

    if options.check:
        prepare_source(options.replace, options.remove, options.revert)

def prepare_source(replace = False, remove = False, revert = False):

    # Directories containing the source code
    print('Preparing source code ...')
    dir_list = [ "Common",
                 "SU2_CFD",
                 "SU2_DEF",
                 "SU2_DOT",
                 "SU2_GEO",
                 "SU2_SOL",
                 "SU2_MSH"]

    file_list = ""

    exclude_dic_lines = {}
    exclude_dic_files = {}

    exclude_file_name = 'preconf.exclude'

#    # Build the dictionaries for line and file excludes that
#    # are defined in the exlude file 'preconf.exclude'.
#    # Syntax:
#    # PathTo/File[:Line1,Line2,...]
#    if os.path.exists(exclude_file_name):
#        print 'Reading \'' + exclude_file_name + '\' ...'
#        with open(exclude_file_name, 'r') as exclude:
#            for line in exclude:
#                exclude_line  = line.split(':')
#                exclude_file = exclude_line[0].rstrip()
#                if len(exclude_line) > 1:
#                    exclude_lines = exclude_line[1].split(',')
#                    for index,item in enumerate(exclude_lines):
#                        exclude_lines[index] = int(item.rstrip())
#                    exclude_dic_lines[exclude_line[0].rstrip()] = exclude_lines
#                else:
#                    exclude_dic_files[exclude_line[0].rstrip()] = [-1]
#    else:
#        print('Exclude file \'' + exclude_file_name + '\' not found. Checking all files.')


    # Hardcoded files that will be skipped
    exclude_dic_files = { 'Common/include/datatype_structure.hpp' : [-1],
                          'Common/include/datatype_structure.inl' : [-1],
                          'Common/include/mpi_structure.hpp' : [-1],
                          'Common/include/mpi_structure.inl' : [-1],
                          'Common/src/datatype_structure.cpp': [-1],
                          'Common/src/mpi_structure.cpp' : [-1] }

    str_double = 'double'

    regex_double = re.compile(r'(^|[^\w])('+str_double+')([^\w]|$)')
    replacement_double = r'\1su2double\3'
    simple_replacements = {'MPI_Reduce'    : 'SU2_MPI::Reduce',
                        'MPI_Allreduce' : 'SU2_MPI::Allreduce',
                        'MPI_Gather'    : 'SU2_MPI::Gather',
                        'MPI_Allgather' : 'SU2_MPI::Allgather',
                        'MPI_Isend'     : 'SU2_MPI::Isend',
                        'MPI_Irecv'     : 'SU2_MPI::Irecv',
                        'MPI_Send'      : 'SU2_MPI::Send',
                        'MPI_Wait'      : 'SU2_MPI::Wait',
                        'MPI_Waitall'   : 'SU2_MPI::Waitall',
                        'MPI_Waitany'   : 'SU2_MPI::Waitany',
                        'MPI_Bsend'     : 'SU2_MPI::Bsend' ,
                        'MPI_Bcast'     : 'SU2_MPI::Bcast',
                        'MPI_Sendrecv'  : 'SU2_MPI::Sendrecv',
                        'MPI_Init'      : 'SU2_MPI::Init',
                        'MPI_Recv'      : 'SU2_MPI::Recv',
                        'MPI_Comm_size' : 'SU2_MPI::Comm_size',
                        'MPI_Comm_rank' : 'SU2_MPI::Comm_rank',
                        'MPI_Init'      : 'SU2_MPI::Init',
                        'MPI_Barrier'   : 'SU2_MPI::Barrier',
                        'MPI_Abort'     : 'SU2_MPI::Abort',
                        'MPI_Request'   : 'SU2_MPI::Request',
                        'MPI_Get_count' : 'SU2_MPI::Get_count',
                        'MPI_Finalize'  : 'SU2_MPI::Finalize',
                        'MPI_Buffer_detach': 'SU2_MPI::Buffer_detach',
                        'MPI_Buffer_attach': 'SU2_MPI::Buffer_attach',
                        'MPI_Status'    : 'SU2_MPI::Status',
                        'sprintf'       : 'SPRINTF'}
    regex_cast_1 = re.compile(r'(^|[^\w|^\\])(int)(\s*\()')
    replacement_cast_1 = r'\1SU2_TYPE::Int\3'
    regex_cast_2 = re.compile(r'\(int\)\s*')

    logfile = open ('preconf.log','w')

    backup_ext = '.orig'
    print('Checking for problems...')
    # Test each source file for the occurrence of missing replacements
    # and print the respective lines.
    for dir in dir_list:
        file_list = glob.glob(dir+os.path.sep+'*[src,include]'+os.path.sep+'*[.cpp,.hpp,.inl]')
        for file in file_list:
            if not file in exclude_dic_files.keys():

                if all([not replace, revert]):
                    # Check if back up file exists
                    if os.path.isfile(file + backup_ext):
                        os.remove(file);
                        shutil.copy(file + backup_ext, file)
                    else:
                        print('Cannot find backup file ' + file + backup_ext)

                # Remove backup files if requested
                if all([not replace, remove]):
                    if os.path.isfile(file + backup_ext):
                        print('Removing' + file + backup_ext)
                        os.remove(file + backup_ext)

                if all([not remove, not revert]):
                    num_found = 0
                    found_line = ""
                    ignore_line = ""
                    new_line = ""
                    for line in fileinput.input(file, inplace = 1, backup = backup_ext):
                        new_line = line.rstrip('\n')
                        if any([re.findall(regex_double, line), find_all(line, simple_replacements), re.findall(regex_cast_1, line)]):
                            if not fileinput.lineno() in exclude_dic_lines.get(file,[]):
                                if replace:
                                    new_line = replace_all(new_line, simple_replacements)
                                    new_line = re.sub(regex_double, replacement_double, new_line)
                                    new_line = re.sub(regex_cast_1, replacement_cast_1, new_line)
                                    found_line = found_line + '\tLine ' + str(fileinput.lineno()) +': ' + line.rstrip() + '\n\t\t => ' + new_line.rstrip() + '\n'
                                else:
                                    found_line = found_line + '\tLine ' + str(fileinput.lineno()) +': ' + line.rstrip() + '\n'
                                num_found = num_found + 1
                            else:
                                ignore_line = ignore_line + 'Ignoring line ' + str(fileinput.lineno()) + ' in ' + file + ' (' + line.rstrip() + ')\n'
                        print(new_line)
                    if num_found > 0:
                        if replace:
                            print('Solved ' + str(num_found) + ' potential problem(s) in ' + file + '.')
                            logfile.write('Solved ' + str(num_found) + ' potential problem(s) in ' + file + ':\n')
                        else:
                            print('Found ' + str(num_found) + ' potential problem(s) in ' + file + '.')
                            logfile.write('Found ' + str(num_found) + ' potential problem(s) in ' + file + ':\n')
                        logfile.write( found_line )
                    else:
                        os.remove(file + backup_ext)

                    if not ignore_line == "":
                        print(ignore_line.rstrip())
            else:
                print('Ignoring file ' + file)
    print('\nPlease check preconf.log to get more information about potential problems.')

def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

def find_all(text, dic):
    for i,j in dic.iteritems():
        if not text.find(i) == -1:
            return True
    return False

def init_codi(argument_dict, modes, mpi_support = False, update = False):

    modules_failed = True
    
    # This information of the modules is used if projects was not cloned using git
    # The sha tag must be maintained manually to point to the correct commit
    sha_version_codi = 'bd4a639c2fe625a80946c8365bd2976a2868cf46'
    github_repo_codi = 'https://github.com/scicompkl/CoDiPack'
    sha_version_medi = '46a97e1d6e8fdd3cb42b06534cff6acad2a49693'
    github_repo_medi = 'https://github.com/SciCompKL/MeDiPack'

    medi_name = 'MeDiPack'
    codi_name = 'CoDiPack'

    alt_name_medi = 'externals/medi'
    alt_name_codi = 'externals/codi'

    # Some log and error files
    log = open( 'preconf.log', 'w' )
    err = open( 'preconf.err', 'w' )
    pkg_environ = os.environ

    codi_status = False
    ampi_status = False

    print("Checking the status of submodules")
    print('=====================================================================')
    # Remove modules if update is requested
    if update:
        if os.path.exists(alt_name_codi):
            print('Removing ' + alt_name_codi)
            shutil.rmtree(alt_name_codi)
        if os.path.exists(alt_name_medi):
            print('Removing ' + alt_name_medi)
            shutil.rmtree(alt_name_medi)

    submodule_check(codi_name, alt_name_codi, github_repo_codi, sha_version_codi, log, err, update)

    if mpi_support:
        submodule_check(medi_name, alt_name_medi, github_repo_medi, sha_version_medi, log, err, update)

    return pkg_environ, True

def submodule_check(name, alt_name, github_rep, sha_tag, log, err, update = False):

    try:
        status = submodule_status(alt_name, update)
        if status:
            print('Found correct version of ' + name + ' in ' + alt_name + '.')

    except RuntimeError:
        if all([os.path.exists(alt_name), not os.path.exists(alt_name + '/' + sha_tag)]):
          print('Found an old or unspecified version of ' + name + ' in ' + alt_name + '.\nUse -u to reset module.')
          sys.exit()
        if not os.path.exists(alt_name):
          print('\ngit command failed (either git is not installed or this is not a git repository).')
          print('\nUsing fall-back method to initialize submodule ' + name)
          download_module(name, alt_name, github_rep, sha_tag, log, err)
        else:
          print('Found correct version of ' + name + ' in ' + alt_name + '.')


def submodule_status(path, update):

    try:
        status = check_output('git submodule status ' + path).decode()
    except RuntimeError:
        raise RuntimeError

    status_indicator = status[0][0]

    if status_indicator == '+':
        sys.stderr.write('WARNING: the currently checked out submodule commit in ' + path + ' does not match the SHA-1 found in the index.\n')
        sys.stderr.write('Use \'git submodule update --init '+ path + '\' to reset the module if necessary.\n')
        return False
    elif any([status_indicator == '-', update]):
        print('Initialize submodule ' + path + ' using git ... ')
        subprocess.check_call('git submodule update --init ' + path, shell = True)

    return True

def download_module(name, alt_name, git_repo, commit_sha, logfile, errorfile):

    print('\nInitializing ' + name + ' \'' + commit_sha + '\'')
    print('=====================================================================')
    # Download package
    try:
        print('Downloading module from ' + git_repo)
        subprocess.check_call('wget -N ' + git_repo + '/archive/' + commit_sha + '.zip', stdout = logfile, stderr = errorfile, shell = True )
    except subprocess.CalledProcessError:
        print('Download of module ' + name + ' failed. See preconf.err for more information.')
        print('To download it manually, perform the following steps:')
        print('\t - Download the zip at \"' + git_repo + '/archive/' + commit_sha + '.zip\"')
        print('\t - Extract the archive to externals/' + alt_name)
        print('\t - Execute command \'touch externals/'+ alt_name + '/' + commit_sha + '\'')
        print('\t - Run preconfigure.py again')
        sys.exit()
    
    # Extract zip archive
    try:
        print('Extracting archive ...')
        subprocess.check_call('unzip -u ' + commit_sha + '.zip', stdout = logfile, stderr = errorfile, shell=True)
    except subprocess.CalledProcessError:
        print('Extraction of module ' + name + ' failed. See preconf.err for more information.')
        sys.exit()

    # Rename folder and create a file to identify the version
    try:
        print('Creating identifier ...')
        subprocess.check_call('mv '+ name + '-' + commit_sha + ' ' + alt_name + ' && touch ' + alt_name + '/' + commit_sha, stdout = logfile, stderr = errorfile, shell = True)
    except subprocess.CalledProcessError:
        print('Renaming of module ' + name + ' failed. See preconf.err for more information.')
        sys.exit()

    # Remove archive
    subprocess.check_call('rm ' + commit_sha + '.zip', shell=True)

def configure(argument_dict,
              conf_environ,
              mpi_support,
              py_wrapper,
              tecio,
              modes,
              made_adolc,
              made_codi):

    # Boostrap to generate Makefile.in
    bootstrap_command = './bootstrap'
    # Set the base command for running configure
    configure_base = '../configure'

    # Add the arguments to the configure command
    for arg in argument_dict:
        configure_base = configure_base + " " + arg + "=" + argument_dict[arg]

    configure_mode = ''
    if mpi_support:
        configure_base = configure_base + ' --enable-mpi'
    if py_wrapper:
        configure_base = configure_base + ' --enable-PY_WRAPPER'
    if tecio:
        configure_base = configure_base + ' --disable-tecio'

    build_dirs = ''
   
    print(  '\nPreparing build environment\n' \
            '=====================================================================')

    run_command(bootstrap_command, 'bootstrap.log', 'bootstrap.err', conf_environ)

    # Create the commands for the different configurations and run configure
    for key in modes:
        if modes[key]:
            print('\nRunning configure in folder ' + key + ' ', end = '')
            if modes[key] == 'CODI':
                if key == 'SU2_DIRECTDIFF':
                    configure_mode = '--enable-codi-forward'
                if key == 'SU2_AD':
                    configure_mode = '--enable-codi-reverse'
                print('using ' + modes[key])
            elif modes[key] == 'ADOLC':
                if key == 'SU2_DIRECTDIFF':
                    configure_mode = '--enable-adolc-forward'
                if key == 'SU2_AD':
                    configure_mode = '--enable-adolc-reverse'
                print('using ' + modes[key])
            elif modes[key] == 'COMPLEX':
                configure_mode = '--enable-complex'
                print('using ' + modes[key])
            else:
                configure_mode = ''
                print('')

            print('=====================================================================')
            log = os.getcwd().rstrip() + '/conf_'+ key+'.log'
            err = os.getcwd().rstrip() + '/conf_'+ key+'.err'

            if not os.path.exists(key):
                os.mkdir(key)

            os.chdir(key)
            run_command(configure_base + ' ' + configure_mode, log, err, conf_environ)
            os.chdir(os.pardir)
            build_dirs += key + ' '

    write_makefile(build_dirs)

    print('\nPre-configuration Summary:\n' \
           '=====================================================================\n'\
          '\tConfiguration sets: '+ build_dirs + '\n')

    print('\tUse "make <install>" to compile (and install) all configured binaries:\n')
    if modes['SU2_BASE']:
        print('\tSU2_CFD            -> General solver for direct, cont. adjoint and linearized equations.\n' \
              '\tSU2_DOT            -> Gradient Projection Code.\n' \
              '\tSU2_DEF            -> Mesh Deformation Code.\n'  \
              '\tSU2_MSH            -> Mesh Adaption Code.\n' \
              '\tSU2_SOL            -> Solution Export Code.\n' \
              '\tSU2_GEO            -> Geometry Definition Code.\n')
    if modes['SU2_AD']:
        print('\tSU2_CFD_AD         -> Discrete Adjoint Solver and general AD support.')
        print('\tSU2_DOT_AD         -> Mesh sensitivity computation and general AD support.')
    if modes['SU2_DIRECTDIFF']:
        print('\tSU2_CFD_DIRECTDIFF -> Direct Differentation Mode.')

    print('\n')
    print('\tPlease be sure to add the $SU2_HOME and $SU2_RUN environment variables,\n' \
           '\tand update your $PATH (and $PYTHONPATH if applicable) with $SU2_RUN.\n' \
           '\n' \
           '\tBased on the input to this configuration, add these lines to your .bashrc file: \n' \
           '\n' \
           '\texport SU2_RUN="'+argument_dict['--prefix']+'/bin"\n' \
           '\texport SU2_HOME="'+os.getcwd().rstrip()+'"\n' \
           '\texport PATH=$PATH:$SU2_RUN\n' \
           '\texport PYTHONPATH=$PYTHONPATH:$SU2_RUN\n')

def run_command(command, log, err, env):

    try:
        logfile = open(log, 'w')
        errfile = open(err, 'w')
        print('Command: ' + command)
        subprocess.check_call(command, env = env, stdout = logfile, stderr = errfile, shell=True)
        print('Logfile written to ' + log)
        logfile.close()
        errfile.close()
    except subprocess.CalledProcessError:
        errfile = open(err, 'r')
        print('\nThere was an error while running command \'' + command + '\'.')
        print('=== Error Log ===')
        print(errfile.read())
        errfile.close()
        sys.exit(1)


def check_output(cmd):
    std, err = subprocess.Popen([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True).communicate()
    if err:
        raise RuntimeError(err)
    return std
def write_makefile(build_dirs):
    print('\nCreating Makefile ...\n')
    makefile = open('Makefile', 'w')

    makefile.writelines(['# This file is auto-generated by preconfigure.py\n',
                         'SUBDIRS = '+ build_dirs + '\n',
                         'INSTALLDIRS = $(SUBDIRS:%=install-%)\n',
                         'CLEANDIRS = $(SUBDIRS:%=clean-%)\n',
                         '\n',
                         'subdirs: $(SUBDIRS)\n',
                         '\n',
                         '$(SUBDIRS):\n',
                         '\t$(MAKE) -C $@\n',
                         '\n',
                         'install: $(INSTALLDIRS)\n',
                         '$(INSTALLDIRS):\n',
                         '\t$(MAKE) -C $(@:install-%=%) install\n',
                         '\n',
                         'clean: $(CLEANDIRS)\n',
                         '$(CLEANDIRS):\n',
                         '\t$(MAKE) -C $(@:clean-%=%) clean\n',
                         '\n',
                         '.PHONY: subdirs $(SUBDIRS)\n',
                         '.PHONY: subdirs $(INSTALLDIRS)\n',
                         '.PHONY: subdirs $(CLEANDIRS)\n',
                         '.PHONY: install\n'])

    makefile.close()

def header():

    print('-------------------------------------------------------------------------\n'\
          '|    ___ _   _ ___                                                      | \n'\
          '|   / __| | | |_  )   Release 6.2.0 \'Falcon\'                            | \n'\
          '|   \__ \ |_| |/ /                                                      | \n'\
          '|   |___/\___//___|   Pre-configuration Script                          | \n'\
          '|                                                                       | \n'\
          '------------------------------------------------------------------------- \n'\
          '| The current SU2 release has been coordinated by the                   | \n'\
          '| SU2 International Developers Society <www.su2devsociety.org>          | \n'\
          '| with selected contributions from the open-source community.           | \n'\
          '------------------------------------------------------------------------- \n'\
          '| The main research teams contributing to the current release are:      | \n'\
          '| - Prof. Juan J. Alonso\'s group at Stanford University.                | \n'\
          '| - Prof. Piero Colonna\'s group at Delft University of Technology.      | \n'\
          '| - Prof. Nicolas R. Gauger\'s group at Kaiserslautern U. of Technology. | \n'\
          '| - Prof. Alberto Guardone\'s group at Polytechnic University of Milan.  | \n'\
          '| - Prof. Rafael Palacios\' group at Imperial College London.            | \n'\
          '| - Prof. Vincent Terrapon\'s group at the University of Liege.          | \n'\
          '| - Prof. Edwin van der Weide\'s group at the University of Twente.      | \n'\
          '| - Lab. of New Concepts in Aeronautics at Tech. Inst. of Aeronautics.  | \n'\
          '------------------------------------------------------------------------- \n'\
          '| Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,       | \n'\
          '|                      Tim Albring, and the SU2 contributors.           | \n'\
          '|                                                                       | \n'\
          '| SU2 is free software; you can redistribute it and/or                  | \n'\
          '| modify it under the terms of the GNU Lesser General Public            | \n'\
          '| License as published by the Free Software Foundation; either          | \n'\
          '| version 2.1 of the License, or (at your option) any later version.    | \n'\
          '|                                                                       | \n'\
          '| SU2 is distributed in the hope that it will be useful,                | \n'\
          '| but WITHOUT ANY WARRANTY; without even the implied warranty of        | \n'\
          '| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      | \n'\
          '| Lesser General Public License for more details.                       | \n'\
          '|                                                                       | \n'\
          '| You should have received a copy of the GNU Lesser General Public      | \n'\
          '| License along with SU2. If not, see <http://www.gnu.org/licenses/>.   | \n'\
          '------------------------------------------------------------------------- \n')



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
