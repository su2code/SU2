#!/usr/bin/env python2

## \file configure.py
#  \brief An extended configuration script.
#  \author T. Albring
#  \version 4.0.0 "Cardinal"
#
# SU2 Lead Developers: Dr. Francisco Palacios (francisco.palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#
# Copyright (C) 2012-2015 SU2, the open-source CFD code.
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

from optparse import OptionParser, BadOptionError
import sys,time, os, subprocess, os.path, glob, re, shutil, fileinput
import commands
from subprocess import call

# "Pass-through" option parsing -- an OptionParser that ignores
# unknown options and lets them pile up in the leftover argument
# list.  Useful to pass unknown arguments to the automake configure.
class PassThroughOptionParser(OptionParser):

    def _process_long_opt(self, rargs, values):
        try:
            OptionParser._process_long_opt(self, rargs, values)
        except BadOptionError, err:
            self.largs.append(err.opt_str)

    def _process_short_opts(self, rargs, values):
        try:
            OptionParser._process_short_opts(self, rargs, values)
        except BadOptionError, err:
            self.largs.append(err.opt_str)

def main():

    # Command Line Options

    parser = PassThroughOptionParser()

    parser.add_option("--enable-directdiff", action="store", type = "string",
                      help="Enable direct differentiation mode support", dest="directdiff_mode", default="")
    parser.add_option("--enable-reverse", action="store", type = "string",
                      help="Enable reverse mode support", dest="reverse_mode", default="")
    parser.add_option("--enable-mpi", action="store_true",
                      help="Enable mpi support", dest="mpi_enabled", default=False)
    parser.add_option("--disable-normal", action="store_true",
                      help="Disable normal mode support", dest="normal_mode", default=False)
    parser.add_option("-c" , "--check", action="store_true",
                      help="Check the source code for potential problems", dest="check", default=False)
    parser.add_option("-r" , "--replace", action="store_true",
                      help="Do a search and replace of necessary symbols. Does a back up source files.", dest="replace", default=False)
    parser.add_option("-d" , "--delete", action="store_true",
                      help="Removes the back up files.", dest="remove", default=False)
    parser.add_option("-v" , "--revert", action="store_true",
                      help="Revert files to original state.", dest="revert", default=False)
    (options, args)=parser.parse_args()

    options.directdiff_mode = options.directdiff_mode.upper()
    options.reverse_mode    = options.reverse_mode.upper()
    conf_environ = os.environ

    made_adolc = False
    made_codi = False

    header()

    modes =  {'NORMAL'     : not options.normal_mode == True,
              'DIRECTDIFF' : options.directdiff_mode,
              'REVERSE'    : options.reverse_mode}

    if not options.check:
        if any([modes["REVERSE"] == 'ADOLC', modes["DIRECTDIFF"] == 'ADOLC']):
            conf_environ, made_adolc = build_adolc(modes,options.mpi_enabled)
        if any([modes["REVERSE"] == 'CODI',  modes["DIRECTDIFF"] == 'CODI']):
            conf_environ, made_codi  = init_codi(modes,options.mpi_enabled)

        configure(args,
                  conf_environ,
                  options.mpi_enabled,
                  modes,
                  made_adolc,
                  made_codi)

    if options.check:
        prepare_source(options.replace, options.remove, options.revert)

def prepare_source(replace = False, remove = False, revert = False):

    # Directories containing the source code
    print 'Preparing source code ...'
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

    # Build the dictionaries for line and file excludes that
    # are defined in the exlude file 'preconf.exclude'.
    # Syntax:
    # PathTo/File[:Line1,Line2,...]
    if os.path.exists(exclude_file_name):
        print 'Reading \'' + exclude_file_name + '\' ...'
        with open(exclude_file_name, 'r') as exclude:
            for line in exclude:
                exclude_line  = line.split(':')
                exclude_file = exclude_line[0].rstrip()
                if len(exclude_line) > 1:
                    exclude_lines = exclude_line[1].split(',')
                    for index,item in enumerate(exclude_lines):
                        exclude_lines[index] = int(item.rstrip())
                    exclude_dic_lines[exclude_line[0].rstrip()] = exclude_lines
                else:
                    exclude_dic_files[exclude_line[0].rstrip()] = [-1]
    else:
        print 'Exclude file \'' + exclude_file_name + '\' not found. Checking all files.'

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
                        'MPI_Bcast'     : 'SU2_MPI::Bcast',
                        'MPI_Sendrecv'  : 'SU2_MPI::Sendrecv',
                        'MPI_Init'      : 'SU2_MPI::Init',
                        'MPI_Recv'      : 'SU2_MPI::Recv',
                        'sprintf'       : 'SPRINTF'}
    regex_cast_1 = re.compile(r'(^|[^\w|^\\])(int)(\s*\()')
    replacement_cast_1 = r'\1SU2_TYPE::Int\3'
    regex_cast_2 = re.compile(r'\(int\)\s*')

    logfile = open ('preconf.log','w')

    backup_ext = '.orig'
    print 'Checking for problems...'
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
                        print 'Cannot find backup file ' + file + backup_ext

                # Remove backup files if requested
                if all([not replace, remove]):
                    if os.path.isfile(file + backup_ext):
                        print 'Removing' + file + backup_ext
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
                        print new_line
                    if num_found > 0:
                        if replace:
                            print 'Solved ' + str(num_found) + ' potential problem(s) in ' + file + '.'
                            logfile.write('Solved ' + str(num_found) + ' potential problem(s) in ' + file + ':\n')
                        else:
                            print 'Found ' + str(num_found) + ' potential problem(s) in ' + file + '.'
                            logfile.write('Found ' + str(num_found) + ' potential problem(s) in ' + file + ':\n')
                        logfile.write( found_line )
                    else:
                        os.remove(file + backup_ext)

                    if not ignore_line == "":
                        print ignore_line.rstrip();
            else:
                print 'Ignoring file ' + file
    print '\nPlease check preconf.log to get more information about potential problems.'

def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

def find_all(text, dic):
    for i,j in dic.iteritems():
        if not text.find(i) == -1:
            return True
    return False

def init_codi(modes, mpi_support = False):

    modules_failed = True
    
    # This information of the modules is used if projects was not cloned using git
    # The sha tag must be maintained manually to point to the correct commit
    sha_version_codi = '3f668b0f5764b11ea592379f1842b41e6f1d489f'
    github_repo_codi = 'https://github.com/scicompkl/CoDiPack'
    sha_version_ampi = '8ab63005c646d5b6a2910521f9c0180ce7076f38'
    github_repo_ampi = 'https://github.com/michel2323/AdjointMPI'

    ampi_name = 'AdjointMPI'
    codi_name = 'CoDiPack'

    alt_name_ampi = 'adjointmpi'
    alt_name_codi = 'CoDi'

    # Some log and error files
    log = open( 'preconf.log', 'w' )
    err = open( 'preconf.err', 'w' )
    pkg_environ = os.environ

    # If project was cloned using git, we can use the 
    # submodule feature of git to initialize the packages
    if all([os.path.exists('.gitmodules')]):
        try:
            print "Updating submodules using git..."
            subprocess.check_call('git submodule update --init', stdout = log, stderr = err, shell = True)
            modules_failed = False
        except subprocess.CalledProcessError: 
            print 'Module initialization using git failed. Using fall-back option ...'
            modules_failed = True
    
    os.chdir('externals')
    
    # If modules still dont exists, use wget to download zip
    if all([not os.path.exists('CoDi/' + sha_version_codi), modules_failed]):
        if all([os.path.exists('CoDi'), not os.path.exists('CoDi/' + sha_version_codi)]):
            print 'Found an old or unspecified version of CoDiPack in externals/CoDi.\nPlease (re)move this folder and try again'
            sys.exit()

        download_module(codi_name, alt_name_codi, github_repo_codi, sha_version_codi, log, err)

         
    if all([not os.path.exists('adjointmpi/' + sha_version_ampi), mpi_support, modules_failed]):
        if all([os.path.exists('adjointmpi'), not os.path.exists('adjointmpi/' + sha_version_ampi)]):
            print 'Found an old or unspecified version of AdjointMPI in externals/adjointmpi.\nPlease (re)move this folder and try again.'
            sys.exit()

        download_module(ampi_name, alt_name_ampi, github_repo_ampi, sha_version_ampi, log, err)

    if mpi_support:
        os.chdir('adjointmpi')
        if not os.path.exists('libAMPI.a'):
             print '\nConfiguring and building AMPI...'
             subprocess.check_call('./configure CFLAGS=-O2 --prefix=$PWD && make', stdout = log, stderr = err, shell=True)
        os.chdir(os.pardir)

    os.chdir(os.pardir)

    return pkg_environ, True

def download_module(name, alt_name, git_repo, commit_sha, logfile, errorfile):

    print 'Initializing ' + name + '...'
   
    # Download package
    try:
        subprocess.check_call('wget ' + git_repo + '/archive/' + commit_sha + '.zip', stdout = logfile, stderr = errorfile, shell = True )
    except subprocess.CalledProcessError:
        print 'Download of module ' + name + ' failed. See preconf.err for more information.'
        sys.exit()
    
    # Extract zip archive
    try:
        subprocess.check_call('unzip -u ' + commit_sha + '.zip', stdout = logfile, stderr = errorfile, shell=True)
    except subprocess.CalledProcessError:
        print 'Extraction of module ' + name + ' failed. See preconf.err for more information.'
        sys.exit()

    # Rename folder and create a file to identify the version
    try:
        subprocess.check_call('mv '+ name + '-' + commit_sha + ' ' + alt_name + ' && touch ' + alt_name + '/' + commit_sha, stdout = logfile, stderr = errorfile, shell = True)
    except subprocess.CalledProcessError:
        print 'Renaming of module ' + name + ' failed. See preconf.err for more information.'
        sys.exit()

def build_adolc(modes, mpi_support = False):

    os.chdir('externals')

    pkg_environ = os.environ
    pkg_environ["PKG_CONFIG_PATH"] = pkg_environ.get("PKG_CONFIG_PATH","") + ":" + subprocess.check_output("pwd").rstrip() + "/adol-c/lib64/pkgconfig"

    configure_command = "./configure --prefix=$PWD"

    # Build adolc
    if any([modes['DIRECTDIFF'], all([modes['REVERSE'], not mpi_support])]):
        pkg_name = "adolc"
        pkg_version = "2.5.3-trunk"

        download_and_compile_adolc(configure_command, False, pkg_name, pkg_version, pkg_environ)

    # If necessary build adolc_ampi
    if all([mpi_support, modes['REVERSE']]):
        print 'MPI currently not supported when using ADOLC.'
        sys.exit()
        #ampi_path = subprocess.check_output("pwd").rstrip() +  "/AdjoinableMPI"
        #pkg_name = "adolc_ampi"
        #pkg_version = "2.5.3-trunk"

        #configure_command = configure_command + " --enable-ampi --with-ampi=" + ampi_path

        #download_and_compile_adolc(configure_command, True, pkg_name, pkg_version, pkg_environ)

    os.chdir(os.pardir)
    return pkg_environ, True

def download_and_compile_adolc(configure_command, ampi_needed, pkg_name, pkg_version, pkg_environ):

    ampi_clone = 'hg clone http://mercurial.mcs.anl.gov/ad/AdjoinableMPI'
    adolc_clone = 'git clone --depth 1 git@gitlab.com:adol-c/adol-c.git'
    adolc_download = 'wget https://gitlab.com/adol-c/adol-c/repository/archive.zip'

    pkg_call = "pkg-config --exists --print-errors ' " + pkg_name + " = " + pkg_version + " ' "

    check_adolc = True
    # pkg-config returns an error if package was not found
    try:
        subprocess.check_call(pkg_call, env = pkg_environ, shell=True)
    except subprocess.CalledProcessError:
        check_adolc = False


    if check_adolc:
        print "Found library " + pkg_name + "-" + pkg_version
    else:
        print "Since library " + pkg_name + "-" + pkg_version + " was not found, it will be downloaded and compiled"
        # Download and build AdjoinableMPI library if mpi support is enabled
        if ampi_needed:
            if not os.path.exists('AdjoinableMPI'):
                try:
                    subprocess.check_call(ampi_clone, shell=True)
                except subprocess.CalledProcessError:
                    print "hg clone failed"
                    pass

            os.chdir('AdjoinableMPI')
            # Generate configure script if it does not exist
            if not os.path.exists('configure'):
                subprocess.call('./autogen.sh', shell=True)
            # Call configure and install it
            subprocess.call('./configure --prefix=$PWD',shell=True)
            subprocess.call('make install',shell=True)
            os.chdir(os.pardir)

        # Download and build ADOL-C
        if not os.path.exists('adol-c'):
            try:
                subprocess.check_call(adolc_clone, shell=True)
            except subprocess.CalledProcessError:
                print "git clone failed, trying wget."
                subprocess.check_call(adolc_download, shell=True)
                subprocess.check_call('unzip archive.zip', shell=True)
                subprocess.check_call('mv adol-c.git adol-c', shell=True)
                os.remove('archive.zip')
                pass

        os.chdir('adol-c')

        # Generate configure script if it does not exist
        if not os.path.exists('configure'):
            subprocess.call('autoreconf -fi',shell=True)

        # Configure and build ADOL-C
        subprocess.call(configure_command ,shell=True)
        subprocess.call('make install',shell=True)

        os.chdir(os.pardir)


def configure(args,
              conf_environ,
              mpi_support,
              modes,
              made_adolc,
              made_codi):

    # Set the base command for running configure
    configure_base = '../configure'

    # Create a dictionary from the arguments
    argument_dict = dict(zip(args[::2],args[1::2]))

    # Set the default installation path (if not set with --prefix)
    argument_dict['--prefix'] = argument_dict.get('--prefix', subprocess.check_output('pwd').rstrip())

    # Add the arguments to the configure command
    for arg in argument_dict:
        configure_base = configure_base + " " + arg + "=" + argument_dict[arg]

    configure_mode = ''
    if mpi_support:
        configure_base = configure_base + ' --enable-mpi'

    build_dirs = ''

    # Create the commands for the different configurations and run configure
    for key in modes:
        if modes[key]:
            print '\nRunning configure in folder ' + key,
            if modes[key] == 'CODI':
                if key == 'DIRECTDIFF':
                    configure_mode = '--enable-codi-forward'
                if key == 'REVERSE':
                    configure_mode = '--enable-codi-reverse'
                print 'using ' + modes[key]
            elif modes[key] == 'ADOLC':
                if key == 'DIRECTDIFF':
                    configure_mode = '--enable-adolc-forward'
                if key == 'REVERSE':
                    configure_mode = '--enable-adolc-reverse'
                print 'using ' + modes[key]
            elif modes[key] == 'COMPLEX':
                configure_mode = '--enable-complex'
                print 'using ' + modes[key]
            else:
                configure_mode = ''
                print ''

            print '====================================================================='
            print  '\tCommand: ' + configure_base + ' ' + configure_mode
            run_command(configure_base + ' ' + configure_mode  , key, conf_environ)
            print '\tLog file written to ' + key +'/conf_'+key+'.log.'
            build_dirs += key + ' '

    write_makefile(build_dirs)

    print '\nPre-configuration Summary:\n' \
           '=====================================================================\n'\
          '\tConfiguration sets: '+ build_dirs + '\n'
    
    print '\tUse "make <install>" to compile (and install) all configured binaries:\n'
    if modes['NORMAL']:
        print '\tSU2_CFD            -> General solver for direct, cont. adjoint and linearized equations.\n' \
              '\tSU2_DOT            -> Gradient Projection Code.\n' \
              '\tSU2_DEF            -> Mesh Deformation Code.\n'  \
              '\tSU2_MSH            -> Mesh Adaption Code.\n' \
              '\tSU2_SOL            -> Solution Export Code.\n' \
              '\tSU2_GEO            -> Geometry Definition Code.\n'
    if modes['REVERSE']:
        print '\tSU2_CFD_REVERSE    -> AD-based Discrete Adjoint Solver.'
    if modes['DIRECTDIFF']:
        print '\tSU2_CFD_DIRECTDIFF -> Direct Differentation Mode.'

    print '\n'
    print  '\tPlease be sure to add the $SU2_HOME and $SU2_RUN environment variables,\n' \
           '\tand update your $PATH (and $PYTHONPATH if applicable) with $SU2_RUN.\n' \
           '\n' \
           '\tBased on the input to this configuration, add these lines to your .bashrc file: \n' \
           '\n' \
           '\texport SU2_RUN="'+argument_dict['--prefix']+'/bin"\n' \
           '\texport SU2_HOME="'+subprocess.check_output('pwd').rstrip()+'"\n' \
           '\texport PATH=$PATH:$SU2_RUN\n' \
           '\texport PYTHONPATH=$PYTHONPATH:$SU2_RUN\n'

def run_command(command, folder, env):
    log = 'conf_'+ folder+'.log'
    err = 'conf_'+ folder+'.err'

    if not os.path.exists(folder):
        os.mkdir(folder)

    os.chdir(folder)
    try:
        logfile = open(log, 'w')
        errfile = open(err, 'w')
        subprocess.check_call(command, env = env, stdout = logfile, stderr = errfile, shell=True)
    except subprocess.CalledProcessError:
        print '\nError in command: ' + command
        print 'See ' + folder +'/conf_'+folder+'.err for details.'
        sys.exit(1)

    os.chdir(os.pardir)

def write_makefile(build_dirs):
    print '\nCreating Makefile ...\n'
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

    print '-------------------------------------------------------------------------\n'\
          '|    ___ _   _ ___                                                      | \n'\
          '|   / __| | | |_  )   Release 4.0.0 \'Cardinal\'                          | \n'\
          '|   \__ \ |_| |/ /                                                      | \n'\
          '|   |___/\___//___|   Pre-configuration Script                          | \n'\
          '|                                                                       | \n'\
          '------------------------------------------------------------------------- \n'\
          '| SU2 Lead Dev.: Dr. Francisco Palacios (francisco.palacios@boeing.com).| \n'\
          '|                Dr. Thomas D. Economon (economon@stanford.edu).        | \n'\
          '------------------------------------------------------------------------- \n'\
          '| SU2 Developers:                                                       | \n'\
          '| - Prof. Juan J. Alonso\'s group at Stanford University.                | \n'\
          '| - Prof. Piero Colonna\'s group at Delft University of Technology.      | \n'\
          '| - Prof. Nicolas R. Gauger\'s group at Kaiserslautern U. of Technology. | \n'\
          '| - Prof. Alberto Guardone\'s group at Polytechnic University of Milan.  | \n'\
          '| - Prof. Rafael Palacios\' group at Imperial College London.            | \n'\
          '------------------------------------------------------------------------- \n'\
          '| Copyright (C) 2012-2015 SU2, the open-source CFD code.                | \n'\
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
          '------------------------------------------------------------------------- \n'



# -------------------------------------------------------------------
#  Run Main Program
# -------------------------------------------------------------------

# this is only accessed if running from command prompt
if __name__ == '__main__':
    main()
