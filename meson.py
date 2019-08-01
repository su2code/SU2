#!/usr/bin/env python3

## \file meson.py
#  \brief An extended meson script for setting up the environment and running meson
#  \author T. Albring
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
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


import sys, os, subprocess, shutil
from pathlib import Path

def setup_environment():

  # This information of the modules is used if projects was not cloned using git
  # The sha tag must be maintained manually to point to the correct commit
  sha_version_codi = '501dcf0305df147481630f20ce37c2e624fb351f'
  github_repo_codi = 'https://github.com/scicompkl/CoDiPack'
  sha_version_medi = 'a95a23ce7585905c3a731b28c1bb512028fc02bb'
  github_repo_medi = 'https://github.com/SciCompKL/MeDiPack'
  sha_version_meson = 'c904d3eefe2a01ca60027e2a5192e1f1c7ca5d9d'
  github_repo_meson = 'https://github.com/mesonbuild/meson'
  sha_version_ninja = 'e0bc2e5fd9036a31d507881e1383adde3672aaef'
  github_repo_ninja = 'https://github.com/ninja-build/ninja'

  medi_name = 'MeDiPack'
  codi_name = 'CoDiPack'
  meson_name = 'meson'
  ninja_name= 'ninja'

  alt_name_medi =   sys.path[0] + '/externals/medi'
  alt_name_codi =   sys.path[0] + '/externals/codi'
  alt_name_meson =  sys.path[0] + '/externals/meson'
  alt_name_ninja =  sys.path[0] + '/externals/ninja'

  # Some log and error files
  log = open( 'config.log', 'w' )
  err = open( 'config.err', 'w' )
  pkg_environ = os.environ

  submodule_check(codi_name, alt_name_codi, github_repo_codi, sha_version_codi, log, err)
  submodule_check(medi_name, alt_name_medi, github_repo_medi, sha_version_medi, log, err)
  submodule_check(meson_name, alt_name_meson, github_repo_meson, sha_version_meson, log, err)
  submodule_check(ninja_name, alt_name_ninja, github_repo_ninja, sha_version_ninja, log, err)
  
  log.close()
  err.close()


def submodule_check(name, alt_name, github_rep, sha_tag, log, err, update = False):
  try:
      status = submodule_status(alt_name, update)
      if status:
          print('Found correct version of ' + name + '.')

  except RuntimeError:
      if all([os.path.exists(alt_name), not os.path.exists(alt_name + '/' + sha_tag)]):
        print('Found an old or unspecified version of ' + name + ' in ' + alt_name)
        sys.exit()
      if not os.path.exists(alt_name):
        print('\ngit command failed (either git is not installed or this is not a git repository).')
        print('\nUsing fall-back method to initialize submodule ' + name)
        download_module(name, alt_name, github_rep, sha_tag, log, err)
      else:
        print('Found correct version of ' + name + ' in ' + alt_name + '.')


def submodule_status(path, update):

  try:
      status   = check_output('git submodule status ' + path).decode()
  except RuntimeError as e:
      print(e)
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
      print('Download of module ' + name + ' failed. See config.err for more information.')
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
      print('Extraction of module ' + name + ' failed. See config.err for more information.')
      sys.exit()

  # Rename folder and create a file to identify the version
  try:
      print('Creating identifier ...')
      subprocess.check_call('mv '+ name + '-' + commit_sha + ' ' + alt_name + ' && touch ' + alt_name + '/' + commit_sha, stdout = logfile, stderr = errorfile, shell = True)
  except subprocess.CalledProcessError:
      print('Renaming of module ' + name + ' failed. See config.err for more information.')
      sys.exit()

  # Remove archive
  subprocess.check_call('rm ' + commit_sha + '.zip', shell=True)

def check_output(cmd):
    std, err = subprocess.Popen([cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell = True).communicate()
    if err:
        raise RuntimeError(err)
    return std

def build_ninja():

  ninjapath = 'externals/ninja'
    
  try:
    p = subprocess.Popen([sys.path[0] + '/ninja', '--version'], stdout=subprocess.PIPE, shell=False)
  except:
    print("ninja executable not found. Building ...")
    p = subprocess.Popen(['./configure.py', '--bootstrap'], stdout=subprocess.PIPE, shell=False, cwd=ninjapath)
    out, err = p.communicate() 
    shutil.copy(ninjapath+'/ninja', '.')

if __name__ == '__main__':
  if sys.version_info[0] < 3:
    raise Exception("Script must be run using Python 3")
   
  # Set up the build environment, i.e. clone or download all submodules
  setup_environment()

  build_ninja()

  # Add paths for meson and ninja to environment
  if os.path.exists('externals/meson/mesonbuild'):
    sys.path.insert(0, str('externals/meson'))
  
  os.environ["NINJA"] = sys.path[0] + "/ninja"

  from mesonbuild import mesonmain

  sys.exit(mesonmain.main())
