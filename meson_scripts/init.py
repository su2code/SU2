#!/usr/bin/env python3

## \file init.py
#  \brief Initializes necessary dependencies for SU2 either using git or it
#         fetches zip files.
#  \author T. Albring
#  \version 7.0.4 "Blackbird"
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

import sys, os, subprocess, shutil, urllib.request, zipfile, time

def remove_file(path, retries=3, sleep=0.1):
  for i in range(retries):
    try:
      os.remove(path)
    except OSError:
      time.sleep(sleep)
    else:
      break


def init_submodules(method = 'auto'):

  cur_dir = sys.path[0]

  # This information of the modules is used if projects was not cloned using git
  # The sha tag must be maintained manually to point to the correct commit
  sha_version_codi = '501dcf0305df147481630f20ce37c2e624fb351f'
  github_repo_codi = 'https://github.com/scicompkl/CoDiPack'
  sha_version_medi = 'edde14f9ac4026b72b1e130f61c0a78e8652afa5'
  github_repo_medi = 'https://github.com/SciCompKL/MeDiPack'
  sha_version_meson = '0435691e83fb7172e2a9635d2eb32d5521089916'
  github_repo_meson = 'https://github.com/mesonbuild/meson'
  sha_version_ninja = '2d15b04e411229cb902332957281622119025e77'
  github_repo_ninja = 'https://github.com/ninja-build/ninja'

  medi_name = 'MeDiPack'
  codi_name = 'CoDiPack'
  meson_name = 'meson'
  ninja_name= 'ninja'
  base_path = cur_dir + os.path.sep + 'externals' + os.path.sep 
  alt_name_medi = base_path + 'medi'
  alt_name_codi = base_path + 'codi'
  alt_name_meson =  base_path + 'meson'
  alt_name_ninja =  base_path + 'ninja'

  if method == 'auto':
    is_git = is_git_directory(cur_dir)
  elif method == 'git':
    is_git = True
  elif method == 'url':
    is_git = False
  else:
    print('Invalid method')
    sys.exit(1)

  # If directory was cloned using git, use submodule feature 
  # to check and initialize submodules if necessary
  if is_git:
    submodule_status(alt_name_codi, sha_version_codi)
    submodule_status(alt_name_medi, sha_version_medi)
    submodule_status(alt_name_meson, sha_version_meson)
    submodule_status(alt_name_ninja, sha_version_ninja)
  # Otherwise download the zip file from git
  else:
    download_module(codi_name, alt_name_codi, github_repo_codi, sha_version_codi)
    download_module(medi_name, alt_name_medi, github_repo_medi, sha_version_medi)
    download_module(meson_name, alt_name_meson, github_repo_meson, sha_version_meson)
    download_module(ninja_name, alt_name_ninja, github_repo_ninja, sha_version_ninja)

def is_git_directory(path = '.'):
  try:
     p = subprocess.call(["git", "branch"], stderr=subprocess.STDOUT, stdout=open(os.devnull, 'w'), cwd=path)
  except FileNotFoundError:
     print("git command not found. Using fall-back method to init submodules")
     return False
  except subprocess.CalledProcessError:
     print("Directory was not cloned using git. Using fall-back method to init submodules")
     return False
  return p == 0

def submodule_status(path, sha_commit):

  if not os.path.exists(path + os.path.sep + sha_commit):

    # Check the status of the submodule
    status = subprocess.run(['git', 'submodule','status', path], stdout=subprocess.PIPE, check = True, cwd = sys.path[0]).stdout.decode('utf-8')

    # The first character of the output indicates the status of the submodule
    # '+' : The submodule does not match the SHA-1 currently in the index of the repository
    # '-' : The submodule is not initialized
    # ' ' : Correct version of submodule is initialized
    status_indicator = status[0][0]

    
    if status_indicator == '+':
      # Write a warning that the sha tags do not match
      sys.stderr.write('WARNING: the currently checked out submodule commit in '
                        + path + ' does not match the SHA-1 found in the index.\n')
      sys.stderr.write('Use \'git submodule update --init '+ path + '\' to reset the module if necessary.\n')
    elif status_indicator == '-':
      # Initialize the submodule if necessary 
      print('Initialize submodule ' + path + ' using git ... ')
      subprocess.run(['git', 'submodule', 'update', '--init', path], check = True, cwd = sys.path[0])

    # Check that the SHA tag stored in this file matches the one stored in the git index
    cur_sha_commit = status[1:].split(' ')[0]
    if (cur_sha_commit != sha_commit):
      print('SHA-1 tag stored in index does not match SHA tag stored in this script.')
  
  

def download_module(name, alt_name, git_repo, commit_sha):
  
  # ZipFile does not preserve file permissions. 
  # This is a workaround for that problem:
  # https://stackoverflow.com/questions/39296101/python-zipfile-removes-execute-permissions-from-binaries
  class MyZipFile(zipfile.ZipFile):
    def _extract_member(self, member, targetpath, pwd):
      if not isinstance(member, zipfile.ZipInfo):
        member = self.getinfo(member)

      targetpath = super()._extract_member(member, targetpath, pwd)

      attr = member.external_attr >> 16
      if attr != 0:
        os.chmod(targetpath, attr)
      return targetpath

  if not os.path.exists(alt_name + os.path.sep + commit_sha):

    if os.path.exists(alt_name) and os.listdir(alt_name): 
      print('Directory ' + alt_name + ' is not empty')
      print('Maybe submodules are already cloned with git?')
      sys.exit(1) 
 
    else:
      print('Downloading ' + name + ' \'' + commit_sha + '\'')
    
      filename = commit_sha + '.zip'

      url = git_repo + '/archive/' + filename

      if not os.path.exists(sys.path[0] + os.path.sep + filename):
        try:
          urllib.request.urlretrieve (url, commit_sha + '.zip')
        except RuntimeError as e:
          print(e)
          print('Download of module ' + name + ' failed.')
          print('Get archive at ' + url)
          print('and place it in the source code root folder')
          print('Run meson.py again')
          sys.exit()
   
      # Unzip file 
      zipf = MyZipFile(sys.path[0] + os.path.sep + filename)
      zipf.extractall(sys.path[0] + os.path.sep + 'externals')
      
      # Remove directory if exists
      if os.path.exists(alt_name):
        os.rmdir(alt_name)

      os.rename(sys.path[0] + os.path.sep + 'externals' + os.path.sep + name + '-' + commit_sha, alt_name)

      # Delete zip file
      remove_file(sys.path[0] + os.path.sep + filename)

      # Create identifier
      f = open(alt_name + os.path.sep + commit_sha, 'w')
      f.close()


   
if __name__ == '__main__':
  if sys.version_info[0] < 3:
    raise Exception("Script must be run using Python 3")
   
  # Set up the build environment, i.e. clone or download all submodules
  init_submodules(sys.argv[1])

  sys.exit(0)

