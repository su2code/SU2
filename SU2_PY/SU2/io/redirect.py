#!/usr/bin/env python

## \file redirect.py
#  \brief python package for file redirection
#  \author T. Lukaczyk, F. Palacios
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

# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import os, sys, shutil, copy, glob
from .tools import add_suffix, make_link, expand_part

# -------------------------------------------------------------------
#  Output Redirection
# -------------------------------------------------------------------
# original source: http://stackoverflow.com/questions/6796492/python-temporarily-redirect-stdout-stderr
class output(object):
    """with SU2.io.redirect_output(stdout,stderr)

    Temporarily redirects sys.stdout and sys.stderr when used in
    a 'with' contextmanager

    Example:
    with SU2.io.redirect_output('stdout.txt','stderr.txt'):
        sys.stdout.write("standard out")
        sys.stderr.write("stanrard error")
        # code
    #: with output redirection

    Inputs:
        stdout - None, a filename, or a file stream
        stderr - None, a filename, or a file stream
    None will not redirect outptu

    """

    def __init__(self, stdout=None, stderr=None):

        _newout = False
        _newerr = False

        if isinstance(stdout, str):
            stdout = open(stdout, "a")
            _newout = True
        if isinstance(stderr, str):
            stderr = open(stderr, "a")
            _newerr = True

        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr
        self._newout = _newout
        self._newerr = _newerr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush()
        self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush()
        self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr

        if self._newout:
            self._stdout.close()
        if self._newerr:
            self._stderr.close()


#: class output()


# -------------------------------------------------------------------
#  Folder Redirection
# -------------------------------------------------------------------
class folder(object):
    """with SU2.io.redirect_folder(folder,pull,link,force) as push

    Temporarily redirects to a working folder, pulling
    and pushing needed files

    Example:

    folder = 'temp'
    pull   = ['file1.txt','file2.txt']
    link   = ['file3.big']
    force  = True

    # original path
    import os
    print(os.getcwd())

    # enter folder
    with SU2.io.redirect_folder(folder,pull,link,force) as push:
        print(os.getcwd())
        # code
        push.append('file4.txt')
    #: with folder redirection

    # returned to original path
    print(os.getcwd())

    Inputs:
        folder - working folder, relative or absolute
        pull   - list of files to pull (copy to working folder)
        link   - list of files to link (symbolic link in working folder)
        force  - True/False overwrite existing files in working folder

    Targets:
        push   - list of files to push (copy to originating path)

    Notes:
        push must be appended or extended, not overwritten
        links in Windows not supported, will simply copy
    """

    def __init__(self, folder, pull=None, link=None, force=True):
        """folder redirection initialization
        see help( folder ) for more info
        """

        if pull is None:
            pull = []
        if link is None:
            link = []

        if not isinstance(pull, list):
            pull = [pull]
        if not isinstance(link, list):
            link = [link]

        origin = os.getcwd()
        origin = os.path.abspath(origin).rstrip("/") + "/"
        folder = os.path.abspath(folder).rstrip("/") + "/"

        self.origin = origin
        self.folder = folder
        self.pull = copy.deepcopy(pull)
        self.push = []
        self.link = copy.deepcopy(link)
        self.force = force

    def __enter__(self):

        origin = self.origin  # absolute path
        folder = self.folder  # absolute path
        pull = self.pull
        push = self.push
        link = self.link
        force = self.force

        # check for no folder change
        if folder == origin:
            return []

        # relative folder path
        # relative = os.path.relpath(folder,origin)

        # check, make folder
        if not os.path.exists(folder):
            os.makedirs(folder)

        # copy pull files
        for name in pull:
            old_name = os.path.abspath(name)
            new_name = os.path.split(name)[-1]
            new_name = os.path.join(folder, new_name)
            if old_name == new_name:
                continue
            if os.path.exists(new_name):
                if force:
                    os.remove(new_name)
                else:
                    continue
            shutil.copy(old_name, new_name)

        # make links
        for name in link:
            old_name = os.path.abspath(name)
            new_name = os.path.split(name)[-1]
            new_name = os.path.join(folder, new_name)
            if old_name == new_name:
                continue
            if os.path.exists(new_name):
                if force:
                    os.remove(new_name)
                else:
                    continue
            make_link(old_name, new_name)

        # change directory
        os.chdir(folder)

        # return empty list to append with files to push to super folder
        return push

    def __exit__(self, exc_type, exc_value, traceback):

        origin = self.origin
        folder = self.folder
        push = self.push
        force = self.force

        # check for no folder change
        if folder == origin:
            return

        # move assets
        for name in push:

            old_name = os.path.abspath(name)
            name = os.path.split(name)[-1]
            new_name = os.path.join(origin, name)

            # links
            if os.path.islink(old_name):
                source = os.path.realpath(old_name)
                if source == new_name:
                    continue
                if os.path.exists(new_name):
                    if force:
                        os.remove(new_name)
                    else:
                        continue
                make_link(source, new_name)

            # moves
            else:
                if old_name == new_name:
                    continue
                if os.path.exists(new_name):
                    if force:
                        os.remove(new_name)
                    else:
                        continue
                shutil.move(old_name, new_name)

        # change directory
        os.chdir(origin)


#: class folder()
