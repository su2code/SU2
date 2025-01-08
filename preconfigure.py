#!/usr/bin/env python3

## \file preconfigure.py
#  \brief An preconfigure script for setting up the build environment
#  \author T. Albring and F. Poli
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

import sys, os, subprocess, shutil, urllib.request, zipfile
from pathlib import Path

sys.path.append(sys.path[0] + os.path.sep + "meson_scripts")
from init import init_submodules
from init import remove_file


def build_ninja():

    # If we are on windows, we don't need to compile ninja, we just download the executable
    if os.name == "nt":
        ninja_exe_url = "https://github.com/ninja-build/ninja/releases/download/v1.9.0/ninja-win.zip"

        # Try to execute ninja, if it fails, download .exe from github
        try:
            subprocess.run(
                [sys.path[0] + os.path.sep + "ninja.exe", "--version"],
                stdout=subprocess.PIPE,
            )
        except OSError:
            print("Downloading ninja ... ")
            try:
                urllib.request.urlretrieve(ninja_exe_url, "ninja-win.zip")
            except:
                print(e)
                print("Download of ninja executable failed.")
                print("Get archive at " + ninja_exe_url)
                print("extract ninja.exe in the source code root folder.")
                print("Run meson.py again.")
                sys.exit(1)

            zipf = zipfile.ZipFile(sys.path[0] + os.path.sep + "ninja-win.zip")
            zipf.extractall(sys.path[0])
            remove_file(sys.path[0] + os.path.sep + "ninja-win.zip")
    else:
        ninjapath = sys.path[0] + os.path.sep + "externals" + os.path.sep + "ninja"
        try:
            subprocess.run(
                [sys.path[0] + os.path.sep + "ninja", "--version"],
                stdout=subprocess.PIPE,
            )
        except OSError:
            print("ninja executable not found. Building ...")
            subprocess.run(["python3", "configure.py", "--bootstrap"], cwd=ninjapath)
            shutil.copy(ninjapath + os.path.sep + "ninja", ".")


def run(
    own_meson=False,
    own_codi=True,
    own_medi=True,
    own_opdi=True,
    own_mpp=True,
    own_cool=True,
    own_mel=True,
    own_mlpcpp=True,
):

    # Set up the build environment, i.e. clone or download submodules
    init_submodules(
        method="auto",
        own_meson=own_meson,
        own_codi=own_codi,
        own_medi=own_medi,
        own_opdi=own_opdi,
        own_mpp=own_mpp,
        own_cool=own_cool,
        own_mel=own_mel,
        own_mlpcpp=own_mlpcpp,
    )

    if own_meson:
        # Build ninja if it cannot be found
        build_ninja()

    # Leave a timestamp
    Path("su2preconfig.timestamp").touch()


if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Script must be run using Python 3")

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--with-own-meson",
        help="download own copies of Meson and Ninja",
        action="store_true",
    )
    parser.add_argument(
        "--no-codi", help="do not download own copy of CoDiPack", action="store_false"
    )
    parser.add_argument(
        "--no-medi", help="do not download own copy of MeDiPack", action="store_false"
    )
    parser.add_argument(
        "--no-opdi", help="do not download own copy of OpDiLib", action="store_false"
    )
    parser.add_argument(
        "--no-mpp", help="do not download own copy of Mutationpp", action="store_false"
    )
    parser.add_argument(
        "--no-coolprop",
        help="do not download own copy of CoolProp",
        action="store_false",
    )
    parser.add_argument(
        "--no-mel", help="do not download own copy of MEL", action="store_false"
    )
    parser.add_argument(
        "--no-mlpcpp",
        help="do not download copy of MLpCpp",
        action="store_false",
    )
    args = parser.parse_args()

    run(
        own_meson=args.with_own_meson,
        own_codi=args.no_codi,
        own_medi=args.no_medi,
        own_opdi=args.no_opdi,
        own_mpp=args.no_mpp,
        own_cool=args.no_coolprop,
        own_mel=args.no_mel,
        own_mlpcpp=args.no_mlpcpp,
    )
