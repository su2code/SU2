#!/usr/bin/env python3

## \file build_cantera.py
#  \brief Builds and installs Cantera as a static C++ library via its own
#         SCons build system, invoked automatically from meson.build so that
#         '-Denable-cantera=true' does not require a separate manual build step.
#  \author SU2 Contributors
#  \version 8.5.0 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2026, SU2 Contributors (cf. AUTHORS.md)
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

import argparse
import multiprocessing
import os
import shutil
import subprocess
import sys

# Only the pieces of Cantera SU2 actually links against (thermo/kinetics/transport)
# are built: no Python/Fortran interfaces, docs, or test suite. Cantera's own SCons
# script fetches its required submodules (fmt, yaml-cpp, sundials) on its own.
#
# Eigen is the one exception: SU2 already vendors its own Eigen copy (externals/eigen)
# and links Eigen-templated code from both SU2 and Cantera into the same executable.
# Letting Cantera vendor a second, different Eigen version is an ODR violation - the
# linker silently merges identically-mangled template instantiations from the two
# versions, which segfaults deep inside Eigen at runtime. Building against SU2's own
# Eigen (system_eigen=y + extra_inc_dirs below) keeps a single, consistent copy.
#
# Same version isn't sufficient on its own, though: Eigen's reduction/QR templates are
# header-only, so each translation unit compiles its own body, and the linker COMDAT-folds
# identically-mangled instantiations from different translation units into one - even
# across different -march targets. SU2 compiles with -march=<cpu-arch> (native by
# default); if Cantera's own compile used a different (or no) -march, the linker could
# keep either compiled body for a shared symbol, mixing SIMD-width assumptions and
# segfaulting. --cpu-arch keeps the two in sync (same fix already applied to CoolProp's
# CMAKE_CXX_FLAGS in meson.build).
CANTERA_SCONS_ARGS = [
    "layout=compact",
    "python_package=n",
    "f90_interface=n",
    "googletest=none",
    "hdf_support=n",
    "doxygen_docs=no",
    "sphinx_docs=no",
    "system_eigen=y",
    "system_fmt=n",
    "system_yamlcpp=n",
    "system_sundials=n",
    "system_highfive=n",
    "logging=warning",
]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", required=True, help="Path to the Cantera submodule checkout"
    )
    parser.add_argument(
        "--prefix", required=True, help="Install prefix for the built library"
    )
    parser.add_argument(
        "--stamp",
        required=True,
        help="File to create on a successful build, for meson tracking",
    )
    parser.add_argument(
        "--boost-inc-dir", default="", help="Explicit Boost header directory (optional)"
    )
    parser.add_argument(
        "--eigen-inc-dir",
        required=True,
        help="SU2's own Eigen include directory (externals/eigen)",
    )
    parser.add_argument(
        "--cpu-arch",
        default="",
        help="Value of SU2's -march= target (cpu-arch option), kept in sync so Eigen template code Cantera and SU2 both instantiate is not ODR-mixed across incompatible vectorization ISAs",
    )
    args = parser.parse_args()

    if not os.path.exists(os.path.join(args.source, "SConstruct")):
        sys.exit(
            "Cantera sources not found in '{}'.\n"
            "Run './preconfigure.py' (or 'git submodule update --init subprojects/cantera') "
            "before configuring with -Denable-cantera=true.".format(args.source)
        )

    scons = shutil.which("scons")
    if scons is None:
        sys.exit(
            "scons was not found on the PATH. Install it (e.g. 'pip install scons' or "
            "'apt install scons') to build Cantera from source."
        )

    # Cantera's SCons script probes for Eigen by searching for "eigen3/Eigen/Dense" first,
    # then plain "Eigen/Dense", across the compiler's whole include search path - not just
    # extra_inc_dirs. If the host has a system libeigen3-dev package, that "eigen3/..."
    # form gets found ahead of SU2's own (unprefixed) copy, defeating the point of pointing
    # this build at SU2's Eigen at all. Force the match onto SU2's copy by presenting it
    # under an "eigen3/" name via a symlink, in a directory searched before system ones.
    eigen_shim_dir = args.prefix.rstrip("/") + "_eigen_shim"
    os.makedirs(eigen_shim_dir, exist_ok=True)
    eigen_shim_link = os.path.join(eigen_shim_dir, "eigen3")
    if os.path.islink(eigen_shim_link) or os.path.exists(eigen_shim_link):
        os.remove(eigen_shim_link)
    os.symlink(args.eigen_inc_dir, eigen_shim_link)

    scons_command = (
        [scons, "build", "install", "prefix=" + args.prefix]
        + CANTERA_SCONS_ARGS
        + ["extra_inc_dirs=" + eigen_shim_dir]
        + ["-j{}".format(multiprocessing.cpu_count())]
    )
    if args.boost_inc_dir:
        scons_command.append("boost_inc_dir=" + args.boost_inc_dir)
    if args.cpu_arch:
        scons_command.append("cc_flags=-march=" + args.cpu_arch)

    result = subprocess.run(scons_command, cwd=args.source)
    if result.returncode != 0:
        sys.exit(
            "Cantera build failed (see output above). If SCons could not find Boost headers, "
            "set -Dcantera-boost-inc-dir=<path-to-boost-headers> and reconfigure, "
            "or install Boost development headers (e.g. the 'libboost-dev' package)."
        )

    # meson tracks this custom_target as always-stale, so re-running scons here is cheap
    # once Cantera itself is up to date; the stamp file just gives meson an output to track.
    open(args.stamp, "w").close()


if __name__ == "__main__":
    main()
