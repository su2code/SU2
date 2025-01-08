#!/usr/bin/env python3

## \file meson.py
#  \brief An extended meson script for setting up the environment and running meson
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

import sys, os

sys.path.append(sys.path[0])
import preconfigure


if __name__ == "__main__":
    if sys.version_info[0] < 3:
        raise Exception("Script must be run using Python 3")

    # Preconfigure
    preconfigure.run(own_meson=True)

    # Add paths for meson and ninja to environment
    os.environ["NINJA"] = sys.path[0] + os.path.sep + "ninja"
    if os.name == "nt":
        os.environ["NINJA"] = os.environ["NINJA"] + ".exe"
    if os.path.exists(
        sys.path[0]
        + os.path.sep
        + "externals"
        + os.path.sep
        + "meson"
        + os.path.sep
        + "mesonbuild"
    ):
        sys.path.insert(
            0, str(sys.path[0] + os.path.sep + "externals" + os.path.sep + "meson")
        )

    from mesonbuild import mesonmain

    sys.exit(mesonmain.main())
