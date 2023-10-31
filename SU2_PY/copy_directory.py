#!/usr/bin/env python
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)

import sys
import os
import shutil

if __name__ == "__main__":
    src = os.path.abspath(sys.argv[1])
    dst = os.path.abspath(sys.argv[2])
    shutil.copytree(src, dst)
