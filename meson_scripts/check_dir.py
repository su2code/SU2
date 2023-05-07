#!/usr/bin/env python3
from pathlib import Path
import argparse
import sys
import os


def check_dir(fn: Path):
    return os.path.exists(fn)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("path", help="path to check")
    P = p.parse_args()
    sys.exit(not check_dir(P.path))
