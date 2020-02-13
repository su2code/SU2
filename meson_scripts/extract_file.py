#!/usr/bin/env python3
from pathlib import Path 
import argparse
import zipfile
import tarfile

def extract_zip(fn: Path, outpath: Path, overwrite: bool = False): 
  outpath = Path(outpath).expanduser().resolve() # need .resolve() in case intermediate relative dir doesn’t exist
  if outpath.is_dir() and not overwrite:
    return
  fn = Path(fn).expanduser().resolve()
  with zipfile.ZipFile(fn) as z:
    z.extractall(str(outpath))


def extract_tar(fn: Path, outpath: Path, overwrite: bool = False):
  outpath = Path(outpath).expanduser().resolve() # need .resolve() in case intermediate relative dir doesn’t exist
  print(outpath)
  if outpath.is_dir() and not overwrite:
    return
  fn = Path(fn).expanduser().resolve()
  if not fn.is_file():
    raise FileNotFoundError(fn)  # keep this, tarfile gives confusing error
  with tarfile.open(fn) as z:
    z.extractall(str(outpath))


if __name__ == "__main__":
  p = argparse.ArgumentParser()
  p.add_argument("infile", help="compressed file to extract")
  p.add_argument("outpath", help="path to extract into")
  P = p.parse_args()

  infile = Path(P.infile)
  if infile.suffix.lower() == ".zip":
    extract_zip(infile, P.outpath, True)
  elif infile.suffix.lower() in (".tar", ".gz", ".bz2", ".xz"):
    extract_tar(infile, P.outpath, True)
  else:
    raise ValueError("Not sure how to decompress {}".format(infile))
