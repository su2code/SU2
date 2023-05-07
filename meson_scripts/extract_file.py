#!/usr/bin/env python3
from pathlib import Path
import argparse
import zipfile
import tarfile


def extract_zip(fn: Path, outpath: Path, overwrite: bool = False):
    outpath = (
        Path(outpath).expanduser().resolve()
    )  # need .resolve() in case intermediate relative dir doesn’t exist
    if outpath.is_dir() and not overwrite:
        return
    fn = Path(fn).expanduser().resolve()
    with zipfile.ZipFile(fn) as z:
        z.extractall(str(outpath))


def extract_tar(fn: Path, outpath: Path, overwrite: bool = False):
    outpath = (
        Path(outpath).expanduser().resolve()
    )  # need .resolve() in case intermediate relative dir doesn’t exist
    print(outpath)
    if outpath.is_dir() and not overwrite:
        return
    fn = Path(fn).expanduser().resolve()
    if not fn.is_file():
        raise FileNotFoundError(fn)  # keep this, tarfile gives confusing error
    with tarfile.open(fn) as z:

        import os

        def is_within_directory(directory, target):

            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)

            prefix = os.path.commonprefix([abs_directory, abs_target])

            return prefix == abs_directory

        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):

            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")

            tar.extractall(path, members, numeric_owner=numeric_owner)

        safe_extract(z, str(outpath))


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
