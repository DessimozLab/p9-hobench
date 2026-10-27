#!/usr/bin/env python3
"""Minimal rsync substitute: copytree.py [--exclude PATTERN ...] SRC/ DST/ (fnmatch on names; merges into DST)."""
import fnmatch, shutil, sys
from pathlib import Path
args, excl = sys.argv[1:], ["__pycache__", "*.pyc"]
while args and args[0] == "--exclude":
    excl.append(args[1]); args = args[2:]
src, dst = Path(args[0]), Path(args[1])
def ign(d, names):
    return {n for n in names if any(fnmatch.fnmatch(n, p) for p in excl)}
shutil.copytree(src, dst, ignore=ign, dirs_exist_ok=True, symlinks=False)
