from __future__ import annotations

import bz2
import gzip
import os
import subprocess
import tempfile
from typing import IO
from zipfile import ZipFile


def resolve_file_for_genome(file_path: str) -> tuple[bool, IO[str] | None, bool]:
    """Return (resolved, handle, iszip). The caller is responsible for closing handle."""
    tempdir = tempfile.mkdtemp()
    resolved = False
    iszip = False
    handle: IO[str] | None = None

    try:
        file_type = subprocess.check_output(["file", file_path]).decode("utf-8").replace("\n", "")
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False, None, False

    if "Zip archive" in file_type:
        iszip = True
        try:
            with ZipFile(file_path, "r") as zf:
                zf.extractall(tempdir)
                namelist = zf.namelist()
        except Exception:
            namelist = []

        candidates = sorted(namelist, key=lambda n: "deCODEme_" in n or "genome" in n, reverse=True)
        for name in candidates:
            f = os.path.join(tempdir, name)
            if not os.path.isfile(f):
                continue
            try:
                ft = subprocess.check_output(["file", f]).decode("utf-8").replace("\n", "")
            except Exception:
                continue
            h = get_handle(f, ft)
            if h:
                handle = h
                resolved = True
                if "deCODEme_" in name or "genome" in name:
                    break
    else:
        handle = get_handle(file_path, file_type)
        if handle:
            resolved = True

    return resolved, handle, iszip


def get_handle(file_path: str, file_type: str) -> IO[str] | None:
    """Open a file for text reading; supports plain text, gzip, and bzip2."""
    try:
        if "ASCII text" in file_type or "RSID sidtune PlaySID compatible" in file_type:
            return open(file_path)
        if "gzip" in file_type:
            return gzip.open(file_path, "rt")
        if "bzip" in file_type:
            return bz2.open(file_path, "rt")
        return open(file_path)
    except Exception:
        return None


def check_for_vcf(handle: IO[str] | None) -> bool:
    """Return True if handle looks like a VCF file (>=10 tab-separated columns, >=5000 data rows).

    Does not close the handle — caller is responsible for closing."""
    if handle is None:
        return False
    vcf_lines = 0
    non_vcf_lines = 0
    columns = 0
    for line in handle:
        if line.startswith("#"):
            columns = len(line.rstrip().split("\t"))
            continue
        if columns == 0:
            columns = len(line.rstrip().split("\t"))
        if columns < 10:
            return False
        if line == "\n":
            continue
        if len(line.rstrip().split("\t")) == columns:
            vcf_lines += 1
        else:
            non_vcf_lines += 1
    total = vcf_lines + non_vcf_lines
    return vcf_lines > 5000 and (total == 0 or non_vcf_lines / total < 0.1)
