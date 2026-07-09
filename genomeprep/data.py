from __future__ import annotations

import os
import pickle
import urllib.request
from dataclasses import dataclass
from typing import Any

import pyfaidx
from platformdirs import user_cache_dir

_DEFAULT_BASE_URL = "https://github.com/changlubio/GenomePrep/releases/download/v1.2.0"
_DATA_FILES = [
    "THE_LIST.dat",
    "badalleles.dat",
    "RS2GRCh37Orien_1.dat",
    "api.23andme.com",
]
_FASTA_FILENAME = "Homo_sapiens.GRCh37.75.dna.toplevel.fa"


@dataclass
class ReferenceData:
    api_23andme: set[str]
    the_list: dict[str, list[str]]
    reverse_snps: dict[str, Any]
    bad_allele_clusters: dict[str, Any]
    faidx: pyfaidx.Fasta


def get_data_dir(data_dir: str | None = None) -> str:
    """Return resolved data directory path.

    Priority: explicit arg > GENOMEPREP_DATA_DIR env var > ~/.cache/genomeprep
    """
    if data_dir is not None:
        return data_dir
    env = os.environ.get("GENOMEPREP_DATA_DIR")
    if env:
        return env
    return user_cache_dir("genomeprep")


def download(dest: str | None = None, base_url: str | None = None) -> str:
    """Download reference data files to *dest* directory. Returns the dest path.

    Downloads: THE_LIST.dat, badalleles.dat, RS2GRCh37Orien_1.dat, api.23andme.com
    Note: The FASTA file (~3 GB) must be supplied manually.
    """
    target = dest or get_data_dir()
    os.makedirs(target, exist_ok=True)
    url_base = (base_url or _DEFAULT_BASE_URL).rstrip("/")
    for filename in _DATA_FILES:
        url = f"{url_base}/{filename}"
        out_path = os.path.join(target, filename)
        if not os.path.exists(out_path):
            print(f"Downloading {filename}...")
            urllib.request.urlretrieve(url, out_path)
        else:
            print(f"Already present: {filename}")
    print(f"Reference data in: {target}")
    print(
        f"Note: {_FASTA_FILENAME} (~3 GB) must be supplied manually. See README."
    )
    return target


def _read_23andme_api(api_path: str) -> set[str]:
    api: set[str] = set()
    with open(api_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split()[2:]
            if len(parts) >= 2:
                api.add(":".join(parts))
    return api


def load_reference_data(datadir: str) -> ReferenceData:
    """Load all reference files from *datadir*.

    Raises FileNotFoundError with helpful messages for missing files.
    The FASTA file is never downloaded automatically — raises FileNotFoundError
    with instructions on where to obtain it.
    """
    list_path = os.path.join(datadir, "THE_LIST.dat")
    bad_path = os.path.join(datadir, "badalleles.dat")
    rev_path = os.path.join(datadir, "RS2GRCh37Orien_1.dat")
    api_path = os.path.join(datadir, "api.23andme.com")
    fasta_path = os.path.join(datadir, _FASTA_FILENAME)

    for p in (list_path, bad_path, rev_path, api_path):
        if not os.path.exists(p):
            raise FileNotFoundError(
                f"Reference file not found: {p}\n"
                "Run 'genomeprep download-data' to fetch reference files."
            )
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(
            f"GRCh37 FASTA not found: {fasta_path}\n"
            f"Download from ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/ "
            f"and place at {fasta_path}"
        )

    with open(list_path, "rb") as fh:
        the_list: dict[str, list[str]] = pickle.load(fh)
    with open(bad_path, "rb") as fh:
        bad_allele_clusters: dict[str, Any] = pickle.load(fh)
    with open(rev_path, "rb") as fh:
        reverse_snps: dict[str, Any] = pickle.load(fh)

    api_23andme = _read_23andme_api(api_path)
    faidx = pyfaidx.Fasta(fasta_path)

    return ReferenceData(
        api_23andme=api_23andme,
        the_list=the_list,
        reverse_snps=reverse_snps,
        bad_allele_clusters=bad_allele_clusters,
        faidx=faidx,
    )
