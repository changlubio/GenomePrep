# GenomePrep — Tested, Typed, Packaged, CI'd Library

**Date:** 2026-07-09  
**Status:** Approved

## Goal

Turn the GenomePrep research codebase into a properly packaged PyPI library with full type annotations, a pytest test suite, and GitHub Actions CI/CD. No new functionality — this is a structural lift.

## Context

The repo currently has three parallel copies of the same code (`bin/`, `bin/utils/`, `utils/`), no `pyproject.toml`, no tests, and no type hints. The Streamlit web app (`app/`) is out of scope for this effort.

---

## 1. Package Structure

All three copies of the code collapse into one authoritative `genomeprep/` package at the repo root. `bin/` and `utils/` are deleted.

```
genomeprep/
├── __init__.py     — process() + public re-exports (DTCparser, DTCrecord)
├── parser.py       — DTCparser, DTCrecord, COMPANY, reasons constants
├── io.py           — resolve_file_for_genome, get_handle, check_for_vcf
├── qc.py           — dtc_read_records, filter_by_badalleles, get_snp_reference, QC constants
├── convert.py      — flip_genotype, output_as_23andme, directconvertSNPtoVCFformat
├── data.py         — DataDir: download(), cache path, env-var override
└── cli.py          — argparse main(), wired up as console entry point
```

Each module has one concern. `__init__.py` orchestrates the pipeline by calling into these modules in the same order as the existing `bin/process.py`.

---

## 2. Packaging & Public API

Build backend: `hatchling`. Version lives in `genomeprep/__init__.py` as `__version__` and is read dynamically.

```toml
[project]
name = "genomeprep"
requires-python = ">=3.10"
dependencies = ["pyfaidx", "pysam", "platformdirs"]

[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "mypy", "ruff"]

[project.scripts]
genomeprep = "genomeprep.cli:main"
```

### CLI subcommands

```
genomeprep process <infile> [-d DATADIR] [-o OUTDIR] [-i INDEX]
genomeprep download-data [--dest DIR]
```

### Python API

```python
from genomeprep import process, DTCparser, DTCrecord
from genomeprep.data import download, get_data_dir
```

### Reference data handling

`data.py` resolves the data directory in priority order:

1. Explicit `data_dir` argument to `process()`
2. `GENOMEPREP_DATA_DIR` environment variable
3. `~/.cache/genomeprep/` (default cache)

`genomeprep download-data` / `genomeprep.data.download()` fetches `THE_LIST.dat`, `badalleles.dat`, `RS2GRCh37Orien_1.dat`, and `api.23andme.com` from GitHub Releases. The GRCh37 FASTA (~3 GB) is documented separately and must be supplied by the user. Reference files are versioned on GitHub Releases and pinned to the library version.

---

## 3. Type Annotations

`mypy --strict` enforced. Key additions:

```python
Company = Literal["23andme", "ancestrydna", "ftdna", "myheritage", "decodeme", "livingdnacsv", "all"]

RejectType = Literal["invalid", "small", "incomplete", "grch36", "grch38", "notgrch37", "vcf", "zip", "ref", "imputation"]

class SNPRecord(TypedDict):
    rsid: str
    chro: str
    pos: str
    call: str
    ref: str
    gav: str
    ba: NotRequired[bool]  # set by filter_by_badalleles; absent before filtering

@dataclass
class ProcessResult:
    status: Literal["Processed", "Invalid", "vcf"] | str
    info: dict[str, int | float | str] | str
```

`DTCparser` and `DTCrecord` get fully annotated attributes. `resolve_file_for_genome` returns `tuple[bool, IO[str] | None, bool]`. No runtime behaviour changes.

---

## 4. Tests

```
tests/
├── conftest.py              — shared fixtures: temp dirs, synthetic genome lines, mock data dir
├── fixtures/
│   ├── 23andme_sample.txt
│   ├── ancestrydna_sample.txt
│   ├── ftdna_sample.csv
│   └── myheritage_sample.csv
├── test_parser.py           — header detection, DTCrecord per company, badsnp logic
├── test_io.py               — resolve_file_for_genome (zip/gzip/plain); uses tutorial/testgenome.zip
├── test_convert.py          — flip_genotype, output_as_23andme, directconvertSNPtoVCFformat
├── test_qc.py               — filter_by_badalleles; dtc_read_records with mocked pyfaidx + pickles
├── test_data.py             — download/cache with mocked HTTP
└── test_cli.py              — subprocess smoke tests for both CLI subcommands
```

`pyfaidx.Fasta` and pickle loads are mocked in `test_qc.py` — no 3 GB genome required to run the suite. Coverage gate: 80% enforced via `--cov-fail-under=80`.

---

## 5. CI/CD

### `ci.yml` — every push and PR

| Job | Command |
|-----|---------|
| lint | `ruff check . && ruff format --check .` |
| typecheck | `mypy genomeprep/ --strict` |
| test | `pytest --cov=genomeprep --cov-fail-under=80` |
| test matrix | Python 3.10, 3.11, 3.12 |

### `publish.yml` — on `v*` tag push

| Step | Detail |
|------|--------|
| build | `python -m build` |
| publish | PyPI via OIDC trusted publisher (no stored API tokens) |

Releasing: `git tag v1.2.0 && git push --tags` triggers the full build and publish pipeline.

---

## Out of Scope

- Streamlit web app (`app/`) — stays as a separate deployment
- VCF input processing (currently rejected with a message, no change)
- New DTC provider support
- Support for Python < 3.10
