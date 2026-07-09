# GenomePrep
To preprocess, quality control and prepare consumer DTC genomes for research

## Tutorial
https://genomeprep.readthedocs.io/en/latest/index.html

## Installation & Usage

### Install

```bash
pip install genomeprep
```

### Download reference data

Before processing, download the required reference files:

```bash
genomeprep download-data
```

Files are saved to `~/.cache/genomeprep/` by default. To use a custom location:

```bash
genomeprep download-data --dest /path/to/data
# or set the environment variable:
export GENOMEPREP_DATA_DIR=/path/to/data
```

> **Note:** The GRCh37 reference FASTA (~3 GB) must be supplied manually. Download
> `Homo_sapiens.GRCh37.75.dna.toplevel.fa` from
> ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/ and place it in
> your data directory.

### Process a genome file

```bash
genomeprep process my_genome.txt -d /path/to/data -o ./outputs
```

### Python API

```python
from genomeprep import process
from genomeprep.data import download, get_data_dir

result = process("my_genome.txt", datadir="/path/to/data")
print(result.status)  # "Processed", "Invalid", etc.
```

## Interactive server and bulk downloads
https://supfam.mrc-lmb.cam.ac.uk/GenomePrep/

## Incentive
The open-source GenomePrep tool-kit, developed on the goodwill of open genome data, addresses the problem of processing raw DTC DNA data in the context of the present: genotype arrays. The output of GenomePrep are DNA datafiles of homogenous formats (23andMe-like or vcf), which enable further research analysis. A single combined data-freeze of genomes that passed checks is also available in official website.

## Feature
1. Developed based on over 7,000 open genetic data
2. Automatic processing any inputs, identify genome from zip files, automatic parsing, converting various DTC genome formats into 23andMe-like format or VCF format. 
3. Automatic recognition of chip array version
4. Supply possibly problematic SNP position filter, stats developed from processing genetic data.

## Publication
C. Lu, B. Greshake Tzovaras, J. Gough, A survey of direct-to-consumer genotype data,and quality control tool (GenomePrep) for research, Computational and Structural Biotechnology Journal(2021), doi: https://doi.org/10.1016/j.csbj.2021.06.040
