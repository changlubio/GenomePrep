# GenomePrep
To preprocess, quality control and prepare consumer DTC genomes for research

## Tutorial
https://genomeprep.readthedocs.io/en/latest/index.html

## Interactive server and bulk downloads
https://supfam.org/GenomePrep/

## Incentive
The open-source GenomePrep tool-kit, developed on the goodwill of open genome data, addresses the problem of processing raw DTC DNA data in the context of the present: genotype arrays. The output of GenomePrep are DNA datafiles of homogenous formats (23andMe-like or vcf), which enable further research analysis. A single combined data-freeze of genomes that passed checks is also available in official website.

## Feature
1. Developed based on over 7,000 open genetic data
2. Automatic processing any inputs, identify genome from zip files, automatic parsing, converting various DTC genome formats into 23andMe-like format or VCF format. 
3. Automatic recognition of chip array version
4. Supply possibly problematic SNP position filter, stats developed from processing genetic data.

## Publication
C. Lu, B. Greshake Tzovaras, J. Gough, A survey of direct-to-consumer genotype data,and quality control tool (GenomePrep) for research, Computational and Structural Biotechnology Journal(2021), doi: https://doi.org/10.1016/j.csbj.2021.06.040

# Run GenomePrep locally

## Download 

Download datadir.tar.gz from Zenodo (https://zenodo.org/records/11408421), which contains dependencies for `bin/process.py`:
* api.23andme.com
* badalleles.dat
* RS2GRCh37Orien_1.dat
* THE_LIST.dat

To download all dependencies, including from public datasets
```{bash}
tar -xvf datadir.tar.gz
cd datadir
wget tp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
```

## Run GenomePrep on a sample genotype array

```
bin/process.py tutorial/testgenome.zip -d ./datadir -o ./outputs -i vcfindex
```

*We analyzed ~5000 OpenSNP genomes in 2020, the number is growing - see how many there are now here[https://opensnp.streamlit.app/]*
