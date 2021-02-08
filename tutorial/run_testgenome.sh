#!/bin/sh

# Follow this to process testgenome.zip 
## 1. Download data files
### Make the data directory and download following files
mkdir datadir
cd datadir
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz
wget https://supfam.mrc-lmb.cam.ac.uk/GenomePrep/datadir/api.23andme.com
wget https://supfam.mrc-lmb.cam.ac.uk/GenomePrep/datadir/badalleles.dat
wget https://supfam.mrc-lmb.cam.ac.uk/GenomePrep/datadir/RS2GRCh37Orien_1.dat
wget https://supfam.mrc-lmb.cam.ac.uk/GenomePrep/datadir/THE_LIST.dat

## 2. Run test genome
### Assuming you are in the tutorial directory now
../process.py testgenome.zip -d datadir

### Sample outputs
'''
Please check your ./outputs dir for a outindex.23andme.checked and outindex.vcf file. 
 Here are some stats
count_snps      638480
count_not_found_in_ref  0
valid_snps      620207
in_ref  0.8892
cluster_id      v5
cluster_mcov    1.0
afterbadallele  613357
final_snps      159259
It took us 52.04032802581787 sec to process. See ya!
'''
