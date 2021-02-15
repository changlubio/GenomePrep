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
../process.py testgenome.zip -d ./datadir -o ./outputs -i outindex

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

# 3. Run non GRCh37 genome
## We use CrossMapy to convert genome build positions, which requires python 3.7
### set up liftover with CrossMap.py 
### ref: http://crossmap.sourceforge.net/
# Note: You don't need this if you are only processing GRCh37 data
wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
#### Note: default position after installation is /usr/local/bin/CrossMap.py; 
#### the path needs to be exported; check by running CrossMap.py --version

../process.py testgenome2.zip -d ./datadir -o ./outputs -i outindex
### Sample outputs
'''
Please check your ./outputs dir for a outindex.23andme.tohg19 file.
Here are some stats for the assembly conversion. (55.95709013938904 sec to convert assembly)
count_lines     959265
valid_snps      959265
count_not_found_in_ref  0
count_in_ref    770008
count_revd_in_ref       3253
final_snps      959265
Now we are processing the GRCh37 file ...
Please check your ./outputs dir for a outindex.23andme.checked and outindex.vcf file.
Here are some stats
badalleles      0
assembly        GRCh37
count_snps      959264
count_not_found_in_ref  0
valid_snps      959264
in_ref  0.8027
cluster_id      c4
cluster_mcov    1.0
afterbadallele  951627
indels  1089
refref  492695
multialleic     12
final_snps      457831
'''
