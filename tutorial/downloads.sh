
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz

# # Fetch genomes from OpenSNP
# get list of all genomes
download_opensnp.py opensnp
# get list of all genomes and download 
download_opensnp.py opensnp --download
# get list of new genomes
download_opensnp.py opensnp_new --old opensnp.tab
# get list and download new genomes
download_opensnp.py opensnp_new --old opensnp.tab --download