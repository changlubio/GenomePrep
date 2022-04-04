#!/usr/bin/env python3

'''
Takes only 1 file as input and do: 
1. Quality Check -> reject or continue
2. Filter Badallele

Need:
THE_LIST.tab
api.23andme.com
'''

MIN_CLUSTER_COVERAGE = 0.8
MIN_COUNT_SNP = 10000

usable_genotypes = set(['A', 'C', 'G', 'T', 'D', 'I'])

refgenomes = {'GRCh37': 'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz'}

import os
import sys
import argparse
import gzip
import time
from collections import Counter
import pickle
import pyfaidx 
from datetime import datetime
import subprocess

from dtc_parser import reasons, DTCparser
from tools import resolve_file_for_genome, flip_genotype, output_as_23andme, get_snp_reference
from snps2vcf import directconvertSNPtoVCFformat

def read_23andme_api(api_23andme_file):
    ## Read 23andme api into memory
    api_23andme = set()
    with open(api_23andme_file, 'r') as inh:
        for line in inh:
            if line.startswith('#'):
                continue
            line = line.rstrip().split()[2:]
            api_23andme.add(':'.join(line))
    return api_23andme

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('-d', '--datadir', default='./datadir')
    parser.add_argument('-o', '--outdir', default='./outputs', help='directory where the quality checked 23andMe or VCF format files are going to be written to')
    args = parser.parse_args()

    subprocess.call('mkdir -p {}'.format(args.outdir), shell=True)

    start = time.time()
    with open(args.infile, 'r') as h:
        infile_list = [x.rstrip() for x in h.readlines()]

    for infile in infile_list:
        status, info = process(infile, args.datadir, args.outdir, os.path.basename(infile))
        print("{} {}sec done.".format(infile, time.time()-start))

def process(infile, datadir, outdir, outindex):

    output_23andme_path = os.path.join(outdir, '.'.join([outindex,'23andme.checked']))
    output_vcf_path = os.path.join(outdir, '.'.join([outindex, 'vcf']))

    badalleles_path = os.path.join(datadir, 'badalleles.dat')
    api_23andme_path = os.path.join(datadir, 'api.23andme.com')
    dna37_path = os.path.join(datadir, 'Homo_sapiens.GRCh37.75.dna.toplevel.fa')
    list_path = os.path.join(datadir, 'THE_LIST.dat')
    reverse_snps_path = os.path.join(datadir, 'RS2GRCh37Orien_1.dat')
    
    status = None
    info = None

    # 1. format parse
    resolved, handle, iszip = resolve_file_for_genome(infile)
    if not resolved:
        status = 'Invalid'
        if iszip:
            info = 'Cannot identify genetic data in your upload'
        else:
            info = 'Cannot identify genetic data in your uploaded ZIP file'
        return status, info
    elif not handle:
        status = 'error!'
        info = "Contact clu@mrc-lmb.cam.ac.uk Error is: we found your data but cannot get a handle"
        return status, info

    # 2. Now it's a flat file, read it with DTC parser first
    dtc_parser = DTCparser(handle, filename=infile)

    # assembly initial checked from file
    if dtc_parser.reject:
        if 'grch' in dtc_parser.reject_type:
            status = 'Invalid'
            info = 'Sorry web server only support processing data from human genome assembly GRCh37.'
            return status, info
        else:
            status = 'Invalid'
            info = 'Error in reading headlines. Please contact us.'
            return status, info
    
    # read the list, 23andme api, fadix, reversed SNPs
    with open(list_path, 'rb') as h:
        theLIST_clust = pickle.load(h)
    
    api_23andme = read_23andme_api(api_23andme_path)

    faidx = pyfaidx.Fasta(dna37_path)

    with open(reverse_snps_path, 'rb') as h:
        RS2GRCh37Orien_1 = pickle.load(h)

    # read all records
    reject, reject_type, store_snps, count_snps, duplicate_snps, count_not_found_in_ref, \
    valid_snps, in_ref, revd_in_ref, \
    in_rev_strand, rev_in_ref, rev_revd_in_ref, \
    not_in_list, cluster_id, cluster_mcov = dtc_read_records(dtc_parser, api_23andme, faidx, theLIST_clust, RS2GRCh37Orien_1)

    excl_duplicates = count_snps - len(duplicate_snps)

    if reject:
        status = reject_type
        info = {
            'badalleles': dtc_parser.total_badsnp,
            'count_snps': count_snps,
            'excl_duplicates': excl_duplicates,
            'count_not_found_in_ref': count_not_found_in_ref,
            'valid_snps': valid_snps,
            'in_ref': in_ref,
            'cluster_id': cluster_id,
            'cluster_mcov': cluster_mcov
        }
        if reject_type == 'invalid':
            if iszip:
                status = 'Invalid'
                info = 'Cannot identify a genetic data in your uploaded ZIP file'
            if count_snps == 0:
                status = 'Invalid'
                info = 'It does not look like a genetic data file. Please try again. Alternatively, contact us if you think it should be right!'
        elif reject_type == 'vcf':
            status = 'vcf'
            info = 'Currently this server is not automatically processing VCF files'
        return status, info
    
    ## Accepted
    # load badalleles
    with open(badalleles_path, 'rb') as h:
        badallelelist=pickle.load(h)
    badalleles = badallelelist[cluster_id]
    store_snps = filter_by_badalleles(store_snps, badalleles)

    genome_info = '''##fileformat=VCFv4.2
##fileDate={0:}
##source=GenomePrepV1.1(dtc_format={2:},cluseterversion={3:},clustermaxcov={4:},from_valid_snps={5:})
##reference={1:}
##phasing=unphased
##INFO=<ID=INDEL,Number=1,Type=String,Description="Insert or Deletion">
##INFO=<ID=GAV,Number=.,Type=String,Description="Genotyping Array Version">
##INFO=<ID=MA,Number=1,Type=String,Description="Genotypes, neither is the same as reference">
##FILTER=<ID=nref,Description="Neigher of the recorded calls is the same as reference">
##FILTER=<ID=bz3,Description="Badallele due to Zscore over 3 compared to 1000G">
##FILTER=<ID=indel,Description="Is INDEL">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'''.format(
        datetime.now().strftime('%y-%m-%d'),
        refgenomes['GRCh37'],
        dtc_parser.company, cluster_id, cluster_mcov, valid_snps
    )

    # filter and output as 23andme
    afterbadallele = output_as_23andme(store_snps, genome_info, output_23andme_path)
    # output to vcf
    final_snps, indels, refref, multialleic = directconvertSNPtoVCFformat(store_snps, genome_info, outindex, output_vcf_path)

    status = 'Processed'
    info = {
        'badalleles': dtc_parser.total_badsnp,
        'assembly': dtc_parser.build,
        'count_snps': count_snps,
        'excl_duplicates': excl_duplicates,
        'count_not_found_in_ref': count_not_found_in_ref,
        'valid_snps': valid_snps,
        'in_ref': in_ref,
        'cluster_id': cluster_id,
        'cluster_mcov': cluster_mcov,
        'afterbadallele':afterbadallele,
        'indels': indels,
        'refref': refref,
        'multialleic': multialleic,
        'final_snps': final_snps
    }
    return status, info

    
def dtc_read_records(dtc_parser, api_23andme, faidx, theLIST_clust, RS2GRCh37Orien_1):

    ## initial reject_type is empty
    reject = False
    reject_type = ""

    count_snps = 0
    count_not_found_in_ref = 0

    vcflines = 0
    nonvcflines = 0
    
    # valid_snps: exclude no calls;
    # also make sure that a reference genotype can be found, otherwise it may be from a different build
    valid_snps = 0
    store_snps = []
    snps_seen = set()
    duplicate_snps = []
    count_genotype_contains_reference = 0
    
    count_not_in_list = 0
    cluster_count = Counter()

    try:
        for record in dtc_parser:
            
            if record.badsnp:
                # print("badsnp encountered file %s, line is %s" % (dtc_parser.path, record.line), file=sys.stderr)
                dtc_parser.total_badsnp += 1
                # check vcf
                if len(record.line.split('\t')) > 9:
                    vcflines += 1
                else:
                    nonvcflines += 1
                continue
            
            ## special treatment for 23andme
            if dtc_parser.company == '23andme' and record.loci not in api_23andme:
                record.badsnp = True
                dtc_parser.total_badsnp += 1
                continue
            
            count_snps += 1

            # resolve duplicates
            snpid = '_'.join([record.chro, record.pos])
            if snpid in snps_seen:
                duplicate_snps.append('_'.join([record.chro, record.pos, record.genotype]))
                continue
            else:
                snps_seen.add(snpid)

            if len(set(record.genotype).intersection(usable_genotypes)) > 0:
                ref = get_snp_reference(record.loci, faidx)

                if ref:
                    if len(ref) == 1:
                        valid_snps += 1
                        # compare to reference
                        if ref in record.genotype:
                            count_genotype_contains_reference += 1
                        
                        # cluster check
                        GAV = []
                        if record.loci in theLIST_clust.keys():
                            for clust in theLIST_clust[record.loci]:
                                cluster_count[clust] += 1
                                GAV.append(clust)
                        else:
                            count_not_in_list += 1
                        
                        store_snps.append({
                            'rsid': record.rsid,
                            'chro': record.chro,
                            'pos': record.pos,
                            'call': record.genotype,
                            'ref': ref,
                            'gav': ','.join(GAV)
                        })

                else:
                    count_not_found_in_ref += 1

    except RecursionError:
        reject = True
        reject_type = 'invalid'

    # invalid if parsed lines < 10000
    if not reject:
        if count_snps < MIN_COUNT_SNP:
            reject = True
            reject_type = 'invalid'

    # reject if <21 chromosomes
    if not reject:
        chromo_set = set([x['chro'] for x in store_snps])
        if len(chromo_set) < 22 or valid_snps < 10000:
            reject = True
            reject_type = 'incomplete'
    
    # assembly check
    if not reject:
        in_ref = round(count_genotype_contains_reference/valid_snps, 4)
        if in_ref < 0.6 or count_not_found_in_ref > 0:
            reject = True
            reject_type = 'notgrch37'
    else:
        in_ref = -1
    
    # imputation check
    if not reject:
        if count_snps > 2000000:
            reject = True
            reject_type = 'imputation'
    
    # reverse strand check
    if not reject:
        in_rev_strand = 0
        count_rev_in_ref = 0
        count_rev_revd_in_ref = 0
        
        count_revd_in_ref = 0

        for snp in store_snps:
            newcall = flip_genotype(snp['call'])
            if snp['ref'] in newcall:
                count_revd_in_ref += 1

            snppos = ':'.join([snp['chro'], snp['pos']])
            if snppos in RS2GRCh37Orien_1.keys():
                in_rev_strand += 1
                
                if snp['ref'] in snp['call']:
                    count_rev_in_ref += 1
                
                if snp['ref'] in newcall:
                    count_rev_revd_in_ref += 1      
            
        rev_revd_in_ref = round(count_rev_revd_in_ref / (in_rev_strand+1), 4)
        rev_in_ref = round(count_rev_in_ref/(in_rev_strand+1),4)
        revd_in_ref = round(count_revd_in_ref / valid_snps, 4)

        if rev_revd_in_ref > 0.1 or in_ref < 0.75 or revd_in_ref > 0.08 or in_ref >0.98:
            reject = True
            reject_type = 'ref'
    else:
        in_rev_strand = rev_in_ref = rev_revd_in_ref = revd_in_ref = -1
    
    # clust info
    if not reject:
        not_in_list = round(count_not_in_list/valid_snps, 3)
        try:
            cluster_id = cluster_count.most_common(1)[0][0]
            cluster_mcov = round(cluster_count.most_common(1)[0][1] / valid_snps, 2)
            if cluster_mcov < MIN_CLUSTER_COVERAGE:
                cluster_id = 'all'
        except:
            cluster_id = cluster_mcov = -1
    else:
        cluster_id = cluster_mcov = not_in_list = -1

    if reject:
        if vcflines > 5000 and nonvcflines/(vcflines+nonvcflines) < 0.1:
            reject_type = 'vcf'
    
    return reject, reject_type, store_snps, count_snps, duplicate_snps, count_not_found_in_ref, \
    valid_snps, in_ref, revd_in_ref, \
    in_rev_strand, rev_in_ref, rev_revd_in_ref, \
    not_in_list, cluster_id, cluster_mcov

def filter_by_badalleles(store_snps, badalleles):
    new_snps = []
    for snp in store_snps:
        snppos = ':'.join([snp['chro'], snp['pos']])
        if snppos in badalleles:
            snp['ba'] = True
        else:
            snp['ba'] = False    
        new_snps.append(snp)
    return new_snps

if __name__ == '__main__':
    main()