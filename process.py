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

from utils.dtc_parser import reasons, DTCparser
from utils.tools import get_handle, flip_genotype, check_for_vcf
from utils.snps2vcf import directconvertSNPtoVCFformat

def get_snp_reference(snp, faidx):
    try:
        chromosome, position = snp.split(':')
    except:
        print("Can't parse snp: %s" % snp, file=sys.stderr)
        return None
    ## see if we can get the reference
    try:
        ref = faidx[chromosome][int(position) - 1].seq
    except KeyError:
        print("No reference found for: %s:%s" % (chromosome, position), file=sys.stderr)
        return None
    return ref

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
    parser.add_argument('-i', '--outindex', default='outindex', help='the index/name for the output files, it will also be used as the index in the VCF-format output')
    args = parser.parse_args()

    subprocess.call('mkdir -p {}'.format(args.outdir), shell=True)

    status, info = process(args.infile, args.datadir, args.outdir, args.outindex)

    if status != 'Processed':
        # if isinstance(info, str):
        # it's rejected
        print('Your file is not a common genotyping data. Please run other scripts specified in the README.MD to process. (It will be automated soon.)')
        print('status:\t{}\ninfo:\t{}'.format(status, info))
    
    else:
        print("Please check your {0:} dir for a {1:}.23andme.checked and {1:}.vcf file. \n Here are some stats".format(
            args.outdir, args.outindex))
        for k, i in info.items():
            print('\t'.join(map(str, [k,i])))
    
    print("See ya!")

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

    start = time.time()

    dtc_parser = DTCparser(infile)

    # no need to do anything if we can't parse
    if dtc_parser.reject:
        status = dtc_parser.reject_type
        info = reasons[dtc_parser.reject_type]
        return status, info

    is_vcf = check_for_vcf(dtc_parser.path)
    if is_vcf:
        status = 'vcf'
        info = 'Processing VCF file'
        return status, info

    # read the list, 23andme api, fadix, reversed SNPs
    with open(list_path, 'rb') as h:
        theLIST_clust = pickle.load(h)
    
    api_23andme = read_23andme_api(api_23andme_path)

    faidx = pyfaidx.Fasta(dna37_path)

    with open(reverse_snps_path, 'rb') as h:
        RS2GRCh37Orien_1 = pickle.load(h)

    # read all records
    reject, reject_type, store_snps, count_snps, count_not_found_in_ref, \
    valid_snps, in_ref, revd_in_ref, \
    in_rev_strand, rev_in_ref, rev_revd_in_ref, \
    not_in_list, cluster_id, cluster_mcov = dtc_read_records(dtc_parser, api_23andme, faidx, theLIST_clust, RS2GRCh37Orien_1)

    if reject:
        status = reject_type
        info = {
            'count_snps': count_snps,
            'count_not_found_in_ref': count_not_found_in_ref,
            'valid_snps': valid_snps,
            'in_ref': in_ref,
            'cluster_id': cluster_id,
            'cluster_mcov': cluster_mcov
        }
        return status, info
    
    ## Accepted
    # load badalleles
    with open(badalleles_path, 'rb') as h:
        badallelelist=pickle.load(h)
    badalleles = badallelelist[cluster_id]
    store_snps = filter_by_badalleles(store_snps, badalleles)

    genome_info = "##GenomePrep at %s\n##company=%s,usebadallele=%s,maxcov=%s,from_valid_snps=%s\n" % (datetime.now().strftime('%d-%m-%y'), dtc_parser.company, cluster_id, cluster_mcov, valid_snps)
    genome_info = genome_info +'\n'.join(dtc_parser.meta_data)

    # filter and output as 23andme
    final_snps = output_as_23andme(store_snps, genome_info, output_23andme_path)
    # output to vcf
    final_snps = directconvertSNPtoVCFformat(store_snps, genome_info, outindex, output_vcf_path)

    status = 'Processed'
    info = {
        'count_snps': count_snps,
        'count_not_found_in_ref': count_not_found_in_ref,
        'valid_snps': valid_snps,
        'in_ref': in_ref,
        'cluster_id': cluster_id,
        'cluster_mcov': cluster_mcov,
        'vcf_snps': final_snps
    }
    return status, info

    
def dtc_read_records(dtc_parser, api_23andme, faidx, theLIST_clust, RS2GRCh37Orien_1):

    ## initial reject_type is empty
    reject = False
    reject_type = ""

    count_snps = 0
    count_not_found_in_ref = 0
    
    # valid_snps: exclude no calls;
    # also make sure that a reference genotype can be found, otherwise it may be from a different build
    valid_snps = 0
    store_snps = []
    count_genotype_contains_reference = 0
    
    count_not_in_list = 0
    cluster_count = Counter()

    try:
        for record in dtc_parser:
            
            if record.badsnp:
                # print("badsnp encountered file %s, line is %s" % (dtc_parser.path, record.line), file=sys.stderr)
                dtc_parser.total_badsnp += 1
                continue
            
            ## special treatment for 23andme
            if dtc_parser.company == '23andme' and record.loci not in api_23andme:
                record.badsnp = True
                dtc_parser.total_badsnp += 1
                continue
            
            count_snps += 1
            if len(set(record.genotype).intersection(usable_genotypes)) > 0:
                ref = get_snp_reference(record.loci, faidx)

                if ref:
                    if len(ref) == 1:
                        valid_snps += 1
                        store_snps.append({
                            'rsid': record.rsid,
                            'chro': record.chro,
                            'pos': record.pos,
                            'call': record.genotype,
                            'ref': ref
                        })
                        # compare to reference
                        if ref in record.genotype:
                            count_genotype_contains_reference += 1
                        
                        # cluster check
                        if record.loci in theLIST_clust.keys():
                            for clust in theLIST_clust[record.loci]:
                                cluster_count[clust] += 1
                        else:
                            count_not_in_list += 1
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
    
    return reject, reject_type, store_snps, count_snps, count_not_found_in_ref, \
    valid_snps, in_ref, revd_in_ref, \
    in_rev_strand, rev_in_ref, rev_revd_in_ref, \
    not_in_list, cluster_id, cluster_mcov

def filter_by_badalleles(store_snps, badalleles):
    new_snps = []
    for snp in store_snps:
        snppos = ':'.join([snp['chro'], snp['pos']])
        if snppos not in badalleles:
            new_snps.append(snp)
    return new_snps

def output_as_23andme(store_snps, genome_info, output_file):
    # Passsed all checks, apply the correct bad snps filter and output in 23andme format

    final_snps = 0
    with open(output_file, 'w') as outh:
        print(genome_info, file=outh)
        for snp in store_snps:
            final_snps += 1
            print('\t'.join([snp['rsid'], snp['chro'], snp['pos'], snp['call']]), file=outh)

    return final_snps

if __name__ == '__main__':
    main()