#!/usr/bin/env python3

'''
GRCh36 outputs from save4check
1. make sure it's GRCh36, remember the reference allele
2. output as temporary vcf file
3. liftover grch36 to grch37
Reject.list: SNPs outside reference, or not similar to the reference genome enough
'''

usable_genotypes = set(['A', 'C', 'G', 'T', 'D', 'I'])

import os
import sys
import argparse
import gzip
import time
from collections import Counter
import pickle
import pyfaidx 
import subprocess
import tempfile
import json

from dtc_parser import DTCparser
from tools import get_snp_reference, flip_genotype, output_as_23andme, resolve_file_for_genome
from process import process

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('-d', '--datadir', default='./datadir')
    parser.add_argument('-o', '--outdir', default='./outputs', help='directory where the quality checked 23andMe or VCF format files are going to be written to')
    parser.add_argument('-i', '--outindex', default='outindex', help='the index/name for the output files, it will also be used as the index in the VCF-format output')
    args = parser.parse_args()

    subprocess.call('mkdir -p {}'.format(args.outdir), shell=True)
    # tmpdir = tempfile.TemporaryDirectory()
    # workingdir = tmpdir.name
    # dna_path = {
    #     '36': os.path.join(datadir, 'Homo_sapiens.NCBI36.54.dna.toplevel.fa'),
    #     '37': os.path.join(datadir, 'Homo_sapiens.GRCh37.75.dna.toplevel.fa'),
    #     '38': os.path.join(datadir, 'Homo_sapiens.NCBI36.54.dna.toplevel.fa')
    # }
    start = time.time()
    status, info = process_assembly(args.infile, args.datadir, args.outdir, args.outindex)

    if status == 'Invalid':
        print('status:\t{}\ninfo:\t{}'.format(status, info))
        
    elif status == 'Processed':
        print("Please check your {0:} dir for a {1:}.23andme.tohg19 file.\nHere are some stats for the assembly conversion. ({2:} sec to convert assembly)".format(
            args.outdir, args.outindex, time.time()-start))
        for k, i in info.items():
            print('\t'.join(map(str, [k,i])))
        print("Now we are processing the GRCh37 file ...", flush=True)
        status, info = process(
            os.path.join(args.outdir, '{}.23andme.tohg19'.format(args.outindex)),
            args.datadir,
            args.outdir,
            args.outindex
        )
        if status == 'Invalid':
            print('status:\t{}\ninfo:\t{}'.format(status, info))
            
        elif status == 'Processed':
            print("Please check your {0:} dir for a {1:}.23andme.checked and {1:}.vcf file.\nHere are some stats".format(
                args.outdir, args.outindex))
            for k, i in info.items():
                print('\t'.join(map(str, [k,i])))
        else:
            # if isinstance(info, str):
            # it's rejected
            print('status:\t{}\ninfo:\t{}'.format(status, info))
            print('TO BE SOLVED')
    else:
        print('status:\t{}\ninfo:\t{}'.format(status, info))
        print('TO BE SOLVED')
        
    print("It took us total {} sec to process.".format(time.time()-start))

def process_assembly(infile, datadir, outdir, outindex, assembly='hg18'):
    output_23andme_path = os.path.join(outdir, '.'.join([outindex,'23andme.tohg19']))
    dna37_path = os.path.join(datadir, 'Homo_sapiens.GRCh37.75.dna.toplevel.fa')

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
    success, temp_info = dtc_records_only(dtc_parser)
    if not success:
        status = 'Invalid'
        info = 'Your data have too many bad records'
        return status, info
    # we are trying to identify build now ... file has badsnps: dtc_parser.badsnps, oksnps: temp_info[1], excludenocal: temp_info[2]
    if assembly == 'hg18':
        liftovermap = os.path.join(datadir, 'hg18ToHg19.over.chain.gz')
    elif assembly == 'hg38':
        liftovermap = os.path.join(datadir, 'hg38ToHg19.over.chain.gz')
    else:
        return 'Invalid', 'We do not support assembly {} yet'.format(assembly)
    faidx = pyfaidx.Fasta(dna37_path)
    store_snps = temp_info[0]
    success, temp_info = liftoverToHg19_snps(liftovermap, store_snps, outdir, outindex, faidx)
    if not success:
        status = 'Error'
        info = temp_info
        return status, info
    if not is_correct_reference(temp_info[1]):
        status = 'Error'
        info = 'Not the right assembly, info is: {}'.format(json.dumps(temp_info[1]))
        return status, info
    store_hg19_snps = temp_info[0]
    header = ['#GenomePrep, liftover assembly done'] + dtc_parser.meta_data
    final_snps = output_as_23andme(store_hg19_snps, '\n'.join(header), output_23andme_path)
    status = 'Processed'
    info = temp_info[1]
    info['final_snps'] = final_snps
    return status, info

def is_correct_reference(counts):
    if counts['count_not_found_in_ref'] ==0 and counts['count_in_ref']/counts['valid_snps'] > 0.75 and counts['count_revd_in_ref']/counts['valid_snps'] < 0.05:
        return True
    return False

def dtc_records_and_ref(dtc_parser, faidx):
    ## initial reject_type is empty    
    count_snps = 0
    count_not_found_in_ref = 0
    count_revd_in_ref = 0
    in_ref = revd_in_ref = -1
    # valid_snps: exclude no calls;
    # also make sure that a reference genotype can be found, otherwise it may be from a different build
    valid_snps = 0
    store_snps = []
    count_genotype_contains_reference = 0
    try:
        for record in dtc_parser:
            if record.badsnp:
                # print("badsnp encountered file %s, line is %s" % (dtc_parser.path, record.line), file=sys.stderr)
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
                        revcall = flip_genotype(record.genotype)
                        if ref in revcall:
                            count_revd_in_ref += 1     
                else:
                    count_not_found_in_ref += 1
    except RecursionError:
        return False, None
        
    in_ref = round(count_genotype_contains_reference/valid_snps, 4)  
    revd_in_ref = round(count_revd_in_ref / valid_snps, 4)
    # if in_ref < 0.6 or count_not_found_in_ref > 0:
    #     reject = True
    #     reject_type = 'notgrch37'
    return True, [store_snps, count_snps, count_not_found_in_ref, valid_snps, in_ref, revd_in_ref]

def dtc_records_only(dtc_parser):
    ## initial reject_type is empty    
    count_snps = 0
    # valid_snps: exclude no calls;
    valid_snps = 0
    store_snps = []
    try:
        for record in dtc_parser:
            if record.badsnp:
                # print("badsnp encountered file %s, line is %s" % (dtc_parser.path, record.line), file=sys.stderr)
                dtc_parser.total_badsnp += 1
                continue
            count_snps += 1
            if len(set(record.genotype).intersection(usable_genotypes)) > 0:
                valid_snps += 1
                store_snps.append({
                    'rsid': record.rsid,
                    'chro': record.chro,
                    'pos': record.pos,
                    'call': record.genotype,
                    'ref': None
                })
    except RecursionError:
        return False, None
    return True, (store_snps, count_snps, valid_snps)

def runCmd(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()
    lines = p.stdout.readlines()
    lines = [x.decode('utf-8').rstrip() for x in lines]
    return p.returncode, lines

def liftoverToHg19_snps(liftovermap, store_snps, output_dir, unique_sampleid, faidx):
    # write store_snps to a temp bed file
    # convert the bed file by CrossMap
    # read the new_snps in temp bed file
    # only do clean up if the conversion is done without error
    temp_hg18_bed = os.path.join(output_dir, unique_sampleid+'.hg18.bed')
    temp_hg19_bed = os.path.join(output_dir, unique_sampleid+'.hg19.bed')
    with open(temp_hg18_bed, 'w') as h:
        for snp in store_snps:
            print('\t'.join([snp['chro'], snp['pos'], str(int(snp['pos'])+1), ':::'.join([snp['rsid'], snp['call']])]), file=h)
    cmd = "CrossMap.py bed {} {} {}".format(liftovermap, temp_hg18_bed, temp_hg19_bed)
    rtc, lines = runCmd(cmd.split())
    if rtc == 0:
        store_hg19_snps, counts = read_bed_file(temp_hg19_bed, faidx)
        runCmd('rm -f {0:}.hg18.bed {0:}.hg19.bed {0:}.hg19.bed.unmap'.format(os.path.join(output_dir, unique_sampleid)).split())
        return True, (store_hg19_snps, counts)
    else:
        return False, 'CrossMap.py existed with error. Error was {}'.format(';'.join(lines))

def read_bed_file(temp_hg19_bed, faidx):
    store_hg19_snps = []
    count_lines = 0
    count_in_ref = 0
    count_revd_in_ref = 0
    count_not_found_in_ref = 0
    with open(temp_hg19_bed, 'r') as h:
        for line in h:
            count_lines += 1
            chro, pos1, _, info = line.rstrip().split('\t')
            rsid, call = info.split(':::')
            if len(set(call).intersection(usable_genotypes)) > 0:
                ref = get_snp_reference(':'.join([chro, pos1]), faidx)
                if ref:
                    if len(ref) == 1:
                        store_hg19_snps.append({
                            'rsid': rsid,
                            'chro': chro,
                            'pos': pos1,
                            'call': call,
                            'ref': ref
                        })
                        # compare to reference
                        if ref in call:
                            count_in_ref += 1
                        revcall = flip_genotype(call)
                        if ref in revcall:
                            count_revd_in_ref += 1     
                else:
                    count_not_found_in_ref += 1
    counts = {
        'count_lines': count_lines,
        'valid_snps': len(store_hg19_snps),
        'count_not_found_in_ref': count_not_found_in_ref,
        'count_in_ref': count_in_ref,
        'count_revd_in_ref': count_revd_in_ref
    }
    return store_hg19_snps, counts


if __name__ == '__main__':
    main()