#!/usr/bin/env python3

'''
Preprocess VCF files for Snowflake, including
0. IGNORE SNPs that a) failed any filter b) in/dels
1. Make sure it's the right assembly, if it's grch36, convert to 37.
2. Fix the chromosome 
3. reset the identifier of the VCF
e.g.
Run save4check_process.py first     ~/git/clu-work/julian_scripts/preprocess_vcf.py -i uid_info
For files that don't require any Filter    ~/git/clu-work/julian_scripts/preprocess_vcf.py --filelist --no_filter -i file.23andme.vcf
'''

import os
import sys
import argparse
import gzip        
import time
from collections import Counter
import pickle
import pyfaidx 
import subprocess

from vcf_parser import *
from tools import get_snp_reference, output_as_23andme
from grch36to37 import liftover36to37_vcf, liftover38to37_vcf

VALID_CHROM = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT'])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filelist', action="store_true", default=None)
    parser.add_argument('--no_filter', action="store_true", default=None, help="Only keep variants that have passed all filters")
    parser.add_argument('-i', '--input')
    parser.add_argument('-o','--output_dir', default='./')
    parser.add_argument('-r', '--rejected_path', default = 'temp_output')
    parser.add_argument('--dna_37', default = '/data/snow/snowflake/Homo_sapiens.GRCh37.75.dna.toplevel.fa')
    parser.add_argument('--daemon', default='1')

    args = parser.parse_args()
    if args.daemon == '1':
        daemon_save4check = '/data/snow/SNV/daemon/save4check/'
    elif args.daemon == '2':
        daemon_save4check = '/data/snow/SNV/daemon2/save4check/'

    output_dir = args.output_dir
    no_filter = args.no_filter

    faidx37 = pyfaidx.Fasta(args.dna_37)

    if not args.filelist:
        # special input from uid_info
        with open(args.input, 'rb') as h:
            uid_info = pickle.load(h)
        fileIDs = [uid for uid in uid_info.keys() if uid_info[uid]['reject_type'] == 'vcf'] 
        fileIDs = sorted(fileIDs)
        # with open('tmp.vcf.uniqueid.queue', 'w') as h:
        #     print('\n'.join(fileIDs), file=h)
        print('Processing %s vcf genomes' % (len(fileIDs)), flush=True)
        infiles = [os.path.join(daemon_save4check, x+'.save4check') for x in fileIDs]
    else:
        infiles = [args.input]

    rejected_path = os.path.join(output_dir, args.rejected_path)
    rh = open(rejected_path, 'w')

    # 1. output grch37 vcf 2. output 23andme 3. output grch36 converted vcf
    outdir_vcf = os.path.join(output_dir, 'vcf_real')
    outdir_23andme = os.path.join(output_dir, 'vcf_23andme')
    outdir_not37 = os.path.join(output_dir, 'vcf_not37')
    subprocess.call('mkdir -p {}'.format(outdir_vcf).split())
    subprocess.call('mkdir -p {}'.format(outdir_23andme).split())
    subprocess.call('mkdir -p {}'.format(outdir_not37).split())

    for idx, infile in enumerate(infiles):
        # if idx < 86:
        #     continue
        # elif idx > 85:
        #     break
        sample_id = os.path.basename(infile).split('.')[0]
        # if sample_id not in ['200096192611867', '314353401426273', '337838520768119', '373878783010182', '406464763350288', '421834574228315', '491797918891557', '492502377927101', '520043128107967', '555257360912882', '612227526888973', '615148754399527', '616898974700572', '687154555290882', '734666682251469', '817057274123376', '828498298801630', '850000636313727', '917792073184305', '946279690527037', '988090253402784']:
        #     continue
        print('processing %s %s' % (idx, infile), file=sys.stdout, flush=True)
        vcfparser = VCFparser(infile)

        if not vcfparser.handle:
            print('ERROR CANNOT OPEN %s' % infile, file=rh, flush=True)
            print('ERROR CANNOT OPEN %s' % infile, file=sys.stdout, flush=True)
            continue
        
        # reset sample id
        if len(vcfparser.sample_ids) != 1:
            if len(vcfparser.sample_ids) == 2 and vcfparser.sample_ids[0]=='unknown':
                # handle special situation: unknown Sample1 
                print('WARNING seperate sample id %s' % infile, file=rh, flush=True)
            else:
                print('ERROR >1 sample %s' % infile, file=rh, flush=True)
                print('ERROR >1 sample %s' % infile, file=sys.stdout, flush=True)
                continue

        start = time.time()

        # reset sample id
        vcfparser.sample_ids = [sample_id]

        ## read records and compare to grch37
        FILTER = 'filter_PASS'
        out_vcf_path = os.path.join(outdir_vcf, sample_id + '.vcf')
        
        try:
            with open(out_vcf_path, 'w') as vcfouthandle:
                store_snps, ref_call, out_call, not_inref37, nocall, invalidchrom, byfilter, indels \
                = read_vcf_records(vcfparser, faidx37, vcfouthandle, no_filter)
        except RecursionError:
            print('ERROR RecursionError %s' % infile, file=rh, flush=True)
            continue

        # #### handle some annoyingly weird veritase genome
        # with open(out_vcf_path, 'w') as vcfouthandle:
        #     store_snps, ref_call, out_call, not_inref37, nocall, invalidchrom, byfilter, indels \
        #     = read_vcf_records_tempspecial(vcfparser, faidx37, vcfouthandle, no_filter)
        # ####
        # it's likely that filter has no PASS values
        if len(store_snps) == 0:
            FILTER = 'filter_NONE'
            vcfparser = VCFparser(infile)
            vcfparser.sample_ids = [sample_id]
            # I don't like VCF files without a filter, put them aside for now
            subprocess.call('rm -f {}'.format(out_vcf_path).split())
            out_vcf_path = os.path.join(outdir_not37, sample_id + '.nofilter.vcf')
            with open(out_vcf_path, 'w') as vcfouthandle:
                store_snps, ref_call, out_call, not_inref37, nocall, invalidchrom, byfilter, indels \
                    = read_vcf_records(vcfparser, faidx37, vcfouthandle, True)
        ######
        
        # discard if incomplete
        chromo_set = set([x['chro'] for x in store_snps])
        if len(chromo_set) < 22 or len(store_snps) < 10000:
            subprocess.call('rm -f {}'.format(out_vcf_path).split())
            print('ERROR incomplete: %s info: %s' % (infile, '\t'.join(map(str, \
            [len(store_snps), not_inref37, nocall, invalidchrom, byfilter, indels]))), file=rh, flush=True)
            print('ERROR incomplete: %s info: %s' % (infile, '\t'.join(map(str, \
            [len(store_snps), not_inref37, nocall, invalidchrom, byfilter, indels]))), file=sys.stdout, flush=True)
            continue

        # if not grch37
        build = 'hg19'
        if vcfparser.build == '36':
            build = 'hg18'
            grch37_vcf_path = os.path.join( outdir_not37, sample_id+'.36to37.vcf')
            final_snps = liftover36to37_vcf(out_vcf_path, grch37_vcf_path)
        elif vcfparser.build == '38':
            build = 'hg38'
            grch37_vcf_path = os.path.join( outdir_not37, sample_id+'.38to37.vcf')
            final_snps = liftover38to37_vcf(out_vcf_path, grch37_vcf_path)
        elif vcfparser.build == '37':
            if not_inref37 > 0:
                build='hg19_prob'
            else:
                build='hg19'
        elif vcfparser.build == 'none':
            if not_inref37 == 0:
                build = 'hg19'
            else:
                build = 'none'
                # move the vcf to new directory
                grch37_vcf_path = os.path.join( outdir_not37, sample_id+'.unknown_build.vcf')
                subprocess.call('mv {0:} {1:}'.format(out_vcf_path, grch37_vcf_path).split())

        # if it's 23andme
        if ref_call/out_call > .1:
            genome_info = '#vcf back to 23andme\n#rsid\tchromosome\tposition\tgenotype'
            out_23andme_path = os.path.join(outdir_23andme, sample_id + '.23andme')
            output_as_23andme(store_snps, genome_info, out_23andme_path)
            if os.path.isfile(out_vcf_path):
                subprocess.call('rm -f {}'.format(out_vcf_path).split())
            vcfTO23andme = '23andme'
        else:
            vcfTO23andme = 'realvcf'
        
        if args.filelist:
            fields = [os.path.basename(infile), sample_id, FILTER, build, str(not_inref37), vcfTO23andme, ref_call, out_call, round(time.time() - start, 2)]
        else:
            fields = [uid_info[sample_id]['realname'], sample_id, FILTER, build, str(not_inref37), vcfTO23andme, ref_call, out_call, round(time.time() - start, 2)]
        print('\t'.join(map(str, fields)), file=rh, flush=True)
        print('\t'.join(map(str, fields)), file=sys.stdout, flush=True)

    rh.close()
    faidx37.close()


def read_vcf_records(vcfparser, faidx37, vcfouthandle, no_filter):

    for line in vcfparser.meta_data:
        print(line, file=vcfouthandle)

    # reset chromosome by reading and rewrite
    # count portion of only reference
    ref_call = out_call = 0
    not_inref37 = 0
    store_snps = []
    nocall = invalidchrom = byfilter = indels = 0
    for record in vcfparser:
        
        call = record.parse_genotype(record.calls[0])

        if not call:
            # Call: does not match expected format for snp
            nocall += 1
            continue
        
        if record.chromosome not in VALID_CHROM:
            # ignore non-standard chromosomes
            invalidchrom += 1
            continue
        
        if not no_filter and record.filter != 'PASS':
            # ignore calls that failed any filter.
            byfilter += 1
            continue
        
        if len(record.wild) != 1 or len(record.mutant) != 1:
            # ignore indels
            indels += 1
            continue
        # elif len(record.wild) < 1 or len(record.mutant) <1:
        #     # ignore if there is no wildtype in VCF -- bad line
        #     nowild += 1
        #     continue

        # we store all SNPs in case we want to output it as 23andme
        genotypes =  record.get_genotype(call)
        if not genotypes:
            nocall += 1
            continue
        
        store_snps.append({
            'rsid': record.id,\
            'chro': record.chromosome,\
            'pos': record.position,\
            'call': genotypes,\
            'ref': record.wild
        })

        OUT_TO_VCF = True
        # percentage of reference calls
        if len(set(call) - set(['0','.'])) == 0:
            # all reference call, we don't output these to vcf
            ref_call += 1
            OUT_TO_VCF = False
        
        # assembly check
        # check_num = 100000
        # if out_call + ref_call < check_num:
        try:
            ref37 = get_snp_reference(':'.join([record.chromosome, record.position]), faidx37)
        except Exception as e:
            print('ref Exception %s' % e, file=sys.stderr)
        if ref37:
            if record.wild[0] != ref37:
                # print('mismatch', record, ref37)
                not_inref37 += 1
                # OUT_TO_VCF = False
        else:
            # output everything, because we will convert them to the correct assembly later
            # OUT_TO_VCF = False 
            not_inref37 += 1
        
        if OUT_TO_VCF:
            print(record, file=vcfouthandle)
            out_call += 1
    
    return store_snps, ref_call, out_call, not_inref37, nocall, invalidchrom, byfilter, indels


def read_vcf_records_tempspecial(vcfparser, faidx37, vcfouthandle, no_filter):

    for line in vcfparser.meta_data:
        print(line, file=vcfouthandle)

    vcfparser.sample_ids = ['xx', 'sample']

    # reset chromosome by reading and rewrite
    # count portion of only reference
    ref_call = out_call = 0
    not_inref37 = 0
    store_snps = []
    nocall = invalidchrom = byfilter = indels = 0
    for record in vcfparser:
        assert(len(record.calls) ==  2)

        record.calls = [record.calls[1]]
        call = record.parse_genotype(record.calls[0])

        if not call:
            # Call: does not match expected format for snp
            nocall += 1
            continue
        
        if record.chromosome not in VALID_CHROM:
            # ignore non-standard chromosomes
            invalidchrom += 1
            continue
        
        if not no_filter and record.filter != 'PASS':
            # ignore calls that failed any filter.
            byfilter += 1
            continue
        
        if len(record.wild) != 1 or len(record.mutant) != 1:
            # ignore indels
            indels += 1
            continue
        # elif len(record.wild) < 1 or len(record.mutant) <1:
        #     # ignore if there is no wildtype in VCF -- bad line
        #     nowild += 1
        #     continue

        # we store all SNPs in case we want to output it as 23andme
        genotypes =  record.get_genotype(call)
        if not genotypes:
            nocall += 1
            continue
        
        store_snps.append({
            'rsid': record.id,\
            'chro': record.chromosome,\
            'pos': record.position,\
            'call': genotypes,\
            'ref': record.wild
        })

        OUT_TO_VCF = True
        # percentage of reference calls
        if len(set(call) - set(['0','.'])) == 0:
            # all reference call, we don't output these to vcf
            ref_call += 1
            OUT_TO_VCF = False
        
        # assembly check
        # check_num = 100000
        # if out_call + ref_call < check_num:
        try:
            ref37 = get_snp_reference(':'.join([record.chromosome, record.position]), faidx37)
        except Exception as e:
            print('ref Exception %s' % e, file=sys.stderr)
        if ref37:
            if record.wild[0] != ref37:
                # print('mismatch', record, ref37)
                not_inref37 += 1
                # OUT_TO_VCF = False
        else:
            # output everything, because we will convert them to the correct assembly later
            # OUT_TO_VCF = False 
            not_inref37 += 1
        
        if OUT_TO_VCF:
            print(record, file=vcfouthandle)
            out_call += 1
    
    return store_snps, ref_call, out_call, not_inref37, nocall, invalidchrom, byfilter, indels

if __name__ == '__main__':
    main()