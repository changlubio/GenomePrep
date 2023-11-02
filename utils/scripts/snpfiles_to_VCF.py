#!/usr/bin/env python3

import argparse
import subprocess
import os
import time
import logging
import multiprocessing
from multiprocessing import Process
import queue
import sys

from snpfile_tools import *

FILE_SIZE_CUTOFF = 400000
 
def produce_file_name_mapping(file_list, dir_badfiles, output_dir, SNPfile_to_vcf_mapping):
    '''
    mapping is a list of input_path and output_path
    '''
    mapping = list("")
    # vcf IDs given for different companies
    index_by_company = {}
    for comp in COMPANY:
        index_by_company[comp] = 0

    with open(SNPfile_to_vcf_mapping, "w") as SNPfile_to_vcf_mapping_handle:
        
        for file_name in file_list:          
            ## 1. determine if we can open the file (ascii or gzip or zip)
            path = test_file(file_name)
            if path is None:
                cmd = 'mv {} {}'.format(file_name, dir_badfiles)
                continue
            ## 2. get the company name, and make sure it's one that we can recognize
            comp_name = read_comments(path)
            if comp_name is None:
                cmd = 'mv {} {}'.format(path, dir_badfiles)
                continue
            ## 3. if we can parse the file, add it to the mapping list
            index_by_company[comp_name] += 1
            vcf_name = "%s_%s" % (comp_name, index_by_company[comp_name])
            
            output_vcf = os.path.join(output_dir,vcf_name + ".vcf")
            print(os.path.basename(file_name) + "\t" + vcf_name, file=SNPfile_to_vcf_mapping_handle)

            mapping.append([file_name, output_vcf])
    
    return mapping

def process_snp_files(dna_file, background_file, mapping, dir_badfiles, badallele_file):
    
    print("Starting pid: %s to process %s files" % (os.getpid(), len(mapping)), file=sys.stderr)
    faidx = pyfaidx.Fasta(dna_file)

    if badallele_file is not None:
        badalleles = read_bad_alleles(badallele_file)

    with pysam.TabixFile(background_file) as tabix_handle:
        
        for [file_name, output_vcf] in mapping:

            start_time = time.time()
            
            ## 1. determine if we can open the file (ascii or gzip or zip)
            path = test_file(file_name)
            if path is None:
                cmd = 'mv {} {}'.format(file_name, dir_badfiles)
                sys.stderr.write("cannot_open: %s" % (file_name))
                continue
            ## 2. get the company name, and make sure it's one that we can recognize
            comp_name = read_comments(path)
            if comp_name is None:
                cmd = 'mv {} {}'.format(path, dir_badfiles)
                sys.stderr.write("cannot_parse: %s\n" % (file_name))
                continue
            ## 3. read the file into sample
            sample = snp_File('',path)
            method = getattr(sample, 'read_' + comp_name)
            method()

            if sample.total < FILE_SIZE_CUTOFF:
                sys.stderr.write("small_file:\t%s\t%s\n" % (sample.total, sample.path))
                cmd = 'mv {} {}'.format(path, dir_badfiles)
                subprocess.call(cmd.split())
                continue                
            ## 4. filter the sample and output as vcf
            sample.format_snps(faidx, tabix_handle)
            tmp = sample.include_total
            if badallele_file is not None:
                sample.filter_badsnp(badalleles)
                sys.stdout.write('%s: before %s\tafter_format_snps: %s\tfiltered_bad_snp: %s\n' % (file_name, sample.total, tmp, sample.include_total))
            else:
                sys.stdout.write('%s: before %s\tafter_format_snps: %s\n' % (file_name, sample.total, tmp))
            sample.print_to_VCF(output_vcf)
        
            sys.stderr.write("Finished:\t%s\ttime: %s\tpid: %s\n" % (file_name, (time.time()-start_time), os.getpid()))

def main():
    """
    Example usage: (on luca) 
    /data/clu/anaconda3/bin/python3 /data/clu/git_working/snowflake_python3/scripts/snpfiles_to_VCF.py \
    -d /home/gough/data/snowflake/Homo_sapiens.GRCh37.75.dna.toplevel.fa \
    -b /home/gough/snowflake/datadir/1.3/GRCh37/background.vcf.gz \
    -i /home/gough/data/rihab/23andmedata/OpenSNP/ \
    -t 50 -o vcf_files -a badallele.tab
    """
    parser = argparse.ArgumentParser(description="Parses raw genotyping files from 23andme, AncestryDNA, deCODEme and FamilyTreeDNA.")
    # There is an unlimited number of files, so nargs is +
    # parser.add_argument("-i","--input_dir", help="Name of the files to be processed.")
    parser.add_argument("-i", "--input_file", help="Name of the files to be processed.", nargs='+')
    parser.add_argument("-o","--output")
    parser.add_argument("-d","--dna_file", help="Homo_sapiens.GRCh37.75.dna.toplevel.fa")
    parser.add_argument("-b","--background_file", help="background.vcf.gz in GRCh37")
    parser.add_argument("-t", "--threads", default=1)
    parser.add_argument("-a","--badallele_file", help="bad allele file", default=None)

    args = parser.parse_args()

    # generate output dirs
    subprocess.call(["mkdir", "-p", args.output])
    dir_output = os.path.join(args.output, 'vcf')
    dir_badfiles = os.path.join(args.output, 'badfiles')
    subprocess.call(["mkdir", "-p", dir_output])
    subprocess.call(["mkdir", "-p", dir_badfiles])
    mapping_file = os.path.join(args.output, 'mapping.txt')

    ## read input files to file_list
    file_list = args.input_file
    
    # produce mapping
    mapping = produce_file_name_mapping(file_list, dir_badfiles, dir_output, mapping_file)
    print("File_mapping complete, total files to process: %s" % (len(mapping)))
    # Start the main-program
    threads = int(args.threads)
    batch_size, _ = divmod(len(mapping),threads)
    
    jobs = []
    for i in range(threads):
        if i<threads-1:
            p=Process(target = process_snp_files, args = (args.dna_file, args.background_file, mapping[i*batch_size:(i+1)*batch_size], dir_badfiles, args.badallele_file))
        else:
            p=Process(target = process_snp_files, args = (args.dna_file, args.background_file, mapping[i*batch_size:], dir_badfiles, args.badallele_file))    
        jobs.append(p)
        p.start()
    for proc in jobs:
        proc.join()
    print("make vcf completed!")

if __name__ == "__main__":
    main()
elif __name__ == "snpfiles_to_VCF":
    print("convert tools imported!")