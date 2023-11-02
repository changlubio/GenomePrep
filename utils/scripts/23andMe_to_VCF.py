#!/usr/bin/env python3

import argparse
import sys

import pyfaidx
import pysam

"""
Converts a 23andMe file to a VCF file. Requires a Snowflake
background file - only variants in the background are kept

Currently only supports the long format version of 23andMe
TODO: Add support for compressed 23andMe file
"""

class InvalidCallException(Exception):
    pass

parser = argparse.ArgumentParser()
parser.add_argument('in_file')
parser.add_argument('dna_file')
parser.add_argument('background_file')
parser.add_argument('out_vcf')
parser.add_argument('ID_number')

args = parser.parse_args()

valid_call = set(['A', 'C', 'G', 'T'])

print("Reading faidx. If no index has been created it will be made now which may take some time...")
faidx = pyfaidx.Fasta(args.dna_file)
print("Done")


with open(args.in_file) as input_handle, \
    open(args.out_vcf, "w") as output_handle, \
    pysam.TabixFile(args.background_file) as tabix_handle :
    # Write VCF Header
    fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', args.ID_number]
    print('#' + '\t'.join(fields), file=output_handle)

    for line in input_handle:
        if line.startswith('#'):
            continue

        rsid, chromosome, position, genotype = line.strip().split()

        genotypes = list(genotype)
        if not (1 <= len(genotypes) <= 2):
            # print("Skipping invalid genotype %s at %s:%s" % (genotype, chromosome, position), file=sys.stderr) 
            continue       

        if not all([g in valid_call for g in genotypes]):
            # print("Skipping non-SNP genotype %s at %s:%s" % (genotype, chromosome, position), file=sys.stderr)     
            continue


        try:
            ref = faidx[chromosome][int(position) - 1].seq
        except KeyError:
            # print("Skipping call at %s:%s as no reference could be found" % (chromosome, position), file=sys.stderr)
            continue

        set_genotypes = set(genotypes)
        if len(set_genotypes) == 1:
            if ref in set_genotypes:
                # All calls are reference
                # We check the snowflake background to find out if there is a known variant there.
                # If not, the background is all reference and we skip. Otherwise, we use the variant (i.e. alt)
                # from the background and set our VCF to reference


                records = list(tabix_handle.fetch(chromosome, int(position) -1, int(position), parser=pysam.asTuple()))
                if len(records) == 0:
                    # print("Skipping call at %s:%s as we are refernece and so is the background" % (chromosome, position), file=sys.stderr)
                    continue            
                if len(records) > 1:
                    # print("Skipping call at %s:%s as multiple entries in bg were found" % (chromosome, position), file=sys.stderr)
                    continue


                record = records[0]

                bg_ref = record[3]
                alt = record[4]

                if bg_ref != ref:
                    # print("Skipping call at %s:%s as ref(%s) != bg ref(%s)" % (chromosome, position, ref, bg_ref), file=sys.stderr)
                    continue

            else:
                # We have one type of call that is not ref - therefore we are all alt
                alt = genotypes[0]


        else:
            # We have two different variants - one of which must be the ref otherwise something has gone wrong
            if ref not in set_genotypes:
                # print("Skipping call at %s:%s as we have distinct alleles %s that are not ref(%s)" % (chromosome, position, genotype, ref), file=sys.stderr)
                continue

            else:
                # One of our calls is ref and one is alt - so alt is whichever is not ref
                non_ref = set_genotypes - set([ref])
                assert(len(non_ref) == 1)
                alt = non_ref.pop()


        calls = []

        for g in genotypes:
            if g == ref:
                calls.append('0')
            elif g == alt:
                calls.append('1')
            else:
                assert(False)

        fields = [chromosome, position, rsid, ref, alt, '.', '.', '.', 'GT', '|'.join(calls)]

        print("\t".join(fields), file=output_handle)

#[clu] is the file not closed?
output_handle.close()
