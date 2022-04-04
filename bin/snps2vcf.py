#!/usr/bin/env python3

import argparse
import sys

import pysam

valid_call = set(['A', 'C', 'G', 'T'])

def directconvertSNPtoVCFformat(store_snps, vcf_info, vcf_name, out_vcf):
    """
    snp:{
            'rsid': record.rsid,
            'chro': record.chro,
            'pos': record.pos,
            'call': record.genotype,
            'ref': ref,
            'ba': True/False
        }
    tabix_handle: pysam.TabixFile(args.background_file) as tabix_handle
    """

    final_snps = 0
    indels = 0
    refref = 0
    multialleic = 0

    with open(out_vcf, 'w') as output_handle:
        # Write VCF Header
        print(vcf_info, file=output_handle)
        fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', vcf_name]
        print('#' + '\t'.join(fields), file=output_handle)

        for snp in store_snps:

            rsid = snp['rsid']
            chromosome = snp['chro']
            position = snp['pos']
            genotype = snp['call']
            ref = snp['ref']
            # INFO
            infox = []
            if snp['gav']=='':
                infox.append('GAV=none')
            else:
                infox.append('GAV={}'.format(snp['gav']))
            # FILTER
            filterx = []
            if snp['ba']:
                filterx.append('bz3')

            genotypes = list(genotype)
            if not (1 <= len(genotypes) <= 2):
                # print("Skipping invalid genotype %s at %s:%s" % (genotype, chromosome, position), file=sys.stderr) 
                continue       

            set_genotypes = set(genotypes)

            if not all([g in valid_call for g in set_genotypes]):
                # skip calls that's not A, T, C, G. namely, skip indels
                indels += 1
                alt = '.'
                infox.append('INDEL={}'.format(genotype))
                filterx.append('indel')

            else:
            
                if len(set_genotypes) == 1:
                    if ref in set_genotypes:
                        # All calls are reference
                        # we ignore it
                        refref += 1
                        alt = '.'
                    else:
                        # All calls are alternative
                        alt = genotypes[0]
                else:
                    if ref in set_genotypes:
                        # one is reference, the other must be alternative
                        non_ref = set_genotypes - set([ref])
                        assert(len(non_ref) == 1)
                        alt = non_ref.pop()
                    else:
                        # Non of the calls are in the reference,
                        # we ignore them
                        multialleic += 1
                        alt = genotypes[0] # set a random alternative of the two calls
                        infox.append('MA={}'.format(str(set_genotypes)))
                        filterx.append('nref')
                        continue
            
            calls = []
            for g in genotypes:
                if g == ref:
                    calls.append('0')
                else:
                    calls.append('1')
            calls = sorted(calls)
            infox = ';'.join(infox)
            if len(filterx) == 0:
                filterx = 'PASS'
            else:
                filterx = ';'.join(filterx)

            fields = [chromosome, position, rsid, ref, alt, '.', filterx, infox, 'GT', '/'.join(calls)]

            print("\t".join(fields), file=output_handle)
            final_snps += 1
    
    return final_snps, indels, refref, multialleic