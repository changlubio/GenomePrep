
'''
Directed to customer -- parser
'''

import subprocess
import os
import re
import time
import logging
import sys

from utils.tools import flip_genotype, resolve_file_for_genome

valid_call = set(['A', 'C', 'G', 'T', 'I', 'D', '_', '-','0'])

valid_chromo = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT'])

COMPANY = ['23andme','ancestrydna', 'ftdna', 'myheritage', 'decodeme', 'livingdnacsv', 'all']

reasons = {'invalid': 'We are unable to properly process your upload. Please check your file and try again.',
'small':'Your genome is too small (<10000 lines).',
'incomplete':'Your genome is incomplete. We only process complete genomes where all chromosomes are sequenced.',
'grch36': 'Your DNA file is generated using human genome assembly GRCh36. ',
'grch38': 'Your DNA file is generated using human genome assembly GRCh38.',
'notgrch37': 'save4check Your DNA file contains positions not from assembly GRCh37.',
'vcf': 'save4check This could be because it is an vcf file.',
'zip': 'save4check This could be because it is compressed/zipped.',
'ref': 'save4check This could be because it has unspecified reverse strand SNPs.',
'imputation': 'save4check This could be because it is imputation.'}


class DTCparser(object):
    """
    Class for parsing a direct-to-customer file. On object creation, the headers will be read and immediately accessible.
    VCF records can then be accessed as an iterator. 
    """

    def __init__(self, handle, filename=''):
        self.path = filename
        self.company = None
        self.build = 'GRCh37'
        self.meta_data = []
        self.total_badsnp = 0

        self.reject = False
        self.reject_type = ''
        self.handle = handle
        
        ## if not rejected, we determine the company and parser upon initialization
        try:
            self.read_headers()
        except:
            raise Exception('Cannot read the file by handle')
            

    def read_headers(self):
        # this will set the company or decide we cannot parse it.

        line = None

        while True:

            try:
                line = self.handle.readline()
            except UnicodeDecodeError:
                self.reject = True
                break

            if isinstance(line, type(bytes())):
                line = line.decode('utf-8')
            
            # Either reacher the end of the file or we found the first non meta line
            if not line or not line.startswith('#'):
                break

            self.meta_data.append(line.rstrip())

        if not line:
            self.reject = True
            self.reject_type = 'invalid'
            self.handle.close()
            return
        
        if len(self.meta_data) == 0:
            # happens when there is a single line header that does not begin with '#'
            self.meta_data.append(line.rstrip())

        # Determine company
        if "23andme" in self.meta_data[0] or "23andMe" in self.meta_data[0]:
            self.company = '23andme'
            for x in self.meta_data:
                if 'build' in x:
                    if '36' in x:
                        self.reject = True
                        self.reject_type = 'grch36'
                        self.build = 'GRCh36'
                        return
                    elif '38' in x:
                        self.reject = True
                        self.reject_type = 'grch38'
                        self.build = 'GRCh38'
                        return
        elif "AncestryDNA" in self.meta_data[0]:
            self.company = 'ancestrydna'
            for x in self.meta_data:
                if 'build' in x:
                    if '36' in x:
                        self.reject = True
                        self.reject_type = 'grch36'
                        self.build = 'GRCh36'
                        return
                    elif '38' in x:
                        self.reject = True
                        self.reject_type = 'grch38'
                        self.build = 'GRCh38'
                        return
            
        elif "RSID,CHROMOSOME,POSITION,RESULT" in self.meta_data[0]:
            self.company = 'ftdna'
            
        elif "MyHeritage" in self.meta_data[0]:
            self.company = 'myheritage'
            
        elif "Name,Variation,Chromosome,Position,Strand,YourCode" in self.meta_data[0]:
            self.company = 'decodeme'
        
        elif "Living DNA customer genotype data download file" in self.meta_data[0]:
            self.company = 'livingdnacsv'
            for x in self.meta_data:
                if 'build' in x:
                    if '36' in x:
                        self.reject = True
                        self.reject_type = 'grch36'
                        self.build = 'GRCh36'
                        return
                    elif '38' in x:
                        self.reject = True
                        self.reject_type = 'grch38'
                        self.build = 'GRCh38'
                        return
        else:
            ## we cannot determine the company from header line, but we still try to parse it with mega parser
            self.company = 'all'


    def __iter__(self):
        return self

    def __next__(self):
        try:
            line = self.handle.readline()
        except:
            # exceptions include lines in bytes etc.
            line = ''

        if not line:
            self.handle.close()
            raise StopIteration

        if isinstance(line, type(bytes())):
            line=line.decode('utf-8')

        record = DTCrecord(line, self.company)

        return record

class DTCrecord(object):
    '''
    create a genotype record of SNP from DTC file line
    '''
    def __init__(self, line, company):
        self.line = line.rstrip()
        self.rsid = self.chro = self.pos = self.genotype = ''
        self.badsnp = False
        method = getattr(self, 'parse_' + company)
        method()
        # 
        self.check_badsnp()
        self.loci = ':'.join([self.chro, self.pos])

    def check_badsnp(self):
        if not self.chro in valid_chromo:
            self.badsnp = True
        
        try:
            int(self.pos)
        except ValueError:
            self.badsnp = True

        genotypes = list(self.genotype)
        if not (1 <= len(genotypes) <= 2):
            self.badsnp = True
        if not all([g in valid_call for g in genotypes]):
            self.badsnp = True

    def parse_23andme(self):
        try:
           (self.rsid, self.chro, self.pos, self.genotype) = self.line.rstrip().split()
        except:
            self.badsnp = True
    
    def parse_ancestrydna(self):
        try:
            (self.rsid, chro, self.pos, allele1, allele2) = self.line.rstrip().split()
        except ValueError:
            try:
                (self.rsid, chro, self.pos, allele1) = self.line.rstrip().split()
                allele2 = '_'
            except:
                self.badsnp = True
                return
        # interpret ancestry dna chromosome
        if (chro == '23' or chro == '25'):
            chro = 'X'
        elif (chro == '24'):
            chro = 'Y'
        elif (chro == '26'):
            chro = 'MT'
        self.chro = chro
        self.genotype = allele1 + allele2
    
    def parse_ftdna(self):
        try:
            (self.rsid, self.chro, self.pos, self.genotype) = self.line.replace('"', '').split(",")
        except ValueError:
            try:
                (self.rsid, self.chro, self.pos, self.genotype) = self.line.split()
            except:
                self.badsnp = True
                return
        if not self.rsid.startswith('rs'):
            self.badsnp = True
    
    def parse_myheritage(self):
        try:
           (self.rsid, self.chro, self.pos, self.genotype) = self.line.replace('"', '').split(",")
        except ValueError:
            try:
                (self.rsid, self.chro, self.pos, self.genotype) = self.line.split()
            except:
                self.badsnp = True
                return
        if not self.rsid.startswith('rs'):
            self.badsnp = True

    def parse_decodeme(self):
        try:
            (self.rsid, variation, self.chro, self.pos, strand, genotype) = self.line.replace(' ', '').split(",")
        except:
            self.badsnp = True
            return
        if strand == '-':
            self.genotype = flip_genotype(genotype)
        else:
            self.genotype = genotype
        if self.chro == 'M':
            self.chro = 'MT'
    
    def parse_livingdnacsv(self):
        try:
            (self.rsid, self.chro, self.pos, self.genotype) = self.line.split()
        except:
            self.badsnp = True
        if not self.rsid.startswith('rs'):
            self.badsnp = True
    
    def parse_all(self):
        # the mega parser
        allele1 = allele2 = None
        try:
            # try 23andme, 4 cols
           (self.rsid, self.chro, self.pos, self.genotype) = self.line.rstrip().split()
        except ValueError:
            try:
                #try ancestrydna, 5 cols
                (self.rsid, self.chro, self.pos, allele1, allele2) = self.line.rstrip().split()
            except ValueError:
                try:
                    # try ftdna, 4 col sep by ,
                    (self.rsid, self.chro, self.pos, self.genotype) = self.line.replace('"', '').split(",")
                except ValueError:
                    try:
                        # decode me
                        (self.rsid, variation, self.chro, self.pos, strand, self.genotype) = self.line.replace(' ', '').split(",")
                        if strand == '-':
                            self.genotype = flip_genotype(self.genotype)
                    except:
                        self.badsnp = True
                        return
        
        # interpret ancestry dna chromosome
        if (self.chro == '23' or self.chro == '25'):
            self.chro = 'X'
        elif (self.chro == '24'):
            self.chro = 'Y'
        elif (self.chro == '26'):
            self.chro = 'MT'
        
        if allele1:
            self.genotype = allele1 + allele2
        
        # interpret decodeme
        if self.chro == 'M':
            self.chro = 'MT'
        
        # if not self.rsid.startswith('rs'):
        #     self.badsnp = True


