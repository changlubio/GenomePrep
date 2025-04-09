
'''
@author: Chang Lu. Created on 30 April 2019

Purpose: compare SNP files
'''

import subprocess
import os
import re
import time
import logging
import sys

import pyfaidx
import pysam

valid_call = set(['A', 'C', 'G', 'T'])
'''
by excluding I, D in valid_call, I have excluded all the indels
'''
valid_chromo = set(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT'])

COMPANY = ['23andme','ancestrydna', 'ftdna', 'myheritage', 'decodeme', 'livingdnacsv', 'livingdnavcf']

# logging.basicConfig(filename=os.path.join(os.getcwd(),'set_ref_alt.log'), level=logging.WARN)

class Snp(object):
    # Each SNP is one line in the original file
    def __init__(self, rsid, chromosome, position, genotype):
        self.rsid = rsid # e.g. RS1234
        self.chromosome = chromosome # 1,2,3...
        self.position = position # e.g. 213412
        self.genotype = genotype # e.g. AC

        ## ref and alt are determined from the GRCh reference genome and 1k background
        self.ref = ""
        self.alt = ""

        self.include = True  # depends on the quality of the snp
    
    @property
    def canonical_id(self):
        return "%s_%s_%s/%s" % (self.chromosome, self.position, self.ref, self.alt)
    
    @property
    def locus(self):
        return "%s:%s" % (self.chromosome, self.position)

    ## define the following function to allow set operations on list of Snps
    def __eq__(self, other):
        return self.chromosome == other.chromosome and self.position == other.position
    def __hash__(self):
        return hash(('chromosome', self.chromosome, 'position', self.position))

    @property
    def is_valid(self):
        
        if not self.chromosome in valid_chromo:
            self.include = False
        try:
            int(self.position)
        except ValueError:
            self.include = False
        
        genotypes = list(self.genotype)
        if not (1 <= len(genotypes) <= 2):
            logging.debug("snp_validity\tinvalid_genotype_length\t%s\t%s" % (self.rsid, self.genotype))
            self.include = False
        if not all([g in valid_call for g in genotypes]):
            logging.debug("snp_validity\tinvalid_genotype\t%s\t%s" % (self.rsid, self.genotype))     
            self.include = False

        return self.include

    def set_snp_ref_alt(self, faidx, tabix_handle):
        """
        set the reference and alternative of SNP, 
        with inputs from the reference GRCh (fasta read as faidx);
                    and the 1k genome background file (tabix_handle)
        Discarded snps are written out in log file
        """

        ## get the reference genotype 
        try:
            ref = faidx[self.chromosome][int(self.position) - 1].seq
        except KeyError:
            logging.debug("set_ref_alt\tno_reference\t%s\t%s\t%s" % (self.rsid, self.chromosome, self.position))
            self.include = False
            return
        ## compare the reference genotype, the SNP, and the background
        set_genotypes = set(self.genotype)
        if len(set_genotypes) == 1:
            # only one variant type
            if ref in set_genotypes:
                # All calls are reference
                # we check the snowflake background to find out if there is a known variant there.
                # If not, the background is all reference and we skip. Otherwise, we use the variant (i.e. alt)
                # from the background and set our VCF to reference
                records = list(tabix_handle.fetch(self.chromosome, int(self.position)-1, int(self.position), parser=pysam.asTuple()))
                if len(records) == 0:
                    logging.debug("set_ref_alt\tno_alternative\t%s\t%s\t%s" % (self.rsid, self.chromosome, self.position))
                    self.include = False
                    return
                elif len(records) > 1:
                    logging.warn("set_ref_alt\tmultiple_entries\t%s\t%s\t%s" % (self.rsid, self.chromosome, self.position))
                    self.include = False
                    return
                
                record = records[0]

                bg_ref = record[3]
                alt = record[4]

                if bg_ref != ref:
                    # len(bg_ref)>1, multiple background references
                    logging.warn("set_ref_alt\tskip\t%s\t%s\t%s\tref(%s)!=bg_ref(%s)" % (self.rsid, self.chromosome, self.position, ref, bg_ref))
                    ##comment below 2lines to yield same results as Ben's legacy code
                    self.include = False
                    return
                
                if record[1] != self.position:
                    logging.warn("set_ref_alt\tskip\t%s\t%s\tposition(%s)!=bg_position(%s)" % (self.rsid, self.chromosome, self.position, record[1]))
                    self.include = False
                    return 
                
            else:
                # We have one type of call that is not ref, therefore we are all alternative
                alt = set_genotypes.pop()
        else:
            # two different variants
            if ref not in set_genotypes:
                logging.debug("set_ref_alt\tskip_triallelic\t%s\t%s\t%s" % (self.rsid, self.chromosome, self.position))
                self.include = False
                return
            else:
                # One of our calls is ref and one is alt - so alt is whichever is not ref
                non_ref = set_genotypes - set([ref])
                try:
                    assert(len(non_ref) == 1)
                except AssertionError:
                    print("ERROR: non_ref is not 1\t%s\t%s\t%s\t%s\t%s" % (self.rsid, self.chromosome, self.position, self.genotype, os.getpid()))
                    self.include = False
                    return
                # if len(non_ref) != 1:
                #     print("non_ref is not 1\t%s\t%s\t%s\t%s\t%s" % (self.rsid, self.chromosome, self.position, self.genotype))
                alt = non_ref.pop()
        
        self.ref = ref
        self.alt = alt
        return

    def is_badallele(self, badalleles):
        if self.locus in badalleles:
            self.include = False


class snp_File:
    def __init__(self, company, path):
        self.company=company
        self.path = path
        self.build = 37
        self.total = 0
        self.snps = []
        self.locus = set()
        
        self.include_total = 0

    def read_23andme(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = '23andme'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (rsid, chro, pos, genotype) = line.split()
                except:
                    sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                    continue
                tmp_snp = Snp(rsid, chro, pos, genotype)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)
        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]

    def read_ancestrydna(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = 'ancestrydna'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (rsid, chro, pos, allele1, allele2) = line.split()
                except ValueError:
                    try:
                        (rsid, chro, pos, allele1) = line.split()
                        allele2 = '_'
                    except:
                        sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                        continue
                # interpret ancestry dna chromosome
                if (chro == '23' or chro == '25'):
                    chro = 'X'
                elif (chro == '24'):
                    chro = 'Y'
                elif (chro == '26'):
                    chro = 'MT'           

                tmp_snp = Snp(rsid, chro, pos, allele1+allele2)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)

        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]

    def read_ftdna(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = 'ftdna'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (rsid, chro, pos, genotype) = line.replace('"', '').split(",")
                except ValueError:
                    try:
                        (rsid, chro, pos, genotype) = line.split()
                    except:
                        sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                        continue
                if not rsid.startswith('rs'):
                    continue
                
                tmp_snp = Snp(rsid, chro, pos, genotype)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)

        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]

    def read_myheritage(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = 'myheritage'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (rsid, chro, pos, genotype) = line.replace('"', '').split(",")
                except:
                    sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                    continue
                if not rsid.startswith('rs'):
                    continue
                
                tmp_snp = Snp(rsid, chro, pos, genotype)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)

        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]
        
    def read_decodeme(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = 'decodeme'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (rsid, aa, chro, pos, bb, genotype) = line.replace(' ', '').split(",")
                except:
                    sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                    continue
                if chro == 'M':
                    chro = 'MT'
                
                tmp_snp = Snp(rsid, chro, pos, genotype)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)

        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]
    
    def read_livingdnacsv(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = 'livingdnacsv'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (rsid, chro, pos, genotype) = line.split()
                except:
                    sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                    continue
                if not rsid.startswith('rs'):
                    continue
                
                tmp_snp = Snp(rsid, chro, pos, genotype)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)

        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]

    def read_livingdnavcf(self):
        try:
            with open(self.path,'r') as input_handle:
                csvfile = input_handle.readlines()
        except:
            sys.stderr.write('Badfile: cannot_read_text\t%s\n' % (self.path))
            self.include = False
            
        ## Add new sample
        self.company = 'livingdnavcf'

        ## Read in the snps
        for line in csvfile:
            if '#' not in line:
                line = line.rstrip()
                try:
                    (chro, pos, rsid, ref, alt, qual, filt, info, form, genotype) = line.split("\t")
                except:
                    sys.stderr.write('bad_snp_format\t%s\t%s\n' % (self.path,line))
                    continue
                
                # interpret living dna contig
                chro = int(chro)
                if chro == 0:
                    continue
                elif chro >= 1 and chro <= 22:
                    chro = str(chro)
                elif chro == 23:
                    chro = 'X'
                elif chro == 24:
                    chro = 'Y'
                # This is actually the Y pseudo autosomal region
                # but just bin into X as a common convention and used by 23andMe
                elif chro == 25:
                    chro = 'X'
                elif chro == 26:
                    chro = 'MT'
                
                if len(genotype) == 1:
                    call = genotype
                elif len(genotype) == 3:
                    call = [genotype[0], genotype[2]]
                else:
                    continue
                genotype=[]
                for c in call:
                    if c == '0':
                        genotype.append(ref)
                    elif c == '1':
                        genotype.append(alt)
                    else:
                        genotype.append('.')

                tmp_snp = Snp(rsid, chro, pos, genotype)
                if not tmp_snp.is_valid:
                    continue
                if tmp_snp.locus not in self.locus:
                    self.locus.add(tmp_snp.locus)
                    self.snps.append(tmp_snp)

        self.total = len(self.snps)
        self.headers = [x for x in csvfile if x.startswith('#')]

    def format_snps(self, faidx, tabix_handle):
        count = 0
        for snp in self.snps:
            snp.set_snp_ref_alt(faidx, tabix_handle)
            if snp.include:
                count += 1
        self.include_total = count
    
    def filter_badsnp(self, badalleles):
        count = 0
        for var in self.snps:
            if var.include:
                var.is_badallele(badalleles)
            if var.include:
                count += 1
        self.include_total = count

    def print_to_VCF(self, output_vcf):
        vcf_name = os.path.basename(output_vcf).split('.')[0]
        with open(output_vcf,"w") as output_handle:
            # Write vcf header
            fields = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', vcf_name]
            print('#' + '\t'.join(fields), file=output_handle)
            # output all snps
            for snp in self.snps:
                if snp.include:
                    calls = []
                    for g in list(snp.genotype):
                        if g == snp.ref:
                            calls.append('0')
                        elif g == snp.alt:
                            calls.append('1')
                        else:
                            print("ERROR!unrecognized_genotype: "+"\t".join([snp.chromosome, snp.position, snp.rsid, snp.genotype]))
                            # assert(False)
                    fields = [snp.chromosome, snp.position, snp.rsid, snp.ref, snp.alt, '.', '.','.','GT','|'.join(calls)]
                    print("\t".join(fields), file=output_handle)

def is_readable_file(file_name):
    try:
        file_type = subprocess.check_output(["file",file_name]).decode("utf-8").replace("\n","")
    except:
        return False
    
    if "ASCII text" in file_type or "RSID sidtune PlaySID compatible" in file_type:
        return True
    else:
        return False   

def test_file(file_name):
    '''
    test if the file is readable, by checking the type or unzip it,
    return the readable file path or None
    '''
    file_type = subprocess.check_output(["file",file_name]).decode("utf-8").replace("\n","")
    if is_readable_file(file_name):
        return file_name
    
    elif "gzip" in file_type:
        file_prename,file_extension = os.path.splitext(file_name)
        # use file extension to tell gunzip, and gunzip would remove the file extension
        if os.path.isfile(file_prename):
            sys.stderr.write("Badfile: duplicates\t%s\n" % file_name)
            return None
        try:
            subprocess.call(["gunzip", "-S", file_extension, file_name])
        except:
            sys.stderr.write("Badfile: cannot_open\t%s\n" % file_name)
            return None
        subprocess.call(["mv", file_prename, file_name])
        new_name = file_prename
        if is_readable_file(new_name):
            return new_name
        else:
            sys.stderr.write("Badfile: cannot_open\t%s\n" % new_name)
            return None
        
    elif "Zip" in file_type:
        cmd = "unzip -o {} -d {}".format(file_name, os.path.dirname(file_name))
        try:
            tmp_out = subprocess.check_output(cmd.split()).decode("utf-8")
        except:
            sys.stderr.write("Badfile: cannot_open\t%s\n" % file_name)
            return None
        new_name = os.path.join( os.path.dirname(file_name), tmp_out.splitlines()[1].split()[1] )
        if is_readable_file(new_name):
            return new_name
        else:
            sys.stderr.write("Badfile: cannot_open\t%s\n" % file_name)
            subprocess.call(['rm','-f',new_name])
            return None
    else:
        sys.stderr.write("Badfile: cannot_open\t%s\n" % file_name)
        return None

def read_comments(path):
    try:
        with open(path, 'r') as input_handle:
            all_lines = input_handle.readlines()
            comments = [x for x in all_lines if '#' in x]
            if len(comments) < 1:
                comments = [all_lines[0]]
    except Exception as ex:
        sys.stderr.write('Badfile: cannot_open_file\t%s, exception is %s\n' % (path, ex))
        return None

    if "23andme" in comments[0] or "23andMe" in comments[0]:
        comp_name = '23andme'
        flag = False
        for x in comments:
            if 'build' in x:
                if '37' in x:
                    flag = True
                elif '36' in x:
                    sys.stderr.write('Badfiles: build_is_36\t%s\n' % (path))
                    return None
        if not flag:
            sys.stderr.write('Badfiles: no_build_info\t%s\n' % (path))
            return None

    elif "AncestryDNA" in comments[0]:
        comp_name = 'ancestrydna'
        flag = False
        for x in comments:
            if 'build' in x:
                if '37' in x:
                    flag = True
                elif '36' in x:
                    sys.stderr.write('Badfiles: build_is_36\t%s\n' % (path))
                    return None
        if not flag:
            sys.stderr.write('Badfiles: no_build_info\t%s\n' % (path))
            return None
        
    elif "RSID,CHROMOSOME,POSITION,RESULT" in comments[0]:
        comp_name = 'ftdna'
        
    elif "MyHeritage" in comments[0]:
        comp_name = 'myheritage'
        
    elif "Name,Variation,Chromosome,Position,Strand,YourCode" in comments[0]:
        comp_name = 'decodeme'
    
    elif "Living DNA customer genotype data download file" in comments[0]:
        comp_name = 'livingdnacsv'
    
    elif re.match('LD[0-9A-Z]+\.total\.vcf', os.path.basename(path)):
        comp_name = 'livingdnavcf'

    else:
        sys.stderr.write('Badfile: unrecognized_company: %s\n' % (path))
        return None
    
    return comp_name

def read_bad_alleles(path):
    bad_alleles = set()
    with open(path, 'r') as bad_alleles_handle:
        for line in bad_alleles_handle:
            temp = line.split()[0]
            # temp is snp position (e.g. 1:23498503)
            bad_alleles.add(temp)

    return bad_alleles