
import os
import subprocess
from zipfile import ZipFile
import gzip
import bz2
import tempfile

# tmpdir = tempfile.TemporaryDirectory()
# resolve_file_for_genome(file_path, tmpdir.name)
# # (don't do) tmpdir.cleanup()
def resolve_file_for_genome(file_path):
    tmpdir = tempfile.TemporaryDirectory()
    tempdir = tmpdir.name
    resolved = False
    resolved_file_path = file_path
    iszip = False
    handle = None

    file_type = subprocess.check_output(["file",file_path]).decode("utf-8").replace("\n","")
    if "Zip archive" in file_type:
        iszip = True
        # extract files
        try:
            with ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(tempdir)
                namelist = zip_ref.namelist()
        except:
            namelist = ''
        # go over the extracted files
        if len(namelist) == 1:
            # if only one file in the zip, reject if we cannot open
            f = os.path.join(tempdir, namelist[0])
            file_type = subprocess.check_output(["file",f]).decode("utf-8").replace("\n","")
            handle = get_handle(f, file_type)
            if handle:
                resolved = True
                resolved_file_path = f

        elif len(namelist) > 1:
            # it will find the last text file, 
            # or the first text file with deCODEme_ or genome in filename
            for f in namelist:
                f = os.path.join(tempdir, f)
                if os.path.isfile(f):
                    file_type = subprocess.check_output(["file",f]).decode("utf-8").replace("\n","")
                    handle = get_handle(f, file_type)
                    if handle:
                        resolved = True
                        resolved_file_path = f
                        if 'deCODEme_' in f or 'genome' in f:
                            # then it must be the one, let's get out
                            break
    else:
        ## reject if cannot open.
        handle = get_handle(file_path, file_type)
        if handle:
            resolved = True

    return resolved, handle, iszip
    

def get_handle(file_path, file_type):
    ## get the read stream handle from a single file 
    # accept txt, gzip， bzip2
    handle = None
    if "ASCII text" in file_type or "RSID sidtune PlaySID compatible" in file_type:
        try:
            handle = open(file_path,'r')
        except:
            pass
    elif "gzip" in file_type:
        try:
            handle = gzip.open(file_path, 'rt')
        except:
            pass
    elif "bzip" in file_type:
        try:
            handle = bz2.open(file_path, 'rt')
        except:
            pass
    else:
        # make last attempt
        try:
            handle = open(file_path, 'r')
        except:
            pass
    return handle

def flip_genotype(genotype):
    new_genotype = ''
    for g in genotype:
        if g == 'A':
            g = 'T'
        elif g == 'T':
            g = 'A'
        elif g == 'C':
            g = 'G'
        elif g == 'G':
            g = 'C'
        new_genotype += g
    return new_genotype

# quick check for vcf, they have at least 10 columns; whereas DTCs have less than 8
def check_for_vcf(handle):
    vcflines = 0
    nonvcflines = 0
    crap_vcf_format = 0
    columns = 0
    if handle:
        for line in handle:
            if line.startswith('#'):
                columns = len(line.rstrip().split('\t'))
                continue
            
            if columns == 0:
                # if there has been no # lines
                columns = len(line.rstrip().split('\t'))
            
            # see if there are >= 10 columns
            if columns < 10:
                return False
            
            if line == '\n':
                crap_vcf_format += 1
                continue

            if len(line.rstrip().split('\t')) == columns:
                vcflines += 1
            else:
                nonvcflines += 1
        handle.close()
    
    # possibly a vcf 
    if vcflines > 5000 and nonvcflines/(vcflines+nonvcflines) < 0.1:
        return True

    return False