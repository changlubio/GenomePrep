import subprocess

def get_handle(file_path, file_type):
    ## get the read stream handle from a single file 
    # accept txt, gzip， bzip2
    handle = None
    reject = False

    if "ASCII text" in file_type or "RSID sidtune PlaySID compatible" in file_type:
        try:
            handle = open(file_path,'r')
        except:
            reject = True
    elif "gzip" in file_type:
        try:
            handle = gzip.open(file_path, 'rt')
        except:
            reject = True
    elif "bzip" in file_type:
        try:
            handle = bz2.open(file_path, 'rt')
        except:
            reject = True
    else:
        # make last attempt
        try:
            handle = open(file_path, 'r')
        except:
            reject = True

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
def check_for_vcf(filename):
    vcflines = 0
    nonvcflines = 0

    crap_vcf_format = 0

    file_type = subprocess.check_output(["file",filename]).decode("utf-8").replace("\n","")
    handle = get_handle(filename, file_type)

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