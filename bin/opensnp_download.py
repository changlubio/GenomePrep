'''
Run periodically to update the database with new data from OpenSNP.

/Volumes/acasis$ 

python git/GenomePrep/bin/opensnp_download.py data/opensnp_uploads.csv data/opensnp_genomes/
'''

import pandas as pd
import numpy as np
import sys
import os
import subprocess

# Read in the data from inputpath
inputpath = sys.argv[1]
uploads = pd.read_csv(inputpath)

# output dir
outputdir = sys.argv[2]
downloaded = os.listdir(outputdir)

# uploads = pd.DataFrame({
#     'nickname': [x['nickname'] for x in all_files], \
#     'userID': [x['userID'] for x in all_files], \
#     'fileID': [x['fileID'] for x in all_files], \
#     'filename': [x['filename'] for x in all_files], \
#     'time': [x['time'] for x in all_files], \
#     'download': [x['download'] for x in all_files]
# })

# Download genomes
for i in range(0, len(uploads)):
    if uploads['filename'][i] in downloaded:
        continue
    print('Downloading ' + uploads['filename'][i])
    subprocess.call(['wget', uploads['download'][i], '-P', outputdir])