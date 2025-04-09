# Path: variable-info/update-outcome-info/get_new_Data_Dictionary.py
# Uses: add_new_fields.sh and the old outcome-info.tsv file
# Date: 2023-10-06
# Author: Chang Lu

import subprocess
import pandas as pd

############### 1. Get the Data Dictionary from UK Biobank ###############
subprocess.call('wget http://biobank.ctsu.ox.ac.uk/%7Ebbdatan/Data_Dictionary_Showcase.csv', shell=True)


############### 2. Convert csv to tsv ###############
# fix encoding format of the csv from Windows to UTF-8
subprocess.call('iconv -f WINDOWS-1252 -t UTF-8 Data_Dictionary_Showcase.csv > Data_Dictionary_Showcase2.csv', shell=True)

# replace \" with '
with open('Data_Dictionary_Showcase2.csv', 'r') as h:
    # read all the text in the file into a single string
    read = h.read()
read = read.replace('\\"', ' ')
with open('Data_Dictionary_Showcase.csv', 'w') as h:
    h.write(read)

# read in the csv, 
dd = pd.read_csv('Data_Dictionary_Showcase.csv', quotechar='"')

# date
from datetime import date
# get today's date
today = date.today()
# format date to be yyyy-mm-dd
today = today.strftime("%Y_%m_%d")
# save as tsv and name the file with the date
dd.to_csv('Data_Dictionary_Showcase_{}.tsv'.format(today), sep='\t', index=False)
dd.to_csv('Data_Dictionary_Showcase.tsv', sep='\t', index=False)

# remove the csv file
subprocess.call('rm Data_Dictionary_Showcase.csv', shell=True)
subprocess.call('rm Data_Dictionary_Showcase2.csv', shell=True)

############### 3. Make updated outcome info file, called outcome-info-new.tsv. ###############
# add new fields
subprocess.call('bash add_new_fields.sh', shell=True)

# remove the tsv file
subprocess.call('rm Data_Dictionary_Showcase.tsv', shell=True)

############## 4. Review fields and manually update PHESANT properties. ###############
