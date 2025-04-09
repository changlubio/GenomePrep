'''
Run periodically to update the database with new data from OpenSNP.

/Volumes/acasis$ python git/GenomePrep/bin/opensnp_update.py
# 2023-09-25
Total number of entries: 6855
Distinct users: 6139
'''
import pandas as pd

import requests
from lxml import html

url = "https://opensnp.org/genotypes?page={}"

def make_updates():

    all_files = []

    # Get table of file information
    req = requests.session()
    for i in range(1, 400):
        # ignore warnings
        requests.packages.urllib3.disable_warnings()
        page = req.get(url.format(i), verify=False)
        tree = html.fromstring(page.text)
        
        # Get file name
        file_names = [x.lstrip('../data/') for x in tree.xpath("//table[@class='table table-hover']//tr//td[5]/a/@href")]
        if len(file_names) == 0:
            break

        # Get file ID
        fileIDs = [x.lstrip() for x in tree.xpath("//table[@class='table table-hover']//td[2]//text()[1]")]
        
        # Get time uplodaded
        times = [x.lstrip() for x in tree.xpath("//table[@class='table table-hover']//td[3]//text()[1]")]

        # Get user ID
        nicknames = [x.lstrip() for x in tree.xpath("//table[@class='table table-hover']//td[1]//text()[1]") if not '\n' in x]

        assert (len(file_names) == len(fileIDs) == len(times) == len(nicknames))

        for j in range(0, len(file_names)):
            assert(fileIDs[j] == file_names[j].split('.')[2])
            all_files.append({
                'nickname': nicknames[j], \
                'userID': file_names[j].split('.')[0], \
                'fileID': fileIDs[j], \
                'filename': file_names[j], \
                'time': times[j], \
                'download': 'https://opensnp.org/data/' + file_names[j]
            })

    uploads = pd.DataFrame({
        'nickname': [x['nickname'] for x in all_files], \
        'userID': [x['userID'] for x in all_files], \
        'fileID': [x['fileID'] for x in all_files], \
        'filename': [x['filename'] for x in all_files], \
        'time': [x['time'] for x in all_files], \
        'download': [x['download'] for x in all_files]
    })
    uploads.to_csv('opensnp_uploads.csv', index=False)

make_updates()