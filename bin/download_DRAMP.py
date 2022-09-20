#!/bin/python3

# TITLE: Download the DRAMP database and format it for mmseqs and diamond

import pandas as pd
import requests
import os

#The ref db directory path
path = '/home/aibrahim/github/Ampcombi_test/amp_combi'

########################################
#  Download_DRAMP
#########################################
#def download_database(path):
#    """
#    This is to download the database (DRAMP)
#    """

final_dir = os.path.join(path, r'DRAMP_db')

#check if the DRAMP folder exists if not create one and change into it 
if not os.path.exists(final_dir):
    os.makedirs(final_dir)
else:
    None

#Download the file and store it in a results directory 
os.chdir(final_dir)
url = 'http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.xlsx'
r = requests.get(url, allow_redirects=True)
with open('general_amps.xlsx', 'wb') as f:
    f.write(r.content)

#Convert excel to tab sep file and write it to a file in the DRAMP_db direct
df = pd.DataFrame(pd.read_excel("general_amps.xlsx"))
df.to_csv(os.path.join(final_dir, 'DRAMP_general_amps'),sep='\t')