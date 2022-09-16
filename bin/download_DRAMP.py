#!/bin/python3

# TITLE: Download the DRAMP database and format it for mmseqs and diamond

import requests
import os

########################################
#  Download_DRAMP
#########################################
#def database(db):    
#   """
#   This is to download the database (DRAMP)
#   """

current_dir = os.getcwd()
final_dir = os.path.join(current_dir, r'DRAMP_db')

#check if the DRAMP folder exists if not create one and change into it 
if os.path.exists(final_dir):
    os.chdir(final_dir)
    if not os.path.exists(final_dir):
        os.makedirs(final_dir) and os.chdir(final_dir)

url = 'http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.xlsx'

r = requests.get(url, allow_redirects=True)

open('general_amps.xlsx', 'wb').write(r.content)