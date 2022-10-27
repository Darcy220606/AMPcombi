#!/bin/python3

# TITLE: Visualise teh complete summary and save it to a HTML file

import subprocess

########################################
#  FUNCTION: DOWNLOAD DRAMP DATABASE AND CLEAN IT
#########################################
def html_generator():
    subprocess.run('HTML.R', text=True)