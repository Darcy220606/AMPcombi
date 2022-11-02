#!/bin/python3

# TITLE: Visualise teh complete summary and save it to a HTML file

import subprocess

########################################
#  FUNCTION: GENERATE AN INTERACTIVE HTML SUMMARY 
#########################################
def html_generator():
    subprocess.run('HTML.R', text=True)