from setuptools import setup, find_packages
   
with open('README.md') as fh:
    long_description = fh.read()

setup(
    name='AMPcombi',
    version='0.1.4',
    author='Anan Ibrahim, Louisa Perelo',
    author_email='ananhamido@hotmail.com, louperelo@gmail.com',
    packages=find_packages(),    
    scripts=['ampcombi/parse/parse.py',
             'ampcombi/parse/amp_database.py',
             'ampcombi/parse/amp_fasta.py',
             'ampcombi/parse/check_input.py',
             'ampcombi/parse/diamond_alignment.sh',
             'ampcombi/parse/diamond_makedb.sh',
             'ampcombi/parse/reformat_tables.py',
             'ampcombi/parse/print_header.py',
             'ampcombi/parse/version.py',
             'ampcombi/summarize/summarize.py'],
    url='http://pypi.python.org/pypi/AMPcombi/',
    license='LICENSE.txt',
    description='A parsing tool for AMP tools.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    keywords=["Proteomics", "Antimicrobial peptides", "Diamond"
              "Standardization", "Formatting","Functional annotation"],
    install_requires=['pandas'],
    entry_points={  
        'console_scripts': [
            'ampcombi-parse = parse:main', 'ampcombi-summarize = summarize:main'
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: POSIX",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Information Analysis"
        ],
)