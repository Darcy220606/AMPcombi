from setuptools import setup

with open('README.md') as fh:
    long_description = fh.read()

setup(
    name='AMPcombi',
    version='0.1.6',
    author='Anan Ibrahim, Louisa Perelo',
    author_email='ananhamido@hotmail.com, louperelo@gmail.com',
    packages=['ampcombi'],
    scripts=['ampcombi/ampcombi.py',
             'ampcombi/amp_database.py',
             'ampcombi/amp_fasta.py',
             'ampcombi/check_input.py',
             'ampcombi/diamond_alignment.sh',
             'ampcombi/diamond_makedb.sh',
             'ampcombi/reformat_tables.py',
             'ampcombi/print_header.py',
             'ampcombi/version.py',
             'ampcombi/visualise_complete_summary.py',
             'ampcombi/HTML.R'],
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
            'ampcombi = ampcombi:main',
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