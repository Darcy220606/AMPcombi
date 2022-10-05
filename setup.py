from distutils.core import setup

setup(
    name='AMPcombi',
    version='0.1.0',
    author='Anan Ibrahim and Louisa Perelo',
    author_email='ananhamido@hotmail.com and louisa.perelo@gmail.com',
    packages=['ampcombi'],
    scripts=['ampcombi/ampcombi.py',
             'ampcombi/amp_database.py',
             'ampcombi/amp_fasta.py',
             'ampcombi/check_input.py',
             'ampcombi/diamond_alignment.sh',
             'ampcombi/diamond_makedb.sh',
             'ampcombi/reformat_tables.py',
             'ampcombi/print_header.py'],
    url='http://pypi.python.org/pypi/AMPcombi/',
    license='LICENSE.txt',
    description='A parsing tool for AMP tools.',
    long_description=open('README.md').read(),
    keywords=["Proteomics", "Antimicrobial peptides", "Diamond"
              "Standardization", "Formatting","Functional annotation"],
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
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        "Topic :: Scientific/Engineering :: Information Analysis"
        ],
)