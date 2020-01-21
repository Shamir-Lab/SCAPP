from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='scapp-dpellow',
    version="0.1",
    author="David Pellow",
    author_email="dpellow@mail.tau.ac.il",
    long_description=long_description,
    long_description_content_type="text/markdown",
    description='SCAPP: an updated algorithm for detecting plasmids from de novo assembly graphs',
    url='https://github.com/dpellow/SCAPP',
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages = find_packages(),#['scapp'], 
    install_requires=[
        'networkx==2.4',
        'pysam==0.15.3',
        'nose==1.3',
        'numpy==1.17'],
    package_data = {'scapp':['data/*/*']},
    include_package_data=True,
    entry_points = {
        "console_scripts": [
            "scapp=scapp.scapp:main",
        ]
    },
    data_files = [('scapp',['scapp/scapp_utils.py', 'config.json', 'params.json','recycle.py', 'make_fasta_from_fastg.py', 'classify_fastg.py','PARAMS.py', 'parse_plasmid_scores.py','create_hits_fasta.py', 'find_plasmid_gene_matches.py'])]
    #data_files = [('',['scapp/scapp_utils.py', 'scapp/config.json', 'scapp/params.json','scapp/recycle.py', 'scapp/make_fasta_from_fastg.py', 'scapp/classify_fastg.py','scapp/PARAMS.py', 'scapp/parse_plasmid_scores.py','scapp/create_hits_fasta.py', 'scapp/find_plasmid_gene_matches.py'])]
)
