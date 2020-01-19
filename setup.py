from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='dpellow-scapp',
    version="1.0",
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
    scripts = ['scapp/bin/scapp.py','scapp/bin/recycle.py', 'scapp/bin/make_fasta_from_fastg.py', 'scapp/bin/classify_fastg.py',\
                'scapp/bin/find_plasmid_gene_matches.py', 'scapp/bin/parse_plasmid_scores.py', \
                'scapp/bin/create_hits_fasta.py'],
    packages = ['scapp'],#,'recyclelib'],
    install_requires=[
        'networkx==2.4',
        'pysam==0.15.3',
        'nose==1.3',
        'numpy==1.17'],
    include_package_data=True,
)
