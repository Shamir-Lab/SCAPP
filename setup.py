from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='scapp',
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
    scripts = ['bin/recycle.py', 'bin/make_fasta_from_fastg.py', 'bin/classify_fastg.py',\
                'bin/find_plasmid_gene_matches.py', 'bin/parse_plasmid_scores.py', \
                'bin/create_hits_fasta.py'],
    packages = ['recyclelib'],
    install_requires=[
        'networkx==2.4',
        'pysam==0.15.3',
        'nose==1.3',
        'numpy==1.17'],
    include_package_data=True,
)
