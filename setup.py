from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='dpellow-scapp',
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
    scripts = ['scapp/scapp.py'],#'scapp/bin/recycle.py', 'scapp/bin/make_fasta_from_fastg.py', 'scapp/bin/classify_fastg.py',\
#                'scapp/bin/find_plasmid_gene_matches.py', 'scapp/bin/parse_plasmid_scores.py', \
#                'scapp/bin/create_hits_fasta.py'],
#    packages = ['scapp'],#find_packages(),
 
    packages=find_packages(where="scapp"),
    package_dir={"": "scapp"}
    install_requires=[
        'networkx==2.4',
        'pysam==0.15.3',
        'nose==1.3',
        'numpy==1.17'],
    package_data = {'scapp':['data/*']},
    include_package_data=True,
#    data_files = ['scapp/bin/PARAMS.py', 'scapp/bin/params.json'],
)
