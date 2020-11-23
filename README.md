# SCAPP

SCAPP assembles plasmids from metagenomic assembly graphs.

- [Installation](#installation)
  * [With Conda](#with-conda)
      - [From Bioconda](#from-bioconda)
      - [From yaml file](#with-install-scapp-yaml-file)
  * [From sources](#from-sources)
      - [Configuring paths to required executables](#configuring-paths-to-required-executables)
      - [Testing your SCAPP installation](#testing-your-scapp-installation)
- [Basic Usage](#basic-usage)
- [Main output files](#main-output-files)
- [Advanced usage](#advanced-usage)
  * [Plasmid-specific genes](#plasmid-specific-genes)

## Installation

### With Conda

#### From Bioconda
You can install directly from Bioconda with `conda install -c bioconda scapp`

#### With install scapp yaml file

Alternatively, you can install SCAPP as a conda package (tested with Miniconda3):
Download the installation file `install_scapp.yaml` in the desired folder. For example:
````
wget https://raw.githubusercontent.com/Shamir-Lab/SCAPP/master/install_scapp.yaml
````

Create and activate the conda environment:
```
conda env create -f install_scapp.yaml
conda activate scapp
```

Now you can run SCAPP by entering the command `scapp`.


### From sources
If not using Conda to install, then download the sources and install according to the following:

SCAPP is written in Python3. SCAPP uses NumPy, NetworkX, pySAM, and nose. The necessary versions of these required dependencies will all be installed by the `setup.py` script.

SCAPP uses [BWA](https://github.com/lh3/bwa) (tested with v0.7.5 and v0.7.17) , [NCBI BLAST+ tools](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (tested with v2.7 and v2.9), and [samtools](https://github.com/samtools/samtools) (tested with v1.9 and v1.10). The executables of these programs should be available on the system on which SCAPP is run.

The [PlasClass classifier](https://github.com/Shamir-Lab/PlasClass) should also be installed in order to use the full functionality of SCAPP.

We recommend using a virtual environment. For example, in Linux, before running `setup.py`:
```
python -m venv scapp-env
source scapp-env/bin/activate
```
To install, download and run setup.py:
```
    git clone https://github.com/Shamir-Lab/SCAPP.git
    cd SCAPP
    python setup.py install
```
It is possible to install as a user without root permissions:
```
python setup.py install --user
```

To install PlasClass, in the folder you would like to install do:
```
git clone https://github.com/Shamir-Lab/PlasClass.git
cd PlasClass
python setup.py install
```

##### Configuring paths to required executables
Note that this step can be skipped if you installed using Conda.

The BWA, samtools, and BLAST+ executables must be available to SCAPP. They can either be added to your `PATH` environment variable, or you can specify the paths to each of them in the file `scapp/config.json`.

For example, if the BWA executable is in `/usr/bin/bwa/` then the line `"BWA_PATH" : "/usr/bin/bwa",` should be completed in the `config.json` file if that location is not in your `PATH`.

##### Testing your SCAPP installation
Once you have completed the above you can run `./run_test.sh` from the outermost SCAPP directory to test your installation.

(If the test fails please ensure you have set up the environment, installed, and configured SCAPP as described. Open an issue on this GitHub page with any problems you run into.)

## Basic Usage
To run the SCAPP pipeline: 
```
scapp -g <fastg graph> -o <output directory> [-k <max k value>] -r1 <reads 1> -r2 <reads 2> [-p <num processes>]
```
If a BAM alignment file of the reads to the assembly graph already exists, then use the following command to avoid re-running the alignment:
```
scapp -g <fastg graph> -o <output directory> [-k <max k value>] -b <BAM file> [-p <num processes>]
```
The common command line options are:

`-g/--graph`: : Assembly graph fastg file.

`-o/--output_dir`: Output directory.

`-k/max_k`: Maximum k value used by the assembler. Default: 55.

`-p/--num_processes`: Number of processes to use. Default: 16.

`-r1/--reads1`: Paired-end reads file 1.

`-r2/--reads2`: Paired-end reads file 2.

`-b/--bam`: BAM alignment file aligning reads to graph nodes `-b` and `-r1`,`-r2` are mutually exclusive.

## Main output files
The output files are written in the directory specified by the user:

`<prefix>.confident_cycs.fasta` is the main output fasta of the plasmid predictions. `<prefix>` is the name of the assembly graph file (minus the `.fastg` suffix).

Under the `intermediate_files` subdirectory `<prefix>_cycs.fasta` is a fasta file of **all** cyclic paths that were considered as potential plasmids before filtering to create the subset that is output.

`intermediate_files/reads_pe_primary.sort.bam(.bai)` is the alignment file for the reads to the assembly graph. It can be re-used with the `-b` option if SCAPP is re-run on the same sample.

Under the `logs` subdirectory, the file `scapp.log` contains all information about the SCAPP run. Other log files in the `logs` subdirectory may be helpful if there is an error or failure in one of the stages that runs BLAST, BWA, or PlasClass.

## Advanced usage
More advanced command line options allow for different stages in the SCAPP pipeline to be modified:

`-sc/--use_scores`: Flag to determine whether to use plasmid scores. Use `False` to turn off plasmid score use. Default: `True`.

`-gh/--use_gene_hits`: Flag to determine whether to use plasmid specific genes. Use `False` to turn off plasmid gene use. Default: `True`.

`-pc/--plasclass`: PlasClass score file. If PlasClass classification of the assembly graph nodes has already been performed, provide the name of the PlasClass output file. (For example: the `intermediate_files/plasclass.out` file from a previous run of SCAPP).

`-pf/--plasflow`: PlasFlow score file. To use PlasFlow scores for the nodes instead of PlasClass, provide the name of the PlasFlow output file.`-pf`,`-pc` are mutually exclusive.

In addition, all of the different thresholds used in the algorithm can be changed by the user:

`-m/--max_CV`: Maximum allowed coefficient of variation for coverage. Default: 0.5.

`-l/--min_length`: Minimum allowed length for potential plasmid. Default: 1000.

`-clft/--classification_thresh`: Threshold for classifying a potential plasmid as a plasmid. Default: 0.5.

`-gm/--gene_match_thresh`: Threshold for % identity and fraction of length covered to determine plasmid gene matches. Default: 0.75.

`-sls/selfloop_score_thresh`: Threshold plasmid score above which a self-loop is considered a potential plasmid. Default: 0.9.

`-slm/--selfloop_mate_thresh`: Threshold fraction of off-loop mate-pairs, below which a self-loop is considered a potential plasmid. Default:0.1.

`-cst/--chromosome_score_thresh`: Threshold score, below which a long node is considered a chromosome node. Default: 0.2.

`-clt/--chromosome_length_thresh`: Threshold length, above which a low scoring node is considered a chromosome node. Default: 10000.

`-pst/--plasmid_score_thresh`: Threshold score, above which a long node is considered a plasmid node. Default: 0.9.

`-plt/--plasmid_length_thresh`: Threshold length, above which a high scoring node is considered a plasmid node. Default: 10000.

`-cd/--good_cyc_dominated_thresh`: Threshold for the maximum fraction of nodes with most mate-pairs off the cycle allowed for the cycleto be considered a potential plasmid. Default: 0.5.

Instead of inputting all of these options on the command-line before each run of SCAPP, the user can change them in the file `bin/params.json`. Set each variable in this file to the desired value and it will be used in SCAPP. Any value passed as a command-line parameter will override the values set in this file.

### Plasmid-specific genes

SCAPP searches for plasmid-specific genes in the assembly graph and potential plasmids. Curated sets of plasmid-specific genes (see the SCAPP manuscript for details) are located in the `scapp/data` directory.

The user can add their own plasmid-specific gene sets in this directory. You may put nucleotide gene sequences in the `data/nt` subdirectory, or amino acid protein sequences in the `data/aa` subdirectory. The sequence files should be in fasta format.

If you wish to remove certain plasmid-specific sets, simply move them out of the `data` directory.


