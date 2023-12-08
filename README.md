# MPI
Microbial Pangenome Interpretation
---------------------------------------------------------
## Function

1) Initialization: The script sets up some initial parameters, such as percent identity, segment length, and the number of threads to be used in subsequent processes.

2) Command Line Options: The script parses command line options to customize the behavior of the pipeline. These options include specifying input data, setting a run identifier, specifying the number of genomes, and more.

3) Parameter Validation: It checks whether the provided parameters are within valid ranges and provides default values for some parameters if they are not provided.

4) Input Sample Handling: It checks if the input data is a directory or a single file. If it's a directory, it combines the contents into a single file.

5) Sequence Partitioning: If the "multiple_chromosomes" flag is enabled, the script performs [sequence partitioning](https://pggb.readthedocs.io/en/latest/rst/tutorials/sequence_partitioning.html). It calculates genetic distances between sequences, identifies communities, and analyzes each community separately.

6) Running the Pipeline: It sets up the required environment variables and runs a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline. This pipeline takes care of various tasks, such as sorting data, running [PGGB](https://github.com/pangenome/pggb), the core pangenome graph building tool used in this workflow, and generating reports.

7) Output Generation: After the pipeline completes, it creates an output directory and stores the results there.

## Setup and quick start

It is recommended to set up a new conda environment for MPI.

1) Clone this repo into your working directory
```bash
git clone https://www.github.com/KyranWissink/MPI 
cd MPI
```

2) Set up the conda environment
```bash
conda env create -f environment.yml
conda activate mpi
```

3) Run bash run.sh for the parameters
```bash
MPI Version: 0.7

Usage:

        bash run.sh [options] -i <input> 

Description:

        Pipeline for pangenome graph creation using pggb

Options:

Mandatory:

        -i --input-sample               input file(s) (fa | fa.gz | dir)

Optional:

        -r --runid                      name for the run. Will also name directories this.

        -mc --multiple-chromosomes      Use this parameter if the sample contains multiple chromosomes.

        -n --number-of-genomes          The number of genomes in the sample

        -p --percent-identity           The lowest similarity between all sequences in percentages

        -poa --poa-parameters           The partial order alignment parameters to use (asm5, asm10, asm20)

        -s --segment-length             Segment length for mapping [default: 10k]

        -t --threads                    Number of threads to use [default: 16]

```

4) View output at output/${runid}/

<br><br>
## Author and affiliation
<img src="https://www.uu.nl/sites/default/files/styles/original_image/public/uu-logo-nl-geenwitruimte.png" height="92" width="291"><br><br>
**Kyran Wissink**<br>Student Bioinformatics and Biocomplexity<br>Utrecht University<br>github.com/KyranWissink<br>k.wissink@students.uu.nl

