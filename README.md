# MPI
Microbial Pangenome Interpretation
---------------------------------------------------------
## Setup and quick start

It is recommended to set up a new conda environment for MPI.

1) Clone this repo into your working directory
```bash
git clone https://www.github.com/KyranWissink/MPI 
cd MPI
```

2) Set up the conda environment
```bash
mamba env create -f environment.yml
mamba activate mpi
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

