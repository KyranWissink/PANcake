# MPI
Microbial Pangenome Interpretation
---------------------------------------------------------
## Setup and quick start

1) Clone this repo into your working directory
```bash
git clone https://www.github.com/KyranWissink/MPI 
cd MPI
```

2) Set up the conda environment
```bash
conda create --name mpi --file spec-file.txt
conda activate mpi
```

3) Run the setup
```bash
bash scripts/setup.sh
```

4) Run bash run.sh for the parameters
```bash
MPI Version: 0.7

Usage:

        bash run.sh [options] -i <input> -o <output_dir>

Description:

        Pipeline for pangenome graph creation using pggb

Options:

Mandatory:

        -i --input-sample               input fasta file or directory. Zipped files are supported.

        -o --output-directory           directory to output all files. Need not be yet created.

Optional:

        -mc --multiple-chromosomes      Use this parameter if the sample contains multiple chromosomes.

        -n --number-of-genomes          The number of genomes in the sample

        -p --percent-identity           The lowest similarity between all sequences in percentages

        -poa --poa-parameters           The partial order alignment parameters to use (asm5, asm10, asm20)

        -s --segment-length             Segment length for mapping [default: 10k]

        -t --threads                    Number of threads to use [default: 16]
Version: 0.7

Usage:

        bash run.sh [options] -i <input> -o <output_dir>

Description:

        Pipeline for pangenome graph creation using pggb

Options:

Mandatory:

        -i --input-sample               input fasta file or directory. Zipped files are supported.

        -o --output-directory           directory to output all files. Need not be yet created.

Optional:

        -mc --multiple-chromosomes      Use this parameter if the sample contains multiple chromosomes.

        -n --number-of-genomes          The number of genomes in the sample

        -p --percent-identity           The lowest similarity between all sequences in percentages

        -poa --poa-parameters           The partial order alignment parameters to use (asm5, asm10, asm20)

        -s --segment-length             Segment length for mapping [default: 10k]

        -t --threads                    Number of threads to use [default: 16]

```
Output will be stored at output/{runid}/

