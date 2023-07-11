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
conda create --name env_name --file pkgs/spec-file.txt
conda activate env_name
```

3) Place your data in the /data folder.

4) Set up the config.yaml correctly. 
This file is in Yaml format, so the layout should not be changed.
```bash
multiple_chromosomes: 0			# MANDATORY 0=no, 1=yes
pggb:					  # PGGB parameters
  haplotypes: 				# OPTIONAL number of haplotypes 
  percent_identity: 			# OPTIONAL lowest similarity percentage between haplotypes
  poa_params: 				# OPTIONAL partial order alignment parameters (based on percent identity)
  segment_length: 3000			# MANDATORY length of the mapped and aligned segment
  threads: 128				# MANDATORY number of threads to use
runid: dep_test				# MANDATORY the runid (the output folder will have this name)
sample: data/testdata/DRB1-3123.fa	# MANDATORY your input data
```
