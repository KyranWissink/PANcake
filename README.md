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

3) Run the setup
```bash
bash scripts/setup.sh
```

4) Set up the config.yaml correctly. 
This file is in Yaml format, so the layout should not be changed.
```bash
multiple_chromosomes: 0/1		# MANDATORY 0=no, 1=yes
pggb:				  	# PGGB parameters
  haplotypes: n				# OPTIONAL number of haplotypes 
  percent_identity: n			# OPTIONAL lowest similarity percentage between haplotypes
  poa_params: asm5/asm10/asm20		# OPTIONAL partial order alignment parameters (based on percent identity)
  segment_length: n			# MANDATORY length of the mapped and aligned segment
  threads: n				# MANDATORY number of threads to use
runid: your_id				# MANDATORY the runid (the output folder will have this name)
sample: path/to/data			# MANDATORY your input data
```
