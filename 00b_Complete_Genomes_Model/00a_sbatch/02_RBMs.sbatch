#!/bin/bash
#SBATCH -J RBMs
#SBATCH -q inferno
#SBATCH -A gts-ktk3-coda20
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH --mem=20G
#SBATCH -t 12:00:00
#SBATCH -o 00_log/02a_RBMs.out
#SBATCH -e 00_log/02a_RBMs.err

# Change to working directory
cd $SLURM_SUBMIT_DIR

# Set python script directory
scripts=/storage/scratch1/9/rconrad6/01_recomb_tests/01c_scripts

# Load Modules Needed
module load parallel/20210922-qv5jj7
module load blast-plus/2.13.0

# check output dir exists or create it
if [ ! -d ${oDir} ]; then mkdir ${oDir}; fi

# run fastANI in parallel
parallel -j 20 --joblog 00_log/02b_parallel.log --colsep ' '  '${scripts}/02b_get_RBMs.py -g1 {1} -g2 {2} -o {3}/{1/.}-{2/.}.rbm' :::: ${infile} ::: ${oDir}

## Run Log:
## sbatch --export oDir=02d_RBM,infile=02e_Allv_gene_file_names.txt ../01b_sbatch/02_RBMs.sbatch