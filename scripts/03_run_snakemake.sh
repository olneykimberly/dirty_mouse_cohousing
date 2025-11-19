#!/bin/bash
#SBATCH --job-name=dirty_CH_processing                                                              
#SBATCH --time=10:00:00                               
#SBATCH --mem=2G
#SBATCH -n 8 # threaded 
#SBATCH --cpus-per-task=1
#SBATCH -o slurm.dirty_CH_processing.job.%j.out
#SBATCH -e slurm.dirty_CH_processing.job.%j.err

# activate conda environment
source $HOME/.bash_profile
module load python3
conda activate dirty_mice

# change directory to where Snakefile is located
CWD="/tgen_labs/jfryer/kolney/dirty_mice/dirty_mouse_cohousing/scripts"
cd $CWD

snakemake --nolock -s Snakefile --jobs 48 --executor slurm --profile slurm_profile --rerun-incomplete --default-resources mem_mb=64000 ntasks=1 threads=8 runtime=550 cpus_per_task=8
