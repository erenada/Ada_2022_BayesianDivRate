#!/bin/sh
#SBATCH --job-name="7_b"
#SBATCH --time=90:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=10   # processor core(s) per node
#SBATCH --mail-user="erenada@uri.edu" #user email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
cd $SLURM_SUBMIT_DIR

module purge

module load MAFFT/7.475-gompi-2020b-with-extensions

python sisrs_07_b_contigs_alignment.py -p 20 -dir /data/schwartzlab/eren/Chapter2
