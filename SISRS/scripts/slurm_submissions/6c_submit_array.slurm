#!/bin/bash
#SBATCH --job-name="6c_step"
#SBATCH --time=120:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=36   # CHANGE processor core(s) per node
#SBATCH --mail-user="erenada@uri.edu" #CHANGE user email address
#SBATCH --mail-type=ALL
#SBATCH --exclusive
#SBATCH --array=[0-31]%4 #CHANGE like 6b - bracketed numbers indicate number of total jobs you need(taxa-1) - 0 based, inclusive; following number = simultaneous
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

cd $SLURM_SUBMIT_DIR

#for advice on array jobs see https://github.com/nreid/using_array_jobs

module purge

#CHANGE if not on a URI system
module load Python/3.7.4-GCCcore-8.3.0
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load Bowtie2/2.3.5.1-GCC-8.3.0
module load FastQC/0.11.8-Java-1.8
module load BBMap/38.81-foss-2019b-Java-1.8
module load Biopython/1.75-foss-2019b-Python-3.7.4
module load Ray/2.3.1-foss-2019b
module load SAMtools/1.10-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

PROCESSORS=36 #CHANGE to the number of processors
OUTFOLDER=/data/schwartzlab/eren/Chapter2 #CHANGE to the analysis directory

SPP=($(cat ${OUTFOLDER}/TaxonList.txt))
echo ${SPP[@]}
echo ${SPP[$SLURM_ARRAY_TASK_ID]}

python3 /data/schwartzlab/eren/Chapter2/scripts/sisrs_06c_align2.py -d $OUTFOLDER -p $PROCESSORS -f ${SPP[$SLURM_ARRAY_TASK_ID]}
