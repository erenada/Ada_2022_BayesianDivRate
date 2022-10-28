#!/bin/bash
#SBATCH --job-name="align_ALL50K"
#SBATCH --time=99:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to processor core(s) per node
#SBATCH --mail-user="erenada@uri.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=END,FAIL
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err
#SBATCH --exclusive

module purge
module load MAFFT/7.475-gompi-2020b-with-extensions

listOfData=$(ls /data/schwartzlab/eren/Chapter2/CONTIGS | tail -n +2 | cut -f 2 | sort | uniq | grep ALL50K)
PTH=(/data/schwartzlab/eren/Chapter2/CONTIGS/ALL50K)
PTH_OUT=(/data/schwartzlab/eren/Chapter2/CONTIGS/ALIGNED/ALL50K)

for type in $listOfData;
do
#  mkdir ${PTH_OUT}/${type}
  FILELIST=$(ls ${PTH})
  for infile in $FILELIST;
  do
    out_file=$(basename ${infile})
    mafft --nwildcard --auto --thread 20 ${infile} > ${PTH_OUT}/${out_file}
  done
done
