#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50M
#SBATCH --account=rosenMRIPrj
#SBATCH --time=12:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
#SBATCH --array=1-903
. ~/.bashrc
conda init bash
conda activate CONDA_ENV
ID=$(($SLURM_ARRAY_TASK_ID - 1))
FILES=merge_out/$ID-out/xmls
mkdir merge_out/$ID-out/xmls-2
python seperateXMLS.py $FILES merge_out/$ID-out/xmls-2
