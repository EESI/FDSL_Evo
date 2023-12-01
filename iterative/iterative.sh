#!/bin/bash
#SBATCH --nodes=6
#SBATCH --ntasks=41
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --account=rosenMRIPrj
#SBATCH --time=3:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
#SBATCH --array=0-100


. ~/.bashrc > /dev/null
conda init bash > /dev/null
conda activate CONDA_ENV

mkdir graphs

AUTODOCKCONFIG=/PATH_TO/AUTODOCK_CONFIG
PDB=/PATH_TO/PROTEIN_PDB
PDBQT=/PATH_TO/PROTEIN_PDBQT
STARTER=/PATH_TO/STARTER_LIGANDS
AADICT=/PATH_TO/AA_DICTIONARY #generated by generateFragmentBank.py
OUTPUT=/PATH_TO/OUTPUT_DIRECTORY

mpiexec -n 41 time python iterative.py -i $SLURM_ARRAY_TASK_ID -c $AUTODOCKCONFIG -p $PDBQT -x $PDB -s $STARTER -a $AADICT -o $OUTPUT