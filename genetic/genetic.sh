#!/bin/bash
#SBATCH --nodes=6
#SBATCH --ntasks=41
#SBATCH --cpus-per-task=6
#SBATCH --mem=10G
#SBATCH --account=rosenMRIPrj
#SBATCH --time=0:30:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
#SBATCH --array=0-50%5

. ~/.bashrc > /dev/null
conda init bash > /dev/null
conda activate CONDA_ENV

mkdir jobResults
mpiexec -n 41 time python genetic.py \
    -i /PATH_TO/INPUT_CSV \
    -c /PATH_TO/VINA_CONFIG \
    -p /PATH_TO/VINA_PROTEIN_PDBQT