#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --account=rosenMRIPrj
#SBATCH --time=2:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def

time python python_scripts/countResiduePairFrequency.py ../main_pickle/combine_pickle.pickle