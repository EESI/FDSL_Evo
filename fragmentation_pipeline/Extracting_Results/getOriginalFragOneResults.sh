#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
. ~/.bashrc
conda init bash
conda activate CONDA_ENV
python python_scripts/getOriginalFragOneResults.py ../main_pickle/combine.csv combinedOrigFragOne.csv /PATH_TO/TOP_FOLDER/ RelA
