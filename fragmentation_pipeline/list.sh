#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5M
#SBATCH --account=rosenMRIPrj
#SBATCH --time=0:01:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def


FILES=improvedLigandsStruct

for filename in $FILES/*; do
	for f in $filename/*; do
		echo $f >> relativePathNames.txt
	done
done
