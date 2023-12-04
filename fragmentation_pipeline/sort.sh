#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M
#SBATCH --time=3:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.

NAME=RelA
mkdir merge_out
for i in {0..5}; do
	mkdir merge_out/$i-out
	mv $i$NAME* merge_out/$i-out	
	echo $i

done
