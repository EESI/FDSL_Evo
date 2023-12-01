#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50M
#SBATCH --account=rosenMRIPrj
#SBATCH --time=2:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
#SBATCH --array=1-5
### the array size should be the number of directories you have as a result of sort.sh OR the size of listOFNames + 1
ID=$(($SLURM_ARRAY_TASK_ID - 1))
FILES=merge_out/$ID-out
cp -R $FILES $TMPDIR
ls -R $TMPDIR

time python singularityCommand.py $TMPDIR/$ID-out/
mkdir $ID-out/xmls
mv $TMPDIR/$ID-out/xmls merge_out/$ID-out
