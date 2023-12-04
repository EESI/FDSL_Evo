#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100M
#SBATCH --time=1:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --array=0-1
FILE=pdb_from_pdbqt/$SLURM_ARRAY_TASK_ID.pdb
echo $TMPDIR
cp -R $FILE $TMPDIR
ls -lR $TMPDIR
FOLDER='pdb'
ext='.pdb'
ID=$(($SLURM_ARRAY_TASK_ID))

time python seperating_v2.py $TMPDIR/$SLURM_ARRAY_TASK_ID$ext $TMPDIR $ID
ls -lR $TMPDIR
cp -r $TMPDIR improvedLigandsStruct/
 