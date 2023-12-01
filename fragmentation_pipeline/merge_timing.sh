#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=90M
#SBATCH --account=rosenMRIPrj
#SBATCH --time=2:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
#SBATCH --array=1-5
. ~/.bashrc
conda init bash
conda activate CONDA_ENV
FILE=RelA.pdb
SCRIPT=merge_processorv2.py
cp $SCRIPT $TMPDIR
cp -r src/ $TMPDIR
cp merge-it.py $TMPDIR
echo $TMPDIR
cp -R $FILE $TMPDIR
ls -lR $TMPDIR
#ls $TMPDIR/improvedLigandsStruct
LIST=($(<relativePathNames.txt))
ID=$(($SLURM_ARRAY_TASK_ID - 1))
TARGET=${LIST[$ID]}
cp -R $TARGET $TMPDIR
echo $TARGET
ls -R $TMPDIR
NAMEB=$ID$FILE
mv $TMPDIR/RelA.pdb $TMPDIR/$NAMEB
mkdir $TMPDIR/Combined
mkdir $TMPDIR/Combined/$ID	
NAMES=($(<dirNames.txt))
NAME=${NAMES[$ID]}
time python $TMPDIR/$SCRIPT $TMPDIR/$NAME  $TMPDIR/$NAMEB $TMPDIR $ID

