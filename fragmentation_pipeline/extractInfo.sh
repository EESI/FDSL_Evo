#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50M
#SBATCH --account=rosenMRIPrj
#SBATCH --time=3:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def
#SBATCH --array=1-5
. ~/.bashrc
conda init bash
conda activate CONDA_ENV
ID=$(($SLURM_ARRAY_TASK_ID - 1))
FILES=merge_out/$ID-out/xmls
LIST=($(<relativePathNames.txt))
TARGET=${LIST[$ID]}
FILES2=$TARGET
LIST2=($(<dirNames.txt))
TARGET2=${LIST2[$ID]}
cp -R $FILES $TMPDIR
mkdir $TMPDIR/Ligands
cp -R $FILES2 $TMPDIR/Ligands
ls -R $TMPDIR
ls $TMPDIR/Ligands

time python extractSmile_RelA_3F4M.py $TMPDIR/xmls/ $TMPDIR/Ligands/$TARGET2/ $TMPDIR
cp $TMPDIR/main.csv merge_out/$ID-out/
cp $TMPDIR/main.pickle merge_out/$ID-out/
