#!/bin/bash

startSETUP=$(date +%s)
. ~/.bashrc > /dev/null
conda init bash > /dev/null
conda activate CONDA_ENV > /dev/null

AUTODOCKCONFIG=$3
PROTEINPDBQT=$4
PROTEINPDB=$5

endSETUP=$(date +%s)

startPREVINA=$(date +%s)
INPUTFILE=($TMPDIR/molFiles/$1_$2.mol) > /dev/null
WD=$(pwd)
mkdir $TMPDIR/input > /dev/null
mkdir $TMPDIR/output > /dev/null
obabel -i mol $INPUTFILE -o pdbqt -O $TMPDIR/input/$1_$2.pdbqt > /dev/null
endPREVINA=$(date +%s)

startVINA=$(date +%s)
#cp $TMPDIR/input/$1_$2.pdbqt /PATH_TO/INPUT_STORAGE
vina --receptor $PROTEINPDBQT --ligand $TMPDIR/input/$1_$2.pdbqt --config $AUTODOCKCONFIG --seed 215 --out $TMPDIR/output/$1_$2_out.pdbqt | awk '/^   1/ {print $2}' 
#cp $TMPDIR/output/$1_$2_out.pdbqt /PATH_TO/OUTPUT_STORAGE

endVINA=$(date +%s)

startMERGE=$(date +%s)
obabel -i pdbqt $TMPDIR/output/$1_$2_out.pdbqt -o pdb -O $TMPDIR/pdbForPLIP/$1_$2_results.pdb > /dev/null
python merge-it.py -ir $PROTEINPDB -if $TMPDIR/pdbForPLIP/$1_$2_results.pdb >> output/$2.txt
endMERGE=$(date +%s)

echo "Elapsed Time SETUP: $(($endSETUP-$startSETUP)) seconds" >> output/$2.txt
echo "Elapsed Time SETUP: $(($endPREVINA-$startPREVINA)) seconds" >> output/$2.txt
echo "Elapsed Time VINA: $(($endVINA-$startVINA)) seconds" >> output/$2.txt
echo "Elapsed Time MERGE: $(($endMERGE-$startMERGE)) seconds" >> output/$2.txt
