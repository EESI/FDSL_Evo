#!/bin/bash

. ~/.bashrc
conda init bash > /dev/null
conda activate drugdesignV3

INPUTFILE=($TMPDIR/molFiles/$1_$2.mol) > /dev/null
echo $1 >> output/$2.txt
mkdir $TMPDIR/input > /dev/null
mkdir $TMPDIR/output > /dev/null
obabel -i mol $INPUTFILE -o pdbqt -O $TMPDIR/input/$1_$2.pdbqt > /dev/null
vina --receptor $4 --ligand $TMPDIR/input/$1_$2.pdbqt --config $3 --seed 215 --cpu 6 --out $TMPDIR/output/$1_$2_out.pdbqt | awk '/^   1/ {print $2}' 
