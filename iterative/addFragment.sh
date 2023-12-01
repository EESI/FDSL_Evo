#!/bin/bash

. ~/.bashrc > /dev/null
conda init bash > /dev/null
conda activate CONDA_ENV > /dev/null

WD=$(pwd)
PROTEINFILENAME=$(basename "$4")
PNAME="${PROTEINFILENAME%.*}"
echo $PNAME
AADICT=$5

startPLIP=$(date +%s)
cp $TMPDIR/pdbForPLIP/$PNAME-dehyd-$3_$2_results-prep.pdb $TMPDIR/pdbForPLIP/$2/$PNAME-dehyd-$1_$2_results-prep.pdb > /dev/null
cd $TMPDIR/pdbForPLIP/$2 > /dev/null
singularity run --bind .:/mnt plip_v2.2.0.simg -f /mnt/$PNAME-dehyd-$1_$2_results-prep.pdb -o /mnt -x > /dev/null
mv report.xml $1_$2.xml > /dev/null
endPLIP=$(date +%s)

echo "FINISHED MERGE"

startADDITION=$(date +%s)
obabel -i pdb $PNAME-dehyd-$1_$2_results-prep.pdb -o mol -O $PNAME-dehyd-$1_$2_results-prep.mol > /dev/null
cd $WD > /dev/null
ls -l $TMPDIR/pdbForPLIP/$2 >> output/$2.txt

python addFragmentToSpecificAtom.py -x $TMPDIR/pdbForPLIP/$2/$1_$2.xml -p $TMPDIR/pdbForPLIP/$2/$PNAME-dehyd-$1_$2_results-prep_protonated.pdb -m $TMPDIR/pdbForPLIP/$2/$PNAME-dehyd-$1_$2_results-prep.mol -g $1 -r $2 -a $AADICT >> output/$2.txt 
ls -l $TMPDIR/pdbForPLIP/$2 >> output/$2.txt
endADDITION=$(date +%s)

echo "Elapsed Time PLIP: $(($endPLIP-$startPLIP)) seconds" >> output/$2.txt
echo "Elapsed Time ADDITION: $(($endADDITION-$startADDITION)) seconds" >> output/$2.txt
