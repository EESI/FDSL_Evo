#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --account=rosenMRIPrj
#SBATCH --time=1:00:00
### Whatever modules you used (e.g. picotte-openmpi/gcc)
### must be loaded to run your code.
### Add them below this line.
#SBATCH --partition=def

#GOES IN FOLDER WITH PDBQT FILES
#SHOULD CHANGE TO USE TMP DIRECTORIES INSTEAD

FILES=pdbqt_files/*.pdbqt
mkdir pdb_from_pdbqt
COUNT=0
for FILE in $FILES
    do echo "Converting $FILE to pdb"
    cut -c-66 "$FILE" > "pdb_from_pdbqt/$COUNT.pdb"
    ((COUNT++))
done
