Below are the steps to run the pipeline. Ensure each step is complete before moving on to the following step. 

Steps:
1. Start with Autodock output files in pdbqt_files

2. Run pdbqt_to_pdb_converter.sh --> pdb_from_pdbqt

3. Run seperatingTiming.sh (seperates models in autodock output) --> improvedLigandsStruct
    - set job array from 0-# of files in pdb_from_pdbqt

4. In improvedLigandsStruct, remove original .pdb files in subdirectories (NOT lowest child directories)
    - if in parent directory, run: "rm improvedLigandsStruct/*/*.pdb"

5. Run list.sh --> relativePathNames.txt

6. Run onlyDirNames.py --> theNames
    - run "time python onlyDirNames.py relativePathNames.txt"

7. Run merge_timing.sh 
    - Set job array from 1- number of lines in dirNames.txt + 1
    - Set the following variables in merge_timing.sh:
        FILE= receptor pdb
        Line 33: change .pdb to receptor pdb
    - In merge_processorv2.py
        Change receptor name to receptor you are currently working on
    - Note:
        All pdb output files will be written into the parent directory, so avoid using an editor that constantly lists the contents of your parent directory for this step. 
        It is recommended to avoid listing the parent directory in this step as it will take too long.
        Run the next script without listing the contents of the parent directory.

8. Run sort.sh --> merge_out
    - Set receptor name and change loop length to number of lines in dirNames.txt + 1

9. Run singularity.sh --> xmls in each directory in merge_out
    - Make sure plip_v2.2.0.simg is in parent directory
    - Set job array length to number of directories in merge_out

10. Run extractInfo.sh --> main.csv and main.pickle in each directory in merge_out
    - Load rdkit module
    - If running receptors RelA and 3F4M, use the python script extractSmile_RelA_3F4M.py on line 29 of extractInfo.sh
    - If running Sprotein receptor:
        Use python script extractSmile_Sprotein.py on line 29 of extractInfo.sh
        Before running extractInfo.sh, run seperateXMLS_Sprotein.sh
        In extractInfo.sh, replace "xmls" with "xmls-2"

11. Run moveMain.sh --> main_pickle
    - Change length in for loop to number of directories in merge_out

12. Run combinePickleExecute.sh --> combine_pickle.pickle and combine.csv in main_pickle
    - Add absolute path of main_pickle as an argument in combinePickleExecute.sh
    - This script's memory usage varries greatly between runs of the pipeline- if not completing, check memory usage and increase accordingly. 

The next set of scripts are located in the Extracting_Results directory.

14. Run countResiduePairFrequency.sh --> residueFrequency.csv and residueFrequency.pickle

15. Run countOrigFragPairFrequency.sh  --> origFragPairFrequency.csv and origFragPairFrequency.pickle

16. Run combineResidues.sh --> combinedResidues.csv

17. Run combineFragOnes.sh --> combinedOrigFragOnes.csv

The following four scripts generate final results

18. Run getDistributiveNoStats.sh --> SpecificToResiduesResults.csv (in parent directory)
    - Change second to last argument in script to absolute path of parent directory 
    - Change last argument to name of receptor

19. Run getOriginalFragmentResults.sh --> GenericToSpecificResults.csv (in parent directory)
    - Change last argument in script to absolute path of parent directory
    - Change last argument to name of receptor

20. Run getOriginalFragOneResults.sh --> CombinedGenericToSpecific.csv (in parent directory)
    - Change last argument in script to absolute path of parent directory
    - Change last argument to name of receptor

21. Run getCombinedResidueResults.sh --> CombinedSpecificToResiduesResults.csv (in parent directory)
    - Change last argument in script to absolute path of parent directory
    - Change last argument to name of receptor

22. To get the top 25% and top 10% of results (sorted by highest binding affinity and then count), run getHighBindingAffinitesByCount.py --> SortedResults
    - Add one argument for receptor name and one for csv file to sort (in that order)
