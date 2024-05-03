# FDSL_Evo
Repository for code from Streamlining Computational Fragment-Based Drug Discovery through Evolutionary Optimization Informed by Ligand-Based Virtual Prescreening

| Folder | Description | 
|--------|-------------|
| genetic | Code to run genetic algorithm |
| iterative | Code to run iterative fine tuning algorithm |
| inputs | Data used as inputs for genetic/iterative algorithms | 
| results | Contains generated ligands and associated scores |
| autodock_config_files | Files to configure Autodock VINA runs | 
| fragmentation_pipeline | Used to generate inputs, from previous work |

## Dependencies

Environment packages can be found in `environment.yml`. In addition, an install of `rdkit`, `vina=1.2.3`, and `openmpi` are required to run the scripts. `Singularity` is also necessary to run `PLIP`, which has documentation that can be found at https://github.com/pharmai/plip?tab=readme-ov-file. Version 2.2.0 was used for this project. 

## Usage

The genetic algorithm is ran using the `genetic.py` script in the genetic folder. The `genetic.sh` script in the genetic folder provides a command line template as well as SLURM parameters to run this script. To merge final results, use the `mergeResults.py` script to generate a single csv of ligands and associated binding affinities with QED scores. Inputs to this algorithm for the proteins presented in the paper are found in the inputs folder. 

The iterative algorithm is ran on the results from the genetic algorithm. The `iterative.py` script within the iterative folder coordinates this. The `iterative.sh` script contains a command line template as well as SLURM parameters to run this script. To run PLIP, be sure to change the PLIP `.smg` file to the name of the installed singularity file. Use `mergeResults.py` script to generate single csv of results.  

Note: Code assumes it is being run on a SLURM cluster, and requires certain environmental variables like $TMPDIR, $SLURM_ARRAY_TASK_ID, and $SLURM_JOB_ID which may not be present on all systems. These are typically used for file handling, and may be adapted based on individual use cases. For greater compatibility, a more portable version is under development. 

Individuals interested in generating novel fragments based off additional ligand libraries or target proteins should view the fragmentation_pipeline folder, which is described at https://doi.org/10.1016/j.jmgm.2023.108669. Note this code is separate from the project described in this paper and has it's own README.txt file. A dataset of prescreened ligands from Autodock VINA is needed to generate fragments using this methodology. 

Feel free to reach out to the authors with any questions. 
