#!/bin/bash -l        
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=25g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=west0883@umn.edu 
#SBATCH -p ram256g
#SBATCH --array=1087,1088,1096


# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID

module load matlab
matlab -nodisplay -nodesktop -r "calculate_ICs_forMSI($SLURM_ARRAY_TASK_ID); exit;"

