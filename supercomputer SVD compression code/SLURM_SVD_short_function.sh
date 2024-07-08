#!/bin/bash -l        
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=200g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=west0883@umn.edu 
#SBATCH -p ram1t
#SBATCH --array=109


# Use '%A' for array-job ID, '%J' for job ID and '%a' for task ID

module load matlab
matlab -nodisplay -nodesktop -r "SVD_forMSI_function($SLURM_ARRAY_TASK_ID); exit;"
