#!/bin/bash -l        
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=200g
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=west0883@umn.edu 
#SBATCH -p ram1t
module load matlab
matlab -nodisplay -r "maxNumCompThreads(1)" < SVD_forMSI_short.m ; exit;