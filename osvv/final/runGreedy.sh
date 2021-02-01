#!/bin/bash                                                         
# 
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kameranis/overlap/graphtools/graphs/social_greedy.out
#SBATCH --error=/scratch/kameranis/overlap/graphtools/graphs/social_greedy.err
#SBATCH --chdir=/scratch/kameranis/overlap/graphtools/osvv/final/
#SBATCH --partition=fast
#SBATCH --exclusive
#SBATCH --job-name=run_greedy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --array=0-1

# bash run_metis.sh social social_results

func=(sweepCut kernighanLin)
name=(SweepCut KernighanLin)
echo "matlab -nodisplay -nosplash -nodesktop -r \"runAlgorithm('../../graphs/social/', '../../graphs/social_results', @${func[$SLURM_ARRAY_TASK_ID]}, '${name[$SLURM_ARRAY_TASK_ID]}'); exit;\" > ../../graphs/SC.out 2> ../../graphs/SC.err;"
matlab -nodisplay -nosplash -nodesktop -r "runAlgorithm('../../graphs/social/', '../../graphs/social_results', @${func[$SLURM_ARRAY_TASK_ID]}, '${name[$SLURM_ARRAY_TASK_ID]}'); exit;" > ../../graphs/SC_$SLURM_ARRAY_TASK_ID.out 2> ../../graphs/SC_$SLURM_ARRAY_TASK_ID.err;
