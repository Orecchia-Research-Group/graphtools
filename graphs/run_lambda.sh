#!/bin/bash                                                  
#
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kameranis/overlap/graphtools/lambda.out
#SBATCH --error=/scratch/kameranis/overlap/graphtools/lambda.err
#SBATCH --chdir=/scratch/kameranis/overlap/graphtools/osvv/final
#SBATCH --partition=fast
#SBATCH --job-name=run_lambda
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=60000
#SBATCH --time=20:00:00
#SBATCH --array=1-10%3


echo "runDirectoryGraphs('$1', '$2', 2, [1/$SLURM_ARRAY_TASK_ID])"
matlab -nodisplay -nosplash -nodesktop -r "runDirectoryGraphs('$1', '$2', 2, [1/$SLURM_ARRAY_TASK_ID]); exit;" > /scratch/kameranis/overlap/graphtools/graphs/lambda_$SLURM_ARRAY_TASK_ID.out 2> /scratch/kameranis/overlap/graphtools/graphs/lambda_$SLURM_ARRAY_TASK_ID.err

