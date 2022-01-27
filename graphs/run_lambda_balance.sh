#!/bin/bash                                                  
#
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kameranis/overlap/graphtools/lambda.out
#SBATCH --error=/scratch/kameranis/overlap/graphtools/lambda.err
#SBATCH --chdir=/scratch/kameranis/overlap/graphtools/osvv/final
#SBATCH --partition=fast
#SBATCH --job-name=run_balance
# #SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 8
#SBATCH --mem=25000
#SBATCH --time=20:00:00
#SBATCH --array=0-440

balances=(0 1 2 3 4 5 10 15 20 25 30 40 50 100 150 200 250 300 350 400 450)
lambda_num=(1 1 2 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1   1   1   1)
lambda_den=(1 2 5 3 4 5 6 7 8 9 10 15 20 30 40 50 75 100 150 200 300)

b_length=21
l_length=21

b_index=$((SLURM_ARRAY_TASK_ID / l_length))
l_index=$((SLURM_ARRAY_TASK_ID % l_length))


LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6:/lib/x86_64-linux-gnu/libexpat.so matlab -nodisplay -nosplash -nodesktop -r "runDirectoryGraphs('$1', '$2', ${lambda_num[$l_index]}, ${lambda_den[$l_index]}, ${balances[$b_index]}); exit;" > /scratch/kameranis/overlap/graphtools/graphs/bl_${lambda_num[$l_index]}_${lambda_den[$l_index]}_${balances[$b_index]}.out 2> /scratch/kameranis/overlap/graphtools/graphs/bl_${lambda_num[$l_index]}_${lambda_den[$l_index]}_${balances[$b_index]}.err

