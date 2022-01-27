#!/bin/bash                                                  
#
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kameranis/overlap/graphtools/lambda_%j_%N.out
#SBATCH --error=/scratch/kameranis/overlap/graphtools/lambda_%j_%N.err
#SBATCH --chdir=/scratch/kameranis/overlap/graphtools/osvv/final
#SBATCH --partition=fast
#SBATCH --job-name=run_lambda
# #SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=8000
#SBATCH --time=20:00:00
#SBATCH --array=0-335

lambda_num=(1 1 2 1 1 1 1 1 1 1 1  1  1  1  1  1  1  1   1   1   1)
lambda_den=(1 2 5 3 4 5 6 7 8 9 10 15 20 30 40 50 75 100 150 200 300)
k=(1 2 3 4)
eta=(10 100 1000 10000)

l_length=21
k_length=4
eta_length=4

l_index=$((SLURM_ARRAY_TASK_ID % l_length))
k_index=$(((SLURM_ARRAY_TASK_ID / l_length) % k_length))
eta_index=$((SLURM_ARRAY_TASK_ID / (l_length * k_length)))

LD_PRELOAD=/lib/x86_64-linux-gnu/libstdc++.so.6:/lib/x86_64-linux-gnu/libexpat.so matlab -nodisplay -nosplash -nodesktop -r "runDirectoryGraphs('$1', '$2', ${lambda_num[$l_index]}, ${lambda_den[$l_index]}, $3, ${k[$k_index]}, ${eta[$eta_index]}); exit;" > /scratch/kameranis/overlap/graphtools/graphs/lambda_${lambda_num[$l_index]}_${lambda_den[$l_index]}_$3_${k[$k_index]}_${eta[$eta_index]}.out 2> /scratch/kameranis/overlap/graphtools/graphs/lambda_${lambda_num[$l_index]}_${lambda_den[$l_index]}_$3_${k[$k_index]}_${eta[$eta_index]}.err

