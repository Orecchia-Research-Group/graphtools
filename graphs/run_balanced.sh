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
#SBATCH --mem=6000
#SBATCH --time=20:00:00
#SBATCH --array=0-20

balances=(0 1 2 3 4 5 10 15 20 25 30 40 50 100 150 200 250 300 350 400 450)

matlab -nodisplay -nosplash -nodesktop -r "runDirectoryGraphs('$1', '$2', 1, 1, ${balances[$SLURM_ARRAY_TASK_ID]}); exit;" > /scratch/kameranis/overlap/graphtools/graphs/balance_${balances[$SLURM_ARRAY_TASK_ID]}.out 2> /scratch/kameranis/overlap/graphtools/graphs/balance_${balances[$SLURM_ARRAY_TASK_ID]}.err

