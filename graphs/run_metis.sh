#!/bin/bash                                                         
# 
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kameranis/overlap/graphtools/graphs/social_metis.out
#SBATCH --error=/scratch/kameranis/overlap/graphtools/graphs/social_metis.err
#SBATCH --chdir=/scratch/kameranis/overlap/graphtools/graphs
#SBATCH --partition=fast
#SBATCH --exclusive
#SBATCH --job-name=run_metis
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=12:00:00

# bash run_metis.sh social social_results

for filename in $1/*.metis
do
    dataset=${filename%.metis}
    set=${dataset##*/}
    echo ${set}
    gpmetis ${filename} 2 -ufactor 100 -ncuts 2000 -niter 1000
    mv ${filename}.part.2 $2/${set}_metis_100.ptn
done
