#!/bin/bash                                                         
# 
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=
#SBATCH --error=
#SBATCH --workdir=
#SBATCH --partition=fast
#SBATCH --exclusive
#SBATCH --job-name=run_metis
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=3:00:00

# bash run_metis.sh social social_results 500

for filename in $1/*.metis
do
    dataset=${filename%.metis}
    set=${dataset##*/}
    echo ${set}
    gpmetis ${filename} 2 -ufactor $3 -ncuts 2000 -niter 1000
    mv ${filename}.part.2 $2/${set}_metis_$3.ptn
done
