#!/bin/bash                                                         
# 
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/home/kameranis/overlap/graphtools/bigclam_run.txt
#SBATCH --error=/home/kameranis/overlap/graphtools/bigclam_error.txt
#SBATCH --workdir=/scratch/kameranis/overlap/graphtools/graphs/walshaw
#SBATCH --partition=fast
#SBATCH --exclusive
#SBATCH --job-name=run_bigclam
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=3:00:00


for filename in $1/*.edgelist
do
    dataset=${filename%.edgelist}
    set=${dataset##*/}
    echo ${set}
    /home/kameranis/Doctoral/snap/examples/Release/bigclam -i:${filename} -o:$2/${set}_ -nt:4 -c:200
    # echo `realpath ${filename}`
done
