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


for filename in $1/*.lemon
do
    dataset=${filename%.lemon}
    set=${dataset##*/}
    # nodes=`head -n 1 ${dataset}.eg2 | awk '{print $1;}'`
    #if [ "$nodes" -gt "20000" ]
    #then
    #    continue
    #fi
    python /home/kameranis/Doctoral/LEMON/LEMON/LEMON.py -d $'\t' -f ${filename} -g ${dataset}.cmty.txt --out $2/${set}_lemon.txt --sd ${dataset}.seed
    # echo `realpath ${filename}`
done
