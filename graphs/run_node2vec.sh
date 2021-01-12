#!/bin/bash
# 
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=/scratch/kameranis/overlap/graphtools/node2vec_run.txt
#SBATCH --error=/scratch/kameranis/overlap/graphtools/node2vec_error.txt
#SBATCH --workdir=/scratch/kameranis/overlap/graphtools/graphs/walshaw
#SBATCH --partition=fast
#SBATCH --job-name=run_node2vec
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=5:00:00

source /scratch/kameranis/overlap/node2vec/n2vpy/bin/activate
for filename in graphs/*.edgelist
do
    dataset=${filename%.edgelist}
    echo $filename
    python /scratch/kameranis/overlap/node2vec/src/main.py --input $filename --output ${dataset}_node2vec --num-walks 20 --walk-length 20 --undirected --dimensions 2 --workers 8 --window-size 10 --q 2
done
