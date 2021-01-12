#!/bin/bash                                                         
# 
#SBATCH --mail-user=kameranis@uchicago.edu
#SBATCH --mail-type=ALL
#SBATCH --output=
#SBATCH --error=
#SBATCH --workdir=
#SBATCH --partition=fast
#SBATCH --exclusive
#SBATCH --job-name=run_kahip
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=1000
#SBATCH --time=3:00:00

# bash run_kahip.sh social social_result

for filename in $1/*.kahip
do
    dataset=${filename%.kahip}
    set=${dataset##*/}
    echo ${set}
    mpirun -n 1 parhip $filename --imbalance 100 --k 2 --preconfiguration=fastsocial --save_partition
    mv tmppartition.txtp $2/${set}_kahip_edge.ptn

    node_separator $filename --imbalance 100 --output_filename=$2/${set}_kahip_node.ptn
    sed -i 's/2/0 1/g' $2/${set}_kahip_node.ptn
done
