for filename in graphs/*.edgelist
do
    dataset=${filename%.edgelist}
    echo $filename
    deepwalk --format edgelist --input $filename --max-memory-data-size 0 --number-walks 50 --representation-size 2 --walk-length 40 --workers 4 --output ${dataset}_deepwalk.txt
done
