for filename in ../synthetic/*.edgelist
do
    dataset=${filename%.edgelist}
    set=${dataset##*/}
    nodes=`head -n 1 ${dataset}.eg2 | awk '{print $1;}'`
    if [ "$nodes" -gt "20000" ]
    then
        continue
    fi
    # /share/Doctoral/snap/examples/Release/agmfitmain -i:${filename} -o:${dataset}_ -c:2 -l:
    echo `realpath ${dataset}`
done
