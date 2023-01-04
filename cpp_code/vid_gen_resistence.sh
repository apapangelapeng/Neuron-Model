#!/bin/bash
current=$1
gamma=$2
for (( r=11; r<=101; r+=10 ))
    do 
    ./propogate_ap_single $current $gamma $r >> /dev/null
    echo "finish csv generation"
    echo $r

    file_name=_$1_$2_$r
    echo $file_name
    echo "starting video gen"

    python3 ../jackson_graphing/vid_gen_single.py -fn $file_name
done