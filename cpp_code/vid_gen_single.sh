current=$1
gamma=$2
diffusion=$3
./propogate_ap_single $current $gamma $diffusion >> /dev/null
echo "finish csv generation"

file_name=_$1_$2_$3
echo $file_name
echo "starting video gen"

python3 ../jackson_graphing/vid_gen_single.py -fn $file_name