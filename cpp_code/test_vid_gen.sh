i=0.06
g=0.995
d=0.98
make propogate_ap
./propogate_ap $i $g $d
python3 ../jackson_graphing/vid_gen.py -i $i  -g $g -d $d