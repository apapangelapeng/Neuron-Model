dif=6/10
for (( i=0; i<=10; i+=2 ))
do
    for (( g=0; g<=30; g+=5 ))
    do 
         for (( d=0; d<=10; d+=1 ))
            do 
                current=$(($i/100))
                gamma=$(($g/10))
                dif=$((dif+1/5/(1/5)**$d))
                echo "running  value varying for current $current, gamma $gamma,diffusion $dif"
                ./propogate_ap $current $g $dif
                python3 ../jackson_graphing/vid_gen.py -i $current  -g $gamma -d $dif
            done
    done
done
./propogate_ap $i $g $d
python3 ../jackson_graphing/vid_gen.py -i $i  -g $g -d $d