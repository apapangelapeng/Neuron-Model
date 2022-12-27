for (( s=0; s<=100; s+=1 ))
do
    ./static_ap $s
done
bash to_python_script.sh
convert -delay 5 -loop 0 $(ls -1 ../graphs/voltage/*.png | sort -V) -quality 95 voltage.mp4
convert -delay 5 -loop 0 $(ls -1 ../graphs/current/*.png | sort -V) -quality 95 current.mp4