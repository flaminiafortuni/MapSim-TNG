#!/bin/bash

#write by hand the snapshot numbers avalaible in the directory "./output" 
declare -a snapshots=("38" "44")


nl=$(cat planes_list.txt | wc -l) #write "$1" instead of planes_list.txt if input from bash
declare -a x
declare -a y
for i in $(seq 1 $nl)
do
    #plane number
    x[i]="$(cat planes_list.txt | awk -v p="$i" '{if(NR==p) print $1}')" #write "$1" instead of planes_list.txt if input from bash
    #snapshot number
    y[i]="$(cat planes_list.txt | awk -v p="$i" '{if(NR==p) print $6}')" #write "$1" instead of planes_list.txt if input from bash
done


for (( i = 0; i < $nl; i++ ))
do
    for (( j=0; j < ${#snapshots[@]}; j++ ))
    do
	if [[ "${snapshots[$j]}" -eq "${y[$i]}" ]];
	   then

	   echo "${x[$i-1]}" > restart_pl.d
	   ./TESTMapSim readTNGxMapSim FreadTNG ${y[$i]}
	   echo "from bash: snap = ${y[$i]}"   	   
	fi
    done        
done
