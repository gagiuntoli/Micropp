#!/bin/bash

files=$(find . -name "micropp-profiling*")


for i in ${files[@]}; do

	rank=$(echo $i | grep -o "[0-9]*")

	#echo $i $rank

done
