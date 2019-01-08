#!/bin/bash

procs=(1 2 4 8 16 32 48)

rm -rf times.dat
t1=$(awk '/time =/{print $3}' res_1.dat)

for i in ${procs[@]}; do

	awk -v t1=$t1 -v p=$i '/time =/{print p "\t" $3 "\t" (t1/$3)}' res_${i}.dat >> times.dat

done
