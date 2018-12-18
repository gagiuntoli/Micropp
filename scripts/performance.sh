#!/bin/bash

source vars.sh

PATH_OUT="./result"
EXEC="../build/test/test_newton_linear"

mkdir -p ${PATH_OUT}

for a in ${factor[@]}; do
	for n in ${sizes[@]}; do
	
		FILE="${PATH_OUT}/out_n${n}_a${a}.txt"
		echo "${EXEC} $n $a"
		${EXEC} $n $a > ${FILE}
	
	done
done
