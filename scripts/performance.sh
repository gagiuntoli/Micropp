#!/bin/bash

source vars.sh

sizes=( 10 20 )

PATH_OUT="./result"
EXEC="../build/test/test_newton_linear"

mkdir -p ${PATH_OUT}

for n in ${sizes[@]}; do

	FILE="${PATH_OUT}/out_${n}.txt"
	echo "${EXEC} $n"
	${EXEC} $n > ${FILE}

done
