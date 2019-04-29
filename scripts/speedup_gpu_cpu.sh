#!/bin/bash

CPU_EXEC="../build_judi_cpu/test/test3d_1"
GPU_EXEC="../build_judi_gpu/test/test3d_1"

sizes=( 25 50 75 100 )

for i in ${sizes[@]}; do

	echo "CPU $i"
	time $CPU_EXEC $i 1 1 &> cpu-${i}.txt

	echo "GPU $i"
	time $GPU_EXEC $i 1 1 &> gpu-${i}.txt

done
