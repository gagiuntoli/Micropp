#!/bin/bash

rm -rf t_ass.dat t_sol.dat sizes.dat

echo "3D"
sizes=(10 20 30 40 50 60)
tstep=5

rm -rf t_ass.dat t_sol.dat sizes.dat

for i in ${sizes[@]}; do

  echo "running case n = " $i

  ./test/test3d_perf_1 $i $i $i 1 $tsteps > out.dat

  echo $i >> sizes.dat
  awk '/ASSEMBLY     ::/{print $3}' out.dat >> t_ass.dat
  awk '/SOLVE        ::/{print $3}' out.dat >> t_sol.dat

done

paste sizes.dat t_ass.dat t_sol.dat > t_vs_n.dat
