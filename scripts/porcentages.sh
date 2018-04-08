#!/bin/bash

echo "2D"
sizes=(10 20 30 40 50 60 70 80 90 100 120 140 160 180 200 220 240)

rm -rf t_ass.dat t_sol.dat sizes.dat

for i in ${sizes[@]}; do

  echo "running case n = " $i
  ../test_1 --nx $i --ny $i --cg_its 2500 > out.dat
  echo $i >> sizes.dat
  awk '/time assembly/{print $4}' out.dat >> t_ass.dat
  awk '/time solve/{print $4}   ' out.dat >> t_sol.dat

done

paste sizes.dat t_ass.dat t_sol.dat > t_vs_n_3d.dat

echo "3D"
sizes=(10 20 30 40 50 60 70 80)

rm -rf t_ass.dat t_sol.dat sizes.dat

for i in ${sizes[@]}; do

  echo "running case n = " $i
  ../test_4 --nx $i --ny $i --nz $i  --cg_its 10000 > out.dat
  echo $i >> sizes.dat
  awk '/time assembly/{print $4}' out.dat >> t_ass.dat
  awk '/time solve/{print $4}   ' out.dat >> t_sol.dat

done

paste sizes.dat t_ass.dat t_sol.dat > t_vs_n_3d.dat
