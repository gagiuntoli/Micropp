#!/bin/bash

EXEC=../build/test/test3d_1 

nx_size=( 2 )
ny_size=( 4 8 12 )
nz_size=( 2 )
dir=1 # sigxx[0] sigyy[1] sigzz[2] sigxy[3] sigxz[4] sigyz[5] 

for nx in ${nx_size[@]}; do 
 for ny in ${ny_size[@]}; do 
  for nz in ${nz_size[@]}; do 
   echo $nx $ny $nz
   $EXEC $nx $ny $nz $dir > out.dat
   awk '/eps/{printf("%lf ", $3)} /sig/{printf("%lf\n", $(3+dir))}' out.dat > s_vs_e_${nx}_${ny}_${nz}.dat
  done
 done
done

echo "set term png ; set output \"curves.png\";" > file
echo "plot \\" >> file

for nx in ${nx_size[@]}; do 
 for ny in ${ny_size[@]}; do 
  for nz in ${nz_size[@]}; do 
   echo "\"s_vs_e_${nx}_${ny}_${nz}.dat\" u 1:2 w lp title '$nx-$ny-$nz',\\" >> file
  done
 done
done

gnuplot file
