#!/bin/bash

EXEC=../../build_release/test/test3d_1 

i=0
out_folder="output_"$i

until [ ! -d $out_folder ]
do
   i=$(($i+1))
   out_folder="output_"$i
done
echo $out_folder
mkdir $out_folder
cd $out_folder

function generic {

# Increases the size independently as nx , ny , nz

nx_size=( 2 4 12 )
ny_size=( 2 4 12 16 )
nz_size=( 2 4 )
dir=1 # sigxx[0] sigyy[1] sigzz[2] sigxy[3] sigxz[4] sigyz[5] 

for nx in ${nx_size[@]}; do 
 for ny in ${ny_size[@]}; do 
  for nz in ${nz_size[@]}; do 
   echo $nx $ny $nz
   $EXEC $nx $ny $nz $dir 130 > out.dat
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

#gnuplot file

}

function uniform {

# Increases the size as nx = ny = nz = nn

nn_size=( 2 4 6 8 )
dires=( 0 1 2 3 4 5 ) # sigxx[0] sigyy[1] sigzz[2] sigxy[3] sigxz[4] sigyz[5] 

for d in ${dires[@]}; do
  for nn in ${nn_size[@]}; do 
     echo "dir = " $d " nn = " $nn
     $EXEC $nn $nn $nn $d 100 > /dev/null
     mv micropp_eps_sig_ctan.dat micropp_epssig_${nn}.dat
  done

  echo "set term png ; set output \"curves_$d.png\";" > gnuplot.in
  echo "plot \\" >> gnuplot.in

  for nn in ${nn_size[@]}; do 
     echo "\"micropp_epssig_${nn}.dat\" u $(($d+1)):$(($d+7)) w lp title '${nn}',\\" >> gnuplot.in
  done

  gnuplot gnuplot.in
done

}

#generic
uniform
cd -
