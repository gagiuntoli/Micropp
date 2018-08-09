#!/bin/bash
#
# Mesh converge test
#
# Uses test3d_1 : 3D 2-Layers microstructure
# with Mat1 ISOLIN and Mat2 Plastic.
#
# Increases the size as nx = ny = nz = nn
#
# Arrangement in micropp output
# epsxx [1] # epsyy [2] # epszz[3] # epsxy[4]  # epsxz[5]  # epsyz[6]
# sigxx [7] # sigyy [8] # sigzz[9] # sigxy[10] # sigxz[11] # sigyz[12]
#

nn_size=( 2 3 4 5 6 7 8 9 10 )
ts=130
dires=( 0 1 2 3 4 5 )

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

for d in ${dires[@]}; do
	
	rm -f max_${d}.dat 
	for nn in ${nn_size[@]}; do 
		echo "dir = " $d " nn = " $nn
		$EXEC $nn $nn $nn $d $ts > /dev/null
		awk -v d_a=$d 'BEGIN{getline;}{print $(d_a + 1) "\t" $(d_a + 7) }' \
			micropp_eps_sig_ctan.dat > sig_vs_eps_${d}_${nn}.dat
		sed '$!d' sig_vs_eps_${d}_${nn}.dat \
			| awk -v n_a=$nn '{print n_a "\t" $2}' >> max_${d}.dat
	done

	echo "set term png;" > gnuplot.in
	echo "set output \"curves_${d}.png\";" >> gnuplot.in
	echo "plot \\" >> gnuplot.in

	for nn in ${nn_size[@]}; do 
		echo "\"sig_vs_eps_${d}_${nn}.dat\" u 1:2 w lp title '${nn}',\\" \
			 >> gnuplot.in
	done

	gnuplot gnuplot.in
done

cd -
