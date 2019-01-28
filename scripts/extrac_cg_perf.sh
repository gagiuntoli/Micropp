#!/bin/bash

source vars.sh

mkdir -p result

PATH_OUT="./result"

for a in ${factor[@]}; do

	echo "#   N     assembly   solve   total" > tvsn_1solve_2rhs_a${a}.txt
	echo "#   N     assembly   solve   total" > tvsn_1solve_2rhs_1mat_a${a}.txt

	for n in ${sizes[@]}; do
	
		file="${PATH_OUT}/out_n${n}_a${a}.txt"
	
		assembly_rhs_time=$(awk '/assembly_rhs/{print $4; exit}' $file)
		assembly_mat_time=$(awk '/assembly_mat/{print $6; exit}' $file)
		solver_time=$(awk '/ell_solve_cgpd/{print $4; exit}' $file)
		total_time=$((solver_time + assembly_rhs_time))
		solver_its=$(awk '/RES/{print $6}' $file)
	
		printf "Extract from ${file}\n"
	
		printf "%-10d\t%-10d\t%-10d\t%-10d\n" \
			$(( n * n * n)) \
			${assembly_rhs_time}\
			${solver_time} \
			$((solver_time + assembly_rhs_time)) \
			>> tvsn_1solve_2rhs_a${a}.txt

		printf "%-10d\t%-10d\t%-10d\t%-10d\n" \
			$(( n * n * n)) \
			$((assembly_rhs_time + assembly_mat_time)) \
			${solver_time} \
			$((solver_time + assembly_rhs_time + assembly_mat_time)) \
			>> tvsn_1solve_2rhs_1mat_a${a}.txt
	
		awk -v tot_time=${solver_time} -v tot_its=${solver_its} '
		/^cgpd/{printf("%-4d\t%e\t%e\n",
		$5, ($5 * tot_time * 1.e-6/ tot_its), $8);
		}' $file > rvsi_n${n}_a${a}.txt
	
	done
done
