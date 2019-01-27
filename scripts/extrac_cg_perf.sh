#!/bin/bash

source vars.sh

mkdir -p result

PATH_OUT="./result"
AWK="/usr/bin/awk"

for a in ${factor[@]}; do

	echo "#   N     assembly   solve   total" > tvsn_linear_a${a}.txt
	/bin/rm -f tvsn_linear_a${a}.txt
	for n in ${sizes[@]}; do
	
		file="${PATH_OUT}/out_n${n}_a${a}.txt"
	
		assembly_time=$(${AWK} '/assembly_rhs/{print $6; exit}' $file)
		assembly_time=$((assembly_time * 2))
		solver_time=$(${AWK} '/ell_solve_cgpd/{print $6; exit}' $file)
		total_time=$((solver_time + assembly_time))
		solver_its=$(${AWK} '/RES/{print $6;}' $file)
	
		printf "Extract from ${file}\n"
	
		printf "$(( n * n * n))\t\
			 ${assembly_time}\t\
			 ${solver_time}\t\
			 ${total_time}\n" >> tvsn_linear_a${a}.txt
	
		${AWK} -v tot_time=${solver_time} -v tot_its=${solver_its} '
		/^cgpd/{
				printf("%-4d\t%e\t%e\n",
				$5, ($5 * tot_time * 1.e-6/ tot_its), $8);
	
		}' $file > rvsi_n${n}_a${a}.txt
		#awk -v tot_time=${solver_time} -v tot_its=${solver_its} '/^cgpd/{print $5 " " $8}' $file > rvsi_n${n}_a${a}.txt
	
	done
done
