#!/bin/bash

source vars.sh

mkdir -p result

PATH_OUT="./result"
AWK="/usr/bin/awk"

echo "#   N     assembly   solve   total" > tvsn_linear.txt
for n in ${sizes[@]}; do

	file="${PATH_OUT}/out_${n}.txt"
	echo "Extract from $file"
	assembly_time=$(${AWK} '/assembly_rhs/{print ($6)}' $file)
	assembly_time=$((assembly_time * 2))
	solver_time=$(${AWK} '/ell_solve_cgpd/{print $6}' $file)
	total_time=$((solver_time + assembly_time))
	total_its=$(${AWK} '/SOLVER_END/{print $4; exit}' $file)

	printf "$(( n * n * n))\t\
		 ${assembly_time}\t\
		 ${solver_time}\t\
		 ${total_time}\n" >> tvsn_linear.txt

	${AWK} -v tot_time=${total_time} -v tot_its=${total_its} '
	BEGIN{flag = 0;}
	{
		if($1 == "SOLVER_START") {
			flag = 1;
			getline;
		}
		if($1 == "SOLVER_END") {
			exit;
		}
		if(flag == 1) {
			printf("%-4d\t%e\t%e\t%e\n",
			$5, ($5 * tot_time * 1.e-6/ tot_its), $8, $11);
		}

	}' $file > rvsi_${n}.txt

done
