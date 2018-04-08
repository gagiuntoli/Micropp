#!/bin/bash

echo "# t_step[1] nr_its[2] nr_tol[3] nr_non_linear[4] LinCriteria[5] epsilon[6] sig11[7] sig22 ... sig13[12]" > aux1.dat

awk '
/Time step/ {printf("%-3d ", $4)}
/NEWTON-R/ {printf("%-2d %e %-2d ", $4, $7, $10)} 
/LinCriteria/ {printf("%-2d ", $3)} 
/e11/ {printf("%e ", $3)}
/Average stress/ {printf("%e %e %e %e %e %e\n", $4, $5, $6, $7, $8, $9)} 
' $1 > aux2.dat

cat aux1.dat aux2.dat > data.dat
rm aux1.dat aux2.dat

#Time step = 58
#NewRap It = 4 Tol = 7.99515
#Average stress = 1926.97 803.158 803.158 5.10311e-12 -7.13371e-12 2.52745e-12 
