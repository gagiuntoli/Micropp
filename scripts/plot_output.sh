#!/bin/bash
awk '
/Time step/ {printf("%d ", $4)}
/NEWTON-R/ {printf("%d %e ", $4, $7)} 
/e11/ {printf("%e ", $3)}
/Average stress/ {printf("%e %e %e %e %e %e\n", $4, $5, $6, $7, $8, $9);getline;} 
' $1 > data.dat



#Time step = 58
#NewRap It = 4 Tol = 7.99515
#Average stress = 1926.97 803.158 803.158 5.10311e-12 -7.13371e-12 2.52745e-12 

