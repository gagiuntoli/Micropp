#!/bin/bash

if [ "$#" -ne 1 ]; then
        echo "Usage ./porcentage.sh <test>"
	exit 1
fi

TEST=$1

rm -rf t_ass.dat t_sol.dat sizes.dat

#sizes=(4 5 6)
sizes=(4 5 6 7 8 9 10 11 12 13 14)
#sizes=(4 5 6 7 8 9 10)
tstep=20

rm -rf t_ass.dat t_sol.dat sizes.dat

for i in ${sizes[@]}; do

  echo "running case n = " $i

  ./${TEST} $i $i $i 1 $tsteps > out.dat

  echo $i " " $(( i * i )) >> sizes.dat
  awk '/^ASSEMBLY/{print $7; exit;}' out.dat >> t_ass.dat
  awk '/^SOLVE/{print $7; exit;}' out.dat >> t_sol.dat

done

echo "#N  N^3   ASSEMBLY   SOLVE" > t_vs_n.dat
paste sizes.dat t_ass.dat t_sol.dat >> t_vs_n.dat

gnuplot -e " \
set terminal postscript eps color font 20; \
set ylabel 'Execution Time [mS]';\
set xlabel 'Problem Size \(Number of Nodes\)';\
titles = 'Assembly Solve';\
set key top left;\

set output 'plot_1.eps'; \
plot for[i=4:3:-1] 't_vs_n.dat' u 1:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \

set output 'plot_2.eps'; \
plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \

set log y;\
set output 'plot_3.eps'; \
plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \

unset log y;\
set log x;\
set output 'plot_4.eps'; \
plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \

set log y;\
set log x;\
set output 'plot_5.eps'; \
plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \
"
