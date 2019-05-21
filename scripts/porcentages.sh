#!/bin/bash

EXEC="../build/test/test3d_1"

sizes=(11 21 31 41 51 61 71 81 91 101)
tstep=1

echo "#N  N^3   ASSEMBLY   SOLVE" > t_vs_n.dat
for i in ${sizes[@]}; do

  echo "running case n = " $i

  $EXEC $i 0 $tsteps > out.dat

  size=$(( i * i * i ))
  tmat=$(awk '/assembly_mat/{print $4; exit;}' out.dat)
  trhs=$(awk '/assembly_rhs/{print $4; exit;}' out.dat)
  tsol=$(awk '/ell_solve_cgpd/{print $4; exit;}' out.dat)

  echo "$size $((tmat + trhs)) $tsol"
  echo "$size $((tmat + trhs)) $tsol" >> t_vs_n.dat

done


#gnuplot -e " \
#set terminal postscript eps color font 20; \
#set ylabel 'Execution Time [mS]';\
#set xlabel 'Problem Size \(Number of Nodes\)';\
#titles = 'Assembly Solve';\
#set key top left;\
#
#set output 'plot_1.eps'; \
#plot for[i=4:3:-1] 't_vs_n.dat' u 1:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \
#
#set output 'plot_2.eps'; \
#plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \
#
#set log y;\
#set output 'plot_3.eps'; \
#plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \
#
#unset log y;\
#set log x;\
#set output 'plot_4.eps'; \
#plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \
#
#set log y;\
#set log x;\
#set output 'plot_5.eps'; \
#plot for[i=4:3:-1] 't_vs_n.dat' u 2:((sum [col=3:i] column(col)) / ${tstep}) with filledcurves x1 title word(titles, i - 2); \
#"
