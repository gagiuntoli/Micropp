#!/usr/bin/gnuplot

set terminal postscript eps color font 20

set key top left

set terminal postscript eps color font 20; \
set ylabel 'Execution Time [mS]';\
set xlabel 'Problem Size \(Number of Nodes\)';\
titles = 'Assembly Solve';\

set output 'tvsn_linear_a1.eps'; \
plot for[i=3:2:-1] 'tvsn_linear_a1.txt' u 1:(sum [col=2:i] column(col) * 1.0e-6) with filledcurves x1 title word(titles, i - 1)

set output 'tvsn_linear_a10.eps'; \
plot for[i=3:2:-1] 'tvsn_linear_a10.txt' u 1:(sum [col=2:i] column(col) * 1.0e-6) with filledcurves x1 title word(titles, i - 1)

set output 'tvsn_linear_a100.eps'; \
plot for[i=3:2:-1] 'tvsn_linear_a100.txt' u 1:(sum [col=2:i] column(col) * 1.0e-6) with filledcurves x1 title word(titles, i - 1)

