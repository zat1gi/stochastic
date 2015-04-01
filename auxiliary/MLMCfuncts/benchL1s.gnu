#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Region Center L1 Flux Values - Benchmark"
set key left

#Graph axis
#set logscale
set autoscale

#Set axes
set xrange [0:10]

#Label axes
set ylabel "Center Valued L1 Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot "../../plots/MLMCfuncts/benchL1s.out" using 3:4:5 pointtype 2 with errorbars

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './benchL1s.ps'
#Screen output
replot
pause -1
