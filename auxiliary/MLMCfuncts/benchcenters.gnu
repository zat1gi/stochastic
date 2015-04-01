#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Cell Center Flux Values - Benchmark"
set key left

#Graph axis
#set logscale
set autoscale

#Set axes
set xrange [0:10]

#Label axes
set ylabel "Center Valued Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot "../../plots/MLMCfuncts/benchcenters.out" using 1:2:3 pointtype 2 with errorbars
                                                           # -1 and 0 good pointtype options as well
#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './benchcenters.ps'
#Screen output
replot
pause -1
