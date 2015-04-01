#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Region Center L2 Norm of Flux Values"
set key bottom right

#Graph axis
#set logscale
set autoscale

#Set axes
set xrange [0:10]

#Label axes
set ylabel "Center Valued L2 Norm of Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot "../../plots/MLMCfuncts/MLMCL2s.out" using 3:4:5 pointtype 2 with errorbars t "MLMC", \
     "../../plots/MLMCfuncts/benchL2s.out" using 3:4:5 pointtype 4 with errorbars t "Benchmark"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './bothL2s.ps'
#Screen output
replot
pause -1
