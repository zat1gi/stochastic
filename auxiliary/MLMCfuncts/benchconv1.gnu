#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Convergence of Functionals at Refinement Levels for Benchmark"
set key right top

#Graph axis
set logscale y
#set autoscale

#Set axes
set xrange [-1:5]
set yrange [0.00000001:1]

#Label axes
set ylabel "Error of Functionals"
set xlabel "Level"

#Enter Plotline

plot "../../plots/MLMCfuncts/benchmarkbias.out" using 1:20:21 w errorbars t "funct"


#Postscript info
set size 1.0,1.0
set view equal xy
set size square
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './benchmarkbias1.ps'
#Screen output
replot
pause -1
