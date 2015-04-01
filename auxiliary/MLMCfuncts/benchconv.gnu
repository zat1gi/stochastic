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

#Label axes
set ylabel "Error of Functionals"
set xlabel "Level"

#Enter Plotline

plot "../../plots/MLMCfuncts/benchmarkbias.out" using 1:2:3 w errorbars t "funct 1", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:4:5 w errorbars t "funct 2", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:6:7 w errorbars t "funct 3", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:8:9 w errorbars t "funct 4", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:10:11 w errorbars t "funct 5", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:12:13 w errorbars t "funct 6", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:14:15 w errorbars t "funct 7", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:16:17 w errorbars t "funct 8", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:18:19 w errorbars t "funct 9", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:20:21 w errorbars t "funct 10", \
     "../../plots/MLMCfuncts/benchmarkbias.out" using 1:22:23 w errorbars t "funct 11"


#Postscript info
set size 1.0,1.0
set view equal xy
set size square
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './benchmarkbias.ps'
#Screen output
replot
pause -1
