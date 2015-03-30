#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Iterative Convergence of Functionals"
set key left bottom

#Graph axis
set logscale y
#set autoscale

#Set axes
#set xrange [1:100000000]

#Label axes
set ylabel "Error of Functionals"
set xlabel "Number of Iterations"

#Enter Plotline

plot "../../plots/MLMCfuncts/iterconv.out" using 1:9 w lines t "functional 1", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:10 w lines t "functional 2", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:11 w lines t "functional 3", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:12 w lines t "functional 4", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:13 w lines t "functional 5", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:14 w lines t "functional 6", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:15 w lines t "functional 7", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:16 w lines t "functional 8", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:17 w lines t "functional 9", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:18 w lines t "functional 10", \
     "../../plots/MLMCfuncts/iterconv.out" using 1:19 w lines t "functional 11"


#Postscript info
set size 1.0,1.0
set view equal xy
set size square
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './iterconv.ps'
#Screen output
replot
pause -1
