#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Iterative Convergence of Functionals"
set key right bottom

#Graph axis
set logscale y
#set autoscale

#Set axes
#set xrange [1:100000000]

#Label axes
set ylabel "Error of Functionals"
set xlabel "Number of Iterations"

#Enter Plotline

plot "../../plots/MLMCfuncts/iterSnconv.out" using 1:9 w lines t "Sn functional 1", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:10 w lines t "Sn functional 2", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:11 w lines t "Sn functional 3", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:12 w lines t "Sn functional 4", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:13 w lines t "Sn functional 5", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:14 w lines t "Sn functional 6", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:15 w lines t "Sn functional 7", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:16 w lines t "Sn functional 8", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:17 w lines t "Sn functional 9", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:18 w lines t "Sn functional 10", \
     "../../plots/MLMCfuncts/iterSnconv.out" using 1:19 w lines t "Sn functional 11", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:9 w lines t "DSA functional 1", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:10 w lines t "DSA functional 2", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:11 w lines t "DSA functional 3", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:12 w lines t "DSA functional 4", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:13 w lines t "DSA functional 5", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:14 w lines t "DSA functional 6", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:15 w lines t "DSA functional 7", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:16 w lines t "DSA functional 8", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:17 w lines t "DSA functional 9", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:18 w lines t "DSA functional 10", \
     "../../plots/MLMCfuncts/iterDSAconv.out" using 1:19 w lines t "DSA functional 11"


#Postscript info
set size 1.0,1.0
set view equal xy
set size square
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './iterconv.ps'
#Screen output
replot
pause -1
