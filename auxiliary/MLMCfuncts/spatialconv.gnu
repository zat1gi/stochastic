#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Spatial Convergence of Functional"
set key left bottom

#Graph axis
set logscale
#set autoscale

#Set axes
set xrange [1:100000000]

#Label axes
set ylabel "Error of Functionals"
set xlabel "Number of Cells"

#Enter Plotline

plot "../../plots/MLMCfuncts/spatialconv.out" using 2:3 w lines t "functional 1", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:4 w lines t "functional 2", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:5 w lines t "functional 3", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:6 w lines t "functional 4", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:7 w lines t "functional 5", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:8 w lines t "functional 6", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:8 w lines t "functional 7", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:10 w lines t "functional 8", \
     "../../plots/MLMCfuncts/spatialconv.out" using 2:11 w lines t "functional 9"

#Postscript info
set size 1.0,1.0
set view equal xy
set size square
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './spatialconv.ps'
#Screen output
replot
pause -1
