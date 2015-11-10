#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Cross Section Distributions"
set key top right

#Graph axis
#set logscale
set autoscale

#Set axes


#Label axes
set ylabel "Average or Variance"
set xlabel "Depth in Slab"

#Enter Plotline
set style line 1 lt 3 lc rgb "blue" pt 1 lw 2
set style line 2 lt 12 lc rgb "blue" pt 2
set style line 3 lt 3 lc rgb "red" pt 1 lw 2
set style line 4 lt 4 lc rgb "red" pt 2


plot "plots/LNxsstats/aveandvar.txt" u 5:1 w linespoints ls 1 t "Exp ave", \
     "plots/LNxsstats/aveandvar.txt" u 5:2 w linespoints ls 2 t "Obs ave", \
     "plots/LNxsstats/aveandvar.txt" u 5:3 w linespoints ls 3 t "Exp var", \
     "plots/LNxsstats/aveandvar.txt" u 5:4 w linespoints ls 4 t "Obs var"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './aveandvar.ps'
#Screen output
replot
pause -1
