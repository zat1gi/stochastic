#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Eigenvalues"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes

#Label axes
set ylabel "Eigenvalue"
set xlabel "Eigenvalue Number"

#Enter Plotline

plot "txt/eigvals.txt" using 1:2 t "Eigenvalues" w linespoints pt 5

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './eigvals.ps'
#Screen output
replot
pause 1