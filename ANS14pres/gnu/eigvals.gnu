#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Eigenvalues" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes

#Label axes
set ylabel "Eigenvalue Magnitude" font "Times-Roman, 19"
set xlabel "Eigenvalue Number" font "Times-Roman, 19"
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key off
#set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "txt/eigvals.txt" using 1:2 t "Eigenvalues" w linespoints pt 5

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './eigvals.ps'
#Screen output
replot
pause 1