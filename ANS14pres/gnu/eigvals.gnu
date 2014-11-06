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
set termopt enhanced

set key top right
set key samplen 5 spacing 3.5 font "Times-Roman, 19"
set style line 2 lt 1 lc rgb "red" lw 1 pt 5 ps 0.5
set style line 2 lt 2 lc rgb "green" lw 1
set style line 3 lt 3 lc rgb "blue" lw 2
set style line 3 lt 3 lc rgb "blue" lw 2
set style line 3 lt 3 lc rgb "blue" lw 2

#Enter Plotline
set style function linespoints
plot "txt/eigvalsgaussian2.txt" using 1:2 t "lam_c/s = 0.05"   w linespoints lt 3 lc rgb "brown" pt 5,\
     "txt/eigvalsgaussian.txt" using 1:2 t "lam_c/s = 0.08"     w linespoints lt 2 lc rgb "brown" pt 5,\
     "txt/eigvals.txt" using 1:2 t "lam_c/s = 0.16" w linespoints lt 1 lc rgb "red"  pt 5,\
     "txt/eigvalsdelta.txt" using 1:2 t "lam_c/s = 0.48"        w linespoints lt 2 lc rgb "blue" pt 5,\
     "txt/eigvalsdelta2.txt" using 1:2 t "lam_c/s = 0.96"       w linespoints lt 3 lc rgb "blue" pt 5

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './eigvals.ps'
#Screen output
replot
pause 1