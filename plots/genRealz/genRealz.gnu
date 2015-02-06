set autoscale
set yrange[   0.28000:   1.12000]

#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Original Realizations"
set key left

#Label axes
set ylabel "Cross Section Value"
set xlabel "x Position"

#Enter Plotline

plot "genRealzplot1.txt" using 1:2 with lines, "genRealzplot2.txt" using 1:2 with lines, "genRealzplot3.txt" using 1:2 with lines

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './genRealzplot.ps'
#Screen output
replot
pause -1
