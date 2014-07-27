#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Generic Plot"
set key left

#Graph axis
#set logscale
set autoscale

#Set axes


#Label axes
set ylabel "y-values"
set xlabel "x-values"

#Enter Plotline

plot "absdiffGam.txt" using 2:3 with lines

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './genericplot.ps'
#Screen output
replot
pause -1
