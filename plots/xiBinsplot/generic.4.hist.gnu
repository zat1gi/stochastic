#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Random Distributions"
set key bottom left

#Graph axis
#set logscale
set autoscale

#Set axes


#Label axes
set ylabel "Probability"
set xlabel "Random Value"

#Enter Plotline

plot "xiBinsAll.txt" using 1:2 with histeps t "number 1", \
     "xiBinsAll.txt" using 1:3 with histeps, \
     "xiBinsAll.txt" using 1:4 with histeps, \
     "xiBinsAll.txt" using 1:5 with histeps

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './customxiBinsplot.ps'
#Screen output
replot
pause -1
