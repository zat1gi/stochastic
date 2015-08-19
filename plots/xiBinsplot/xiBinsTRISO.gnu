#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Random Distributions for TRISO Fuel in Pebble"
set key top left

#Graph axis
#set logscale
set autoscale

#Set axes


#Label axes
set ylabel "Probability"
set xlabel "Random Value"

#Enter Plotline

plot "plots/xiBinsplot/xiBinsTRISOplot1-3.txt" using 1:2 with histeps t "Eig 1", \
     "plots/xiBinsplot/xiBinsTRISOplot1-3.txt" using 1:3 with histeps t "Eig 2", \
     "plots/xiBinsplot/xiBinsTRISOplot1-3.txt" using 1:4 with histeps t "Eig 3", \
     "plots/xiBinsplot/xiBinsTRISOplot4-6.txt" using 1:2 with histeps t "Eig 4", \
     "plots/xiBinsplot/xiBinsTRISOplot4-6.txt" using 1:3 with histeps t "Eig 5", \
     "plots/xiBinsplot/xiBinsTRISOplot4-6.txt" using 1:4 with histeps t "Eig 6", \
     "plots/xiBinsplot/xiBinsTRISOplot4-6.txt" using 1:5 with histeps t "unit Gauss"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './customTRISOxiBinsplot.ps'
#Screen output
replot
pause -1
