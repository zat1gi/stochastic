#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Random Parameter Distributions"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Probability"
set xlabel "Value"

#Enter Plotline

plot "txt/xiBins.txt" using 1:2 with histeps t "Eigenmode 1", \
     "txt/xiBins.txt" using 1:3 with histeps t "Eigenmode 2", \
     "txt/xiBins.txt" using 1:4 with histeps t "Eigenmode 3", \
     "txt/xiBins.txt" using 1:5 with histeps t "Eigenmode 4"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './xiBins.ps'
#Screen output
replot
pause 1