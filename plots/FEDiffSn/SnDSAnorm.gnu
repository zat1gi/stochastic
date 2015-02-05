#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Diff Between Sn Flux and DSA Flux in Slab"
set key left

#Label axes
set ylabel "Difference"
set xlabel "x Position"

#Enter Plotline

plot "SnDSAnorm.txt" using 1:2 with lines

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './SnDSAnorm.ps'
#Screen output
pause 2
replot

