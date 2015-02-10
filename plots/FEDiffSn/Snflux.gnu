#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux in Slab - Sn Solve"
set key left

#Label axes
set ylabel "Flux"
set xlabel "x Position"

#Enter Plotline

plot "Snflux.txt" using 1:2 with lines

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './Snflux.ps'
#Screen output
pause 2
replot