#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux in Slab, Method of Manufactured Solutions"
set key left

#Label axes
set ylabel "Flux"
set xlabel "x Position"

#Enter Plotline

plot "fluxmoms.txt" using 1:2 with lines, "fluxmoms.txt" using 1:3 with circles

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './fluxmoms.ps'
#Screen output
pause 4
replot
