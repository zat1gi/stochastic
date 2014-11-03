#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Eigenfunctions" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes

#Label axes
set ylabel "Eigenfunction Value" font "Times-Roman, 19"
set xlabel "Position in Slab" font "Times-Roman, 19"
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key bottom right
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "txt/Eigfunc.txt" using 1:2 with lines t "Eigenfunction 1", \
     "txt/Eigfunc.txt" using 1:3 with lines t "Eigenfunction 2", \
     "txt/Eigfunc.txt" using 1:4 with lines t "Eigenfunction 3", \
     "txt/Eigfunc.txt" using 1:5 with lines t "Eigenfunction 4"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './Eigfunc.ps'
#Screen output
replot
pause 1