#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Spatial Co Maintained" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Percent Co Maintained" font "Times-Roman, 19"
set xlabel "Position in Slab" font "Times-Roman, 19"
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key bottom center
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline
plot "txt/CoEff.txt" using 1:2 with lines t "Eigenmode 1", \
     "txt/CoEff.txt" using 1:3 with lines t "Eigenmode 2", \
     "txt/CoEff.txt" using 1:4 with lines t "Eigenmode 3", \
     "txt/CoEff.txt" using 1:5 with lines t "Eigenmode 4"
#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './CoEff.ps'
#Screen output
replot
pause 1