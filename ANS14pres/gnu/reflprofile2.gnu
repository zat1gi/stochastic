#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Reflection PDF" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Probability" font "Times-Roman, 19"
set xlabel "Current on Reflection Boundary" font "Times-Roman, 19"
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key top right
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "txt/radMCtranreflprofile.txt"   u 3:4 w histeps t "TMC",\
     "txt/radWoodtranreflprofile.txt" u 3:4 w histeps t "WMC"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './reflprofile2.ps'
#Screen output
replot
pause 1
