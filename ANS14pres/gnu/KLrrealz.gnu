set termoption dash
set ytics border
set ytics mirror
set bars small
set title "KL Reconstructed Cross Sections" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Cross Section Value" font "Times-Roman, 19"
set xlabel "Position in Slab" font "Times-Roman, 19"
set style line 1 lt 1 lc rgb "magenta" lw 1
set style line 2 lt 1 lc rgb "cyan" lw 1
set style line 3 lt 1 lc rgb "brown" lw 1
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key bottom right
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline
plot "txt/KLrrealzplot.txt" using 1:2 with lines t "Realization 1" ls 1, \
     "txt/KLrrealzplot.txt" using 1:3 with lines t "Realization 2" ls 2, \
     "txt/KLrrealzplot.txt" using 1:4 with lines t "Realization 3" ls 3

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './KLrrealz.ps'
#Screen output
replot
pause 1