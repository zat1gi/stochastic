set termoption dash
set ytics border
set ytics mirror
set bars small
set title "Eigenfunction/Realization Integration" font "Times-Roman, 21"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Eigenfunction Value/Cross Section Value" font "Times-Roman, 21"
set xlabel "Position in Slab" font "Times-Roman, 21"
set style line 1 lt 1 lc rgb "red"     lw 1
set style line 2 lt 2 lc rgb "green"   lw 1
set style line 3 lt 3 lc rgb "blue"    lw 2
set style line 4 lt 1 lc rgb "magenta" lw 1
set style line 5 lt 1 lc rgb "cyan"    lw 1
set style line 6 lt 1 lc rgb "brown"   lw 1
set ytics font "Times-Roman, 22"
set xtics font "Times-Roman, 22"

set key bottom right
set key samplen 5 spacing 3.5 font "Times-Roman, 23"


#Enter Plotline
plot "txt/genRealzplot11.txt" using 1:2 with lines t "Realization 1" ls 1, \
     "txt/Eigfunc.txt" using 1:5 with lines t "Eigenfunction 4" ls 4

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './genRealz_Eigfunc.ps'
#Screen output
replot
pause 1