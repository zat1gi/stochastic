set termoption dash
set ytics border
set ytics mirror
set bars small
set title "KL Reconstructed Cross Sections"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Cross Section Value"
set xlabel "Position in Slab"
set style line 1 lt 1 lc rgb "red"     lw 1
set style line 2 lt 2 lc rgb "green"   lw 1
set style line 3 lt 3 lc rgb "blue"    lw 2
set style line 4 lt 1 lc rgb "magenta" lw 1
set style line 5 lt 1 lc rgb "cyan"    lw 1
set style line 6 lt 1 lc rgb "brown"   lw 1

#Enter Plotline
plot "txt/genRealzplot11.txt" using 1:2 with lines t "Realization 1" ls 1, \
     "txt/genRealzplot19.txt" using 1:2 with lines t "Realization 2" ls 2, \
     "txt/genRealzplot21.txt" using 1:2 with lines t "Realization 3" ls 3, \
     "txt/KLrrealzplot.txt" using 1:2 with lines t "KL Realization 1" ls 4, \
     "txt/KLrrealzplot.txt" using 1:3 with lines t "KL Realization 2" ls 5, \
     "txt/KLrrealzplot.txt" using 1:4 with lines t "KL Realization 3" ls 6

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './genRealz_KLrrealz.ps'
#Screen output
replot
pause 1