set termoption dash
set autoscale
set yrange[   0.29404:   0.63632]
#Main graph frame
set ytics border
set ytics mirror
set bars small
#set title "Original Realizations" font "Times-Roman, 19"
set key left

#Label axes
#set ylabel "Cross Section Value" font "Times-Roman, 19"
#set xlabel "Position in Slab" font "Times-Roman, 19"
set ytics font "Times-Roman, 30"
set xtics font "Times-Roman, 30"
set key top right
set key samplen 5 spacing 3.5 font "Times-Roman, 30"
unset colorbox
set style line 1 lt 1 lc rgb "red" lw 1
set style line 2 lt 2 lc rgb "green" lw 1
set style line 3 lt 3 lc rgb "blue" lw 2
#Enter Plotline

plot "txt/genRealzplot11.txt" using 1:2 with lines t "Realization 1" ls 1

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 19
set output './genRealz1.ps'
#Screen output
replot
pause 1
