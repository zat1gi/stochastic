#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux Over Domain: Material 1" font "Times-Roman, 21"
set key left bottom
set autoscale
set logscale y
#Label axes
set ylabel "Flux" font "Times-Roman, 21"
set xlabel "Position in Slab" font "Times-Roman, 21"
set ytics font "Times-Roman, 21"
set xtics font "Times-Roman, 21"
set key top right
set key samplen 5 spacing 3.5 font "Times-Roman, 21"
#Enter Plotline

plot \
  "txt/radMC_fluxmat.out"   u 1:2 w lines t "Markov Mix, TMC",\
  "txt/radWood_fluxmat.out" u 1:2 w lines t "Markov Mix, WMC"
#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './fluxmat1_2.ps'
#Screen output
replot
pause 1
