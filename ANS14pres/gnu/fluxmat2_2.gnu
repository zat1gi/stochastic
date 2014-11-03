#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux Over Domain, Material 2, Point on Tracklength Binning"
set key left bottom
set autoscale
set logscale y
#Label axes
set ylabel "Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot \
  "txt/radMC_fluxmat.out"   u 1:4 w lines t "Markov Mix, TMC",\
  "txt/radWood_fluxmat.out" u 1:4 w lines t "Markov Mix, WMC"
#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './fluxmat2_2.ps'
#Screen output
replot
pause 1
