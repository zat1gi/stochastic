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
  "plots/fluxplots/radMC_fluxmat.out" u 1:4 t "Markov Mix, TMC",\
  "plots/fluxplots/radWood_fluxmat.out" u 1:4 t "Markov Mix, WMC",\
  "plots/fluxplots/LPMC_fluxmat.out" u 1:3 t "Levermore-Pomraning MC",\
  "plots/fluxplots/atmixMC_fluxmat.out" u 1:2 t "(Mat Irr) Atomic Mix, TMC",\
#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './flux.ps'
#Screen output
replot

pause -1
