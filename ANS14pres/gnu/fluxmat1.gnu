#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux Over Domain, Material 1, Point on Tracklength Binning"
set key left bottom
set autoscale
set logscale y
#Label axes
set ylabel "Flux"
set xlabel "Position in Slab"
#Enter Plotline

plot \
  "txt/radMC_fluxmat.out"   u 1:2 t "Markov Mix, TMC",\
  "txt/radWood_fluxmat.out" u 1:2 t "Markov Mix, WMC",\
  "txt/KLWood_fluxmat.out"  u 1:2 t "(Mat Irr) KL Reconstructions, WMC",\
  "txt/LPMC_fluxmat.out"    u 1:2 t "Levermore-Pomraning MC",\
  "txt/atmixMC_fluxmat.out" u 1:2 t "(Mat Irr) Atomic Mix, TMC",\
#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './fluxmat1_5.ps'
#Screen output
replot
pause 5
