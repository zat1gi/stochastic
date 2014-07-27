
#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux (Cell) Over Stochastic and Reconstructed Realizations"
set key left
set autoscale

#Label axes
set ylabel "Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot "plots/radWood_cellflux.txt" w lines t "Original Realz, Woodcock Sampling",\
       "plots/radWood_cellflux.txt" u 1:2:3 t "Or Rz, rWS Error" w yerrorbars

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './cellflux.ps'
#Screen output
replot
pause 1
