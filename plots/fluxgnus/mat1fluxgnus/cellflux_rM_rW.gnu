
#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux in Material 1 (Cell) Over Stochastic Realizations"
set key left
set autoscale

#Label axes
set ylabel "Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot \
     "plots/radMC_fbcellflux.txt" w lines t "Original Realz, Traditional MC",\
       "plots/radMC_fbcellflux.txt" u 1:2:3 t "Or Rz, TMC Error" w yerrorbars, \
     "plots/radWood_fbcellflux.txt" w lines t "Original Realz, Woodcock Sampling",\
       "plots/radWood_fbcellflux.txt" u 1:2:3 t "Or Rz, rWS Error" w yerrorbars

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './mat1cellflux.ps'
#Screen output
replot
pause -1
