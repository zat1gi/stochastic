
#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Flux in Material 2 (Cell) Over Stochastic Realizations"
set key left
set autoscale

#Label axes
set ylabel "Flux"
set xlabel "Position in Slab"

#Enter Plotline

plot "plots/radMC_fbcellflux.txt" u 1:4 w lines t "Original Realz, Traditional MC",\
       "plots/radMC_fbcellflux.txt" u 1:4:5 t "Or Rz, TMC Error" w yerrorbars

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './mat2cellflux.ps'
#Screen output
replot
pause -1
