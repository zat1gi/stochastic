
#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Ensemble Material Distributions in Original Realziations"
set key left
set autoscale

#Label axes
set ylabel "Ensemble Mat"
set xlabel "Position in Slab"

#Enter Plotline

plot "plots/ensemble_mat.txt" u 1:2:3 t "Mat 1" w yerrorbars, \
     "plots/ensemble_mat.txt" u 1:4:5 t "Mat 2" w yerrorbars, \
     "plots/ensemble_mat.txt" u 1:6:7 t "Both Mats" w yerrorbars

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './ensemble_mat.ps'
#Screen output
replot
pause -1
