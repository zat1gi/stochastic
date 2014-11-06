#Main graph frame
set termoption dash
set ytics border
set ytics mirror
set bars small
set title "Flux Over Domain" font "Times-Roman, 21"
set key left bottom
set autoscale
set logscale y
#Label axes
set ylabel "Flux" font "Times-Roman, 21"
set xlabel "Position in Slab" font "Times-Roman, 21"
set ytics font "Times-Roman, 21"
set xtics font "Times-Roman, 21"
set termopt enhanced

set key top right
set key samplen 5 spacing 3.5 font "Times-Roman, 21"
#Enter Plotline

plot \
     'txt/radWood_fluxmat.out' u 1:2 t 'phi_1 WMC' w lines lt 1 lc rgb "brown" lw 2,\
     'txt/radWood_fluxmat.out' u 1:4 t 'phi_2 WMC' w lines lt 1 lc rgb "brown" lw 2,\
     'txt/radWood_fluxmat.out' u 1:($2*0.25+$4*0.75) t '<phi>_{} WMC' w lines lt 1 lc rgb "blue" lw 3,\
     'txt/KLWood_fluxmat.out'  u 1:2 t '<phi>_{KL} WMC' w lines lt 2 lc rgb "red" lw 3

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './fluxmatall_4.ps'
#Screen output
replot
pause 1
