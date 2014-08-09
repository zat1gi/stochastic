#Main graph frame
set ytics border
set ytics mirror
set bars small
#set title "Transmission PDF" font "Times-Roman, 19"
set key left

#Graph axis
#set logscale
set autoscale

#Set axes


#Label axes
set ylabel "Probability" font "Times-Roman, 19"
set xlabel "Slab Reflection and Transmission Current" font "Times-Roman, 19"
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key top center
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "radMCtranreflprofile.ANS2014.txt"   u 3:4 w histeps t "TMC Reflectance",\
     "radWoodtranreflprofile.ANS2014.txt" u 3:4 w histeps t "WMC Reflectance",\
     "radMCtranreflprofile.ANS2014.txt"   u 1:2 w histeps t "TMC Transmittance",\
     "radWoodtranreflprofile.ANS2014.txt" u 1:2 w histeps t "WMC Transmittance"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './tranreflprofile.ps'
#Screen output
replot
pause -1
