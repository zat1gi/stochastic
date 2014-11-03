#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Reflection PDF" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Probability" font "Times-Roman, 19"
set xlabel "Current on Reflection Boundary" font "Times-Roman, 19"
set style line 6 lt 1 lc rgb "brown" lw 1
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key top right
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "txt/radMCtranreflprofile.txt"   u 1:2 w histeps t "TMC transmittance",\
     "txt/radWoodtranreflprofile.txt" u 1:2 w histeps t "WMC transmittance",\
     "txt/KLWoodtranreflprofile.txt"  u 1:2 w histeps t "KLWMC transmittance",\
     "txt/radMCtranreflprofile.txt"   u 3:4 w histeps t "TMC reflectance",\
     "txt/radWoodtranreflprofile.txt" u 3:4 w histeps t "WMC reflectance",\
     "txt/KLWoodtranreflprofile.txt"  u 3:4 w histeps t "KLWMC reflectance" ls 6

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './tranreflprofile3.ps'
#Screen output
replot
pause 1
