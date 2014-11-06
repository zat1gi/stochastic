#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Transmittance and Reflectance PDFs" font "Times-Roman, 19"
set key left
#Graph axis
#set logscale
set autoscale
#Set axes
#Label axes
set ylabel "Probability" font "Times-Roman, 19"
set xlabel "Current on Boundary" font "Times-Roman, 19"
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"
set xrange[0.0:0.6]
set key top center#at 0,27 left
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "txt/radMCtranreflprofile.txt"   u 1:2 w histeps t "TMC Transmittance",\
     "txt/radWoodtranreflprofile.txt" u 1:2 w histeps t "WMC Transmittance",\
     "txt/radMCtranreflprofile.txt"   u 3:4 w histeps t "TMC Reflectance",\
     "txt/radWoodtranreflprofile.txt" u 3:4 w histeps t "WMC Reflectance"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './tranreflprofile2.ps'
#Screen output
replot
pause 1
