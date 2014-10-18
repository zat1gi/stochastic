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
set ytics font "Times-Roman, 19"
set xtics font "Times-Roman, 19"

set key top center
set key samplen 5 spacing 3.5 font "Times-Roman, 19"

#Enter Plotline

plot "plots/tranreflprofile/radMCtranreflprofile.txt"   u 3:4 w histeps t "TMC",\
     "plots/tranreflprofile/radWoodtranreflprofile.txt" u 3:4 w histeps t "WMC",\
     "plots/tranreflprofile/KLWoodtranreflprofile.txt"  u 3:4 w histeps t "KLWMC",\
     "plots/tranreflprofile/LPtranreflprofile.txt"      u 1:2 w histeps t "LPMC",\
     "plots/tranreflprofile/atmixtranreflprofile.txt"   u 1:2 w histeps t "atmixTMC"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './reflprofile.ps'
#Screen output
replot
pause 1
