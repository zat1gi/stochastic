
#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "<-reflection      Transport Results     transmission->"
set key left
set autoscale

#Label axes
set ylabel "Percent Reflected or Transmitted"
set xlabel "<-reflection     Result and Case Type      transmission->"

#Enter Plotline

plot "plots/transAdams.txt" using 1:2:3  title "Adams refl"     with yerrorbars, \
     "plots/transAdams.txt" using 1:7:8  title "Adams refl2"    with yerrorbars, \
     "plots/transAdams.txt" using 6:4:5  title "Adams trans"    with yerrorbars, \
     "plots/transAdams.txt" using 6:9:10 title "Adams trans2"   with yerrorbars, \
     "plots/radtrans.txt"   using 1:2:3  title "radtrans refl"  with yerrorbars, \
     "plots/radtrans.txt"   using 6:4:5  title "radtrans trans" with yerrorbars, \
     "plots/radWood.txt"    using 1:2:3  title "radWood refl"   with yerrorbars, \
     "plots/radWood.txt"    using 6:4:5  title "radWood trans"  with yerrorbars

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './transport.ps'
#Screen output
replot
pause 5