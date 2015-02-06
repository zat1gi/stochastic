set dgrid3d         40,        40
set view 60,130                         
                                        
set title "Positional Correlation Plots"
set ylabel "\"y\" or \"x2\" position"   
set xlabel "\"x\" or \"x1\" position"   
                                        
splot "Correlation.txt" using 1:2:3 with lines title "Correlation Expected", "Correlation.txt" using 1:2:4 with lines title "Correlation Yielded "

set size 1.0,0.6                                        
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './Correlation.ps'
replot
pause -1
