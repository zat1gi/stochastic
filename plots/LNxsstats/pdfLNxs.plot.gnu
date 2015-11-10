#Main graph frame
set ytics border
set ytics mirror
set bars small
set title "Cross Section Distributions"
set key top right

#Graph axis
#set logscale
set autoscale

#Set axes


#Label axes
set ylabel "Probability"
set xlabel "Cross Section Value"

#Enter Plotline

plot "plots/LNxsstats/pdfLNxs.txt" using 1:2 with histeps t "Depth 1", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:3 with histeps t "Depth 2", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:4 with histeps t "Depth 3", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:5 with histeps t "Depth 4", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:6 with histeps t "Depth 5", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:7 with histeps t "Depth 6", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:8 with histeps t "Depth 7", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:9 with histeps t "Depth 8", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:10 with histeps t "Depth 9", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:11 with histeps t "Depth 10", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:12 with histeps t "Depth 11", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:13 with histeps t "Depth 12", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:14 with histeps t "Depth 13", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:15 with histeps t "Depth 14", \
     "plots/LNxsstats/pdfLNxs.txt" using 1:16 with histeps t "Depth 15"

#Postscript info
set size 1.0,0.6
set terminal postscript portrait enhanced color dashed lw 2 "Times-Roman" 11
set output './pdfLNxs.ps'
#Screen output
replot
pause 1
