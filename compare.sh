#This bash script is used to compare output of two versions of stochastic
#for debugging purposes
#This needs to be run from the main stochastic directory, with the copy with
#which to compare in an adjacent directory named 'stochastic2'.
#That which you wish to compare should include the string 'flag'.

./astochastic | grep 'flag' > ../stochastic2/print1.out

cd ../stochastic2/stochastic
./astochastic | grep 'flag' > ../print2.out

cd ../
diff print1.out print2.out > diff.out

emacs print1.out &
emacs print2.out &
emacs diff.out &
