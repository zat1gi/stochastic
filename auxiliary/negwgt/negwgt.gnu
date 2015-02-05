set title "Weight Adjustment Values as Function of Arbitrary XS"
set xlabel 'sigt*'
set ylabel 'weight adjustment'

sigs=-20.0
sigt=0.5

f(x) = abs(sigs)/x
 g(x) = abs((x-sigt)/x)
 g(x) = (x+abs(sigt))/x

set xrange[0:10]
set ytics 1

plot f(x) t 'sigs/sigt*', g(x) t 'abs((x-sigt)/x)', f(x)+g(x) t 'both, sigs=-1.5, sigt=0.5'
pause -1


