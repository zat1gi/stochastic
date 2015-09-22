mv gmon.out auxiliary/profiletools/gmon.out
gprof astochastic auxiliary/profiletools/gmon.out -p > auxiliary/profiletools/gmon.txt
head -30 auxiliary/profiletools/gmon.txt
