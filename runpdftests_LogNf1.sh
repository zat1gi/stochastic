./manipinput.sh 7 2000000, 20000000,
./manipinput.sh 11 Gaus, LogN,
./manipinput.sh 11 f2, f1,
./manipinput.sh 6 1100,11, 450,5,
./astochastic
mv plots/tranreflprofile/GaussKLtranreflprofile.txt plots/tranreflprofile/GaussKLtranreflprofile_LogNf1.txt
mv plots/tranreflprofile/reflprofile.pdf plots/tranreflprofile/reflprofile_LogNf1.pdf
mv plots/tranreflprofile/tranprofile.pdf plots/tranreflprofile/tranprofile_LogNf1.pdf

mv plots/fluxplots/GaussKL_fluxall.out plots/fluxplots/GaussKL_fluxall_LogNf1.out
mv plots/fluxplots/fluxall.pdf plots/fluxplots/fluxall_LogNf1.pdf
