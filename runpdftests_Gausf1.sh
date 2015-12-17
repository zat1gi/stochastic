./manipinput.sh 7 20000, 2000,
./manipinput.sh 11 LogN, Gaus,
./manipinput.sh 11 f2, f1,
./manipinput.sh 6 1100,11, 450,5,
./astochastic
mv plots/tranreflprofile/GaussKLtranreflprofile.txt plots/tranreflprofile/GaussKLtranreflprofile_Gausf1.txt
mv plots/tranreflprofile/reflprofile.pdf plots/tranreflprofile/reflprofile_Gausf1.pdf
mv plots/tranreflprofile/tranprofile.pdf plots/tranreflprofile/tranprofile_Gausf1.pdf

mv plots/fluxplots/GaussKL_fluxall.out plots/fluxplots/GaussKL_fluxall_Gausf1.out
mv plots/fluxplots/fluxall.pdf plots/fluxplots/fluxall_Gausf1.pdf
