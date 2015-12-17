./manipinput.sh 7 20000000, 2000000,
./manipinput.sh 11 LogN, Gaus,
./manipinput.sh 11 f1, f2,
./manipinput.sh 6 450,5, 1100,11,
./astochastic
mv plots/tranreflprofile/GaussKLtranreflprofile.txt plots/tranreflprofile/GaussKLtranreflprofile_Gausf2.txt
mv plots/tranreflprofile/reflprofile.pdf plots/tranreflprofile/reflprofile_Gausf2.pdf
mv plots/tranreflprofile/tranprofile.pdf plots/tranreflprofile/tranprofile_Gausf2.pdf

mv plots/fluxplots/GaussKL_fluxall.out plots/fluxplots/GaussKL_fluxall_Gausf2.out
mv plots/fluxplots/fluxall.pdf plots/fluxplots/fluxall_Gausf2.pdf
