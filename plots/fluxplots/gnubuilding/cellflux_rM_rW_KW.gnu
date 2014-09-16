
"plots/radMC_cellflux.txt" w lines t "Original Realz, Traditional MC",\
       "plots/radMC_cellflux.txt" u 1:2:3 t "Or Rz, TMC Error" w yerrorbars, \
     "plots/radWood_cellflux.txt" w lines t "Original Realz, Woodcock Sampling",\
       "plots/radWood_cellflux.txt" u 1:2:3 t "Or Rz, rWS Error" w yerrorbars, \
     "plots/KLWood_cellflux.txt" w lines t "Reconstructed Realz, Woodcock Sampling",\
       "plots/KLWood_cellflux.txt" u 1:2:3 t "Re Rz, KWS Error" w yerrorbars


