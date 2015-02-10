#Makefile
FC = gfortran

OBJ = mods/stochastic.f90
MODS = mods/mcnp_random.f90 mods/variables.f90 mods/utilities.f90 mods/timeman.f90 mods/genRealz.f90 mods/Loadcase.f90 mods/KLreconstruct.f90 mods/KLresearch.f90 mods/radtransMC.f90 mods/MLMC.f90 mods/FEDiffSn.f90

#Builds Targets
mods: $(MODS)
	$(FC) -c $(MODS)
	rm -f *.o

exec: $(OBJ)
	make mods
	$(FC) -fbounds-check -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

debug: $(OBJ)
	make mods
	$(FC) -g -fbounds-check -Wall -fbacktrace -finit-real=nan -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/


clean:
	rm -f test mods/*.mod mods/*~
	rm -f *mod *mod~
	rm -f *.o

run: $(OBJ)
	make mods
	$(FC) -fbounds-check -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/
	./astochastic
