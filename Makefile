#Makefile
FC = mpif90

OBJ = mods/stochastic.F90
MODS = mods/mcnp_random.f90 mods/utilities.F90 mods/variables.F90 mods/timeman.f90 mods/genRealz.f90 mods/Loadcase.f90 mods/KLconstruct.f90 mods/KLresearch.f90 mods/radtransMC.F90

#Builds Targets
mods: $(MODS)
	$(FC) -c $(MODS)
	rm -f *.o

exec: $(OBJ)
	make mods
	$(FC) -fbounds-check -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

mpiexec: $(OBJ)
	make mods
	$(FC) -fbounds-check -DUSE_MPI -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

profile: $(OBJ)
	make mods
	$(FC) -pg -fbounds-check -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/
	#to use profile tools: 'make profile'; './astochastic'; './auxiliary/profiletools/profile.sh'

debug: $(OBJ)
	make mods
	$(FC) -g -fbounds-check -Wall -fbacktrace -finit-real=nan -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

mpidebug: $(OBJ)
	make mods
	$(FC) -g -DUSE_MPI -fbounds-check -Wall -fbacktrace -finit-real=nan -o astochastic  $(OBJ) $(MODS)
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
