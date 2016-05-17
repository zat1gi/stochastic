#Makefile
FC = mpif90

OBJ = mods/stochastic.F90
MODS = mods/mt_stream_f90-1.11/mt_kind_defs.F90 mods/mt_stream_f90-1.11/mt_stream.F90 mods/mt_stream_f90-1.11/f_jump_ahead_coeff/gf2xe.F90 mods/mt_stream_f90-1.11/f_jump_ahead_coeff/f_get_coeff.F90 mods/utilities.F90 mods/variables.F90 mods/timeman.f90 mods/genRealz.f90 mods/Loadcase.f90 mods/KLconstruct.f90 mods/KLresearch.f90 mods/radtransMC.F90

#Builds Targets
mods: $(MODS)
	$(FC) -c $(MODS)
	rm -f *.o

exec: $(OBJ)
	make mods
	$(FC) -fno-range-check -fbounds-check -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

mpiexec: $(OBJ)
	make mods
	$(FC) -fno-range-check -fbounds-check -DUSE_MPI -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

profile: $(OBJ)
	make mods
	$(FC) -pg -fno-range-check -fbounds-check -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/
	#to use profile tools: 'make profile'; './astochastic'; './auxiliary/profiletools/profile.sh'

debug: $(OBJ)
	make mods
	$(FC) -g -fno-range-check -fbounds-check -Wall -fbacktrace -finit-real=nan -o astochastic  $(OBJ) $(MODS)
	mv *.mod mods/

mpidebug: $(OBJ)
	make mods
	$(FC) -g -fno-range-check -DUSE_MPI -fbounds-check -Wall -fbacktrace -finit-real=nan -o astochastic  $(OBJ) $(MODS)
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
