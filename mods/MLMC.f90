module MLMC
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine UQ_MLMC( icase )
  !This subroutine is the main driver for the MLMC method.
  !It's first implementation is planned to be using FEDiffSn as the deterministic solver
  !to demonstrate MLMC.  It's second implementation is planned to be using MC transport
  !as the solver.  Potentially the third is with MC transport over random binary media.
  use MLMCvars, only: MLMC_TOL, numMLMCcells
  use FEDiffSn, only: setflvarspassedtrue
  integer :: icase

  integer :: Level
  real(8) :: MLMCerrest
  logical :: flMLMC
  logical :: flread = .false.

  call MLMCallocate                                    !allocate variables
  call MLMCinitialize( flMLMC, Level )                 !initialize variables
  call setflvarspassedtrue                             !tell FE mod to accept input from here

  !main MLMC loop
  do while(flMLMC)
print *,"Level:",Level

    !1 Determine number of cells, make array which holds number of cells as a function of level L
    !1 aka add new Level and necessary cells
    if(Level /= 0) print *,"--reallocate to accomodate new Level--"
if(Level /= 0 .and. flread) read *
    if(Level /= 0) call MLMCaddLevel( Level )
    if(Level /= 0) print *,"--reallocate to accomodate new Level--"
if(Level /= 0 .and. flread) read *
if(Level /= 0) print *

    !2 Solve initial samples for Level, solve V~(estimated variance) for each level L
    print *,"--evaluate baseline samples, get variance estimate--"
if(flread) read *
    call MLMCevalNewSamps( Level )
    print *,"--evaluate baseline samples, get variance estimate--"
if(flread) read *
print *

    !3 Using V~s(estimated variance), compute optimal M~s(estimated opt num of samples)
    !3.2 expand arrays to hold new M~s(number of samples)
    print *,"--compute optimal number of samples--"
if(flread) read *
    call MLMCcomputeOptSamps( Level )
    print *,"--compute optimal number of samples--"
if(flread) read *
print *
    print *,"--reallocate to accomodate new samples--"
if(flread) read *
    call MLMCaddSamples( Level )
    print *,"--reallocate to accomodate new samples--"
if(flread) read *
print *

    !4 Evaluate any new samples needed, Gave and Gvar calculated here
    print *,"--evaluate extra samples needed--"
if(flread) read *
    call MLMCevalNewSamps( Level )
    print *,"--evaluate extra samples needed--"
if(flread) read *
print *

    !5 test total error, set flMLMC==.false.?
    if(Level>1) print *,"--calculate estimated error--"
if(Level>1 .and. flread) read *
    if(Level>1) MLMCerrest = MLMCcalcErrEst( Level )
    if(Level>1) print *,"MLMCerrest:",MLMCerrest
    if(Level>1) print *,"MLMC_TOL  :",MLMC_TOL
    if(Level>1) print *,"--calculate estimated error--"
if(Level>1 .and. flread) read *
    !MLMCerrest = 0.1d0
    if(Level>1 .and. MLMCerrest<=MLMC_TOL) then
      flMLMC = .false.
    else
      Level = Level + 1
!if(Level==4) flMLMC = .false.
    endif
  enddo !loops over realizations

  call MLMCdeallocate

  end subroutine UQ_MLMC




  subroutine MLMCallocate
  !Allocate arrays here.  Try to only allocate in order to make paralellization easier later.
  use MLMCvars, only: numMLMCcells, M_optsamps, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, bnumMLMCsamps, numcellsLevel0, ncellwidth

  if(.not.allocated(numMLMCcells)) allocate(numMLMCcells(0:0))
  if(.not.allocated(M_optsamps)) allocate(M_optsamps(3,0:0))

  if(.not.allocated(Q_ufunctional)) allocate(Q_ufunctional(bnumMLMCsamps,0:0))
  if(.not.allocated(G_ufunctional)) allocate(G_ufunctional(bnumMLMCsamps,0:0))
  if(.not.allocated(Gave)) allocate(Gave(0:0))
  if(.not.allocated(Gvar)) allocate(Gvar(0:0))
  if(.not.allocated(ncellwidth)) allocate(ncellwidth(0:0))

  end subroutine MLMCallocate


  subroutine MLMCinitialize( flMLMC, Level )
  !Initialize/load values to variables here.
  use MLMCvars, only: numMLMCcells, numcellsLevel0, MLMC_failprob, C_alpha, &
                      M_optsamps, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, bnumMLMCsamps, ncellwidth
  use utilities, only: erfi
  integer :: Level
  logical :: flMLMC

  Level  = 0
  flMLMC = .true.  

  C_alpha = sqrt(2.0d0)*erfi(1.0d0-MLMC_failprob)

  numMLMCcells    = numcellsLevel0
  M_optsamps(1,0) = bnumMLMCsamps
  M_optsamps(2,0) = 0
  M_optsamps(3,0) = 0

  Q_ufunctional = 0.0d0
  G_ufunctional = 0.0d0
  Gave          = 0.0d0
  Gvar          = 0.0d0

  ncellwidth    = real(1,8)/real(numcellsLevel0,8)

  end subroutine MLMCinitialize


  subroutine MLMCdeallocate
  !Deallocate arrays here.
  use MLMCvars, only: numMLMCcells, M_optsamps, Q_ufunctional, G_ufunctional, &
                      Gave, Gvar, ncellwidth

  if(allocated(numMLMCcells)) deallocate(numMLMCcells)
  if(allocated(M_optsamps)) deallocate(M_optsamps)

  if(allocated(Q_ufunctional)) deallocate(Q_ufunctional)
  if(allocated(G_ufunctional)) deallocate(G_ufunctional)
  if(allocated(Gave)) deallocate(Gave)
  if(allocated(Gvar)) deallocate(Gvar)
  if(allocated(ncellwidth)) deallocate(ncellwidth)

  end subroutine MLMCdeallocate





  subroutine MLMCaddLevel( Level )
  !Increase size of arrays to handle new Level and in case of response function
  !the associated number of cells as well.
  use MLMCvars, only: numMLMCcells, nextLevelFactor, Q_ufunctional, &
                      G_ufunctional, Gave, Gvar, M_optsamps, bnumMLMCsamps, &
                      ncellwidth
  integer :: Level

  integer, allocatable :: tiarray1(:)
  integer, allocatable :: tiarray2(:,:)

  real(8), allocatable :: trarray1(:)
  real(8), allocatable :: trarray2(:,:)

  !add new Level to numMLMCcells and populate it
  call move_alloc(numMLMCcells,tiarray1)
  allocate(numMLMCcells(0:size(tiarray1)))
  numMLMCcells = 0
  numMLMCcells(0:size(tiarray1)-1) = tiarray1
  numMLMCcells(Level) = numMLMCcells(Level-1) * nextLevelFactor
  deallocate(tiarray1)

  !add new Level to M_optsamps and populate it
  call move_alloc(M_optsamps,tiarray2)
  allocate(M_optsamps(3,0:size(tiarray2(1,:))))
  M_optsamps = 0
  M_optsamps(:,0:size(tiarray2(1,:))-1) = tiarray2
  M_optsamps(1,Level) = bnumMLMCsamps
  deallocate(tiarray2)
print *,"just reallocated, M_optsamps(1,:):",M_optsamps(1,:)
print *,"just reallocated, M_optsamps(2,:):",M_optsamps(2,:)

  !add new Level to Q_ufunctional and initialize it
  call move_alloc(Q_ufunctional,trarray2)
  allocate(Q_ufunctional(size(trarray2(:,1)),0:size(trarray2(1,:))))
  Q_ufunctional = 0.0d0
  Q_ufunctional(:,0:size(trarray2(1,:))-1) = trarray2
  deallocate(trarray2)

  !add new Level to G_ufunctional and initialize it
  call move_alloc(G_ufunctional,trarray2)
  allocate(G_ufunctional(size(trarray2(:,1)),0:size(trarray2(1,:))))
  G_ufunctional = 0.0d0
  G_ufunctional(:,0:size(trarray2(1,:))-1) = trarray2
  deallocate(trarray2)

  !add new Level to Gave and initialize it
  call move_alloc(Gave,trarray1)
  allocate(Gave(0:size(trarray1(:))))
  Gave = 0.0d0
  Gave(0:size(trarray1(:))-1) = trarray1
  deallocate(trarray1)

  !add new Level to Gvar and initialize it
  call move_alloc(Gvar,trarray1)
  allocate(Gvar(0:size(trarray1(:))))
  Gvar = 0.0d0
  Gvar(0:size(trarray1(:))-1) = trarray1
  deallocate(trarray1)

  !add new Level to ncellwidth and initialize it
  call move_alloc(ncellwidth,trarray1)
  allocate(ncellwidth(0:size(trarray1(:))))
  ncellwidth = 0.0d0
  ncellwidth(0:size(trarray1(:))-1) = trarray1
  ncellwidth(Level) = ncellwidth(Level-1) / real(nextLevelFactor,8)
  deallocate(trarray1)

  end subroutine MLMCaddLevel


  subroutine MLMCaddSamples( Level )
  !Increase size of arrays to handle new optimal number of samples.
  use MLMCvars, only: Q_ufunctional, G_ufunctional, M_optsamps
  integer :: Level

  real(8), allocatable :: trarray2(:,:)

  !add space for new samples to Q_functional and initialize the space
  call move_alloc(Q_ufunctional,trarray2)
  allocate(Q_ufunctional(maxval(M_optsamps),0:size(trarray2(1,:))-1))
  Q_ufunctional = 0.0d0
  Q_ufunctional(1:size(trarray2(:,1)),:) = trarray2
  deallocate(trarray2)

  !add space for new samples to G_ufunctional and initialize the space
  call move_alloc(G_ufunctional,trarray2)
  allocate(G_ufunctional(maxval(M_optsamps),0:size(trarray2(1,:))-1))
  G_ufunctional = 0.0d0
  G_ufunctional(1:size(trarray2(:,1)),:) = trarray2
  deallocate(trarray2)

  end subroutine MLMCaddSamples



  !3 Using V~s(estimated variance), compute optimal M~s(estimated opt num of samples)
  subroutine MLMCcomputeOptSamps( Level )
  use MLMCvars, only: MLMC_TOLsplit, MLMC_TOL, C_alpha, Gave, Gvar, M_optsamps, &
                      ncellwidth,linsolveEff,numDimensions
  integer :: Level

  integer :: ilevel
  real(8) :: workterm
  real(8) :: accterm !last term in estimate, sum over levels, sum as go

  workterm = 0.0d0
  accterm  = 0.0d0

print *,"M_optsamps(1,:):",M_optsamps(1,:)
  do ilevel=0,Level
    !solve for work term
    workterm = ncellwidth(ilevel)**(-linsolveEff*numDimensions)
    !accumlate next part of accumulating term
    accterm = accterm + sqrt(Gvar(ilevel)*workterm)
print *,"ilevel:",ilevel,"  accterm:",accterm
print *,"term1:",(MLMC_TOLsplit*MLMC_TOL/C_alpha)**(-2.0d0)
print *,"term2:",sqrt(abs(Gvar(ilevel)/workterm))
    !calculate new optimal samples estimate
    M_optsamps(3,ilevel) = ceiling(  (MLMC_TOLsplit*MLMC_TOL/C_alpha)**(-2.0d0) * &
                                     sqrt(abs(Gvar(ilevel)/workterm)) * &
                                     accterm                             )                         

print *,"M_optsamps(3,ilevel):",M_optsamps(3,ilevel)
    !use larger of new value and old value
    M_optsamps(1,ilevel) = max(M_optsamps(2,ilevel),M_optsamps(3,ilevel))
    !if new level's estimate higher, make backwards compatible
    if(ilevel/=0) then
      if(M_optsamps(1,ilevel) > M_optsamps(1,ilevel-1)) &
         M_optsamps(1,0:ilevel-1) = M_optsamps(1,ilevel)
    endif

print *,"theoretical optimal M for level",ilevel,":",M_optsamps(3,ilevel)
print *,"practical   optimal M for level",ilevel,":",M_optsamps(1,ilevel)
  enddo

  end subroutine MLMCcomputeOptSamps




  !4 Evaluate any new samples needed
  subroutine MLMCevalNewSamps( Level )
  !Evaluate any new samples needed, recompute functional values
  use MLMCvars, only: M_optsamps, Gave, Gvar, G_ufunctional
  use rngvars, only: setrngappnum, rngappnum, rngstride
  use utilities, only: mean_and_var_p
  use mcnp_random, only: RN_init_particle
  integer :: Level

  integer :: ilevel,isamp, isamplow

print *,"new samps M_optsamps(1,:):",M_optsamps(1,:)
print *,"new samps M_optsamps(2,:):",M_optsamps(2,:)
print *,"new samps M_optsamps(3,:):",M_optsamps(3,:)
  do ilevel = 0,Level                                          !search each level
print *,"ilevel:",ilevel
    if( M_optsamps(1,ilevel)>M_optsamps(2,ilevel) ) then       !if new samps to compute
      isamplow = max(M_optsamps(2,ilevel)+1,1)                 !set lowest samp num
print *,"isamplow:",isamplow
      do isamp = isamplow,M_optsamps(1,ilevel)                 !cycle through samps to compute
        call setrngappnum('MLMCsamp')                          !set rng unique to sample
        call RN_init_particle( int(rngappnum*rngstride+isamp,8) )
        call sampleInput( ilevel )                             !update solver input info
        call solveSamples( ilevel,isamp )                      !solves QoI with and applies norm
if(mod(isamp,1000)==0) print *,"ilevel:",ilevel,"  isamp:",isamp
if(mod(isamp,1000)==0) print *,"G_ufunctional(",isamp,",",ilevel,"):",G_ufunctional(isamp,ilevel)
      enddo
    endif
    call mean_and_var_p( G_ufunctional(:,ilevel),&             !solve ave and var of functionals
                         size(G_ufunctional(:,ilevel)),Gave(ilevel),Gvar(ilevel) )
print *,"just formed Gave(",ilevel,"):",Gave(ilevel)
print *,"just formed Gvar(",ilevel,"):",Gvar(ilevel)
    M_optsamps(2,ilevel) = M_optsamps(1,ilevel)                !save old # of opt samps
  enddo
  end subroutine MLMCevalNewSamps



  subroutine sampleInput( ilevel )
  !This subroutine samples input parameters for this solve and passes them as needed
  use genSampvars, only: specialprob, nummat, param1, param2, param1_mean, param1_uncert, &
                         param2_mean, param2_uncert
  use MLMCvars, only: numMLMCcells
  use FEDiffSn, only: sigt, c, numcells,      flvarspassed
  use mcnp_random, only: rang

  integer :: ilevel
  real(8) :: loc_sigs, loc_siga
  logical :: flsiga            !have I already sampled siga?

  !The following preparation of sigt and c is designed for use with FEDiffSn, 
  !and is written on the assumption that there is only 1 material present.
  !When similar parsing is written later, code with that in mind.
  flsiga = .false.

  !sample and prepare sigt
  if( param1(1)=='sigt' ) then
    if( param1(2)=='sigt1-abs' ) then
      sigt = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_uncert(1)
    elseif( param1(2)=='sigt1-frac' ) then
      sigt = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_mean(1)*param1_uncert(1)
    endif
  elseif( param1(1)=='sigs' ) then
    if( param1(2)=='sigs1-abs' ) then
      loc_sigs = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_uncert(1)
    elseif( param1(2)=='sigs1-frac' ) then
      loc_sigs = param1_mean(1) + (rang()*2.0d0-1.0d0) * param1_mean(1)*param1_uncert(1)
    endif
    if( param2(2)=='siga1-abs' ) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='siga1-frac' ) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
    flsiga = .true.
    sigt = loc_sigs + loc_siga
  endif

  !sample and prepare c
  if( param2(1)=='c' ) then
    if( param2(2)=='c1-abs' ) then
      c = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='c1-frac' ) then
      c = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
  elseif( param2(1)=='siga' ) then
    if( param2(2)=='siga1-abs' .and. .not.flsiga) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='siga1-frac' .and. .not.flsiga) then
      loc_siga = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
    c = merge(1.0d0-loc_siga/sigt,loc_sigs/(loc_sigs+loc_siga),.not.flsiga)
  elseif( param2(1)=='sigs' ) then
    if( param2(2)=='sigs1-abs' ) then
      loc_sigs = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_uncert(1)
    elseif( param2(2)=='sigs1-frac' ) then
      loc_sigs = param2_mean(1) + (rang()*2.0d0-1.0d0) * param2_mean(1)*param2_uncert(1)
    endif
    c = loc_sigs/sigt
  endif

  if( c>1.0d0 ) then
    stop "--'c' sampled as greater than 1, check your input!"
  endif

  numcells = numMLMCcells(ilevel)

  end subroutine sampleInput



  subroutine solveSamples( ilevel,isamp )
  !This subroutine drives the solver for a new response function flux for a new sample
  !collects solution, and deallocated module variables from FEDiffSn.
  use MLMCvars, only: Q_ufunctional, G_ufunctional, nextLevelFactor, numcellsLevel0, &
                      ncellwidth
  use FEDiffSn, only: FEmain,FEDiffSn_externaldeallocate, &
                      solve,phidiff,phiSnl,phiSnr,phiDSAl,phiDSAr

  integer :: ilevel, isamp, i, icell, first, last
  real(8), allocatable :: flux(:),cellwidth

  call FEmain

  !collect flux from FEmain results
  if(allocated(flux)) deallocate(flux)
  if(solve(3)==1) then
    allocate(flux(size(phiDSAl)))
    flux = 0.0d0
    flux = (phiDSAl + phiDSAr) / 2.0d0
  elseif(solve(2)==1) then
    allocate(flux(size(phiSnl)))
    flux = 0.0d0
    flux = (phiSnl + phiSnr) / 2.0d0
  elseif(solve(1)==1) then
    allocate(flux(size(phidiff)))
    flux = 0.0d0
    flux = phidiff
  endif

  !L2 norm over all domain for Q_ufunctional
  Q_ufunctional(isamp,ilevel) = 0.0d0
  do i=1,numcellsLevel0*nextLevelFactor**ilevel
    Q_ufunctional(isamp,ilevel) = Q_ufunctional(isamp,ilevel) + flux(i)**2*ncellwidth(ilevel)
  enddo
  Q_ufunctional(isamp,ilevel) = sqrt(Q_ufunctional(isamp,ilevel))

  !!Here lies my first attempt at an L2 norm.  I think this was wrong, I think this is the 
  !!L2 of the L1 norm in each original cell...
  !do i=1,numcellsLevel0
  !  first = 1+(i-1)*nextLevelFactor**ilevel
  !  last  = i*nextLevelFactor**ilevel
  !  Q_ufunctional(isamp,ilevel) = Q_ufunctional(isamp,ilevel) + &
  !            ( sum(flux(first:last))/(last-first+1) )**2
  !enddo
  !Q_ufunctional(isamp,ilevel) = sqrt(Q_ufunctional(isamp,ilevel))

  !solve G_ufunctional
  if( ilevel==0 ) then
    G_ufunctional(isamp,ilevel) = Q_ufunctional(isamp,ilevel)
  else
    G_ufunctional(isamp,ilevel) = Q_ufunctional(isamp,ilevel) - Q_ufunctional(isamp,ilevel-1)
  endif

  call FEDiffSn_externaldeallocate

  end subroutine solveSamples



  !5 test total error, set flMLMC==.false.?
  function MLMCcalcErrEst( Level )
  use MLMCvars, only: Gave, Gvar, C_alpha, ncellwidth, nextLevelFactor, M_optsamps
             
  real(8) :: MLMCcalcErrEst, err1, C_w, err2, Vsum
  integer :: Level, ilevel
  real(8) :: qq = 2.0d0  !Lect 10, pg 3, >=1, order of error, example of =2

  C_w  = max( abs(Gave(Level  ))/(ncellwidth(Level  )**qq*(nextLevelFactor**qq-1)), &
              abs(Gave(Level-1))/(ncellwidth(Level-1)**qq*(nextLevelFactor**qq-1))    )
  err1 = C_w * ncellwidth(Level)**qq
print *,"err1, disc err:",err1
  Vsum = 0.0d0
  do ilevel=0,Level
    Vsum = Vsum + abs(Gvar(ilevel)/M_optsamps(1,ilevel))
  enddo
  err2 = C_alpha * sqrt(Vsum)
print *,"err contribution Gvar       :",Gvar
print *,"err contribution M_optsamps :",M_optsamps
print *,"err contribution Vsum       :",Vsum,"   sqrt(Vsum):",sqrt(Vsum)
print *,"err2, MC   err:",err2
  MLMCcalcErrEst = err1 + err2

  end function MLMCcalcErrEst


end module MLMC
