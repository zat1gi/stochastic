module rngvars
  implicit none
  integer(8)           :: rngseed              ! user input which seeds rng stream

  integer              :: rngappnum            ! application number
  integer(8)           :: rngstride=1668163541 ! num of 'part histories' to skip between apps
contains
subroutine setrngappnum( rngapp )
  !this subroutine sets rngappnum according to what application, to standardize
  !rngs used for that application
  character(*) :: rngapp
  select case (rngapp)
    case ("genRealz")      !allow reproducible binary geometries
      rngappnum = 1
    case ("KLRealz")       !allow reproducible KL geometries
      rngappnum = 2
    case ("radtrans")      !allow reproducible transport/correlated transport
      rngappnum = 3
  end select
end subroutine setrngappnum
#ifdef USE_MPI
subroutine bcast_rngvars_vars
  use mpi
  implicit none
  integer :: ierr

  call MPI_Bcast(rngseed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(rngstride, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_rngvars_vars
#endif
end module rngvars


module timevars
  implicit none
  real(8)              :: t1                   ! initial computer time for run
end module timevars


module genRealzvars
  implicit none
  !inputs
  real(8)              :: Adamscase            ! load special cases or no
  real(8),dimension(2) :: sig                  ! cross sections
  real(8),dimension(2) :: scatrat              ! scattering ratios
  real(8),dimension(2) :: lam                  ! ave path lengths
  real(8)              :: s                    ! slab thickness
  integer              :: numRealz             ! number of realz to create
  real(8)              :: GBsigave, GBsigvar   ! Gauss-based sig average and variance
  real(8)              :: GBscatrat            ! Gauss-based scattering ratio
  real(8)              :: GBlamc               ! Gauss-based correlation length
  real(8)              :: GBs                  ! Gauss-based slab thickness
  character(6)         :: chgeomtype           ! 'contin'uous or 'binary'

  integer              :: pltgenrealznumof     !
  character(7)         :: pltgenrealz(4)       !
  integer, allocatable :: pltgenrealzwhich(:)  !
  !non inputs
  real(8)              :: P(2)                 ! probability of mat1 or 2 for binary mixtures
  real(8)              :: lamc                 ! "correlation length" for binary mixtures  
  integer, allocatable :: matType(:)           ! material type in cell
  real(8), allocatable :: matLength(:)         ! material boundaries
  real(8)              :: atmixsig             ! atomically mixed cross section value
  real(8)              :: atmixscatrat         ! atomically mixed scattering ratio
  real(8)              :: sigave_              ! maintains input sigave when sigave is translated to LN
  real(8)              :: sigvar               ! variance of total xs in KL expansion
  real(8)              :: sigvar_              ! maintains input sigvar when sigvar is translated to LN
  real(8)              :: scatvar              ! variance of scat xs in KL expansion, 'material' mode
  real(8)              :: absvar               ! variance of abs  xs in KL expansion, 'material' mode
  real(8)              :: sigave               ! weighted average sigma value for binary mixtures
  real(8)              :: sigscatave           ! average scattering xs
  real(8)              :: sigabsave            ! average absorption xs

  integer              :: largesti             !
  integer              :: numPath(2)           !
  real(8)              :: totLength(2)         ! summed length of each material in stoc family
  integer              :: nummatSegs           ! number of material segments in realz (rid, use size?)
  real(8)              :: perFirstTally(2)     ! percentage of first material as mat1 or 2, bin mix
  real(8)              :: devFirstTally(2)     ! stdev of first material as mat1 or 2, bin mix
  real(8)              :: matFirstTally(2)=0   ! tally of first material as mat1 or 2, bin mix
  real(8)              :: sumPath(2)           ! (consider making these guys local)
  real(8)              :: sqrPath(2)           !
  real(8), allocatable :: matdxs(:,:,:)        ! delta x of different materials for flux tallying
  integer              :: numPosRealz          ! tally of num of realz that are always positive
  integer              :: numNegRealz          ! tally of num of realz that are at some point neg
contains
#ifdef USE_MPI
subroutine bcast_genRealzvars_vars
  use mpi
  implicit none
  integer :: ierr

  call MPI_Bcast(Adamscase, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numRealz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBsigave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBscatrat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBlamc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chgeomtype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(pltgenrealznumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(lamc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(atmixsig, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(atmixscatrat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigave_, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigvar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigvar_, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(scatrat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(absvar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigscatave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigabsave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(largesti, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nummatSegs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numPosRealz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numNegRealz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_genRealzvars_vars


subroutine bcast_genRealzvars_alloc_de(flalloc)
  use mpi
  implicit none
  logical :: flalloc
  integer :: ierr

  return
end subroutine bcast_genRealzvars_alloc_de


subroutine bcast_genRealzvars_arrays
  use mpi
  implicit none
  integer :: ierr

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_genRealzvars_arrays
#endif

end module genRealzvars



!--- KL ---!
module KLvars  !"KLresearch" and "KLconstruct"
  implicit none
  !inputs
  integer              :: binNumof             ! number of bins (for xi?)
  integer              :: numEigs              ! number of KL eigenmodes, if==0, defer to next two
  integer              :: snumEigs             ! number of KL eigenmodes for scattering xs
  integer              :: anumEigs             ! number of KL eigenmodes for absorption xs
  integer              :: numSlice             ! number of points to plot eigenfunction at
  integer              :: levsrefEig           !
  real(8)              :: binSmallBound        ! smallest xi value for xi bins
  real(8)              :: binLargeBound        ! largest xi value for xi bins
  logical              :: flmeanadjust         ! flag, perform mean adjust or not
  real(8)              :: meanadjust_tol       ! tolerance for new mean adjustment
  integer              :: numrefinesameiter    ! num of iters for extrema determ w/ no change to stop
  character(3)         :: Gaussrandtype        ! 'BM'-Box-Muller sampling or 'inv' inverse sampling
  character(4)         :: chGBcase             ! special GB cases, 'none', 'f1', 'f2' (Fichtl 1, 2)

  character(7)         :: pltxiBins(4)         !
  character(7)         :: pltxiBinsgauss       !
  character(7)         :: pltEigf(4)           !
  character(7)         :: pltCo(4)             !
  character(7)         :: Corropts(2)          !
  integer, allocatable :: pltEigfwhich(:)      ! 
  integer, allocatable :: pltxiBinswhich(:,:)  ! 
  integer, allocatable :: pltCowhich(:)        ! number of eigenmodes to plot with
  integer              :: pltxiBinsnumof       ! 
  integer              :: pltEigfnumof         !
  integer              :: pltConumof           !
  integer              :: KLrnumpoints         ! Number of points to plot in KL reconstructions
  character(7)         :: pltKLrealz(4)        !
  integer              :: pltKLrealznumof      ! number or realizations to plot
  integer, allocatable :: pltKLrealzwhich(:,:) ! realization number and num of eigenmodes to plot

  logical              :: flGaussdiffrand=.true.! want GB meths to use different rand vars for xi samps?
  character(4)         :: chGausstype          ! Gauss-Based mode: 'Gaus','LogN','ChiS'
  character(7)         :: chLNmode = 'Glamc'   ! LN cov and lamc:'Glamc'-Gausslamc,'fitlamc'-expfit,'numeric'
  character(7)         :: chLNxschecktype      ! type of cross section to analyze
  integer              :: numLNxspts           ! number of x-values to perform checks at
  integer              :: numLNxsbins          ! number of bins to use in creating pdf of values
  character(7)         :: chLNxsplottype       ! 'noplot','preview', or 'plot'
  !non-inputs
  real(8), allocatable :: gam(:)               ! solutions to eigenvalue transcendental
  real(8), allocatable :: alpha(:)             ! other form of gam
  real(8), allocatable :: Ak(:)                ! normalization coefficients in KL expansion
  real(8), allocatable :: Eig(:)               ! eigenvalues ok KL expansion
  real(8), allocatable :: xi(:,:)              ! array of chosen xi values for reusing reconstructions
  real(8)              :: sigsmeanadjust=0.0d0 ! positive translation of sigs mean xs (mat-based mode)
  real(8)              :: sigameanadjust=0.0d0 ! positive translation of siga mean xs (mat-based mode)
  integer              :: mostinBin            !
  integer              :: Corrnumpoints        ! Number of points used when confirming covariance func
  real(8), allocatable :: binPDF(:,:)          ! 
  real(8), allocatable :: binBounds(:)         !
  real(8)              :: binSize              !
  real(8), allocatable :: KLrx(:)              !
  real(8), allocatable :: KLrxi(:)             !
  real(8), allocatable :: KLrxivalsa(:,:)      ! abs psuedo-random numbers for KL media, if same use this
  real(8), allocatable :: KLrxivalss(:,:)      ! scat psuedo-random numbers for KL media
  real(8), allocatable :: pltKLrealzarray(:,:) !
  real(8), allocatable :: KLrxisig(:)          !

  integer              :: KLrmaxnumzeros       ! number of zeros in realization with most, for allocating
  real(8), allocatable :: KLr_maxima(:,:)      ! location of local maxima in KL reconstructions
  real(8), allocatable :: KLzerostot(:,:)      ! location of zeros in tot xs KL reconstructions
  real(8), allocatable :: KLzerosabs(:,:)      ! location of zeros in abs xs KL reconstructions
  real(8), allocatable :: KLzerosscat(:,:)     ! location of zeros in scat xs KL reconstructions
  real(8), allocatable :: KLzerostotn(:,:)     ! location of zeros in native total xs KL reconstructions
contains
#ifdef USE_MPI
subroutine bcast_KLvars_vars
  use mpi
  implicit none
  integer :: ierr

  call MPI_Bcast(binNumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numEigs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(snumEigs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(anumEigs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numSlice, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(levsrefEig, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(binSmallBound, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(binLargeBound, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(flmeanadjust, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(meanadjust_tol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numrefinesameiter, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Gaussrandtype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chGBcase, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(pltxiBinsgauss, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltxiBinsnumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltEigfnumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltConumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(KLrnumpoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltKLrealznumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(flGaussdiffrand, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chGausstype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chLNmode, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chLNxschecktype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numLNxspts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numLNxsbins, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chLNxsplottype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(sigsmeanadjust, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigameanadjust, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mostinBin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Corrnumpoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(binSize, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(KLrmaxnumzeros, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_KLvars_vars


subroutine bcast_KLvars_alloc_de(flalloc)
  use mpi
  implicit none
  logical :: flalloc
  integer :: ierr

  return
end subroutine bcast_KLvars_alloc_de


subroutine bcast_KLvars_arrays
  use mpi
  implicit none
  integer :: ierr

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_KLvars_arrays
#endif

end module KLvars



!--- Transport ---!
module MCvars
  implicit none
  !inputs
  character(7)         :: chTrantype           ! type of transport to perform
  integer              :: numParts             ! number of particles
  integer              :: maxnumParts          ! maximum number of particles to run per realization
  real(8)              :: reflrelSEMtol        ! reflection relative SEM tolerance must converge to
  real(8)              :: tranrelSEMtol        ! transmission relative SEM tolerance must converge to
  integer              :: mindatapts           ! minimum number of data points to trust statistics
  integer              :: LPamnumParts         ! number of particles for LP or atomic mix
  integer              :: trannprt             ! how often to print to screen
  integer              :: fluxnumcells         ! number of cells for flux profile

  character(6)         :: rodOrplanar          ! transport geometry mode
  character(8)         :: sourceType           ! 'leftbeam','leftiso', or 'intern' - souce type
  character(7)         :: pltflux(4)           ! plot material irrespective flux
  character(7)         :: pltmatflux           ! plot material respective fluxes
  character(7)         :: pltfluxtype          ! full 'track' length flux plot or 'point' flux binning? 
  integer              :: trprofile_binnum     ! number of bins for tran and refl profiles
  character(7)         :: binplot              ! 'plot', 'noplot', 'preview', pdf of leakage
  logical              :: flCR_MCSC=.false.    ! correlated random number seeds at particle histories

  logical              :: flnegxs=.false.      ! allow trans on neg xs? 'yes', or throw out realz 'no'

  !non inputs
  integer, allocatable :: numPartsperj(:)      ! number of particles actually ran per realization
  real(8), allocatable :: ABreflection(:,:)    ! Adams/Brantley Reflection Values
  real(8), allocatable :: ABtransmission(:,:)  ! Adams/Brantley Transmission Values
                                               ! First rank: 1 average, 2 stdev
                                               ! Second rank: 1 AdMC, 2 AdLP, 3 BrMC, 4 BrLP, 5 Bratmix

  integer              :: radtrans_int         !

  real(8), allocatable :: fluxfaces(:)         ! mesh for flux tallies
  logical              :: flfluxplot=.false.   ! run fluxplot stuff or not?
  logical              :: flfluxplotall=.false.! run fluxplot for material irrespective
  logical              :: flfluxplotmat=.false.! run fluxplot for material respective
  real(8), allocatable :: fluxall(:,:)         ! flux of a method, all mats, (fluxnumbins,numRealz)
  real(8), allocatable :: fluxmat1(:,:)        ! flux of a method, mat1, (fluxnumbins,numRealz)
  real(8), allocatable :: fluxmat2(:,:)        ! flux of a method, mat2, (fluxnumbins,numRealz)
  real(8), allocatable :: fluxmatnorm(:,:,:)   ! amount of material in bins, (fluxnumbins,numRealz,2) 
  real(8), allocatable :: stocMC_fluxall(:,:)  ! flux in cells, (fluxnumbins,nummeths,2)
  real(8), allocatable :: stocMC_fluxmat1(:,:) ! flux in mat1, '2' is mean and var
  real(8), allocatable :: stocMC_fluxmat2(:,:) ! flux in mat2

  real(8), allocatable :: reflect(:)           ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: transmit(:)          ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: absorb(:)            ! slab absorption tally, TMC on binary mixtures

  real(8), allocatable :: LPamMCsums(:)        ! for LP or atomic mix, sum of tallies, refl, tran, abs

                                               ! stoc means mean and variance accross stochastic domain
                                               ! MC means for MC transport solves in spatial domain
  real(8), allocatable :: stocMC_reflection(:) ! reflection, mean & var in stoc space
  real(8), allocatable :: stocMC_transmission(:)! trans, mean & var in stoc space
  real(8), allocatable :: stocMC_absorption(:) ! absorption, mean & var in stoc space


  integer              :: Wood_rej(2)          ! generic Woodcock rejection tally
  integer              :: numpnSamp(2)         ! number of pos, neg samples (in KL reconstuctions)
  real(8)              :: areapnSamp(4)        ! area pos,neg sampled: totpos, totneg, maxpos, maxneg
  integer              :: numcSamp(2)          ! number of neg/tot KL scatrat samples

  real(8),allocatable  :: binmaxind(:)         ! These are used in ceilings for WMC
  real(8),allocatable  :: binmaxes(:)          ! max value in each bin
  real(8),allocatable  :: fbinmax(:)           ! forward  facing bin max values
  real(8),allocatable  :: bbinmax(:)           ! backward facing bin max values

  real(8)              :: position             ! 'current' position of particle
  real(8)              :: oldposition          ! most recent position of particle
  real(8)              :: mu                   ! 'current' direction of particle
  integer              :: nceilbin             ! number of bins in ceiling calcs for WMC
  contains

#ifdef USE_MPI
subroutine reduceMCresults
  use mpi
  use mpiaccess
  implicit none
  integer :: ierr

  if(jobid==0) then
    if(allocated(fluxall)) &
    call MPI_Reduce(MPI_IN_PLACE, fluxall, size(fluxall), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(fluxmat1)) &
    call MPI_Reduce(MPI_IN_PLACE, fluxmat1, size(fluxmat1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(fluxmat2)) &
    call MPI_Reduce(MPI_IN_PLACE, fluxmat2, size(fluxmat2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(fluxmatnorm)) &
    call MPI_Reduce(MPI_IN_PLACE, fluxmatnorm, size(fluxmatnorm), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(reflect)) &
    call MPI_Reduce(MPI_IN_PLACE, reflect, size(reflect), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(transmit)) &
    call MPI_Reduce(MPI_IN_PLACE, transmit, size(transmit), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(absorb)) &
    call MPI_Reduce(MPI_IN_PLACE, absorb, size(absorb), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(LPamMCsums)) &
    call MPI_Reduce(MPI_IN_PLACE, LPamMCsums, size(LPamMCsums), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, Wood_rej, size(Wood_rej), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, numpnSamp, size(numpnSamp), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, areapnSamp, size(areapnSamp), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(MPI_IN_PLACE, numcSamp, size(numcSamp), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  else
    if(allocated(fluxall)) &
    call MPI_Reduce(fluxall, fluxall, size(fluxall), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(fluxmat1)) &
    call MPI_Reduce(fluxmat1, fluxmat1, size(fluxmat1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(fluxmat2)) &
    call MPI_Reduce(fluxmat2, fluxmat2, size(fluxmat2), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(fluxmatnorm)) &
    call MPI_Reduce(fluxmatnorm, fluxmatnorm, size(fluxmatnorm), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(reflect)) &
    call MPI_Reduce(reflect, reflect, size(reflect), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(transmit)) &
    call MPI_Reduce(transmit, transmit, size(transmit), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(absorb)) &
    call MPI_Reduce(absorb, absorb, size(absorb), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if(allocated(LPamMCsums)) &
    call MPI_Reduce(LPamMCsums, LPamMCsums, size(LPamMCsums), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(Wood_rej, Wood_rej, size(Wood_rej), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(numpnSamp, numpnSamp, size(numpnSamp), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(areapnSamp, areapnSamp, size(areapnSamp), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_Reduce(numcSamp, numcSamp, size(numcSamp), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  endif
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
end subroutine reduceMCresults

subroutine bcast_MCvars_vars
  use mpi
  implicit none
  integer :: ierr
  call MPI_Bcast(chTrantype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numParts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(maxnumParts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(reflrelSEMtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(tranrelSEMtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mindatapts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(LPamnumParts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(trannprt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fluxnumcells, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(rodOrplanar, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sourceType, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltmatflux, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltfluxtype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(trprofile_binnum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(binplot, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(flCR_MCSC, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(flnegxs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(radtrans_int, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(flfluxplot, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(flfluxplotall, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(flfluxplotmat, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(position, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(oldposition, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mu, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nceilbin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_MCvars_vars


subroutine bcast_MCvars_alloc_de(flalloc)
  use mpi
  implicit none
  logical :: flalloc
  integer :: ierr

  return
end subroutine bcast_MCvars_alloc_de


subroutine bcast_MCvars_arrays
  use mpi
  implicit none
  integer :: ierr

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_MCvars_arrays
#endif

end module MCvars



!--- UQ space (MC, SC) ---!
module UQvars
  implicit none
  !inputs
  character(5)         :: chUQtype             ! 'MC', 'LagSC', 'PCESC', UQ approach
  integer, allocatable :: Qs(:)                ! for SC, order in each dimension
  !non inputs
  real(8), allocatable :: UQwgts(:)            ! for 'MC', 1/numRealz, for 'xxxSC', cubature wgts
contains
#ifdef USE_MPI
subroutine bcast_UQvars_vars
  use mpi
  implicit none
  integer :: ierr

  call MPI_Bcast(chUQtype, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_UQvars_vars


subroutine bcast_UQvars_alloc()
  use mpi
  use genRealzvars, only: numRealz
  use KLvars, only: snumEigs,anumEigs,flGaussdiffrand
  implicit none
  integer :: ierr

  if(.not.allocated(Qs)) then
    if(flGaussdiffrand) then
      allocate(Qs(snumEigs+anumEigs))
    else
      allocate(Qs(snumEigs))
    endif
    Qs = 0
  endif
  if(.not.allocated(UQwgts)) then
    allocate(UQwgts(numRealz))
    UQwgts = 0d0
  endif
  return
end subroutine bcast_UQvars_alloc


subroutine bcast_UQvars_dealloc()
  use mpi
  implicit none
  integer :: ierr

  if(allocated(Qs)) deallocate(Qs)
  if(allocated(UQwgts)) deallocate(UQwgts)
  return
end subroutine bcast_UQvars_dealloc


subroutine bcast_UQvars_arrays
  use mpi
  implicit none
  integer :: ierr

  if(allocated(Qs)) call MPI_Bcast(Qs, size(Qs), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(allocated(UQwgts)) call MPI_Bcast(UQwgts, size(UQwgts), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_UQvars_arrays
#endif

end module UQvars 

