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
  real(8)              :: GBsigaave,GBsigsave  ! Gauss-based abs and scat xs ave
  real(8)              :: GBsigavar,GBsigsvar  ! Gauss-based abs and scat xs var
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
  real(8)              :: scatvar              ! variance of scat xs in KL expansion, 'material' mode
  real(8)              :: absvar               ! variance of abs  xs in KL expansion, 'material' mode
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
  integer              :: numPosRealz          ! tally of num of realz that are always positive
  integer              :: numNegRealz          ! tally of num of realz that are at some point neg
contains
#ifdef USE_MPI
subroutine bcast_genRealzvars_vars
  use mpi
  implicit none
  integer :: ierr

  call MPI_Bcast(sig, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(scatrat, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(lam, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltgenrealz, 28, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(P, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numPath, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(totLength, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(perFirstTally, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(devFirstTally, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(matFirstTally, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sumPath, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sqrPath, 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(Adamscase, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(s, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numRealz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBsigaave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBsigavar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBsigsave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBsigsvar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBlamc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(GBs, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chgeomtype, 6, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(pltgenrealznumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(lamc, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(atmixsig, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(atmixscatrat, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(scatvar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(absvar, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigscatave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigabsave, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(largesti, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nummatSegs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numPosRealz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numNegRealz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_genRealzvars_vars


subroutine bcast_genRealzvars_alloc()
  use mpi
  implicit none

  if(.not.allocated(pltgenrealzwhich)) then
    allocate(pltgenrealzwhich(pltgenrealznumof))
    pltgenrealzwhich = 0
  endif
  return
end subroutine bcast_genRealzvars_alloc


subroutine bcast_genRealzvars_dealloc()
  use mpi
  implicit none

  if(allocated(pltgenrealzwhich)) deallocate(pltgenrealzwhich)
  return
end subroutine bcast_genRealzvars_dealloc


subroutine bcast_genRealzvars_arrays
  use mpi
  implicit none
  integer :: ierr

  if(allocated(pltgenrealzwhich)) call MPI_Bcast(pltgenrealzwhich, size(pltgenrealzwhich), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
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
  real(8), allocatable :: KLrxi(:)             !
  real(8), allocatable :: KLrxivalsa(:,:)      ! abs psuedo-random numbers for KL media, if same use this
  real(8), allocatable :: KLrxivalss(:,:)      ! scat psuedo-random numbers for KL media
  real(8), allocatable :: pltKLrealzarray(:,:) !
  real(8), allocatable :: KLrxisig(:)          !

  integer              :: KLrmaxnumzeros       ! number of zeros in realization with most, for allocating
  real(8), allocatable :: KLzerosabs(:,:)      ! location of zeros in abs xs KL reconstructions
  real(8), allocatable :: KLzerosscat(:,:)     ! location of zeros in scat xs KL reconstructions
  real(8), allocatable :: KLzerostotn(:,:)     ! location of zeros in native total xs KL reconstructions
contains
#ifdef USE_MPI
subroutine bcast_KLvars_vars
  use mpi
  implicit none
  integer :: ierr

  call MPI_Bcast(pltxiBins, 28, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltEigf, 28, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltCo, 28, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Corropts, 14, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltKLrealz, 28, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

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
  call MPI_Bcast(Gaussrandtype, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chGBcase, 4, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(pltxiBinsgauss, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltxiBinsnumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltEigfnumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltConumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(KLrnumpoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltKLrealznumof, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(flGaussdiffrand, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chGausstype, 4, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chLNmode, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chLNxschecktype, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numLNxspts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numLNxsbins, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(chLNxsplottype, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(sigsmeanadjust, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sigameanadjust, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mostinBin, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Corrnumpoints, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(binSize, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(KLrmaxnumzeros, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_KLvars_vars



subroutine bcast_KLvars_alloc()
  use mpi
  use genRealzvars, only: numRealz
  implicit none

  if(.not.allocated(pltEigfwhich)) then
    allocate(pltEigfwhich(pltEigfnumof))
    pltEigfwhich = 0
  endif
  if(.not.allocated(pltxiBinswhich)) then
    allocate(pltxiBinswhich(2,pltxiBinsnumof))
    pltxiBinswhich = 0
  endif
  if(.not.allocated(pltCowhich)) then
    allocate(pltCowhich(pltConumof))
    pltCowhich = 0
  endif
  if(.not.allocated(pltKLrealzwhich)) then
    allocate(pltKLrealzwhich(3,pltKLrealznumof))
    pltKLrealzwhich = 0
  endif

  if(.not.allocated(gam)) then
    allocate(gam(numEigs))
    gam = 0d0
  endif
  if(.not.allocated(alpha)) then
    allocate(alpha(numEigs))
    alpha = 0d0
  endif
  if(.not.allocated(Ak)) then
    allocate(Ak(numEigs))
    Ak = 0d0
  endif
  if(.not.allocated(Eig)) then
    allocate(Eig(numEigs))
    Eig = 0d0
  endif
  if(.not.allocated(xi)) then
    allocate(xi(numRealz,numEigs))
    xi = 0d0
  endif

  if(.not.allocated(binPDF)) then
    allocate(binPDF(binNumof,numEigs+1))
    binPDF = 0d0
  endif
  if(.not.allocated(binBounds)) then
    allocate(binBounds(binNumof+1))
    binBounds = 0d0
  endif
  if(.not.allocated(KLrxi)) then
    allocate(KLrxi(KLrnumpoints))
    KLrxi = 0d0
  endif
  if(.not.allocated(KLrxivalsa)) then
    allocate(KLrxivalsa(numRealz,anumEigs))
    KLrxivalsa = 0d0
  endif
  if(.not.allocated(KLrxivalss)) then
    allocate(KLrxivalss(numRealz,snumEigs))
    KLrxivalss = 0d0
  endif
  if(.not.allocated(pltKLrealzarray)) then
    allocate(pltKLrealzarray(KLrnumpoints,pltKLrealznumof+1))
    pltKLrealzarray = 0d0
  endif
  if(.not.allocated(KLrxisig)) then
    allocate(KLrxisig(KLrnumpoints))
    KLrxisig = 0d0
  endif

  if(.not.allocated(KLzerosabs)) then
    allocate(KLzerosabs(KLrmaxnumzeros,numRealz))
    KLzerosabs = 0d0
  endif
  if(.not.allocated(KLzerosscat)) then
    allocate(KLzerosscat(KLrmaxnumzeros,numRealz))
    KLzerosscat = 0d0
  endif
  if(.not.allocated(KLzerostotn)) then
    allocate(KLzerostotn(KLrmaxnumzeros,numRealz))
    KLzerostotn = 0d0
  endif

  return
end subroutine bcast_KLvars_alloc


subroutine bcast_KLvars_dealloc()
  use mpi
  implicit none

  if(allocated(pltEigfwhich)) deallocate(pltEigfwhich)
  if(allocated(pltxiBinswhich)) deallocate(pltxiBinswhich)
  if(allocated(pltCowhich)) deallocate(pltCowhich)
  if(allocated(pltKLrealzwhich)) deallocate(pltKLrealzwhich)

  if(allocated(gam)) deallocate(gam)
  if(allocated(alpha)) deallocate(alpha)
  if(allocated(Ak)) deallocate(Ak)
  if(allocated(Eig)) deallocate(Eig)
  if(allocated(xi)) deallocate(xi)

  if(allocated(binPDF)) deallocate(binPDF)
  if(allocated(binBounds)) deallocate(binBounds)
  if(allocated(KLrxi)) deallocate(KLrxi)
  if(allocated(KLrxivalsa)) deallocate(KLrxivalsa)
  if(allocated(KLrxivalss)) deallocate(KLrxivalss)
  if(allocated(pltKLrealzarray)) deallocate(pltKLrealzarray)
  if(allocated(KLrxisig)) deallocate(KLrxisig)

  if(allocated(KLzerosabs)) deallocate(KLzerosabs)
  if(allocated(KLzerosscat)) deallocate(KLzerosscat)
  if(allocated(KLzerostotn)) deallocate(KLzerostotn)
  return
end subroutine bcast_KLvars_dealloc



subroutine bcast_KLvars_arrays
  use mpi
  implicit none
  integer :: ierr

  if(allocated(pltEigfwhich)) call MPI_Bcast(pltEigfwhich, size(pltEigfwhich), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(allocated(pltxiBinswhich)) call MPI_Bcast(pltxiBinswhich, size(pltxiBinswhich), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(allocated(pltCowhich)) call MPI_Bcast(pltCowhich, size(pltCowhich), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(allocated(pltKLrealzwhich)) call MPI_Bcast(pltKLrealzwhich, size(pltKLrealzwhich), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  if(allocated(gam)) call MPI_Bcast(gam, size(gam), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(alpha)) call MPI_Bcast(alpha, size(alpha), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(Ak)) call MPI_Bcast(Ak, size(Ak), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(Eig)) call MPI_Bcast(Eig, size(Eig), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(xi)) call MPI_Bcast(xi, size(xi), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if(allocated(binPDF)) call MPI_Bcast(binPDF, size(binPDF), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(binBounds)) call MPI_Bcast(binBounds, size(binBounds), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(KLrxi)) call MPI_Bcast(KLrxi, size(KLrxi), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(KLrxivalsa)) call MPI_Bcast(KLrxivalsa, size(KLrxivalsa), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(KLrxivalss)) call MPI_Bcast(KLrxivalss, size(KLrxivalss), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(pltKLrealzarray)) &
                      call MPI_Bcast(pltKLrealzarray, size(pltKLrealzarray), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(KLrxisig)) call MPI_Bcast(KLrxisig, size(KLrxisig), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if(allocated(KLzerosabs)) call MPI_Bcast(KLzerosabs, size(KLzerosabs), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(KLzerosscat)) call MPI_Bcast(KLzerosscat, size(KLzerosscat), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(KLzerostotn)) call MPI_Bcast(KLzerostotn, size(KLzerostotn), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
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
    if(allocated(numPartsperj)) &
    call MPI_Reduce(MPI_IN_PLACE, numPartsperj, size(numPartsperj), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
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
    if(allocated(numPartsperj)) &
    call MPI_Reduce(numPartsperj, numPartsperj, size(numPartsperj), MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
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

  call MPI_Bcast(pltflux, 28, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(Wood_rej, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numpnSamp, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(areapnSamp, 4, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numcSamp, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(chTrantype, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(numParts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(maxnumParts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(reflrelSEMtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(tranrelSEMtol, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mindatapts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(LPamnumParts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(trannprt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(fluxnumcells, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  call MPI_Bcast(rodOrplanar, 6, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(sourceType, 8, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltmatflux, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(pltfluxtype, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(trprofile_binnum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(binplot, 7, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
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


subroutine bcast_MCvars_alloc()
  use mpi
  use genRealzvars, only: numRealz
  implicit none

  if(.not.allocated(numPartsperj)) then
    allocate(numPartsperj(numRealz))
    numPartsperj = 0
  endif
  if(.not.allocated(ABreflection)) then
    allocate(ABreflection(2,5))
    ABreflection = 0d0
  endif
  if(.not.allocated(ABtransmission)) then
    allocate(ABtransmission(2,5))
    ABtransmission = 0d0
  endif

  if(.not.allocated(fluxfaces)) then
    allocate(fluxfaces(fluxnumcells+1))
    fluxfaces = 0d0
  endif
  if(.not.allocated(fluxall)) then
    allocate(fluxall(fluxnumcells,numRealz))
    fluxall = 0d0
  endif
!  if(.not.allocated(fluxmat1)) then
!    allocate(fluxmat1(,))
!    fluxmat1 = 0d0
!  endif
!  if(.not.allocated(fluxmat2)) then
!    allocate(fluxmat2(,))
!    fluxmat2 = 0d0
!  endif
!  if(.not.allocated(fluxmatnorm)) then
!    allocate(fluxmatnorm(,,))
!    fluxmatnorm = 0d0
!  endif

  if(.not.allocated(reflect)) then
    allocate(reflect(numRealz))
    reflect = 0d0
  endif
  if(.not.allocated(transmit)) then
    allocate(transmit(numRealz))
    transmit = 0d0
  endif
  if(.not.allocated(absorb)) then
    allocate(absorb(numRealz))
    absorb = 0d0
  endif
  if(.not.allocated(LPamMCsums)) then
    allocate(LPamMCsums(3))
    LPamMCsums = 0d0
  endif
  return
end subroutine bcast_MCvars_alloc


subroutine bcast_MCvars_dealloc()
  use mpi
  implicit none

  if(allocated(numPartsperj)) deallocate(numPartsperj)
  if(allocated(ABreflection)) deallocate(ABreflection)
  if(allocated(ABtransmission)) deallocate(ABtransmission)

  if(allocated(fluxfaces)) deallocate(fluxfaces)
  if(allocated(fluxall)) deallocate(fluxall)
  if(allocated(fluxmat1)) deallocate(fluxmat1)
  if(allocated(fluxmat2)) deallocate(fluxmat2)
  if(allocated(fluxmatnorm)) deallocate(fluxmatnorm)

  if(allocated(reflect)) deallocate(reflect)
  if(allocated(transmit)) deallocate(transmit)
  if(allocated(absorb)) deallocate(absorb)
  if(allocated(LPamMCsums)) deallocate(LPamMCsums)

  if(allocated(binmaxind)) deallocate(binmaxind)
  if(allocated(binmaxes)) deallocate(binmaxes)
  if(allocated(fbinmax)) deallocate(fbinmax)
  if(allocated(bbinmax)) deallocate(bbinmax)
  return
end subroutine bcast_MCvars_dealloc


subroutine bcast_MCvars_arrays
  use mpi
  use mpiaccess
  implicit none
  integer :: ierr
  if(allocated(numPartsperj)) call MPI_Bcast(numPartsperj, size(numPartsperj), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if(allocated(ABreflection)) call MPI_Bcast(ABreflection, size(ABreflection), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(ABtransmission)) call MPI_Bcast(ABtransmission, size(ABtransmission), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if(allocated(fluxfaces)) call MPI_Bcast(fluxfaces, size(fluxfaces), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(fluxall)) call MPI_Bcast(fluxall, size(fluxall), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(fluxmat1)) call MPI_Bcast(fluxmat1, size(fluxmat1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(fluxmat2)) call MPI_Bcast(fluxmat2, size(fluxmat2), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(fluxmatnorm)) call MPI_Bcast(fluxmatnorm, size(fluxmatnorm), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if(allocated(reflect)) call MPI_Bcast(reflect, size(reflect), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(transmit)) call MPI_Bcast(transmit, size(transmit), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(absorb)) call MPI_Bcast(absorb, size(absorb), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(LPamMCsums)) call MPI_Bcast(LPamMCsums, size(LPamMCsums), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if(allocated(binmaxind)) call MPI_Bcast(binmaxind, size(binmaxind), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(binmaxes)) call MPI_Bcast(binmaxes, size(binmaxes), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(fbinmax)) call MPI_Bcast(fbinmax, size(fbinmax), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if(allocated(bbinmax)) call MPI_Bcast(bbinmax, size(bbinmax), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
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

  call MPI_Bcast(chUQtype, 5, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_UQvars_vars


subroutine bcast_UQvars_alloc()
  use mpi
  use genRealzvars, only: numRealz
  use KLvars, only: snumEigs,anumEigs,flGaussdiffrand
  implicit none

print *,"flGaussdiffrand:",flGaussdiffrand,"s/anumEigs",snumEigs,anumEigs
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

  if(allocated(Qs)) deallocate(Qs)
  if(allocated(UQwgts)) deallocate(UQwgts)
  return
end subroutine bcast_UQvars_dealloc


subroutine bcast_UQvars_arrays
  use mpi

  use mpiaccess
  implicit none
  integer :: ierr
print *,"inn jobid-1:",jobid,"  size(Qs):",size(Qs)
  if(allocated(Qs)) call MPI_Bcast(Qs, size(Qs), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
print *,"inn jobid-2:",jobid
  if(allocated(UQwgts)) call MPI_Bcast(UQwgts, size(UQwgts), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
print *,"inn jobid-3:",jobid
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  return
end subroutine bcast_UQvars_arrays
#endif

end module UQvars 

