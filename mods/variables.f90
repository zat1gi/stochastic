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
    case ("radMC")         !transport           !possible correlation 1
      rngappnum = 1
    case ("radWood")                            !possible correlation 1
      rngappnum = 2
    case ("KLWood")                             !possible correlation 2
      rngappnum = 3
    case ("LPMC")
      rngappnum = 4
    case ("atmixMC")
      rngappnum = 5
    case ("GaussKL")                            !possible correlation 2
      rngappnum = 6
    case ("genRealzKLres") !material generation !possible correlation 3
      rngappnum = 10
    case ("genRealzTMC")                        !possible correlation 3
      rngappnum = 11
    case ("genRealzWMC")                        !possible correlation 3
      rngappnum = 12
    case ("KLRealzMarkov")                      !possible correlation 4
      rngappnum = 13
    case ("KLRealzGaussB")                      !possible correlation 4
      rngappnum = 14
  end select
  end subroutine setrngappnum
end module rngvars


module timevars
  implicit none
  real(8), allocatable :: time(:)              ! keeps cumulative time data for everything
                                               ! 1)genR, 2)radMC, 3)radWood, 4)KLnoise, 5)KLcol 
                                               ! 6)KLrec, 7)KLWood, 8)LPMC, 9)atmixMC, 10)GaussKL
  integer, parameter   :: ntime = 10           ! size of 'time' array
  real(8), allocatable :: FOM(:,:)             ! Figure of Merit = (1/(Var*cputime))
                                               ! (1==radMC/2==radWood/3==KLWood, 1==refl/2==tran)     
  real(8)              :: t1                   ! initial computer time for run
  real(8)              :: runtime              ! total cpu time code has been running
  integer, allocatable :: totparts             ! total particles to be run
  integer, allocatable :: cumparts             ! cumulative particles run

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

  integer              :: pltgenrealznumof     !
  character(7)         :: pltgenrealz(4)       !
  integer, allocatable :: pltgenrealzwhich(:)  !
  !non inputs
  logical              :: flprint              ! flag to print if printing selected
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
  logical              :: flCorrMarkov=.false. ! correlated sampling for Markov realizations?
  logical              :: flCorrRealz =.false. ! correlated KLres (for KL recon) and Markov realz?
  logical              :: flGBgeom    =.true.  ! Gauss-based geom? (or defer to Markov-based)

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

end module genRealzvars



!--- KL ---!
module KLvars  !"KLresearch" and "KLreconstruct"
  implicit none
  !inputs
  character(3)         :: KLres                ! perform "KL research" operations?
  character(3)         :: KLrec                ! perform "KL reconstruct" operations?
  character(3)         :: KLnoise              ! perform "KL noise" operations?
  character(3)         :: KLvarcalc            ! calculate amount of variance kept in eigs? 'yes' 'no
  real(8)              :: KLvarkept_tol        ! tolerance of rel eig size for var calcs
  integer              :: binNumof             ! number of bins (for xi?)
  integer              :: numEigs              ! number of eigenmodes to calculate
  integer              :: numSlice             ! number of points to plot eigenfunction at
  integer              :: levsrefEig           !
  real(8)              :: binSmallBound        ! smallest xi value for xi bins
  real(8)              :: binLargeBound        ! largest xi value for xi bins
  logical              :: flmeanadjust         ! flag, perform mean adjust or not
  real(8)              :: meanadjust_tol       ! tolerance for new mean adjustment
  integer              :: numrefinesameiter    ! num of iters for extrema determ w/ no change to stop
  logical              :: flMarkov=.false.     ! run Markov-based KL routines?
  logical              :: flGauss =.false.     ! run Gauss-random-based KL routines?
  character(3)         :: Gaussrandtype        ! 'BM'-Box-Muller sampling or 'inv' inverse sampling

  character(7)         :: pltxiBins(4)         !
  character(7)         :: pltxiBinsgauss       !
  character(7)         :: pltEigf(4)           !
  character(7)         :: pltCo(4)             !
  character(7)         :: Corropts(2)          !
  integer, allocatable :: pltEigfwhich(:)      ! 
  integer, allocatable :: pltxiBinswhich(:,:)  ! 
  integer, allocatable :: pltCowhich(:,:)      ! 
  integer              :: pltxiBinsnumof       ! 
  integer              :: pltEigfnumof         !
  integer              :: pltConumof           !
  integer              :: KLrnumpoints         ! Number of points to plot in KL reconstructions
  character(7)         :: pltKLrealz(4)        !
  integer              :: pltKLrealznumof      !
  integer, allocatable :: pltKLrealzwhich(:,:) !

  logical              :: flglGaussdiffrand=.true.! global - want GB meths to use different rand vars for xi samps?
  logical              :: flglLN = .false.     ! global flag to denote 'GaussKL' is Log-normal
  character(7)         :: chLNmode = 'Glamc'   ! Log-normal cov and lamc:'Glamc'-Gausslamc,'fitlamc'-expfit,'numeric'
  logical              :: flLNxscheck = .false.! perform mean, var, and pdf check of LN xs distributions?
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
  real(8)              :: meanadjust = 0.0d0   ! positive translation of mean xs (totxs mode)
  real(8)              :: sigsmeanadjust=0.0d0 ! positive translation of sigs mean xs (mat-based mode)
  real(8)              :: sigameanadjust=0.0d0 ! positive translation of siga mean xs (mat-based mode)
  logical              :: flCorrKL=.false.     ! correlated random numbers for KL Markov and GaussB realz?
  logical              :: flLN=.false.         ! local flag to denote we are using Log-normal modeling
  integer              :: mostinBin            !
  integer              :: Corrnumpoints        ! Number of points used when confirming covariance func
  real(8), allocatable :: binPDF(:,:)          ! 
  real(8), allocatable :: binBounds(:)         !
  real(8), allocatable :: AllEig(:)            ! many more Eigs than used, for variance kept calcs
  real(8), allocatable :: Allgam(:)            ! many more gam vals than used, for variance kept calcs
  real(8), allocatable :: varmain(:)           ! variance maintained at certain eigenmode
  real(8)              :: binSize              !
  real(8), allocatable :: KLrx(:)              !
  real(8), allocatable :: KLrxi(:)             !
  real(8), allocatable :: KLrxivals(:,:)       ! psuedo-random numbers for KL media
  real(8), allocatable :: KLrxivalss(:,:)      ! same, but for scattering xs when Gauss-Based methods (diff xis)
  logical              :: flGaussdiffrand      ! local Gauss-Based methods to use different rand vars for xi samps?
  real(8), allocatable :: pltKLrealzarray(:,:) !
  real(8), allocatable :: KLrxisig(:)          !

  integer              :: KLrmaxnumzeros       ! number of zeros in realization with most, for allocating
  real(8), allocatable :: KLr_maxima(:,:)      ! location of local maxima in KL reconstructions
  real(8), allocatable :: KLzerostot(:,:)      ! location of zeros in tot xs KL reconstructions
  real(8), allocatable :: KLzerosabs(:,:)      ! location of zeros in abs xs KL reconstructions
  real(8), allocatable :: KLzerosscat(:,:)     ! location of zeros in scat xs KL reconstructions
  real(8), allocatable :: KLzerostotn(:,:)     ! location of zeros in native total xs KL reconstructions

end module KLvars



!--- Transport ---!
module MCvars
  implicit none
  !inputs
  character(7)         :: chTrantype           ! type of transport to perform
  integer              :: numParts             ! number of particles
  integer              :: LPamnumParts         ! number of particles for LP or atomic mix
  integer              :: trannprt             ! how often to print to screen
  integer              :: fluxnumcells         ! number of cells for flux profile

  character(6)         :: rodOrplanar          ! transport geometry mode
  character(6)         :: sourceType           ! 'intern' or 'left', distributed or beam source
  character(7)         :: pltflux(4)           ! plot material irrespective flux
  character(7)         :: pltmatflux           ! plot material respective fluxes
  character(7)         :: pltfluxtype          ! full 'track' length flux plot or 'point' flux binning? 
  integer              :: trprofile_binnum     ! number of bins for tran and refl profiles
  character(7)         :: binplot              ! 'plot', 'noplot', 'preview', pdf of leakage

  logical              :: flnegxs=.false.      ! allow trans on neg xs? 'yes', or throw out realz 'no'
  logical              :: fldistneg=.false.    ! allow on the fly smoothing of negs?

  !non inputs
  logical              :: flCorrMC=.false.     ! correlated random number MC transport

  real(8), allocatable :: ABreflection(:,:)    ! Adams/Brantley Reflection Values
  real(8), allocatable :: ABtransmission(:,:)  ! Adams/Brantley Transmission Values
                                               ! First rank: 1 average, 2 stdev
                                               ! Second rank: 1 AdMC, 2 AdLP, 3 BrMC, 4 BrLP, 5 Bratmix

  integer              :: radtrans_int         !

  real(8), allocatable :: fluxfaces(:)         ! mesh for flux tallies
  logical              :: flfluxplot           ! run fluxplot stuff or not?
  logical              :: flfluxplotall        ! run fluxplot for material irrespective
  logical              :: flfluxplotmat        ! run fluxplot for material respective
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
  real(8)              :: disthold             ! used in a workaround for neg vals in KLWood

  real(8)              :: position             ! 'current' position of particle
  real(8)              :: oldposition          ! most recent position of particle
  real(8)              :: mu                   ! 'current' direction of particle
  integer              :: nceilbin             ! number of bins in ceiling calcs for WMC
end module MCvars

