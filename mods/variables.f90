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
    case ("radMC")
      rngappnum = 1
    case ("radWood")
      rngappnum = 2
    case ("KLWood")
      rngappnum = 3
    case ("LPMC")
      rngappnum = 4
    case ("atmixMC")
      rngappnum = 5
    case ("genRealz")
      rngappnum = 6
    case ("KLRealz")
      rngappnum = 7
  end select
  end subroutine setrngappnum
end module rngvars


module timevars
  implicit none
  real(8), allocatable :: time(:)              ! keeps cumulative time data for everything
                                               ! 1)genR, 2)radMC, 3)radWood, 4)KLnoise, 
                                               ! 5)KLcol, 6)KLrec, 7)KLWood, 8)LPMC, 9)atmixMC
  integer, parameter   :: ntime = 9            ! size of 'time' array
  real(8), allocatable :: FOM(:,:)             ! Figure of Merit = (1/(Var*cputime))
                                               ! (1==radMC/2==radWood/3==KLWood, 1==refl/2==tran)     
  integer, parameter   :: nFOM = 3             ! size of 'FOM' array
  real(8)              :: t1                   ! initial computer time for run
  real(8)              :: runtime              ! total cpu time code has been running
  integer, allocatable :: totparts(:)          ! total particles to be run in each method
  integer, allocatable :: cumparts(:)          ! cumulative particles run in each method

end module


module genRealzvars
  implicit none
  !inputs
  real(8)              :: Adamscase            ! load special cases or no
  real(8),dimension(2) :: sig                  ! cross sections
  real(8),dimension(2) :: scatrat              ! scattering ratios
  real(8),dimension(2) :: lam                  ! ave path lengths
  real(8)              :: s                    ! slab thickness
  integer              :: numRealz             ! number of realz to create
  
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
  real(8)              :: CoExp                ! variance of total xs in KL expansion
  real(8)              :: Coscat               ! variance of scat xs in KL expansion, 'material' mode
  real(8)              :: Coabs                ! variance of abs  xs in KL expansion, 'material' mode
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

end module



!--- KL ---!
module KLvars  !"KLresearch" and "KLreconstruct"
  implicit none
  !inputs
  character(3)         :: KLres                ! perform "KL research" operations?
  character(3)         :: KLrec                ! perform "KL reconstruct" operations?
  character(3)         :: KLnoise              ! perform "KL noise" operations?
  integer              :: KLrnumRealz          ! num of realizations to reconstruct
  integer              :: KLrprintat           ! print at this many realizations
  character(8)         :: KLxigentype          ! 'totxs' 'material', xi gen method, allows dissimilar cs?
  character(3)         :: KLvarcalc            ! calculate amount of variance kept in eigs? 'yes' 'no
  real(8)              :: KLvarkept_tol        ! tolerance of rel eig size for var calcs
  integer              :: binNumof             ! number of bins (for xi?)
  integer              :: numEigs              ! number of eigenmodes to calculate
  integer              :: numSlice             ! number of points to plot eigenfunction at
  integer              :: levsrefEig           !
  real(8)              :: binSmallBound        ! smallest xi value for xi bins
  real(8)              :: binLargeBound        ! largest xi value for xi bins
  character(3)         :: KLadjust             ! flag, perform or not
  real(8)              :: meanadjust_tol       ! tolerance for new mean adjustment


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
  integer              :: KLrnumpoints(2)      !
  character(7)         :: pltKLrrealz(4)       !
  integer              :: pltKLrrealznumof     !
  integer, allocatable :: pltKLrrealzwhich(:,:)!
  character(7),allocatable :: pltKLrrealzPointorXi(:) ! which type of reconstruction, pt or xi based?
  !non-inputs
  real(8), allocatable :: gam(:)               ! solutions to eigenvalue transcendental
  real(8), allocatable :: alpha(:)             ! other form of gam
  real(8), allocatable :: Ak(:)                ! normalization coefficients in KL expansion
  real(8), allocatable :: Eig(:)               ! eigenvalues ok KL expansion
  real(8), allocatable :: xi(:,:)              ! array of chosen xi values for reusing reconstructions
  real(8)              :: meanadjust = 0d0     ! positive translation of mean xs

  integer              :: mostinBin            !
  integer              :: Corrnumpoints        !
  real(8), allocatable :: binPDF(:,:)          ! 
  real(8), allocatable :: binBounds(:)         !
  real(8), allocatable :: AllEig(:)            ! many more Eigs than used, for variance kept calcs
  real(8), allocatable :: Allgam(:)            ! many more gam vals than used, for variance kept calcs
  real(8), allocatable :: varmain(:)           ! variance maintained at certain eigenmode
  real(8)              :: binSize              !
  integer              :: negcnt               ! tally of number of realizations that have gone neg
  real(8), allocatable :: KLrx(:)              !
  real(8), allocatable :: KLrxi(:)             !
  real(8), allocatable :: KLrxivals(:,:)       !
  real(8), allocatable :: pltKLrrealzarray(:,:)!
  real(8), allocatable :: KLrrandarray(:,:,:)  !
  real(8), allocatable :: KLrsig(:)            !
  real(8), allocatable :: KLrxisig(:)          !

end module KLvars



!--- Transport ---!
module MCvars
  implicit none
  !inputs
  character(8)         :: probtype             ! 'material' or 'coefs', what is random?
  character(3)         :: radMC                ! perform TMC on binary mixtures operations?
  character(3)         :: radWood              ! perform WMC on binary mixtures operations?
  character(3)         :: KLWood               ! perform WMC on KL reconstructions operations?
  character(3)         :: LPMC                 ! perform MC in the LP sense for binary mixtures?
  character(3)         :: atmixMC              ! perform TMC over atomic mix of binary mixtures?
  character(3)         :: WAMC                 ! perform weight adjust MC over KL resonstructions?
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
  character(7)         :: radMCbinplot         ! 'plot', 'noplot', 'preview', pdf of leakage
  character(7)         :: radWoodbinplot       ! 'plot', 'noplot', 'preview', pdf of leakage
  character(7)         :: KLWoodbinplot        ! 'plot', 'noplot', 'preview', pdf of leakage

  character(3)         :: allowneg             ! allow tranport on neg xs? 'yes', 'no'
  character(3)         :: distneg              ! allow on the fly smoothing of negs? 'yes', 'no'

  integer              :: refsigMode           ! 1-3, 1) userinp, 2) adaptive, 3) adaptive woodcock
  real(8)              :: userrefsig           ! option 7, manually set refsig value as this
  real(8)              :: wgtmax               ! limit max weights tallied with WAMC method
  real(8)              :: wgtmin               ! limit min weights tallied with WAMC method
  character(3)         :: wgtmaxmin            ! chop wgt tallies off at max and min values?
  integer              :: negwgtbinnum         ! number of bins for adaptive negwgts
  integer              :: nwvalsperbin         ! number of values in bin (and one edge value)

  !non inputs
  integer, parameter   :: numPosMCmeths = 6    ! total number of MC transport methods available

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
  real(8), allocatable :: stocMC_fluxall(:,:,:) ! flux of MC meths in cells, (fluxnumbins,nummeths,2)
  real(8), allocatable :: stocMC_fluxmat1(:,:,:)! flux of MC meths in mat1, '2' is mean and var
  real(8), allocatable :: stocMC_fluxmat2(:,:,:)! flux of MC meths in mat2

  real(8), allocatable :: reflect(:)           ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: transmit(:)          ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: absorb(:)            ! slab absorption tally, TMC on binary mixtures

  real(8), allocatable :: LPamMCsums(:)        ! for LP or atomic mix, sum of tallies, refl, tran, abs

                                               ! stoc means mean and variance accross stochastic domain
                                               ! MC means for MC transport solves in spatial domain
  real(8), allocatable :: stocMC_reflection(:,:) ! reflection of MC meths, mean & var in stoc space
  real(8), allocatable :: stocMC_transmission(:,:)! trans of MC meths, mean & var in stoc space
  real(8), allocatable :: stocMC_absorption(:,:) ! absorption of MC meths, mean & var in stoc space


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
  integer, allocatable :: MCcaseson(:)         ! reference of on or not on, cases selected or not
  character(7), allocatable :: MCcases(:)      ! library of MC transport cases
  integer              :: nceilbin             ! number of bins in ceiling calcs for WMC

  real(8)              :: weight               ! weight of particle history for WAMC
  real(8)              :: refsig               ! arbitrary reference xs for WAMC distance and weight calcs
  real(8)              :: maxratio             ! for adaptive negwgts, max acceptable increase ratio
  real(8), allocatable :: negwgtsigs(:,:)      ! sigs for ad negwgts, rank 2, 1=sigs, 2=sigt, 3=refsig
end module MCvars



!--- UQ_MLMC   ---!
module MLMCvars
  implicit none
  !inputs
  character(3)         :: detMLMC              ! perform deterministic MLMC on varying coefs? yes or no
  real(8)              :: MLMC_TOL             ! Tolerance to iterate until achieved
  real(8)              :: MLMC_TOLsplit        ! Phi - Tol split param, how much to type I or II err
  real(8)              :: MLMC_failprob        ! alpha - prob of fail to converge for MC in stoc dim
  integer              :: numcellsLevel0       ! number of cells in initial Level of MLMC 
  integer              :: nextLevelFactor      ! factor to increase number of cells in each Level by
  integer              :: bnumMLMCsamps        ! baseline number of samples for new Level

  !non inputs
  integer, parameter   :: numPosMLMCmeths = 1  ! total number of MLMC transport methods available
  integer, allocatable :: MLMCcaseson(:)       ! reference of on or not on, cases selected or not
  character(7), allocatable :: MLMCcases(:)    ! library of MLMC transport cases
  integer, allocatable :: numMLMCcells(:)      ! number of cells in each Level of MLMC
  real(8)              :: C_alpha              ! coefficient based on conv fail prob for MC in stoc dim
  integer, allocatable :: M_optsamps(:,:)      ! optimal # of samps (1-new est/2-old est,Level)

  real(8), allocatable :: uflux(:,:,:)         ! response function, flux here (samp#,Level,cell)
  real(8), allocatable :: Q_ufunctional(:,:)   ! function of u, (samp#,Level)
  real(8), allocatable :: G_ufunctional(:,:)   ! function of u MLMC form, (samp#,Level)
  real(8), allocatable :: Gave(:)              ! ave of funct of u MLMC form for each Level
  real(8), allocatable :: Gvar(:)              ! var of funct of u MLMC form for each Level

end module MLMCvars

