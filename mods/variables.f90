module timevars
  implicit none
  real(8), allocatable :: time(:)              ! keeps cumulative time data for everything
                                               ! genR,radMC,radWood,KLnoise,KLcol,KLrec,KLWood,LPMC,atmixMC
  integer, parameter   :: ntime = 9            ! size of 'time' array
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
  real(8)              :: P(2)                 ! probability of mat1 or 2 for binary mixtures
  real(8)              :: lamc                 ! "correlation length" for binary mixtures  
  integer, allocatable :: matType(:)           ! material type in cell
  real(8), allocatable :: matLength(:)         ! material boundaries
  real(8)              :: atmixsig             ! atomically mixed cross section value
  real(8)              :: atmixscatrat         ! atomically mixed scattering ratio

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
  character(3)         :: KLvarcalc            ! calculate amount of variance kept in eigs? 'yes' 'no
  real(8)              :: KLvarkept_tol        ! tolerance of rel eig size for var calcs
  integer              :: binNumof             ! number of bins (for xi?)
  integer              :: numEigs              ! number of eigenmodes to calculate
  integer              :: numSlice             ! number of points to plot eigenfunction at
  integer              :: levsrefEig           !
  real(8)              :: binSmallBound        ! smallest xi value for xi bins
  real(8)              :: binLargeBound        ! largest xi value for xi bins

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
  real(8)              :: sigave               ! weighted average sigma value for binary mixtures

  integer              :: mostinBin            !
  integer              :: Corrnumpoints        !
  real(8), allocatable :: binPDF(:,:)          ! 
  real(8), allocatable :: binBounds(:)         !
  real(8), allocatable :: AllEig(:)            ! many more Eigs than used, for variance kept calcs
  real(8), allocatable :: Allgam(:)            ! many more gam vals than used, for variance kept calcs
  real(8), allocatable :: varmain(:)           ! variance maintained at certain eigenmode
  integer              :: tnumEigs             ! temp numEigs, local use "KLreval", nice to get rid of
  real(8)              :: binSize              !
  real(8)              :: CoExp                !
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
  character(3)         :: radMC                ! perform TMC on binary mixtures operations?
  character(3)         :: radWood              ! perform WMC on binary mixtures operations?
  character(3)         :: KLWood               ! perform WMC on KL reconstructions operations?
  character(3)         :: LPMC                 ! perform MC in the LP sense for binary mixtures?
  character(3)         :: atmixMC              ! perform TMC over atomic mix of binary mixtures?
  integer              :: numParts             ! number of particles
  integer              :: LPamnumParts         ! number of particles for LP or atomic mix
  integer              :: trannprt             ! how often to print to screen
  integer              :: fluxnumcells         ! number of cells for flux profile
  integer              :: pfnumcells           ! number of cells for flux profile (?) !get rid of me

  character(6)         :: rodOrplanar          ! transport geometry mode
  character(6)         :: sourceType           ! 'intern' or 'left', distributed or beam source
  character(6)         :: plotflux(2)          ! 
  character(6)         :: results              !
  character(7)         :: pltflux(4)           ! plot material irrespective flux
  character(7)         :: pltmatflux           ! plot material respective fluxes
  character(7)         :: pltfluxtype          ! full 'track' length flux plot or 'point' flux binning? 
  integer              :: trprofile_binnum     ! number of bins for tran and refl profiles
  character(7)         :: radMCbinplot         ! 'plot', 'noplot', 'preview'
  character(7)         :: radWoodbinplot       ! 'plot', 'noplot', 'preview'
  character(7)         :: KLWoodbinplot        ! 'plot', 'noplot', 'preview'

  character(3)         :: allowneg             ! allow tranport on neg xs? 'yes', 'no'
  character(3)         :: distneg              ! allow on the fly smoothing of negs? 'yes', 'no'

  !non inputs
  integer, parameter   :: numPosMCmeths = 5    ! total number of MC transport methods available

  real(8), allocatable :: ABreflection(:,:)    ! Adams/Brantley Reflection Values
  real(8), allocatable :: ABtransmission(:,:)  ! Adams/Brantley Transmission Values
                                               ! First rank: 1 average, 2 stdev
                                               ! Second rank: 1 AdMC, 2 AdLP, 3 BrMC, 4 BrLP, 5 Bratmix

  integer              :: radtrans_int         !

  real(8), allocatable :: fluxfaces(:)         ! mesh for flux tallies
  real(8), allocatable :: flux(:,:)            ! flux tally both mats, TMC over binary mixtures !rid of me!replace me with "fluxo", and change it to just "flux"!!!
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


  real(8), allocatable :: Woodf(:,:)           ! flux tally both mats, generic Woodcock
  real(8), allocatable :: radWoodf(:,:)        ! flux tally both mats, WMC over binary mixtures
  real(8), allocatable :: KLWoodf(:,:)         ! flux tally both mats, WMC over KL reconstruction
  real(8), allocatable :: fWoodf(:,:)          ! flux tally in 1st mat, generic Woodcock
  real(8), allocatable :: fradWoodf(:,:)       ! flux tally in 1st mat, WMC over binary mixtures
  real(8), allocatable :: fKLWoodf(:,:)        ! flux tally in 1st mat, WMC over KL reconstructions
  real(8), allocatable :: bWoodf(:,:)          ! flux tally in 2nd mat, generic Woodcock
  real(8), allocatable :: bradWoodf(:,:)       ! flux tally in 2nd mat, WMC over binary mixtures
  real(8), allocatable :: bKLWoodf(:,:)        ! flux tally in 2nd mat, WMC over KL reconstructions

  real(8), allocatable :: initcur(:) !here for politics only, finish new MC driver, get rid of me

  real(8), allocatable :: reflect(:)           ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: Woodr(:)             ! generic reflection tally for Woodcock routines
  real(8), allocatable :: radWoodr(:)          ! slab reflection tally, WMC on binary mixtures
  real(8), allocatable :: KLWoodr(:)           ! slab reflection tally, WMC on KL reconstructions

  character(3)         :: Wood                 ! toggle flag between 'rad' and 'KL'

  real(8), allocatable :: transmit(:)          ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: Woodt(:)             ! generic transmission tally for Woodcock routines
  real(8), allocatable :: radWoodt(:)          ! slab transmission tally, WMC on binary mixtures
  real(8), allocatable :: KLWoodt(:)           ! slab transmission tally, WMC on KL reconstructions

  real(8), allocatable :: absorb(:)            ! slab absorption tally, TMC on binary mixtures
  real(8), allocatable :: Wooda(:)             ! generic absorption tally for Woodcock routines
  real(8), allocatable :: radWooda(:)          ! slab absorption tally, WMC on binary mixtures
  real(8), allocatable :: KLWooda(:)           ! slab absorption tally, WMC on KL reconstructions

  real(8), allocatable :: LPamMCsums(:)        ! for LP or atomic mix, sum of tallies, refl, tran, abs

                                               ! stoc means mean and variance accross stochastic domain
                                               ! MC means for MC transport solves in spatial domain
  real(8), allocatable :: stocMC_reflection(:,:) ! reflection of MC meths, mean & var in stoc space
  real(8), allocatable :: stocMC_transmission(:,:)! trans of MC meths, mean & var in stoc space
  real(8), allocatable :: stocMC_absorption(:,:) ! absorption of MC meths, mean & var in stoc space


  integer              :: Wood_rej(2)          ! generic Woodcock rejection tally
  integer              :: radWood_rej(2)       ! Woodcock rejection tally, WMC on binary mixtures
  integer              :: KLWood_rej(2)        ! Woodcock rejection tally, WMC on KL reconstructions
  integer              :: numpnSamp(2)         ! tally of positive and negative KL reconstructions (?)
  real(8)              :: areapnSamp(4)        ! tals pos&neg area: tot pos, tot neg, max pos, max neg

  real(8),allocatable  :: binmaxind(:)         ! These are used in ceilings for WMC
  real(8),allocatable  :: binmaxes(:)          ! 
  real(8),allocatable  :: fbinmax(:)           ! 
  real(8),allocatable  :: bbinmax(:)           ! 
  real(8)              :: disthold             ! used in a workaround for neg vals in KLWood

  real(8)              :: position             ! 'current' position of particle
  real(8)              :: oldposition          ! most recent position of particle
  real(8)              :: mu                   ! 'current' direction of particle
  integer, allocatable :: MCcaseson(:)         ! reference of on or not on, cases selected or not
  character(7), allocatable :: MCcases(:)      ! library of MC transport cases
  integer              :: nceilbin             ! number of bins in ceiling calcs for WMC

end module MCvars

