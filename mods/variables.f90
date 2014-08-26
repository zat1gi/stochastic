module timevars
  implicit none
  real(8), allocatable :: time(:)              ! keeps cumulative time data for everything
                                               ! genR,radMC,radWood,KLnoise,KLcol,KLrec,KLWood
  integer, parameter   :: ntime = 7            ! size of 'time' array

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
  character(7)         :: plotmatdxs           !
  integer, allocatable :: pltgenrealzwhich(:)  !
  !non inputs
  real(8)              :: P(2)                 ! probability of mat1 or 2 for binary mixtures
  real(8)              :: lamc                 ! "correlation length" for binary mixtures  
  integer, allocatable :: matType(:)           ! material type in cell
  real(8), allocatable :: matLength(:)         ! material boundaries

  integer              :: largesti             !
  integer              :: numPath(2)           !
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
  integer              :: KLrnumRealz          ! num of realizations to reconstruct
  integer              :: KLrprintat           ! print at this many realizations
  character(3)         :: KLvarcalc            ! calculate amount of variance kept in eigs? 'yes' 'no
  real(8)              :: KLvarkept_tol        ! tolerance of rel eig size for var calcs
  integer              :: binNumof             ! number of bins (for xi?)
  integer              :: numEigs              ! number of eigenmodes to calculate
  integer              :: numSlice             ! 
  integer              :: levsrefEig           !
  real(8)              :: binSmallBound        !
  real(8)              :: binLargeBound        !

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
  real(8)              :: totLength(2)         !
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
  integer              :: numParts             ! number of particles
  integer              :: trannprt             ! how often to print to screen
  integer              :: pfnumcells           ! number of cells for flux profile (?)

  character(6)         :: rodOrplanar          ! transport geometry mode
  character(6)         :: sourceType           ! 'intern' or 'left', distributed or beam source
  character(6)         :: plotflux(2)          ! 
  character(6)         :: results              !
  character(7)         :: pltflux(4)           ! 
  integer              :: trprofile_binnum     ! number of bins for tran and refl profiles
  character(7)         :: radMCbinplot         ! 'plot', 'noplot', 'preview'
  character(7)         :: radWoodbinplot       ! 'plot', 'noplot', 'preview'
  character(7)         :: KLWoodbinplot        ! 'plot', 'noplot', 'preview'

  character(3)         :: allowneg             ! allow tranport on neg xs? 'yes', 'no'
  character(3)         :: distneg              ! allow on the fly smoothing of negs? 'yes', 'no'

  !non inputs
  integer              :: radtrans_int         !

  real(8), allocatable :: initcur(:)           ! initial current into the slab

  real(8), allocatable :: reflect(:)           ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: Woodr(:)             ! generic reflection tally for Woodcock routines
  real(8), allocatable :: radWoodr(:)          ! slab reflection tally, WMC on binary mixtures
  real(8), allocatable :: KLWoodr(:)           ! slab reflection tally, WMC on KL reconstructions

  real(8), allocatable :: transmit(:)          ! slab reflection tally, TMC on binary mixtures
  real(8), allocatable :: Woodt(:)             ! generic transmission tally for Woodcock routines
  real(8), allocatable :: radWoodt(:)          ! slab transmission tally, WMC on binary mixtures
  real(8), allocatable :: KLWoodt(:)           ! slab transmission tally, WMC on KL reconstructions

  real(8), allocatable :: absorb(:)            ! slab absorption tally, TMC on binary mixtures
  real(8), allocatable :: Wooda(:)             ! generic absorption tally for Woodcock routines
  real(8), allocatable :: radWooda(:)          ! slab absorption tally, WMC on binary mixtures
  real(8), allocatable :: KLWooda(:)           ! slab absorption tally, WMC on KL reconstructions

  character(3)         :: Wood                 ! toggle flag between 'rad' and 'KL'

  integer              :: Wood_rej(2)          ! generic Woodcock rejection tally
  integer              :: radWood_rej(2)       ! Woodcock rejection tally, WMC on binary mixtures
  integer              :: KLWood_rej(2)        ! Woodcock rejection tally, WMC on KL reconstructions
  integer              :: numpnSamp(2)         ! tally of positive and negative KL reconstructions (?)
  real(8)              :: areapnSamp(4)        ! tals pos&neg area: tot pos, tot neg, max pos, max neg

end module MCvars

