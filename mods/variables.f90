!--- genRealz ---!
module genRealzvars
  implicit none
  !inputs
  real(8) :: Adamscase                         ! load special cases or no
  real(8), dimension(2) :: sig                 ! cross sections
  real(8), dimension(2) :: scatrat             ! scattering ratios
  real(8), dimension(2) :: lam                 ! ave path lengths
  real(8) :: s                                 ! slab thickness
  integer :: numRealz                          ! number of realz to create
end module



!--- KL ---!
module KLvars  !"KLresearch" and "KLreconstruct"
  implicit none
  !inputs
  character(3)         :: KLvarcalc            ! calculate amount of variance kept in eigs? 'yes' 'no
  real(8)              :: KLvarkept_tol        ! tolerance of rel eig size for var calcs
  integer              :: binNumof             ! number of bins (for xi?)
  integer              :: numEigs              ! number of eigenmodes to calculate
  integer              :: numSlice             ! 
  integer              :: levsrefEig           !

  integer, allocatable :: pltEigfwhich(:)      ! 
  integer, allocatable :: pltxiBinswhich(:,:)  ! 
  integer, allocatable :: pltCowhich(:,:)      ! 
  integer              :: pltxiBinsnumof       ! 
  integer              :: pltEigfnumof         !
  integer              :: pltConumof           !
  !non-inputs
  real(8), allocatable :: gam(:)               ! solutions to eigenvalue transcendental
  real(8), allocatable :: alpha(:)             ! other form of gam
  real(8), allocatable :: Ak(:)                ! normalization coefficients in KL expansion
  real(8), allocatable :: Eig(:)               ! eigenvalues ok KL expansion
  real(8), allocatable :: xi(:,:)              ! array of chosen xi values for reusing reconstructions

  real(8), allocatable :: binPDF(:,:)          ! 
  real(8), allocatable :: binBounds(:)         !
  real(8), allocatable :: AllEig(:)            ! many more Eigs than used, for variance kept calcs
  real(8), allocatable :: Allgam(:)            ! many more gam vals than used, for variance kept calcs
  real(8), allocatable :: varmain(:)           ! variance maintained at certain eigenmode
  integer              :: tnumEigs             ! temp numEigs, local use "KLreval", nice to get rid of

end module KLvars



!--- Transport ---!
module MCvars
  implicit none
  !inputs
  integer :: trprofile_binnum                  ! number of bins for tran and refl profiles
  character(7) :: radMCbinplot                 ! 'plot', 'noplot', 'preview'
  character(7) :: radWoodbinplot               ! 'plot', 'noplot', 'preview'
  character(7) :: KLWoodbinplot                ! 'plot', 'noplot', 'preview'
end module MCvars

