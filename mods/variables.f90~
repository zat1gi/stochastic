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
  character(3) :: KLvarcalc                    ! calculate amount of variance kept in eigs? 'yes' 'no
  real(8) :: KLvarkept_tol                     ! tolerance of rel eig size for var calcs
  !non-inputs
  real(8), allocatable :: AllEig(:)            ! many more Eigs than used, for variance kept calcs
  real(8), allocatable :: Allgam(:)            ! many more gam vals than used, for variance kept calcs
  real(8), allocatable :: varmain(:)           ! variance maintained at certain eigenmode
end module KLvars



!--- Transport ---!
module transport_sharedvars
  implicit none
  !inputs
  integer :: trprofile_binnum                  ! number of bins for tran and refl profiles
end module transport_sharedvars

module radMCvars !Traditional MC on original Markov realizations
  implicit none
  !inputs
  character(7) :: radMCbinplot                 ! 'plot', 'noplot', 'preview'
end module radMCvars

module radWoodvars !Woodcock sampling on original Markov realizations
  implicit none
  !inputs
  character(7) :: radWoodbinplot               ! 'plot', 'noplot', 'preview'
end module radWoodvars

module KLWoodvars  !Woodcock sampling on KL reconstructed realizations
  implicit none
  !inputs
  character(7) :: KLWoodbinplot                ! 'plot', 'noplot', 'preview'
end module KLWoodvars
