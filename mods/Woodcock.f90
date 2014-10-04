module Woodcock
  use mcnp_random
  use timeman
  use radtransMC
  use utilities
  use KLmeanadjust
  implicit none

CONTAINS
  ! print statements in this module use 600-699

  subroutine Woodnegstats
  use genRealzvars, only: numRealz
  use KLvars, only: negcnt
  use MCvars, only: numpnSamp, areapnSamp, distneg, KLWood, allowneg

  real(8) :: pos,neg

  open(unit=100,file="Woodnegstats.out")

  if(KLWood=='yes' .and. allowneg=='yes') then
    if(distneg=='no')  write(100,*) "--Negative smoothing stats, neg smoothing off--"
    if(distneg=='yes') write(100,*) "--Negative smoothing stats, neg smoothing on--"
    600 format("  Neg realz   : ",f8.5,"%, ",i21," /",i21)
    601 format("  Neg samples : ",f8.5,"%, ",i21," /",i21)
    602 format("  Neg area    : ",f8.5,"%")
    603 format("  Ave neg samp: ",f11.4,"   Ave pos samp: ",f11.4)
    604 format("  Max neg samp: ",f11.4,"   Max pos samp: ",f11.4)

    write(100,600) real(negcnt,8)/real(numRealz,8),negcnt,numRealz
    pos = real(numpnSamp(1),8)
    neg = real(numpnSamp(2),8)
    write(100,601) neg/(pos+neg),numpnSamp(2),numpnSamp(1)+numpnSamp(2)
    write(100,602) -areapnSamp(2)/(areapnSamp(1)-areapnSamp(2))
    write(100,603) areapnSamp(2)/numpnSamp(2),areapnSamp(1)/numpnSamp(1)
    write(100,604) areapnSamp(4),areapnSamp(3)
    write(100,*)
  else
    write(100,*)
  endif

  close(unit=100)
  call system("mv Woodnegstats.out texts")
  end subroutine Woodnegstats


end module
