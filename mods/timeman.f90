module timeman
  implicit none

CONTAINS
  ! print statements in this module use 1100-1200


!! time tracking funcs and subs

  subroutine initialize_t1
  !sets t1, the beginning of a new repetative operation
  use timevars, only: t1
  call cpu_time(t1)
  end subroutine initialize_t1



  subroutine timeupdate( chtype,ndone,ntotal )
  !uses t1 and new t2 to give update on time for new repetative operation
  use genRealzvars, only: numRealz
  use MCvars, only: numPartsperj
  use timevars, only: t1
  integer :: ndone,ntotal,      orda,ordm
  real(8) :: t2,                ave,maxv,eps=0.0000000001d0
  character(*) :: chtype

  call cpu_time(t2)

  if(chtype=='radMC' .or. chtype=='radWood' .or. chtype=='KLWood' .or. chtype=='GaussKL') then
    1101 format(A9,"   ",f6.1,"% of ",i7,"   time/est:",f7.2,"/",f7.2," min   h/r-ave/max: ",f3.1,"E",i1,"/",f3.1,"E",i1)
    ave = real(sum(numPartsperj),8)/numRealz+eps
    orda= floor(log(ave)/log(10d0))
    maxv= real(maxval(numPartsperj),8)+eps 
    ordm= floor(log(maxv)/log(10d0))
    write(*,1101) chtype,real(ndone,8)/ntotal*100d0,ntotal,(t2-t1)/60.0d0,(t2-t1)*ntotal/ndone/60.0d0, &
                  ave/(10.0**orda),orda,maxv/(10.0**ordm),ordm
  else
    1100 format(A9,"   ",f6.1,"% of ",i7,"   time/est:",f7.2,"/",f7.2," min")
    write(*,1100) chtype,real(ndone,8)/ntotal*100d0,ntotal,(t2-t1)/60.0d0,(t2-t1)*ntotal/ndone/60.0d0
  endif

  end subroutine


end module timeman

