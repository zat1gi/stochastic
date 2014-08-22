module KLmeanadjust
  implicit none

  real(8) :: meanadjust =0d0           ! positive translation of mean xs
  character(3) :: KLadjust             ! flag, perform or not
  real(8) :: meanadjust_tol            ! tolerance for new mean adjustment

  real(8) :: xl,xr                     ! endpoints to integrate on
  real(8) :: step                      ! size of largest step when searching for next point
  real(8) :: stol = 0.00000001         ! refinestep cross section tolerance

  real(8) :: aveposarea                ! average positive area (after ignore neg)
  real(8) :: avenegarea                ! average negative area (amount ignored)
  real(8) :: perposdomain              ! percent of domain positive (in percent)
  real(8) :: pernegdomain              ! percent of domain negative (in percent)
CONTAINS
  ! print statements in this module use # 500-599

  subroutine KLadjustmean( lamc )
  !This subroutine is the master for setting the value of "meanadjust", which will conserve
  !the mean of the reconstructions after ignoring negative values in transport within the 
  !chosen tolerance
  use genRealzvars, only: s
  use KLvars,       only: alpha, Ak, Eig, numEigs, sigave, KLrnumRealz
  real(8) :: lamc

  integer :: j,adjustiter
  real(8) :: intsigave,areacont,xmid

  !integrate on all as check
  intsigave = 0d0
  do j=1,KLrnumRealz
    intsigave = intsigave + &
                KLrxi_integral(j,lamc,0d0,s)/KLrnumRealz/s
  enddo
  500 format("  Integrator/reconstruction check - sigave: ",f8.5,"  intsigave: ",f8.5,"  relerr: ",es10.2)
  write(*,500) sigave,intsigave,abs(sigave-intsigave)/sigave


  !meanadjust solver
  step = s/(numEigs*32d0) !hard parameter, relative step size
  adjustiter=0
  do
    adjustiter=adjustiter+1
    aveposarea = 0d0
    perposdomain = 0d0
    avenegarea = 0d0
    pernegdomain = 0d0

    print *,"Beginning mean adjustment iteration ",adjustiter
    do j=1,KLrnumRealz
!print *,"adjusting realization: ",j
      xr = KLrxi_point(j,lamc,0d0)
      xl = 0d0
      do
!print *,"about to find next point"
        xr = findnextpoint(j,lamc)
!print *,"found next point: ",xr
        areacont = KLrxi_integral(j,lamc,xl,xr)/KLrnumRealz/s

        xmid = KLrxi_point(j,lamc,(xr+xl)/2d0)
        if(xmid>0d0) then
          aveposarea = aveposarea + areacont
          perposdomain = perposdomain + (xr - xl)/KLrnumRealz/s * 100
        else
          avenegarea = avenegarea + areacont
          pernegdomain = pernegdomain + (xr - xl)/KLrnumRealz/s * 100
        endif

        xl = xr
        if(xr>=s) exit
      enddo !sum within realization
    enddo !loop over realizations

    if(abs(aveposarea-sigave)/sigave<meanadjust_tol) exit
    print *,"avenegarea: ",avenegarea,"  aveposarea: ",aveposarea
    meanadjust = meanadjust + (sigave - aveposarea)!/s
    if(abs(sigave-aveposarea)<meanadjust_tol) exit
    print *,"meanadjust: ",meanadjust
    !if(adjustiter==6) exit
  enddo !loop over calculating adjustment
    
  end subroutine KLadjustmean



  function findnextpoint(j,lamc)
  !This function searches ahead and finds either 1) next time a reconstructed realization
  !changes signs, or 2) the end of the slab
  use genRealzvars, only: s
  use KLvars,       only: alpha, Ak, Eig, numEigs, sigave, KLrnumRealz
  integer :: j
  real(8) :: lamc

  real(8) :: findnextpoint
  real(8) :: curx,oldx,curs,olds !position, then sigma value

  curx = xl
  curs = KLrxi_point(j,lamc,curx)
!print *,"outside: curx: ",curx,"  step: ",step
  do 
    oldx = curx
    olds = curs

    curx = curx + step
    curs = KLrxi_point(j,lamc,curx)
!print *,"inside1: curx: ",curx,"  step: ",step
    if(curs*olds<0d0) then
      curx = refinenextpoint(j,lamc,oldx,curx)
!print *,"inside2: curx: ",curx,"  step: ",step
      exit
    elseif(curx>=s) then
      curx = s
!print *,"inside3: curx: ",curx,"  step: ",step
      exit
    endif
!print *,"inside4: curx: ",curx,"  step: ",step

  enddo
  findnextpoint = curx

  end function findnextpoint



  function refinenextpoint(j,lamc,oldx,curx)
  !This function takes a range and zeroes in on transision in sign of cross section within tolerance
  use genRealzvars, only: s
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave, KLrnumRealz
  integer :: j
  real(8) :: lamc,oldx,curx

  real(8) :: refinenextpoint,stepsign,curs,olds,curstep

  stepsign = -1d0
  curstep = step
  curs = KLrxi_point(j,lamc,curx)
  do
    curstep = curstep/2d0
    oldx = curx
    olds = curs
    curx = curx + curstep*stepsign
    curs = KLrxi_point(j,lamc,curx)

    if(abs(curs)<stol) then
      refinenextpoint = curx
      exit
    endif
    stepsign = stepsign * merge(1d0,-1d0,curs*olds>0d0)
  enddo

  end function refinenextpoint








  function Eigfunc(Ak,alpha,lamc,xpos)
  ! Used in KLrxi_point to call reconstructed realization at given point
  real(8) :: Ak,alpha,lamc,xpos,Eigfunc

  Eigfunc = Ak * ( sin(alpha*xpos) + lamc*alpha*cos(alpha*xpos) )

  end function Eigfunc




  function KLrxi_point(j,lamc,xpos)
  ! Evaluates KL reconstructed realizations at a given point
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave, KLrxivals
  integer :: j
  real(8) :: lamc,xpos
  real(8) :: KLrxi_point

  integer :: curEig
  real(8) :: Eigfterm

  KLrxi_point = sigave + meanadjust
  do curEig=1,numEigs
    Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,xpos)
    KLrxi_point = KLrxi_point + sqrt(Eig(curEig)) * Eigfterm * KLrxivals(j,curEig)
  enddo

  end function KLrxi_point




  function Eigfuncint(Ak,alpha,lamc,xl,xr)
  ! Used in KLrxi_integral to integrate on KL reconstructed realizations
  real(8) :: Ak,alpha,lamc,xl,xr,Eigfuncint

  Eigfuncint = Ak * (         (-cos(alpha*xr)+cos(alpha*xl))/alpha &
                       + lamc*( sin(alpha*xr)-sin(alpha*xl))           )

  end function Eigfuncint



  function KLrxi_integral(j,lamc,xl,xr)
  ! This function integrates on KL reconstructed realizations from xl to xr
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave, KLrxivals
  integer :: j
  real(8) :: lamc,xl,xr
  real(8) :: KLrxi_integral

  integer :: curEig
  real(8) :: Eigfintterm

  KLrxi_integral = (sigave + meanadjust) * (xr - xl)
  do curEig=1,numEigs
    Eigfintterm = Eigfuncint(Ak(curEig),alpha(curEig),lamc,xl,xr)
    KLrxi_integral = KLrxi_integral + sqrt(Eig(curEig)) * Eigfintterm * KLrxivals(j,curEig)
  enddo

  end function KLrxi_integral




end module KLmeanadjust
