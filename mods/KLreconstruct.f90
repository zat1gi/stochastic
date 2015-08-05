module KLreconstruct
  use utilities
  use mcnp_random
  implicit none

CONTAINS
  ! print statemtns in this module use # 500-599


  subroutine KLreconstructions(icase)
  !Master subroutine for those which create and plot realizations for Markov KL or 
  !Gauss-based KL.  Placed here to declutter multiple instances in 'stochastic.f90'.
  use KLmeanadjust, only: KLadjustmean
  use KLvars, only: flmeanadjust
  integer :: icase

  call KLrmeshgen         !creates mesh for fixed x and xi material constructions
  call KLrgenrealz(icase) !selects array of random variables xi and tests for negativity
  if(flmeanadjust) call KLadjustmean !adjusts mean after lopping neg cross sections
  call KLrplotrealz       !plots reconstructed realizations

  end subroutine KLreconstructions



  subroutine KLrmeshgen
  !This subroutine creates a mesh based on selected frequency of 
  !sampling in x for a fixed point reconstruction, and then for a fixed xi
  !construction.  The fixed xi construction is the only one we care about.
  use genRealzvars , only: s
  use KLvars, only: KLrnumpoints, KLrx, KLrxi, pltKLrrealzPointorXi

  integer :: i
  real(8) :: KLrxstepsize

  if(pltKLrrealzPointorXi(1)=='fpoint') then
    if(allocated(KLrx)) deallocate(KLrx)
    allocate(KLrx(KLrnumpoints(1)))
    KLrx = 0                             !create mesh for fixed point KL reconstruction
    KLrxstepsize = s / KLrnumpoints(1)
    do i=1,KLrnumpoints(1)
      KLrx(i) = KLrxstepsize*i - KLrxstepsize/2
    enddo
  endif

  if(pltKLrrealzPointorXi(1)=='fxi') then
    if(allocated(KLrxi)) deallocate(KLrxi)
    allocate(KLrxi(KLrnumpoints(2)))
    KLrxi = 0                            !create mesh for fixed xi KL reconstruction
    KLrxstepsize = s / KLrnumpoints(2)
    do i=1,KLrnumpoints(2)
      KLrxi(i) = KLrxstepsize*i - KLrxstepsize/2
    enddo
  endif

  end subroutine KLrmeshgen




  subroutine KLrgenrealz(icase)
  !This subroutine reconstructs realizations based upon the KL expansion.
  !It reconstructs based upon the fixed point and fixed xi methods.
  !It also passes an array of selected random variables xi to be plotted in KLreval.
  use rngvars, only: rngappnum, rngstride, setrngappnum
  use timevars, only: time
  use utilities, only: TwoGaussrandnums, erfi
  use genRealzvars, only: s, lamc, sigave, numPosRealz, posRealz
  use KLvars,       only: gam, alpha, Ak, Eig, binPDF, binNumof, numEigs, &
                          KLrnumpoints, KLrnumRealz, KLrprintat, pltKLrrealz, &
                          pltKLrrealznumof, pltKLrrealzwhich, KLrx, KLrxi, KLrxivals, &
                          pltKLrrealzarray, KLrrandarray, KLrsig, KLrxisig, &
                          pltKLrrealzPointorXi, Gaussrandtype, flCorrKL
  use MCvars, only: MCcases
  use timeman, only: KL_timeupdate
  use mcnp_random, only: RN_init_particle
  integer :: i,j,curEig,w,u,icase
  real(8) :: KLsigtemp,Eigfterm,xiterm,rand,rand1,tt1,tt2,xiterms(2)
  logical :: flrealzneg
  character(5) :: flKLtype = 'KLrec'

  call cpu_time(tt1)

  write(*,*) "Starting method: ",flKLtype  
  do j=1,KLrnumRealz
    !set random number application
    if(MCcases(icase)=='KLWood' .or. (MCcases(icase)=='GaussKL' .and. flCorrKL)) then
      call setrngappnum('KLRealzMarkov')
    elseif(MCcases(icase)=='GaussKL') then
      call setrngappnum('KLRealzGaussB')
    endif
    !set random number based on application
    call RN_init_particle( int(rngappnum*rngstride+j,8) )

    if(pltKLrrealzPointorXi(1)=='fpoint') then !create a realization, fixed point
      KLrsig = 0
      do i=1,KLrnumpoints(1)
        !This process not collapsed because KLrrandarray is needed, and does not fit
        !the form of the function.
        KLrsig(i) = sigave
        do curEig=1,numEigs
          Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,KLrx(i))
          rand = rang()
          do u=1,pltKLrrealznumof   !capture rand if useful to plot later
            if( pltKLrrealzwhich(1,u)==j ) then
              KLrrandarray(i,curEig,u+1) = rand
            endif
          enddo
          call select_from_PDF( binPDF,binNumof,numEigs,xiterm,rand )
          KLrsig(i) = KLrsig(i) + sqrt(Eig(curEig)) * Eigfterm * xiterm
        enddo
      enddo
      612 format("  ",f14.8)     !print sigma values to text file, fixed point
      open(unit=10,file="KLrsig.txt")
      do i=1,KLrnumpoints(1)
        write(10,612,advance="no") KLrsig(i)
      enddo
      write(10,*)
    endif

    if(pltKLrrealzPointorXi(1)=='fxi') then !create a realization, fixed xi
      KLrxisig = 0
      do curEig=1,numEigs + mod(numEigs,2)  !select xi values
          rand = rang()
          if((MCcases(icase)=='KLWood' .or. MCcases(icase)=='WAMC') .and. curEig<=numEigs) then
            call select_from_PDF( binPDF,binNumof,curEig,xiterm,rand )
          elseif(MCcases(icase)=='GaussKL' .and. Gaussrandtype=='BM') then
            if(mod(curEig,2)==1) rand1 = rand
            if(mod(curEig,2)==0) then
              call TwoGaussrandnums(rand1,rand,xiterms)
              KLrxivals(j,curEig-1) = xiterms(1)
              xiterm = xiterms(2)
            endif
          elseif(MCcases(icase)=='GaussKL' .and. Gaussrandtype=='inv') then
            xiterm = sqrt(2.0d0)*erfi(2.0d0*rand-1.0d0)
          endif
        if(curEig<=numEigs) KLrxivals(j,curEig) = xiterm
      enddo

      flrealzneg=.false.
      call KLr_negsearch( j,flrealzneg )
      if(.not.flrealzneg) then  !counts the number of realz that contain a positive value
        numPosRealz=numPosRealz+1
        posRealz(j) = 1
        print *,"numNegRealz  : ",j-numPosRealz,"          realz#: ",j
      endif

      do i=1,KLrnumpoints(2)  !create realization
        KLrxisig(i) = KLrxi_point(j,KLrxi(i),chxstype='total')
      enddo
      open(unit=11,file="KLrxisig.txt") !print sigma values to text file, fixed xi
      do i=1,KLrnumpoints(2)
        write(11,612,advance="no") KLrxisig(i)
      enddo
      write(11,*)
    endif


    if(mod(j,KLrprintat)==0) call KL_timeupdate( j,tt1,flKLtype )
  enddo

  end subroutine KLrgenrealz




  subroutine KLrplotrealz
  !This subroutine uses the stored array of pseudo-random numbers used in KLrgenrealz
  !to plot the selected reconstructed realizations.
  use genRealzvars, only: lamc, sigave, numRealz, numPosRealz
  use KLvars,      only: gam, alpha, Ak, Eig, binPDF, binNumof, numEigs, &
                         KLrnumpoints, pltKLrrealz, pltKLrrealznumof, &
                         pltKLrrealzwhich, KLrx, KLrxi, pltKLrrealzarray, KLrrandarray, &
                         KLrsig, KLrxisig, pltKLrrealzPointorXi

  integer :: i,curEig,m,KLrnumpts,tnumEigs
  real(8) :: KLsigtemp,Eigfterm,xiterm,rand

  call system("mv KLrsig.txt plots/KLsigvals")
  call system("mv KLrxisig.txt plots/KLsigvals")

  if( pltKLrrealz(1) .NE. 'noplot' ) then  !plot using generic plotter
    do m=1,pltKLrrealznumof
      tnumEigs=pltKLrrealzwhich(2,m)


      if( pltKLrrealzPointorXi(m) .EQ. 'fpoint' ) then  !create a realz, fixed point
        KLrnumpts=KLrnumpoints(1)
        KLrsig = 0
        do i=1,KLrnumpoints(1)
          !This process not collapsed because KLrrandarray is needed, and does not fit
          !the form of the function.
          KLrsig(i) = sigave
          do curEig=1,tnumEigs
            Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,KLrx(i))
            rand = KLrrandarray(i,curEig,m+1)
            call select_from_PDF( binPDF,binNumof,numEigs,xiterm,rand )
            KLrsig(i) = KLrsig(i) + sqrt(Eig(curEig)) * Eigfterm * xiterm
          enddo
          pltKLrrealzarray(i,1)   = KLrx(i)    !record x values
          pltKLrrealzarray(i,m+1) = KLrsig(i)  !record that realization
        enddo
      endif



      if( pltKLrrealzPointorXi(m) .EQ. 'fxi' ) then  !create a realz, fixed xi
        KLrnumpts=KLrnumpoints(2)
        KLrxisig = 0
        do i=1,KLrnumpoints(2)
          KLrxisig(i) = KLrxi_point(pltKLrrealzwhich(1,m),KLrxi(i),chxstype='total',tnumEigsin=tnumEigs)
          pltKLrrealzarray(i,1)   = KLrxi(i)     !record x values
          pltKLrrealzarray(i,m+1) = KLrxisig(i)  !record that realization
        enddo
      endif
    enddo



    call generic_plotter( KLrnumpts,pltKLrrealznumof,pltKLrrealzarray,&
                          pltKLrrealz )

    call system("mv genericplot.txt plots/KLrrealzplot/KLrrealzplot.txt")
    call system("mv genericplot.ps  plots/KLrrealzplot/KLrrealzplot.ps")
    call system("mv genericplot.pdf plots/KLrrealzplot/KLrrealzplot.pdf")
  endif

  print *," Total num reconstructed realz w/ neg value: ",numRealz-numPosRealz
  print *,

  end subroutine KLrplotrealz





  subroutine KLr_negsearch( j,flrealzneg )
  !This subroutine searches for negative values in KL realizations.
  !If negative value found, set flrealzneg=.true., otherwise remain .false..
  use genRealzvars, only: s
  use KLvars, only: alpha, Ak, Eig, numEigs
  integer :: j
  logical :: flrealzneg

  integer :: i,k,l
  integer :: nminnersteps = 12
  integer :: nmrefine     = 7
  real(8) :: innerstep,outerstep,refinestep
  real(8) :: minsig,minpos,xsig,xpos,minpos_o,minsig_o

  outerstep  = s/numEigs
  innerstep  = s/numEigs/nminnersteps
  refinestep =innerstep*0.24d0

  do i=1,numEigs
    minpos=(outerstep*(i-1))
    minsig=KLrxi_point(j,minpos,chxstype='total')
    do k=2,nminnersteps
      xpos=(outerstep*(i-1)+innerstep*(k-1))
      xsig= KLrxi_point(j,xpos,chxstype='total')
      if(xsig<minsig) then
        minsig=xsig
        minpos=xpos
      endif
    enddo

    do k=1,nmrefine
    minpos_o=minpos
    minsig_o=minsig
      do l=1,5
        xpos=minpos_o-2*refinestep+((k-1)*refinestep)
        if(xpos<0) xpos=0.0d0
        if(xpos>s) xpos=s
        xsig= KLrxi_point(j,xpos,chxstype='total')
        if(xsig<minsig) then
          minsig=xsig
          minpos=xpos
        endif
      enddo
    enddo
    if(minsig<0) then
      flrealzneg=.true.
print *,"minpos",minpos,"minsig",minsig
      exit
    endif

  enddo

  end subroutine KLr_negsearch



  function KLrxi_integral(j,xl,xr)
  !This function integrates on KL reconstructed realizations from xl to xr.
  !Integration is always on total cross section, but 'material'- or 'totxs'-
  !cross sections are adjusted with or without meanadjust.
  use genRealzvars, only: lamc, sigave, Coscat, Coabs
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrxivals, meanadjust, flmatbasedxs

  integer :: j
  real(8) :: xl,xr
  real(8) :: KLrxi_integral

  integer :: curEig
  real(8) :: Eigfintterm, Coterm

  if(flmatbasedxs) then
    Coterm = (sqrt(Coscat)+sqrt(Coabs))**2
  elseif(.not.flmatbasedxs) then
    Coterm = 1.0d0
  endif

  KLrxi_integral = 0d0
  do curEig=1,numEigs
    Eigfintterm = Eigfuncint(Ak(curEig),alpha(curEig),lamc,xl,xr)
    KLrxi_integral = KLrxi_integral + sqrt(Eig(curEig)) * &
                                           Eigfintterm * KLrxivals(j,curEig)
  enddo
  KLrxi_integral = (sigave + meanadjust) * (xr - xl) + (sqrt(Coterm) * KLrxi_integral)

  end function KLrxi_integral



  function Eigfuncint(Ak,alpha,lamc,xl,xr)
  ! Used in KLrxi_integral to integrate on KL reconstructed realizations
  real(8) :: Ak,alpha,lamc,xl,xr,Eigfuncint

  Eigfuncint = Ak * (         (-cos(alpha*xr)+cos(alpha*xl))/alpha &
                       + lamc*( sin(alpha*xr)-sin(alpha*xl))           )
  end function Eigfuncint



  recursive function KLrxi_point(j,xpos,chxstype,tnumEigsin) result(KL_point)
  !Evaluates KL reconstructed realizations at a given point.
  !It has options for total, scattering only, or absorption only cross sectional values.
  !It has an option to solve for less than the available number of eigenvalues.
  !It can function when in 'material'-based or 'totxs'-based mode.
  !It can function when adjusting mean or not adjusting mean.
  !In some of these settings this function calls itself.
  use genRealzvars, only: lamc, P, sig, scatrat, Coscat, Coabs, sigscatave, sigabsave, &
                          sigave, CoExp, scatrat
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrxivals, meanadjust, flmatbasedxs, flmeanadjust

  integer :: j
  real(8) :: xpos
  real(8) :: KL_point
  character(*) :: chxstype
  integer, optional :: tnumEigsin

  integer :: curEig,tnumEigs
  real(8) :: Eigfterm, Coterm, avesigval
  real(8) :: totxsmean, scatxsnomean, totxsnomean, absxsnomean
  logical :: flsolve !solve do-loop (or iteratively use this routine)

  flsolve=.true.
  tnumEigs = merge(tnumEigsin,numEigs,present(tnumEigsin))

  select case (chxstype)
    case ("total")
      if(flmatbasedxs) then
        avesigval = sigave
        Coterm    = (sqrt(Coscat)+sqrt(Coabs))**2
      elseif(.not.flmatbasedxs) then
        avesigval = sigave
        Coterm    = 1d0
      endif
    case ("scatter")
      if(flmatbasedxs) then
        if(.not.flmeanadjust) then
          avesigval = sigscatave
          Coterm    = Coscat
        elseif(flmeanadjust) then
          totxsmean = KLrxi_point(j,xpos,chxstype='total',tnumEigsin=tnumEigs)
          flmeanadjust = .false. !temporary turn off to get scattering ratio
          scatxsnomean = KLrxi_point(j,xpos,chxstype='scatter',tnumEigsin=tnumEigs)
          totxsnomean = KLrxi_point(j,xpos,chxstype='total',tnumEigsin=tnumEigs)
          flmeanadjust = .true.
          KL_point = totxsmean * (scatxsnomean / totxsnomean)
          flsolve = .false.
        endif
      elseif(.not.flmatbasedxs) then
        totxsnomean = KLrxi_point(j,xpos,chxstype='total',tnumEigsin=tnumEigs)
        KL_point = totxsnomean * scatrat(1)
        flsolve = .false.
      endif
    case ("absorb")
      if(flmatbasedxs) then
        if(.not.flmeanadjust) then
          avesigval = sigabsave
          Coterm    = Coabs
        elseif(flmeanadjust) then
          totxsmean = KLrxi_point(j,xpos,chxstype='total',tnumEigsin=tnumEigs)
          flmeanadjust = .false. !temporary turn off to get scattering ratio
          absxsnomean = KLrxi_point(j,xpos,chxstype='absorb',tnumEigsin=tnumEigs)
          totxsnomean = KLrxi_point(j,xpos,chxstype='total',tnumEigsin=tnumEigs)
          flmeanadjust = .true.
          KL_point = totxsmean * (absxsnomean / totxsnomean)
          flsolve = .false.
        endif
      elseif(.not.flmatbasedxs) then
        totxsnomean = KLrxi_point(j,xpos,chxstype='total',tnumEigsin=tnumEigs)
        KL_point = totxsnomean * scatrat(1)
        flsolve = .false.
      endif
  end select

  if(flsolve) then
    KL_point = 0d0
    do curEig=1,tnumEigs
      Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,xpos)
      KL_point = KL_point + sqrt(Eig(curEig)) * Eigfterm * KLrxivals(j,curEig)
    enddo
    KL_point = (avesigval + meanadjust) + (sqrt(Coterm) * KL_point)
  endif

  end function KLrxi_point




  function Eigfunc(Ak,alpha,lamc,xpos)
  ! Used in KLrxi_point to call reconstructed realization at given point
  real(8) :: Ak,alpha,lamc,xpos,Eigfunc

  Eigfunc = Ak * ( sin(alpha*xpos) + lamc*alpha*cos(alpha*xpos) )

  end function Eigfunc


end module KLreconstruct










module KLmeanadjust
  implicit none

  real(8) :: xl,xr                     ! endpoints to integrate on
  real(8) :: step                      ! size of largest step when searching for next point
  real(8) :: stol = 0.00000001         ! refinestep cross section tolerance

  real(8) :: aveposarea                ! average positive area (after ignore neg)
  real(8) :: avenegarea                ! average negative area (amount ignored)
  real(8) :: perposdomain              ! percent of domain positive (in percent)
  real(8) :: pernegdomain              ! percent of domain negative (in percent)
CONTAINS
  ! print statements in this module use # 500-599

  subroutine KLadjustmean
  !This subroutine is the master for setting the value of "meanadjust", which will conserve
  !the mean of the reconstructions after ignoring negative values in transport within the 
  !chosen tolerance
  use genRealzvars, only: s, sigave
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrnumRealz, meanadjust, meanadjust_tol
  use KLreconstruct, only: KLrxi_point, KLrxi_integral

  integer :: j,adjustiter
  real(8) :: intsigave,areacont,xmid

  !initialize meanadjust
  meanadjust = 0.0d0

  !integrate on all as check
  intsigave = 0d0
  do j=1,KLrnumRealz
    intsigave = intsigave + &
                KLrxi_integral(j,0d0,s)/KLrnumRealz/s
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
      xr = KLrxi_point(j,0d0,chxstype='total')
      xl = 0d0
      do
        !find next point and area between these two
        xr = findnextpoint(j)
        areacont = KLrxi_integral(j,xl,xr)/KLrnumRealz/s

        !calc aveposarea for tol check, also negstat tallies
        xmid = KLrxi_point(j,(xr+xl)/2d0,chxstype='total')
        if(xmid>0d0) then
          aveposarea = aveposarea + areacont
          perposdomain = perposdomain + (xr - xl)/KLrnumRealz/s * 100
        else
          avenegarea = avenegarea + areacont
          pernegdomain = pernegdomain + (xr - xl)/KLrnumRealz/s * 100
        endif

        !advance lagging point
        xl = xr
        if(xr>=s) exit
      enddo !sum within realization
    enddo !loop over realizations

    if(abs(aveposarea-sigave)/sigave<meanadjust_tol) exit
    print *,"avenegarea: ",avenegarea,"  aveposarea: ",aveposarea
    meanadjust = meanadjust + (sigave - aveposarea)
    print *,"meanadjust: ",meanadjust
  enddo !loop over calculating adjustment
    
  end subroutine KLadjustmean



  function findnextpoint(j)
  !This function searches ahead and finds either 1) next time a reconstructed realization
  !changes signs, or 2) the end of the slab
  use genRealzvars, only: s, sigave
  use KLvars,       only: alpha, Ak, Eig, numEigs, KLrnumRealz
  use KLreconstruct, only: KLrxi_point
  integer :: j

  real(8) :: findnextpoint
  real(8) :: curx,oldx,curs,olds !position, then sigma value

  curx = xl
  curs = KLrxi_point(j,curx,chxstype='total')
  do 
    oldx = curx
    olds = curs

    curx = curx + step
    curs = KLrxi_point(j,curx,chxstype='total')
    if(curs*olds<0d0) then
      curx = refinenextpoint(j,oldx,curx)
      exit
    elseif(curx>=s) then
      curx = s
      exit
    endif

  enddo
  findnextpoint = curx

  end function findnextpoint



  function refinenextpoint(j,oldx,curx)
  !This function takes a range and zeroes in on transition in sign of cross section within tolerance
  use genRealzvars, only: s, sigave
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrnumRealz
  use KLreconstruct, only: KLrxi_point
  integer :: j
  real(8) :: oldx,curx

  real(8) :: refinenextpoint,stepsign,curs,olds,curstep

  stepsign = -1d0
  curstep = step
  curs = KLrxi_point(j,curx,chxstype='total')
  do
    curstep = curstep/2d0
    oldx = curx
    olds = curs
    curx = curx + curstep*stepsign
    curs = KLrxi_point(j,curx,chxstype='total')

    if(abs(curs)<stol) then
      refinenextpoint = curx
      exit
    endif
    stepsign = stepsign * merge(1d0,-1d0,curs*olds>0d0)
  enddo

  end function refinenextpoint


end module KLmeanadjust
