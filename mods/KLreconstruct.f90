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
  !It tests for negative realizations, rejecting and replacing them is specified.
  !It also passes an array of selected random variables xi to be plotted in KLreval.
  use rngvars, only: rngappnum, rngstride, setrngappnum
  use timevars, only: time
  use utilities, only: TwoGaussrandnums, erfi
  use genRealzvars, only: s, lamc, sigave, numPosRealz, numNegRealz
  use KLvars,       only: gam, alpha, Ak, Eig, binPDF, binNumof, numEigs, &
                          KLrnumpoints, KLrnumRealz, KLrprintat, pltKLrrealz, &
                          pltKLrrealznumof, pltKLrrealzwhich, KLrx, KLrxi, KLrxivals, &
                          pltKLrrealzarray, KLrrandarray, KLrsig, KLrxisig, &
                          pltKLrrealzPointorXi, Gaussrandtype, flCorrKL, flmeanadjust
  use MCvars, only: MCcases, flnegxs, KLWood, GaussKL
  use timeman, only: KL_timeupdate
  use mcnp_random, only: RN_init_particle
  integer :: i,tentj,realj,curEig,w,u,icase
  real(8) :: KLsigtemp,Eigfterm,xiterm,rand,rand1,tt1,tt2,xiterms(2)
  logical :: flrealzneg, flacceptrealz
  logical :: flpurpose(3)=.false. !1)neg or not, 2)max vals, 3)zeros

  call cpu_time(tt1)

  write(*,*) "Starting method: KLrec"
  !set purposes for KLr_realzanalysis (purpose,[planned purpos])
  flpurpose(1)=.true.                                       !neg or not (no neg, [neg analysis])
  if(KLWood=='yes' .or. GaussKL=='yes') flpurpose(2)=.true. !max vals (KL WMC)
  if(flmeanadjust) flpurpose(3)=.true.                      !zeros (mean adjust, [neg analysis])

  tentj=0
  realj=1
  do
    tentj=tentj+1
    flacceptrealz=.true.

    !set random number application
    if(MCcases(icase)=='KLWood' .or. (MCcases(icase)=='GaussKL' .and. flCorrKL)) then
      call setrngappnum('KLRealzMarkov')
    elseif(MCcases(icase)=='GaussKL') then
      call setrngappnum('KLRealzGaussB')
    endif
    !set random number based on application
    call RN_init_particle( int(rngappnum*rngstride+tentj,8) )

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
            if( pltKLrrealzwhich(1,u)==realj ) then
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
              KLrxivals(realj,curEig-1) = xiterms(1)
              xiterm = xiterms(2)
            endif
          elseif(MCcases(icase)=='GaussKL' .and. Gaussrandtype=='inv') then
            xiterm = sqrt(2.0d0)*erfi(2.0d0*rand-1.0d0)
          endif
        if(curEig<=numEigs) KLrxivals(realj,curEig) = xiterm
      enddo

      !count num of realz w/ neg xs, set flag to accept or reject realz
      flrealzneg=.false.
!      call KLr_negsearch( realj,flrealzneg )
      call KLr_analyzerealz( realj,flrealzneg,flpurpose )
      if(.not.flrealzneg) numPosRealz=numPosRealz+1
      if(     flrealzneg) then
        if(.not.flnegxs) flacceptrealz=.false.
        numNegRealz=numNegRealz+1
        print *,"numNegRealz  : ",numNegRealz," tentative realz#: ",tentj
      endif

      do i=1,KLrnumpoints(2)  !create realization
        KLrxisig(i) = KLrxi_point(realj,KLrxi(i),chxstype='total')
      enddo
      open(unit=11,file="KLrxisig.txt") !print sigma values to text file, fixed xi
      do i=1,KLrnumpoints(2)
        write(11,612,advance="no") KLrxisig(i)
      enddo
      write(11,*)
    endif

    if(flnegxs) then
      if(mod(realj,KLrprintat)==0)       call KL_timeupdate( realj,tt1,'KLrec' )
      if(KLrnumRealz==realj      ) exit
    else
      if(mod(numPosRealz,KLrprintat)==0) call KL_timeupdate( realj,tt1,'KLrec' )
      if(KLrnumRealz==numPosRealz) exit
    endif
    if(flacceptrealz) realj=realj+1
  enddo

  end subroutine KLrgenrealz




  subroutine KLrplotrealz
  !This subroutine uses the stored array of pseudo-random numbers used in KLrgenrealz
  !to plot the selected reconstructed realizations.
  use genRealzvars, only: lamc, sigave, numRealz, numPosRealz, numNegRealz
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

  print *," Total num reconstructed realz w/ neg value: ",numNegRealz,"/",numNegRealz+numPosRealz
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
    if(minsig<0.0d0) then
      flrealzneg=.true.
      print *,"minpos",minpos,"minsig",minsig
      exit
    endif

  enddo

  end subroutine KLr_negsearch



  subroutine KLr_analyzerealz( j,flrealzneg,flpurpose )
  !This subroutine can 1) discern if realization contains negativity
  !2) produce all local max values for use in WMC
  !3) produce all zeros for use in negativity analysis and mean adjustment
  !It first finds bounds for all extrema, then locates and labels extrema.
  !It then uses labeled extrema as bounds for zeros, which are found and labeled.
  !Routine exits when all chosen purposes accomplished.
  use genRealzvars, only: s, numRealz
  use KLvars, only: KLr_maxima, KLr_zeros, numEigs
  logical :: flpurpose(3) !1)neg or not, 2)max vals, 3)zeros
  logical :: flrealzneg
  integer :: j

  logical :: flaccomplished(3), flfinished
  integer :: i, exi, imax, izer, exsize
  logical, allocatable :: flextremamax(:), flextremapos(:)
  real(8), allocatable :: derivloc(:), derivval(:)
  real(8), allocatable :: extrema(:)
  real(8) :: rtemp, KLpoint1, KLpoint2

  flaccomplished = .false.
  flfinished     = .false.

  call KLr_findextremabounds( j,derivloc,derivval )              !find extrema bounds

  exsize = arithmaticsum(1,numEigs,1,numEigs)
  allocate(extrema(exsize+1))
  allocate(flextremamax(exsize+1))
  allocate(flextremapos(exsize+1))
  extrema = 0.0d0

                                                                 !find extream/negativities
  !label extrema (left boundary)
  flextremamax(1) = merge(.true.,.false.,KLrxi_point(j,extrema(1),'total') > &
                                         KLrxi_point(j,extrema(1)+0.000000000001d0,'total'))
  flextremapos(1) = merge(.true.,.false.,KLrxi_point(j,extrema(1),'total')>0.0d0)

  exi = 2           !index of next extrema to search for
  do i=1,size(derivloc)-1
    if(derivval(i)*derivval(i+1)<0.0d0) then  !points are extrema bounds
      !find extrema
      extrema(exi) = KLr_findzeros( j,derivloc(i),derivloc(i+1),derivval(i),derivval(i+1),flderiv=.true. )
      !label extrema (local max or min)
      flextremamax(exi) = merge(.true.,.false.,KLr_concavity( j,extrema(exi)         )< 0.0d0)
      flextremapos(exi) = merge(.true.,.false.,KLrxi_point(   j,extrema(exi),'total' )>=0.0d0)
      !test if extrema is negative
      if(.not.flextremapos(exi)) then
        flrealzneg = .true.
        if(flpurpose(1)) flaccomplished(1)=.true.       !a negativitiy found
      endif
      !advance extrema index
      exi = exi + 1

      flfinished = flfinishedtest(flpurpose,flaccomplished)
      if(flfinished) exit
    endif
  enddo
  !label extrema (right boundary)
  extrema(exi) = s
  flextremamax(exi) = merge(.true.,.false.,KLrxi_point(j,extrema(exi),'total') > &
                                           KLrxi_point(j,extrema(exi)-0.000000000001d0,'total'))
  flextremapos(exi) = merge(.true.,.false.,KLrxi_point(   j,extrema(exi),'total' )>0.0d0)
  if(flpurpose(1)) flaccomplished(1)=.true.             !any negativities are found
  flfinished = flfinishedtest(flpurpose,flaccomplished)


  if(flpurpose(2) .and. .not.flfinished) then                        !collect maxima
    if(j==1 .and. allocated(KLr_maxima)) deallocate(KLr_maxima)
    if(.not.allocated(KLr_maxima)) then
      allocate(KLr_maxima(numRealz,exsize))
      KLr_maxima = 0.0d0
    endif
    imax = 1
    do i=1,exi
      if(flextremamax(i)) then
        KLr_maxima(j,imax) = extrema(i)
        imax = imax+1
      endif
    enddo
  endif
  if(flpurpose(2)) flaccomplished(2)=.true.             !all maxes found
  flfinished = flfinishedtest(flpurpose,flaccomplished)


  if(.not.flfinished) then                                           !find all zeros
    if(j==1 .and. allocated(KLr_zeros)) deallocate(KLr_zeros)
    if(.not.allocated(KLr_zeros)) then
      allocate(KLr_zeros(numRealz,exsize+numEigs))
       KLr_zeros = 0.0d0
    endif
    izer = 1
    do i=1,exi-1
      !if extrema are bounds of zero then find zero
      if(flextremapos(i) .neqv. flextremapos(i+1)) then
        KLpoint1 = KLrxi_point(j,extrema(i),'total')
        KLpoint2 = KLrxi_point(j,extrema(i+1),'total')
        KLr_zeros(j,izer) = KLr_findzeros( j,extrema(i),extrema(i+1),KLpoint1,KLpoint2,flderiv=.false. )
        izer = izer + 1
      endif
    enddo
    if(flpurpose(3)) flaccomplished(3)=.true.           !all zeros found
  endif

  deallocate(derivloc)
  deallocate(derivval)
  deallocate(extrema)
  deallocate(flextremamax)
  deallocate(flextremapos)


  end subroutine KLr_analyzerealz


  function flfinishedtest( flpurpose,flaccomplished )
  !This function tests if flpurpose and flaccomplished are the same and returns true of false.
  logical :: flpurpose(:),flaccomplished(:)
  logical :: flfinishedtest
  integer :: i

  flfinishedtest = .false.
  do i=1,size(flpurpose)
    if(flpurpose(i) .neqv. flaccomplished(i)) exit
    if(i==size(flpurpose)) flfinishedtest = .true.
  enddo
  end function flfinishedtest




  subroutine KLr_findextremabounds( j,derivloc,derivval )
  !This takes the derivative of the KL expansion at points on a successively refined grid 
  !searching for bounds on all zeros of the derivative of the KL expansion (extrema of KL).
  use genRealzvars, only: s
  use KLvars, only: numEigs, numextremasameiter
  integer :: j
  real(8), allocatable :: derivloc(:), derivval(:), derivnew(:)
  real(8), allocatable :: derivloc_(:),derivval_(:),derivnew_(:)

  integer :: i, realzwochange, arsum, tallychanges, tallychanges_, itersametally
  real(8) :: step

  arsum = arithmaticsum(1,numEigs,1,numEigs)
  allocate(derivloc(arsum))
  allocate(derivval(arsum))
  allocate(derivnew(arsum))
  tallychanges = 0
  tallychanges_= 0
  itersametally= 0

  do
    call move_alloc(derivloc,derivloc_)    !keep old values, but create room for new ones
    call move_alloc(derivval,derivval_)
    call move_alloc(derivnew,derivnew_)
    allocate(derivloc(2*size(derivloc_)-1))
    allocate(derivval(2*size(derivval_)-1))
    allocate(derivnew(2*size(derivnew_)-1))
    derivnew = 0
    do i=1,size(derivloc_)
      derivloc(2*i-1) = derivloc_(i)
      derivval(2*i-1) = derivval_(i)
      derivnew(2*i-1) = derivnew_(i)
    enddo
    if(size(derivnew_)==arsum) derivnew = 0
    deallocate(derivloc_)
    deallocate(derivval_)
    deallocate(derivnew_)

    step = real(s,8)/real(size(derivloc)-1,8)     !setup location values
    derivloc(1) = 0.0d0
    do i=2,size(derivloc)
      derivloc(i) = derivloc(i-1) + step
    enddo

    do i=1,size(derivloc)                         !solve all new needed values
     if(derivnew(i)==0) then
       derivval(i) = KLr_derivative( j,derivloc(i) )
       derivnew(i) = 1
     endif
    enddo
    
    tallychanges = 0                              !count number of ranges found
    do i=1,size(derivloc)-1
      if(derivval(i)*derivval(i+1)<0.0d0) tallychanges = tallychanges + 1
    enddo

    if(tallychanges==tallychanges_) then          !tally iters with same number of changes
      itersametally = itersametally+1
    else
      itersametally = 0
    endif
    tallychanges_ = tallychanges

    if(itersametally>=numextremasameiter) exit                      !exit loop if not changed for several iters
  enddo

  deallocate(derivnew)
  end subroutine KLr_findextremabounds




  function KLr_findzeros( j,loc1,loc2,val1,val2,flderiv )
  !This function accepts bounds on an zero in either the derivative or value of the KL
  !expansion and returns the location of the zero.
  !The algorithm assumes val1 or val2 is negative, the other is positive, and there
  !is only one transition between positive and negative between the two.
  integer :: j
  real(8) :: loc1,loc2
  real(8) :: KLr_findzeros
  logical :: flderiv  !.true. means deriv, .false. means value

  real(8) :: loc3 !loc1 (left), loc2 (right), loc3 (middle)
  real(8) :: val1,val3,val2 !evaluation of loc#, either 'deriv' or 'value'

  do
    !check for valid bounds
    if(val1*val2>0.0d0) then
      print *,"bounds for zero both sign!:",val1,val2
      stop
    endif

    !new location and value
    loc3 = (loc1 + loc2) / 2.0d0
    val3 = merge(KLr_derivative(j,loc3),KLrxi_point(j,loc3,'total'),flderiv)

    !if found zero, exit loop
    if(abs(val3)<0.000000000001d0) exit

    !advance location (and value) towards zero
    if(val1*val3>=0.0d0) then
      loc1=loc3
      val1=val3
    else
      loc2=loc3
    endif
  enddo

  KLr_findzeros = loc3
  end function KLr_findzeros




  function KLr_concavity( j,xloc )
  !This function evaluates the value of the second derivative of a total cross section
  !constructed with the KL expansion at a location in the domain.
  !Positive values mean concave up, negative values mean concave down.
  !Originally added to test whether extrema is maxima or minima.
  use genRealzvars, only: lamc, Coscat, Coabs
  use KLvars, only: numEigs, alpha, Eig, Ak, KLrxivals, flmatbasedxs
  real(8) :: xloc, KLr_concavity
  integer :: j

  integer :: curEig
  real(8) :: Coterm

  if(flmatbasedxs) then
    Coterm = (sqrt(Coscat)+sqrt(Coabs))**2
  elseif(.not.flmatbasedxs) then
    Coterm = 1.0d0
  endif

  KLr_concavity = 0d0
  do curEig=1,numEigs
    KLr_concavity = KLr_concavity   + sqrt(Eig(curEig)) * KLrxivals(j,curEig) * &
                                            Ak(curEig)  *    alpha(curEig)**2 * &
                    (-sin( alpha(curEig)*xloc ) - lamc*alpha(curEig)*cos( alpha(curEig)*xloc ))
  enddo
  KLr_concavity = sqrt(Coterm) * KLr_concavity

  end function KLr_concavity





  function KLr_derivative( j,xloc )
  !This function evaluates the value of the derivative of a total cross section
  !constructed with the KL expansion at a location in the domain.
  use genRealzvars, only: lamc, Coscat, Coabs
  use KLvars, only: numEigs, alpha, Eig, Ak, KLrxivals, flmatbasedxs
  real(8) :: xloc, KLr_derivative
  integer :: j

  integer :: curEig
  real(8) :: Coterm

  if(flmatbasedxs) then
    Coterm = (sqrt(Coscat)+sqrt(Coabs))**2
  elseif(.not.flmatbasedxs) then
    Coterm = 1.0d0
  endif

  KLr_derivative = 0d0
  do curEig=1,numEigs
    KLr_derivative = KLr_derivative + sqrt(Eig(curEig)) * KLrxivals(j,curEig) * &
                                            Ak(curEig)  *       alpha(curEig) * &
                     (cos( alpha(curEig)*xloc ) - lamc*alpha(curEig)*sin( alpha(curEig)*xloc ))
  enddo
  KLr_derivative = sqrt(Coterm) * KLr_derivative

  end function KLr_derivative




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



  function KLrxi_point(j,xpos,chxstype,tnumEigsin) result(KL_point)
  !!Evaluates KL reconstructed realizations at a given point.
  !!It has options for total, scattering only, or absorption only cross sectional values.
  !!It has an option to solve for less than the available number of eigenvalues.
  !!It can function when in 'material'-based or 'totxs'-based mode.
  !!It can function when adjusting mean or not adjusting mean.
  !!In some of these settings this function calls itself.
  use genRealzvars, only: lamc, scatrat, Coscat, Coabs, sigscatave, sigabsave, &
                          sigave
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrxivals, meanadjust, flmatbasedxs, &
                    sigsmeanadjust, sigameanadjust
  use utilities, only: Heavi

  integer :: j
  real(8) :: xpos
  real(8) :: KL_point
  character(*) :: chxstype
  integer, optional :: tnumEigsin

  real(8) :: sigt, siga, sigs, KL_sum
  integer :: curEig,tnumEigs
  real(8) :: Eigfterm

  tnumEigs = merge(tnumEigsin,numEigs,present(tnumEigsin))

  !solve summation of KL terms to use below
  KL_sum = 0d0
  do curEig=1,tnumEigs
    Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,xpos)
    KL_sum   = KL_sum + sqrt(Eig(curEig)) * Eigfterm * KLrxivals(j,curEig)
  enddo

  if(.not.flmatbasedxs) then
    !determine underlying value
    sigt = sigave + meanadjust + KL_sum
    !determine point value
    select case (chxstype)
      case ("total")
        KL_point = sigt
      case ("scatter")
        KL_point = sigt * scatrat(1)
      case ("absorb")
        KL_point = sigt * (1.0d0 - scatrat(1))
      case ("scatrat")
        KL_point = scatrat(1)
        print *,"why did you call me, you already know this info!"
    end select
  elseif(flmatbasedxs) then
    !cross section values
    if(chxstype .ne. 'scatter') &
      siga = sigabsave  + sigameanadjust + sqrt(Coabs)  * KL_sum
    if(chxstype .ne. 'absorb') &
      sigs = sigscatave + sigsmeanadjust + sqrt(Coscat) * KL_sum

    !determine point value
    select case (chxstype)
      case ("total")
        KL_point = Heavi(sigs)*sigs + Heavi(siga)*siga
      case ("scatter")
        KL_point = sigs
      case ("absorb")
        KL_point = siga
      case ("scatrat")
        KL_point = Heavi(sigs)*sigs / ( Heavi(sigs)*sigs + Heavi(siga)*siga )
    end select
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
