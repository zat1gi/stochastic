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
  use KLvars, only: flmeanadjust,flmatbasedxs
  integer :: icase

  call KLrmeshgen         !creates mesh for fixed x and xi material constructions
  call KLrgenrealz(icase) !selects array of random variables xi and tests for negativity
  if(flmeanadjust .and. .not.flmatbasedxs) call KLadjustmean('totaln') !adjusts tot mean after lopping neg cross sections
  if(flmeanadjust .and.      flmatbasedxs) call KLadjustmean('scatter') !adjusts scat mean after lopping neg xss
  if(flmeanadjust .and.      flmatbasedxs) call KLadjustmean('absorb') !adjusts abs mean after lopping neg xss
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
  logical :: flrealzneg, flacceptrealz, flfindzeros
  logical :: flpurpose(3)=.false. !1)neg or not, 2)max vals, 3)zeros

  call cpu_time(tt1)

  write(*,*) "Starting method: KLrec"
  if(flmeanadjust) flfindzeros=.true.                      !zeros (mean adjust, [neg analysis])

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
      call KLr_negsearch( realj, 'scatter', flrealzneg )
      if(.not.flrealzneg) call KLr_negsearch( realj, 'absorb' , flrealzneg )
      !call KLr_negsearch( realj, 'totaln' , flrealzneg ) !don't need, only need if abs or scat is neg
      !!flrealzneg=.false.  !can use this routine, but it is much slower, more thourough, but mothballed!
      !!call KLr_realznegandzeros( realj, 'scatter', flrealzneg, flfindzeros )
      !!if(flfindzeros .or. .not.flrealzneg) call KLr_realznegandzeros( realj, 'absorb', flrealzneg, flfindzeros )

      if(.not.flrealzneg) numPosRealz=numPosRealz+1
      if(     flrealzneg) then
        if(.not.flnegxs) flacceptrealz=.false.
        numNegRealz=numNegRealz+1
        print *,"numNegRealz  : ",numNegRealz," tentative realz#: ",tentj
      endif

      do i=1,KLrnumpoints(2)  !create realization
        KLrxisig(i) = KLr_point(realj,KLrxi(i),'totale')
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
          KLrxisig(i) = KLr_point(pltKLrrealzwhich(1,m),KLrxi(i),'totale',tnumEigsin=tnumEigs)
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





  subroutine KLr_negsearch( j,chxstype,flrealzneg )
  !This subroutine searches for negative values in KL realizations.
  !If negative value found, set flrealzneg=.true., otherwise remain .false..
  use genRealzvars, only: s
  use KLvars, only: alpha, Ak, Eig, numEigs
  integer :: j
  character(*) :: chxstype
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
    minsig=KLr_point(j,minpos,chxstype)
    do k=2,nminnersteps
      xpos=(outerstep*(i-1)+innerstep*(k-1))
      xsig= KLr_point(j,xpos,chxstype)
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
        xsig= KLr_point(j,xpos,chxstype)
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



  subroutine KLr_realznegandzeros( j,chxstype,flrealzneg,flfindzeros )
  !This subroutine discerns if a realization contains negativity for the cross section of
  !interest, and if desired the zeros.  First is searches in each segment of a realization for
  !bounds on zeros.  If determining whether the realization contains negativity is the
  !only goal (flfindzeros==.false.), the subroutine will exit with this information at the
  !first sighting of negativity.  Otherwise it will then cycle through each set of bounds on
  !a zero and find and store the zeros.
  use genRealzvars, only: s, numRealz
  use KLvars, only: KLzerostot, KLzerosabs, KLzerosscat, KLzerostotn, numEigs, numrefinesameiter, KLrmaxnumzeros
  use utilities, only: arithmaticsum, geometricsum
  integer :: j
  character(*) :: chxstype
  logical :: flrealzneg
  logical :: flfindzeros !find zeros (or only test if neg or not)?
  integer :: numslabsecs = 10

  integer :: arsum, secpts, isec, ipt, minfinalsize, izero
  real(8) :: zl,zr
  real(8), allocatable :: zloc(:), zval(:) !zloc and zval for finding zeros in a segment of the domain
  real(8), allocatable :: zlocmaster(:,:), zvalmaster(:,:) !these hold zloc and zval for each material segment zval
  real(8), allocatable :: zlocmaster_(:,:),zvalmaster_(:,:)!temporary arrays for above when enlarging
  integer, allocatable :: zlocsizes(:) !size of zloc(:) held in zlocmaster(:,#)
  real(8), allocatable :: KLzeros(:) !zeros of realization
  real(8), allocatable :: KLzerosabs_(:,:),KLzerosscat_(:,:),KLzerostotn_(:,:) !temporary arrays for increasing size

  KLrmaxnumzeros = 0

  !find bounds on zeros, if only care if negative abort when negative sampled
  arsum = arithmaticsum(1,numEigs,1,numEigs)
  arsum = merge(arsum,60,arsum>60)           !limit max # estimated pts per slab
  secpts = ceiling(real(arsum,8)/real(numslabsecs))
  secpts = merge(secpts,3,secpts<3)          !limit min # pts per segment
  minfinalsize = secpts*2**(numrefinesameiter-1)-geometricsum(1,2,numrefinesameiter-1) !initial size plus min adjustments
  allocate(zlocmaster(minfinalsize,numslabsecs))
  allocate(zvalmaster(minfinalsize,numslabsecs))
  allocate(zlocsizes(numslabsecs))
  do isec=1,numslabsecs                       !cycle through slab segments
    if(allocated(zloc)) deallocate(zloc)      !set up variables
    if(allocated(zval)) deallocate(zval)
    allocate(zloc(minfinalsize))
    allocate(zval(minfinalsize))
    zloc = 0.0d0
    zval = 0.0d0
    zl = s/real(numslabsecs,8)*real(isec-1,8)
    zr = s/real(numslabsecs,8)*real(isec  ,8)
    call KLr_refinezerogrid( j,zloc,zval,zl,zr,secpts,flrealzneg,flfindzeros,chxstype,order=0 ) !find zero bounds in seg
    if(.not.flfindzeros .and. flrealzneg) exit             !if only want whether realz is neg and is, exit
!print *,"zval:",zval
    zlocsizes(isec)=size(zloc)
    if(size(zloc)>size(zlocmaster(:,1))) then !enlarge master arrays if needed
      call move_alloc(zlocmaster,zlocmaster_)
      call move_alloc(zvalmaster,zvalmaster_)
      allocate(zlocmaster(size(zloc),numslabsecs))
      allocate(zvalmaster(size(zval),numslabsecs))
      zlocmaster(1:size(zlocmaster_(:,1)),1:numslabsecs) = zlocmaster_
      zvalmaster(1:size(zvalmaster_(:,1)),1:numslabsecs) = zvalmaster_
      deallocate(zlocmaster_)
      deallocate(zvalmaster_)
    endif
    zlocmaster(1:size(zloc),isec) = zloc      !store values
    zvalmaster(1:size(zloc),isec) = zval
  enddo
  deallocate(zloc)
  deallocate(zval)

  if(flfindzeros) then
    !search in each cell which contains a zero and find the zero
    allocate(KLzeros(KLrmaxnumzeros))
    izero = 1
    do isec = 1,numslabsecs
      do ipt = 1,zlocsizes(isec)-1
        if(zvalmaster(ipt,isec)*zvalmaster(ipt+1,isec)<0) then
          !call zero finder and store
          KLzeros(izero) = KLr_findzeros(j,zlocmaster(ipt,isec),zlocmaster(ipt+1,isec), &
                                           zvalmaster(ipt,isec),zvalmaster(ipt+1,isec),chxstype,order=0) 
          izero = izero + 1
        endif
      enddo
    enddo

    !store zeros in absorb module array
    if(chxstype .eq. 'absorb') then
      if(.not.allocated(KLzerosabs)) allocate(KLzerosabs(KLrmaxnumzeros,numRealz))
      if(size(KLzerosabs(:,1))<KLrmaxnumzeros) then
        call move_alloc(KLzerosabs,KLzerosabs_)
        allocate(KLzerosabs(KLrmaxnumzeros,numRealz))
        KLzerosabs( 1:size(KLzerosabs_(:,1)) , 1:size(KLzerosabs_(1,:)) ) = KLzerosabs_
        deallocate(KLzerosabs_)
      endif
      KLzerosabs(1:size(KLzeros),j) = KLzeros
    endif
    !store zeros in scatter module array
    if(chxstype .eq. 'scatter') then
      if(.not.allocated(KLzerosscat)) allocate(KLzerosscat(KLrmaxnumzeros,numRealz))
      if(size(KLzerosscat(:,1))<KLrmaxnumzeros) then
        call move_alloc(KLzerosscat,KLzerosscat_)
        allocate(KLzerosscat(KLrmaxnumzeros,numRealz))
        KLzerosscat( 1:size(KLzerosscat_(:,1)) , 1:size(KLzerosscat_(1,:)) ) = KLzerosscat_
        deallocate(KLzerosscat_)
      endif
      KLzerosscat(1:size(KLzeros),j) = KLzeros
    endif
    !store zeros in totaln module array
    if(chxstype .eq. 'totaln') then
      if(.not.allocated(KLzerostotn)) allocate(KLzerostotn(KLrmaxnumzeros,numRealz))
      if(size(KLzerostotn(:,1))<KLrmaxnumzeros) then
        call move_alloc(KLzerostotn,KLzerostotn_)
        allocate(KLzerostotn(KLrmaxnumzeros,numRealz))
        KLzerostotn( 1:size(KLzerostotn_(:,1)) , 1:size(KLzerostotn_(1,:)) ) = KLzerostotn_
        deallocate(KLzerostotn_)
      endif
      KLzerostotn(1:size(KLzeros),j) = KLzeros
    endif


    deallocate(zlocmaster)
    deallocate(zvalmaster)
    deallocate(zlocsizes)
    if(allocated(KLzeros)) deallocate(KLzeros)
  endif !if(flfindzeros)

  end subroutine KLr_realznegandzeros




  subroutine KLr_refinezerogrid( j,zloc,zval,zl,zr,secpts,flrealzneg,flfindzeros,chxstype,order )
  !This subroutine searches for zeros in the 'order' derivative of the KL expansion for the 'chxstype'
  !cross section, for the segment of the realization between 'zl' and 'zr', starting the search with
  !'secpts' number of locations.  Locations and values are returned in 'zloc' and 'zval'.
  !'j' is realization number.  If only interested in whether realizations are negative or not
  !(.not.'flfindzeros') then exit with 'flrealzneg'=.true. .
  use KLvars, only: numrefinesameiter, KLrmaxnumzeros
  integer :: j, order
  integer, intent(in) :: secpts
  character(*) :: chxstype
  real(8) :: zl,zr
  real(8), allocatable :: zloc(:), zval(:)
  real(8), allocatable :: zloc_(:),zval_(:)
  integer, allocatable :: znew(:), znew_(:)
  logical :: flrealzneg,flfindzeros

  integer :: i, l, numpts, numpts_, realzwochange, tallychanges, tallychanges_, itersametally
  real(8) :: step

  allocate(znew(size(zloc)))
  znew = 0
  tallychanges = 0
  tallychanges_= 0
  itersametally= 0
  numpts       = secpts
  numpts_      = numpts

  outerdoloop: do
    if(numpts>size(zloc)) then         !enlarge array sizes, keep old values and flags
      deallocate(zloc)
      call move_alloc(zval,zval_)
      call move_alloc(znew,znew_)
      allocate(zloc(2*size(zval_)-1))
      allocate(zval(2*size(zval_)-1))
      allocate(znew(2*size(znew_)-1))
      zloc = 0.0d0
      zval = 0.0d0
      znew = 0
      do i=1,size(zval_)
        zval(2*i-1) = zval_(i)
        znew(2*i-1) = znew_(i)
      enddo
      deallocate(zval_)
      deallocate(znew_)
    elseif(numpts_ .ne. numpts) then   !same array sized, keep old values and flags (set old slots to zero)
      do i=1,numpts_
        l = numpts_-i+1
        zval(2*l-1) = zval(l)
        znew(2*l-1) = znew(l)
        zval(l)     = 0.0d0
        znew(l)     = 0
      enddo
    endif

    step = (zr-zl)/real(numpts-1,8)     !setup location values
    zloc(1) = zl
    do i=2,numpts
      zloc(i) = zloc(i-1) + step
    enddo
!print *,"zloc:",zloc
    do i=1,numpts                       !solve all new needed values
     if(znew(i)==0) then
       zval(i) = KLr_point(j,zloc(i),chxstype,orderin=order)
!print *,"zloc(i):",zloc(i)
!print *,"zval(i):",zval(i),"  chxstype:",chxstype
       if(order==0 .and. .not.flrealzneg  .and. zval(i)<0.0d0) then  !set neg realz flag
         flrealzneg = .true.
         if(.not.flfindzeros) exit outerdoloop        !if only want whether realz is neg ans is, exit
       endif
       znew(i) = 1
     endif
    enddo

    tallychanges = 0                              !count number of ranges found
    do i=1,numpts-1
      if(zval(i)*zval(i+1)<0.0d0) tallychanges = tallychanges + 1
    enddo

    if(tallychanges==tallychanges_) then          !tally iters with same number of changes
      itersametally = itersametally+1
    else
      itersametally = 0
    endif
    tallychanges_ = tallychanges
    numpts_       = numpts
    numpts        = numpts * 2 - 1

    if(itersametally>=numrefinesameiter) exit                      !exit loop if not changed for several iters
  enddo outerdoloop

  KLrmaxnumzeros = KLrmaxnumzeros + tallychanges

  deallocate(znew)
  end subroutine KLr_refinezerogrid




  function KLr_findzeros( j,loc1,loc2,val1,val2,chxstype,order )
  !This function accepts bounds on an zero in any type of cross section and any
  !order of derivative of that cross section and returns the location of the zero.
  !The algorithm assumes val1 or val2 is negative, the other is positive, and there
  !is only one transition between positive and negative between the two.
  integer :: j, order
  real(8) :: loc1,loc2
  real(8) :: KLr_findzeros
  character(*) :: chxstype

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
    val3 = KLr_point(j,loc3,chxstype,orderin=order)

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




  function KLrxi_integral(j,xl,xr,chxstype,tnumEigsin) result(KL_int)
  !This function integrates on KL reconstructed realizations from xl to xr.
  !Integration is on either total, scattering, or absorbing cross section.
  !Routine included mean adjust for any of these.
  use genRealzvars, only: lamc, sigave, Coscat, Coabs, scatrat, sigscatave, sigabsave
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrxivals, meanadjust, flmatbasedxs, &
                    sigsmeanadjust, sigameanadjust
  use utilities, only: Heavi
  integer :: j
  real(8) :: xl,xr,KL_int
  character(*) :: chxstype
  integer, optional :: tnumEigsin

  integer :: curEig, tnumEigs
  real(8) :: Eigfintterm, KL_sum, sigs, siga, sigt

  tnumEigs = merge(tnumEigsin,numEigs,present(tnumEigsin))

  !solve summation of KL terms to use below
  KL_sum = 0d0
  do curEig=1,tnumEigs
    Eigfintterm = Eigfuncint(Ak(curEig),alpha(curEig),lamc,xl,xr)
    KL_sum   = KL_sum + sqrt(Eig(curEig)) * Eigfintterm * KLrxivals(j,curEig)
  enddo

  if(.not.flmatbasedxs) then
    !determine underlying value
    sigt = (sigave + meanadjust)*(xr - xl) + KL_sum
    !determine point value
    select case (chxstype)
      case ("totale")
        KL_int = Heavi(sigt)*sigt
      case ("totaln")
        KL_int = sigt
      case ("scatter")
        KL_int = sigt * scatrat(1)
      case ("absorb")
        KL_int = sigt * (1.0d0 - scatrat(1))
    end select
  elseif(flmatbasedxs) then
    !cross section values
    if(chxstype .ne. 'scatter') &
      siga = (sigabsave  + sigameanadjust)*(xr - xl) + sqrt(Coabs)  * KL_sum
    if(chxstype .ne. 'absorb') &
      sigs = (sigscatave + sigsmeanadjust)*(xr - xl) + sqrt(Coscat) * KL_sum

    !determine point value
    select case (chxstype)
      case ("totale")
        KL_int = Heavi(sigs)*sigs + Heavi(siga)*siga
      case ("totaln")
        KL_int = sigs + siga
      case ("scatter")
        KL_int = sigs
      case ("absorb")
        KL_int = siga
    end select
  endif

  end function KLrxi_integral



  function Eigfuncint(Ak,alpha,lamc,xl,xr)
  ! Used in KLrxi_integral to integrate on KL reconstructed realizations
  real(8) :: Ak,alpha,lamc,xl,xr,Eigfuncint

  Eigfuncint = Ak * (         (-cos(alpha*xr)+cos(alpha*xl))/alpha &
                       + lamc*( sin(alpha*xr)-sin(alpha*xl))           )
  end function Eigfuncint



  function KLr_point(j,xpos,chxstype,tnumEigsin,orderin) result(KL_point)
  !Evaluates KL reconstructed realizations at a given point.
  !It has options for total, scattering only, absorption only, or scattering ratio.
  !It has an option to solve for less than the available number of eigenvalues.
  !It can solve any derivative order of the KL process with optional argument 'orderin'.
  !It can function when in 'material'-based or 'totxs'-based mode.
  !It can function when adjusting mean or not adjusting mean.
  !'totaln', total-native is xs w/o setting to 0, 'totale', total-effective is w/ 0 setting.
  use genRealzvars, only: lamc, scatrat, Coscat, Coabs
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrxivals, flmatbasedxs
  use utilities, only: Heavi

  integer :: j
  real(8) :: xpos
  real(8) :: KL_point
  character(*) :: chxstype
  integer, optional :: tnumEigsin
  integer, optional :: orderin

  real(8) :: sigt, siga, sigs, KL_sum, Eigfterm
  integer :: curEig,tnumEigs,order
  real(8) :: tmeanadjust,tsigsmeanadjust,tsigameanadjust,tsigave,tsigscatave,tsigabsave

  !load any optional values or their default
  tnumEigs = merge(tnumEigsin,numEigs,present(tnumEigsin))
  order    = merge(orderin   ,0      ,present(orderin)   )

  !solve summation of KL terms to use below
  KL_sum = 0d0
  do curEig=1,tnumEigs
    Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,xpos,order)
    KL_sum   = KL_sum + sqrt(Eig(curEig)) * Eigfterm * KLrxivals(j,curEig)
  enddo

  !set non-x-dependent values based order
  call KLr_setmeans(order,tmeanadjust,tsigsmeanadjust,tsigameanadjust,tsigave,tsigscatave,tsigabsave)
  !solve value at point
  if(.not.flmatbasedxs) then
    !determine underlying value
    sigt = tsigave + tmeanadjust + KL_sum
    !determine point value
    select case (chxstype)
      case ("totale") !if deriv, no Heaviside
        KL_point = merge(Heavi(sigt)*sigt,sigt,order==0)
      case ("totaln")
        KL_point = sigt
      case ("scatter")
        KL_point = sigt * scatrat(1)
      case ("absorb")
        KL_point = sigt * (1.0d0 - scatrat(1))
      case ("scatrat")
        KL_point = scatrat(1)
    end select
  elseif(flmatbasedxs) then
    !cross section values
    if(chxstype .ne. 'scatter') &
      siga = tsigabsave  + tsigameanadjust + sqrt(Coabs)  * KL_sum
    if(chxstype .ne. 'absorb') &
      sigs = tsigscatave + tsigsmeanadjust + sqrt(Coscat) * KL_sum
    !determine point value
    select case (chxstype)
      case ("totale") !if deriv, no Heaviside
        KL_point = merge(Heavi(sigs)*sigs + Heavi(siga)*siga, sigs+siga, order==0)
      case ("totaln")
        KL_point = sigs + siga
      case ("scatter")
        KL_point = sigs
      case ("absorb")
        KL_point = siga
      case ("scatrat")
        KL_point = Heavi(sigs)*sigs / ( Heavi(sigs)*sigs + Heavi(siga)*siga )
    end select
  endif
  end function KLr_point




  function Eigfunc(Ak,alpha,lamc,xpos,orderin)
  ! Evaluates Eigenfunction term for any order (0+) derivative of the eigenfunction
  real(8) :: Ak,alpha,lamc,xpos,Eigfunc
  integer, optional :: orderin
  integer :: order

  order = merge(orderin,0,present(orderin))
  select case (mod(order,4))
    case (0)
      Eigfunc = Ak * ( sin(alpha*xpos) + lamc*alpha*cos(alpha*xpos) )
    case (1)
      Eigfunc = Ak * (-cos(alpha*xpos) + lamc*alpha*sin(alpha*xpos) )
    case (2)
      Eigfunc = Ak * (-sin(alpha*xpos) - lamc*alpha*cos(alpha*xpos) )
    case (3)
      Eigfunc = Ak * ( cos(alpha*xpos) - lamc*alpha*sin(alpha*xpos) )
  end select
  Eigfunc = Eigfunc/(alpha**order)

  end function Eigfunc



  subroutine KLr_setmeans(order,tmeanadjust,tsigsmeanadjust,tsigameanadjust,tsigave,tsigscatave,tsigabsave)
  !This subroutine sets values for non-x-dependent terms based on derivative order
  use genRealzvars, only: sigave, sigscatave, sigabsave
  use KLvars, only: meanadjust, sigsmeanadjust, sigameanadjust

  integer :: order
  real(8) :: tmeanadjust,tsigsmeanadjust,tsigameanadjust,tsigave,tsigscatave,tsigabsave
  if(order==0) then
    tmeanadjust     = meanadjust
    tsigsmeanadjust = sigsmeanadjust
    tsigameanadjust = sigameanadjust
    tsigave         = sigave
    tsigscatave     = sigscatave
    tsigabsave      = sigabsave
  else
    tmeanadjust     = 0.0d0
    tsigsmeanadjust = 0.0d0
    tsigameanadjust = 0.0d0
    tsigave         = 0.0d0
    tsigscatave     = 0.0d0
    tsigabsave      = 0.0d0
  endif

  end subroutine KLr_setmeans



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

  subroutine KLadjustmean(chxstype)
  !This subroutine is the master for setting the value of "meanadjust", which will conserve
  !the mean of the reconstructions after ignoring negative values in transport within the 
  !chosen tolerance
  use genRealzvars, only: s, sigave, sigscatave, sigabsave
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrnumRealz, meanadjust, meanadjust_tol, sigsmeanadjust, &
                    sigameanadjust
  use KLreconstruct, only: KLr_point, KLrxi_integral

  character(*) :: chxstype

  integer :: j,adjustiter
  real(8) :: intsigave,areacont,xmid

  !initialize meanadjust
  if(chxstype=='scatter') sigsmeanadjust = 0.0d0
  sigameanadjust = 0.0d0
  meanadjust = 0.0d0

  !integrate on all as check
  intsigave = 0d0
  do j=1,KLrnumRealz
    intsigave = intsigave + &
                KLrxi_integral(j,0d0,s,chxstype='totaln')/KLrnumRealz/s
  enddo
  write(*,*)
  500 format("  Integrator/reconstruction check - sigave: ",f8.5,"  intsigave: ",f8.5,"  relerr: ",es10.2)
  write(*,500) sigave,intsigave,abs(sigave-intsigave)/sigave


  !meanadjust solver
  step = s/(numEigs*32d0) !hard parameter, relative step size

  adjustiter=0
  do
    adjustiter=adjustiter+1
    aveposarea   = 0d0
    perposdomain = 0d0
    avenegarea   = 0d0
    pernegdomain = 0d0

    print *,"Beginning mean adjustment iteration ",adjustiter," for cross section:",chxstype
    do j=1,KLrnumRealz
      xr = KLr_point(j,0d0,chxstype)
      xl = 0d0
      do
        !find next point and area between these two
        xr = findnextpoint(j,chxstype)
        areacont = KLrxi_integral(j,xl,xr,chxstype=chxstype)/KLrnumRealz/s

        !calc aveposarea for tol check, also negstat tallies
        xmid = KLr_point(j,(xr+xl)/2d0,chxstype)
        if(xmid>0d0) then
          aveposarea   = aveposarea + areacont
          perposdomain = perposdomain + (xr - xl)/KLrnumRealz/s * 100
        else
          avenegarea   = avenegarea + areacont
          pernegdomain = pernegdomain + (xr - xl)/KLrnumRealz/s * 100
        endif

        !advance lagging point
        xl = xr
        if(xr>=s) exit
      enddo !sum within realization
    enddo !loop over realizations

    print *,"avenegarea: ",avenegarea,"  aveposarea: ",aveposarea
    if(chxstype=='totaln') then
      meanadjust     = meanadjust     + (sigave - aveposarea)
      if(abs(aveposarea-sigave)/sigave<meanadjust_tol) exit
    elseif(chxstype=='scatter') then
      sigsmeanadjust = sigsmeanadjust + (sigscatave - aveposarea)
      meanadjust     = sigsmeanadjust +    sigameanadjust
      print *,"sigsmeanadjust: ",sigsmeanadjust
      if(abs(aveposarea-sigscatave)/sigave<meanadjust_tol) exit
    elseif(chxstype=='absorb') then
      sigameanadjust = sigameanadjust + (sigabsave - aveposarea)
      meanadjust     = sigsmeanadjust +    sigameanadjust
      print *,"sigameanadjust: ",sigameanadjust
      if(abs(aveposarea-sigabsave)/sigave<meanadjust_tol) exit
    endif
    print *,"meanadjust    : ",meanadjust
  enddo !loop over calculating adjustment
    
  end subroutine KLadjustmean



  function findnextpoint(j,chxstype)
  !This function searches ahead and finds either 1) next time a reconstructed realization
  !changes signs, or 2) the end of the slab
  use genRealzvars, only: s, sigave
  use KLvars,       only: alpha, Ak, Eig, numEigs, KLrnumRealz
  use KLreconstruct, only: KLr_point
  integer :: j
  character(*) :: chxstype

  real(8) :: findnextpoint
  real(8) :: curx,oldx,curs,olds !position, then sigma value

  curx = xl
  curs = KLr_point(j,curx,chxstype)
  do 
    oldx = curx
    olds = curs

    curx = curx + step
    curs = KLr_point(j,curx,chxstype)
    if(curs*olds<0d0) then
      curx = refinenextpoint(j,oldx,curx,chxstype)
      exit
    elseif(curx>=s) then
      curx = s
      exit
    endif

  enddo
  findnextpoint = curx

  end function findnextpoint



  function refinenextpoint(j,oldx,curx,chxstype)
  !This function takes a range and zeroes in on transition in sign of cross section within tolerance
  use genRealzvars, only: s, sigave
  use KLvars, only: alpha, Ak, Eig, numEigs, KLrnumRealz
  use KLreconstruct, only: KLr_point
  integer :: j
  real(8) :: oldx,curx
  character(*) :: chxstype

  real(8) :: refinenextpoint,stepsign,curs,olds,curstep

  stepsign = -1d0
  curstep = step
  curs = KLr_point(j,curx,chxstype)
  do
    curstep = curstep/2d0
    oldx = curx
    olds = curs
    curx = curx + curstep*stepsign
    curs = KLr_point(j,curx,chxstype)

    if(abs(curs)<stol) then
      refinenextpoint = curx
      exit
    endif
    stepsign = stepsign * merge(1d0,-1d0,curs*olds>0d0)
  enddo

  end function refinenextpoint


end module KLmeanadjust
