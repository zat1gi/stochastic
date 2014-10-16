module KLreconstruct
  use utilities
  use timeman
  use KLmeanadjust
  use mcnp_random
  implicit none

CONTAINS
  ! print statemtns in this module use # 500-599


  subroutine KLrmeshgen
  !This subroutine creates a mesh based on selected frequency of 
  !sampling in x for a fixed point reconstruction, and then for a fixed xi
  !construction.  The fixed xi construction is the only one we care about.
  use genRealzvars , only: s
  use KLvars, only: KLrnumpoints, KLrx, KLrxi

  integer :: i
  real(8) :: KLrxstepsize

  allocate(KLrx(KLrnumpoints(1)))
  KLrx = 0                             !create mesh for fixed point KL reconstruction
  KLrxstepsize = s / KLrnumpoints(1)
  do i=1,KLrnumpoints(1)
    KLrx(i) = KLrxstepsize*i - KLrxstepsize/2
  enddo

  allocate(KLrxi(KLrnumpoints(2)))
  KLrxi = 0                            !create mesh for fixed xi KL reconstruction
  KLrxstepsize = s / KLrnumpoints(2)
  do i=1,KLrnumpoints(2)
    KLrxi(i) = KLrxstepsize*i - KLrxstepsize/2
  enddo

  end subroutine KLrmeshgen







  subroutine KLrgenrealz
  !This subroutine reconstructs realizations based upon the KL expansion.
  !It reconstructs based upon the fixed point and fixed xi methods.
  !It also passes an array of selected random variables xi to be plotted in KLreval.
  use timevars, only: time
  use genRealzvars, only: s, lamc, sigave
  use KLvars,       only: gam, alpha, Ak, Eig, binPDF, binNumof, numEigs, &
                          KLrnumpoints, KLrnumRealz, KLrprintat, negcnt, pltKLrrealz, &
                          pltKLrrealznumof, pltKLrrealzwhich, KLrx, KLrxi, KLrxivals, &
                          pltKLrrealzarray, KLrrandarray, KLrsig, KLrxisig, CoExp
  integer :: i,j,curEig,w,u
  real(8) :: KLsigtemp,Eigfterm,xiterm,rand,tt1,tt2
  character(3) :: neg
  character(5) :: flKLtype = 'KLrec'

  call cpu_time(tt1)

  write(*,*) "Starting method: ",flKLtype  
  do j=1,KLrnumRealz

    KLrsig = 0          !create a realization, fixed point
    do i=1,KLrnumpoints(1)
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
                               !^ sqrt(CoExp) to remove Co
      enddo
    enddo
    612 format("  ",f14.8)     !print sigma values to text file, fixed point
    open(unit=10,file="KLrsig.txt")
    do i=1,KLrnumpoints(1)
      write(10,612,advance="no") KLrsig(i)
    enddo
    write(10,*)



    KLrxisig = 0      !create a realization, fixed xi
    do curEig=1,numEigs  !select xi values
      rand = rang()
      call select_from_PDF( binPDF,binNumof,curEig,xiterm,rand )
      KLrxivals(j,curEig) = xiterm
    enddo

    neg='no'
    call KLr_negsearch( j,neg )
    if(neg=='yes') then  !counts the number of realz that contain a negative value
      negcnt=negcnt+1
      print *,"negcnt  : ",negcnt,"          realz#: ",j
    endif

    do i=1,KLrnumpoints(2)  !create realization
      KLrxisig(i) = KLrxi_point(j,KLrxi(i))
    enddo
    open(unit=11,file="KLrxisig.txt") !print sigma values to text file, fixed xi
    do i=1,KLrnumpoints(2)
      write(11,612,advance="no") KLrxisig(i)
    enddo
    write(11,*)


    if(mod(j,KLrprintat)==0) call KL_timeupdate( j,tt1,flKLtype )
  enddo

  end subroutine KLrgenrealz







  subroutine KLrplotrealz
  !This subroutine uses the stored array of "random" numbers used in KLrgenrealz
  !to plot the selected reconstructed realizations.
  use genRealzvars, only: lamc, sigave
  use KLvars,      only: gam, alpha, Ak, Eig, binPDF, binNumof, numEigs, &
                         KLrnumpoints, negcnt, pltKLrrealz, pltKLrrealznumof, &
                         pltKLrrealzwhich, KLrx, KLrxi, pltKLrrealzarray, KLrrandarray, &
                         KLrsig, KLrxisig, pltKLrrealzPointorXi, CoExp
  use radtransMC, only: KLrxi_point2

  integer :: i,curEig,m,KLrnumpts,tnumEigs
  real(8) :: KLsigtemp,Eigfterm,xiterm,rand

  call system("mv KLrsig.txt plots")
  call system("mv KLrxisig.txt plots")

  if( pltKLrrealz(1) .NE. 'noplot' ) then  !plot using generic plotter
    do m=1,pltKLrrealznumof
      tnumEigs=pltKLrrealzwhich(2,m)


      if( pltKLrrealzPointorXi(m) .EQ. 'fpoint' ) then  !create a realz, fixed point
        KLrnumpts=KLrnumpoints(1)
        KLrsig = 0
        do i=1,KLrnumpoints(1)
          KLrsig(i) = sigave
          do curEig=1,tnumEigs
            Eigfterm = Eigfunc(Ak(curEig),alpha(curEig),lamc,KLrx(i))
            rand = KLrrandarray(i,curEig,m+1)
            call select_from_PDF( binPDF,binNumof,numEigs,xiterm,rand )
            KLrsig(i) = KLrsig(i) + sqrt(Eig(curEig)) * Eigfterm * xiterm
                                   !^ sqrt(CoExp) to remove Co
          enddo
          pltKLrrealzarray(i,1)   = KLrx(i)    !record x values
          pltKLrrealzarray(i,m+1) = KLrsig(i)  !record that realization
        enddo
      endif



      if( pltKLrrealzPointorXi(m) .EQ. 'fxi' ) then  !create a realz, fixed xi
        KLrnumpts=KLrnumpoints(2)
        KLrxisig = 0
        do i=1,KLrnumpoints(2)
          KLrxisig(i) = KLrxi_point2(pltKLrrealzwhich(1,m),KLrxi(i),'total  ',tnumEigsin=tnumEigs)
          pltKLrrealzarray(i,1)   = KLrxi(i)     !record x values
          pltKLrrealzarray(i,m+1) = KLrxisig(i)  !record that realization
        enddo
      endif
    enddo



    call generic_plotter( KLrnumpts,pltKLrrealznumof,pltKLrrealzarray,&
                          pltKLrrealz )

    call system("mv genericplot.txt plots/KLrrealzplot.txt")
    call system("mv genericplot.ps  plots/KLrrealzplot.ps")
    call system("mv genericplot.pdf plots/KLrrealzplot.pdf")
    call system("rm plots/KLrrealzplot.eps")
    call system("mv genericplot.eps plots/KLrrealzplot.eps")
  endif

  print *," Total num reconstructed realz w/ neg value: ",negcnt
  print *,

  end subroutine KLrplotrealz





  subroutine KLr_negsearch( j,neg )
  use genRealzvars, only: s
  use KLvars, only: alpha, Ak, Eig, numEigs
  integer :: j
  character(3) :: neg

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
    minsig=KLrxi_point(j,minpos)
    do k=2,nminnersteps
      xpos=(outerstep*(i-1)+innerstep*(k-1))
      xsig= KLrxi_point(j,xpos)
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
        xsig= KLrxi_point(j,xpos)
        if(xsig<minsig) then
          minsig=xsig
          minpos=xpos
        endif
      enddo
    enddo
    if(minsig<0) then
      neg='yes'
print *,"minpos",minpos,"minsig",minsig
      exit
    endif

  enddo

  end subroutine KLr_negsearch



end module KLreconstruct
