module FEDiffSn
  implicit none
  !true if vars after next comment are passed in (as opposed to read)
  logical              :: flvarspassed = .false.

  !variables that can be passed to this module
  real(8)              :: sigt
  real(8)              :: c
  integer              :: numcells

  !variables that are passed from this module
  integer, allocatable :: solve(:) !solve(1)==1, diffyes, solve(2)==1, Snyes, solve(3)==1, DSAyes
  real(8), allocatable :: phidiff(:)
  real(8), allocatable :: phiSnl(:) ,phiSnr(:)
  real(8), allocatable :: phiDSAl(:),phiDSAr(:)

CONTAINS

  !---------------------------------------------------------------------------------
  !
  !  Here begins what was before the main program
  !
  !---------------------------------------------------------------------------------

  subroutine FEmain
  use utilities, only: gauss_leg_quad

  !Input Vars
  integer :: numangs
  real(8) :: a,curRight,curLeft,constq,phistart,tol
  character(5) :: qtype !'evend'const src,'func1'quad src,'moms'methodofmansolns src
  character(7) :: SnBCs !'vacuum' for vacuum
  logical :: flplot

  !Derived Vars
  integer :: iter
  real(8) :: siga,dx,error,tic,toc,Sntime,Sntimeperiter,DSAtime,DSAtimeperiter,timehere
  real(8),allocatable :: phipos(:),q(:),                x(:)
  real(8),allocatable :: phiSnlold(:),phiSnrold(:)
  real(8),allocatable :: psil(:,:),psir(:,:),psiBCl(:),psiBCr(:),qr(:),ql(:),mu(:),wgts(:)
  real(8),allocatable :: phiHsl(:),phiHsr(:),phiresl(:),phiresr(:)
  real(8),allocatable :: phiDSAlold(:),phiDSArold(:),phiDSAdiff(:)
  real(8) :: alpha,beta

  !Read and Prepare parameters
  call input(            a,curRight,curLeft,constq,qtype,phistart,&
                         solve,numangs,SnBCs,tol,flplot )
       siga=             sigt*(1.0-c)


  !----------------- Solve diffusion problem (solve(1)==1) -----------------!
  if(solve(1)==1) then
    call meshgen(        numcells,a,x,dx ) !create mesh in x
    call fluxinitcont(   numcells,phidiff,phistart ) !create phi=phistart
    if(qtype=='moms') call FEMomsPosited( numcells,sigt,siga,a,x,phipos,qtype ) !Set moms pos soln
    call diffsourceset(  numcells,a,x,dx,curRight,curLeft,constq,qtype,q,sigt,siga ) !Set q w/ BCs
    alpha = 0.25
    beta  = 0.5
    call FEDiff(         numcells,sigt,c,a,siga,dx,phidiff,q,alpha,beta ) !Solve diff Prob
    if(flplot) call difffluxplot(   numcells,x,phidiff,phipos,qtype ) !Plot diff flux
    if(solve(2)==1 .or. solve(3)==1) deallocate(q)
    if(solve(2)==1 .or. solve(3)==1) deallocate(x)
  endif


  !----------------- Solve Sn problem (solve(2)==1) --------------------!
  if(solve(2)==1) then
!    print *
!    print *,"Beginning Sn solve"
    call cpu_time(tic)
    error = 1.0
    call meshgen(        numcells,a,x,dx ) !create mesh in x
    call gauss_leg_quad( numangs,mu,wgts,real(-1.0,8),real(1.0,8) ) !create mesh in mu
    call fluxinitdisc(   numcells,phiSnl,phistart ) !create phil=phistart
    call fluxinitdisc(   numcells,phiSnr,phistart ) !create phir=phistart
    call fluxinitdisc(   numcells,phiSnlold,phistart ) !create old phil for error calc
    call fluxinitdisc(   numcells,phiSnrold,phistart ) !create old phir for error calc
    call psiinit(        numcells,psil,numangs ) !create psi=0.0
    call psiinit(        numcells,psir,numangs ) !create psi=0.0
    call Snsourceset(    numcells,a,x,dx,constq,qtype,qr,ql,sigt,siga,SnBCs,&
                         psiBCl,psiBCr,numangs)!Set qs & psiBCs
    iter=0
    do while (error>tol)
      iter=iter+1
      call FESn(         numcells,phiSnl,phiSnr,psil,psir,x,dx,qr,ql,sigt,siga,mu,&
                         psiBCl,psiBCr,numangs ) !Solve an Sn sweep
      call FESncalcphi(  numcells,phiSnl,phiSnr,psil,psir,wgts,numangs ) !Calc phi from psi
      if(iter>1) call FESnerror( numcells,phiSnl,phiSnr,phiSnlold,phiSnrold,error ) !Check converge?
      call FESnNewToOld( numcells,phiSnl,phiSnr,phiSnlold,phiSnrold )
    enddo
!    print *,"iteration: ",iter
!    call cpu_time(toc)
!    Sntime = toc-tic
!    Sntimeperiter = Sntime/iter
!    1000 format(" time     : ",f12.7," sec       Time/iter: ",f12.7," sec")
!    write(*,1000) Sntime,Sntimeperiter
!    print *
    if(flplot) call Snfluxplot(     numcells,x,phiSnl,phiSnr,qtype ) !Plot Sn flux
    if(solve(3)==1) deallocate(x)
    if(solve(3)==1) deallocate(qr)
    if(solve(3)==1) deallocate(ql)
    if(solve(3)==1) deallocate(psil)
    if(solve(3)==1) deallocate(psir)
    if(solve(3)==1) deallocate(mu)
    if(solve(3)==1) deallocate(wgts)
    if(solve(3)==1) deallocate(psiBCl)
    if(solve(3)==1) deallocate(psiBCr)
  endif



  !-------------------- Solve Sn with DSA (solve(3)==1) ------------------!
  if(solve(3)==1) then
!    print *
!    print *,"Beginning DSA solve"
    call cpu_time(tic)
    error = 1.0
    !init for Sn
    call meshgen(        numcells,a,x,dx ) !create mesh in x
    call gauss_leg_quad( numangs,mu,wgts,real(-1.0,8),real(1.0,8) ) !create mesh in mu
    call alpha_beta(     numangs,mu,wgts,alpha,beta ) !solve for alpha and beta
    call fluxinitdisc(   numcells,phiDSAl,phistart ) !create phil=phistart
    call fluxinitdisc(   numcells,phiDSAr,phistart ) !create phir=phistart
    call fluxinitdisc(   numcells,phiHsl,phistart ) !create half-step phil to update
    call fluxinitdisc(   numcells,phiHsr,phistart ) !create half-step phil to update
    call fluxinitdisc(   numcells,phiresl,phistart ) !create residual phil for update
    call fluxinitdisc(   numcells,phiresr,phistart ) !create residual phil for update
    call fluxinitdisc(   numcells,phiDSAlold,phistart ) !create old phil for error calc
    call fluxinitdisc(   numcells,phiDSArold,phistart ) !create old phir for error calc
    call psiinit(        numcells,psil,numangs ) !create psi=0.0
    call psiinit(        numcells,psir,numangs ) !create psi=0.0
    call Snsourceset(    numcells,a,x,dx,constq,qtype,qr,ql,sigt,siga,SnBCs,&
                         psiBCl,psiBCr,numangs)!Set qs & psiBCs
    !init for diff
    call fluxinitcont(   numcells,phiDSAdiff,phistart ) !create phi=phistart
    iter=0
    do while (error>tol)
      iter=iter+1
      !Step 1, Sn solve for half step solution of flux
      call FESn(         numcells,phiDSAlold,phiDSArold,psil,psir,x,dx,qr,ql,sigt,siga,mu,&
                         psiBCl,psiBCr,numangs ) !Solve an Sn sweep
      call FESncalcphi(  numcells,phiHsl,phiHsr,psil,psir,wgts,numangs ) !Calc phi from psi
      !Step 2, Diffusion solve with half step informed source
      call DSAdiffsource(numcells,a,x,dx,q,sigt,siga,phiHsl,phiHsr,&
                         phiDSAlold,phiDSArold,iter ) !Set q w/ BCs
      call FEDiff(       numcells,sigt,c,a,siga,dx,phiDSAdiff,q,alpha,beta ) !Solve diff Prob
      !Step 3, Simplified WLA to relate cont & discont flux to get residual, add to half step
      call SWLA(         numcells,phiDSAdiff,phiDSAlold,phiDSArold,phiHsl,phiHsr,&
                         phiresl,phiresr,sigt,siga,dx,alpha,beta ) !Solve SWLA for residual
      call HalfResSum(   numcells,phiDSAl,phiDSAr,phiresl,phiresr,phiHsl,phiHsr ) !hs + res

      if(iter>1) call FESnerror( numcells,phiDSAl,phiDSAr,phiDSAlold,phiDSArold,error ) !Converge?
      call FESnNewToOld( numcells,phiDSAl,phiDSAr,phiDSAlold,phiDSArold )
    enddo
!    print *,"iteration: ",iter
!    call cpu_time(toc)
!    DSAtime = toc-tic
!    DSAtimeperiter = DSAtime/iter
!    1001 format(" time     : ",f12.7," sec       Time/iter: ",f12.7," sec")
!    write(*,1001) DSAtime,DSAtimeperiter
!    print *
    if(flplot) call DSAfluxplot(     numcells,x,phiDSAl,phiDSAr,qtype ) !Plot DSA flux
  endif



  if(solve(2)==1 .and. solve(3)==1) then
    call compareSnDSA(  numcells,x,phiSnl,phiSnr,phiDSAl,phiDSAr,qtype,flplot ) !soln(Sn=DSA)?
  endif

  call FEDiffSn_internaldeallocate( phipos,q,x,phiSnlold, &
                                 phiSnrold,psil,psir,psiBCl,psiBCr,qr,ql,mu,wgts, &
                                 phiHsl,phiHsr,phiresl,phiresr, &
                                 phiDSAlold,phiDSArold,phiDSAdiff )

  end subroutine FEmain



  !---------------------------------------------------------------------------------
  !
  !  Here begins what was before in the io module
  !
  !---------------------------------------------------------------------------------



  subroutine input( a,curRight,curLeft,constq,qtype,phistart,solve,numangs,&
                    SnBCs,tol,flplot )
  use genRealzvars, only: s
  !Reads in values from input file
  integer :: numangs
  integer,allocatable :: solve(:) !solve(1)==1, diffyes, solve(2)==1, Snyes
  real(8) :: a,curRight,curLeft,constq,phistart,tol
  character(5) :: qtype !'evend' for constant source, 'func1' for quadratic source
  character(7) :: SnBCs !'vacuum' for vacuum
  logical      :: flplot       !plot or not

  character(3) :: plotprofile  !plot or not, store info in logical form
  character(3) :: dumchar !lets me skip a line of character input so I can label
  real(8)      :: dumreal !lets me skip a line of real input

  allocate(solve(3))

  open(unit=2,file="auxiliary/FEDiffSn.inp")
  read(2,*) dumchar    !Problem Type
  read(2,*) solve(1),solve(2),solve(3)
  read(2,*) plotprofile

  read(2,*) dumchar    !Material and Mesh
  if(flvarspassed) then
    read(2,*) dumreal,dumreal
    read(2,*) dumreal
  else
    read(2,*) sigt,c
    read(2,*) numcells
  endif
  read(2,*) a
  a=s !override value with value from other input file
  read(2,*) numangs

  read(2,*) dumchar    !Source
  read(2,*) qtype
  read(2,*) constq
  read(2,*) curRight,curLeft
  read(2,*) SnBCs

  read(2,*) dumchar    !Initflux
  read(2,*) phistart

  read(2,*) dumchar    !Iteration Inputs
  read(2,*) tol

  close(unit=2)

  if(plotprofile=='yes') then
    flplot = .true.
  elseif(plotprofile=='no') then
    flplot = .false.
  else
    stop "FEDiffSn for plotting enter 'yes' or 'no'"
  endif

  end subroutine input




  subroutine difffluxplot( numcells,x,phi,phipos,qtype )
  !Prints x and phi center-cell values and plots
  integer :: numcells
  real(8) :: x(:),phi(:),phipos(:)
  character(5) :: qtype

  integer :: i

  if(qtype=='evend' .or. qtype=='func1') then
    open(unit=3,file="diffflux.txt")
    100 format("#  xpos   diffphi")
    101 format(f10.6,"  ",f10.6)
    write(3,100)
    do i=1,numcells
      write(3,101) (x(i)+x(i+1))/2,(phi(i)+phi(i+1))/2
    enddo
    close(unit=3)

    call system("gnuplot plots/FEDiffSn/diffflux.gnu")
    call system("ps2pdf diffflux.ps")
    call system("mv diffflux.txt diffflux.ps diffflux.pdf plots/FEDiffSn")
  elseif(qtype=='moms') then
    open(unit=3,file="fluxmoms.txt")
    102 format("#  xpos   diffphi   phipos")
    103 format(f10.6,"  ",f10.6,"  ",f10.6)
    write(3,102)
    do i=1,numcells
      write(3,103) (x(i)+x(i+1))/2,(phi(i)+phi(i+1))/2,(phipos(i)+phipos(i+1))/2
    enddo
    close(unit=3)

    call system("gnuplot plots/FEDiffSn/fluxmoms.gnu")
    call system("ps2pdf fluxmoms.ps")
    call system("mv fluxmoms.txt fluxmoms.ps fluxmoms.pdf plots/FEDiffSn")
  endif

  end subroutine difffluxplot


  subroutine Snfluxplot( numcells,x,phir,phil,qtype )
  !Prints x and phi center-cell values and plots
  integer :: numcells
  real(8) :: x(:),phir(:),phil(:)
  character(5) :: qtype

  integer :: i

  if(qtype=='evend' .or. qtype=='func1') then
    open(unit=3,file="Snflux.txt")
    110 format("#  xpos   Snphi")
    111 format(f10.6,"  ",f10.6)
    write(3,110)
    do i=1,numcells
      write(3,111) (x(i)+x(i+1))/2,(phil(i)+phir(i))/2
    enddo
    close(unit=3)

    call system("gnuplot plots/FEDiffSn/Snflux.gnu")
    call system("ps2pdf Snflux.ps")
    call system("mv Snflux.txt Snflux.ps Snflux.pdf plots/FEDiffSn")
  endif

  end subroutine Snfluxplot



  subroutine DSAfluxplot( numcells,x,phir,phil,qtype )
  !Prints x and phi center-cell values and plots
  integer :: numcells
  real(8) :: x(:),phir(:),phil(:)
  character(5) :: qtype

  integer :: i

  if(qtype=='evend' .or. qtype=='func1') then
    open(unit=4,file="DSAflux.txt")
    120 format("#  xpos   DSAphi")
    121 format(f10.6,"  ",f10.6)
    write(4,120)
    do i=1,numcells
      write(4,121) (x(i)+x(i+1))/2,(phil(i)+phir(i))/2
    enddo
    close(unit=4)

    call system("gnuplot plots/FEDiffSn/DSAflux.gnu")
    call system("ps2pdf DSAflux.ps")
    call system("mv DSAflux.txt DSAflux.ps DSAflux.pdf plots/FEDiffSn")
  endif

  end subroutine DSAfluxplot




  subroutine compareSnDSA( numcells,x,phiSnl,phiSnr,phiDSAl,phiDSAr,qtype,flplot )
  !Prints norm of center-cell values for Sn and DSA soltions
  integer :: numcells
  real(8) :: x(:),phiSnr(:),phiSnl(:),phiDSAl(:),phiDSAr(:)
  character(5) :: qtype
  logical :: flplot

  integer :: i

  if(qtype=='evend' .or. qtype=='func1') then
    open(unit=5,file="SnDSAnorm.txt")
    130 format("#  xpos   SnDSAerr")
    131 format(f10.6,"  ",f20.16)
    write(5,130)
    do i=1,numcells
      write(5,131) (x(i)+x(i+1))/2,(phiSnl(i)+phiSnr(i))/2-(phiDSAl(i)+phiDSAr(i))/2
    enddo
    close(unit=5)

    if(flplot) then
      call system("gnuplot plots/FEDiffSn/SnDSAnorm.gnu")
      call system("ps2pdf SnDSAnorm.ps")
      call system("mv SnDSAnorm.txt SnDSAnorm.ps SnDSAnorm.pdf plots/FEDiffSn")
    endif
  endif

  end subroutine compareSnDSA



  !---------------------------------------------------------------------------------
  !
  !  Here begins what was before in the FEsetup module
  !
  !---------------------------------------------------------------------------------

!----------------- Used Solely in Diffusion -----------------!

  subroutine diffsourceset(  numcells,a,x,dx,curRight,curLeft,constq,qtype,q,sigt,siga )
  !Sets source for diffusion options
  integer :: numcells
  real(8) :: a,x(:),dx,curRight,curLeft,constq,sigt,siga
  real(8),allocatable :: q(:)
  character(5) :: qtype

  integer :: i

  allocate(q(numcells+1))
  q=0.0
  if(qtype=='evend') then !constant source
    q(1) = constq*dx/2.0 + 2*curLeft
    do i=2,numcells
      q(i) = 2*constq*dx/2.0
    enddo
    q(numcells+1) = constq*dx/2.0 - 2*curRight
  elseif(qtype=='func1') then !quadratic source
    q(1) = quadql(a,x(1),x(2)) + 2*curLeft
    do i=2,numcells
      q(i) = quadqr(a,x(i-1),x(i))+quadql(a,x(i),x(i+1))
    enddo
    q(numcells+1) = quadqr(a,x(numcells),x(numcells+1)) - 2*curRight
  elseif(qtype=='moms') then !method of man solns with quad source
    q(1) = momsl(a,x(1),x(2),sigt,siga) + 2*curLeft
    do i=2,numcells
      q(i) = momsr(a,x(i-1),x(i),sigt,siga)+momsl(a,x(i),x(i+1),sigt,siga)
    enddo
    q(numcells+1) = momsr(a,x(numcells+1),x(numcells),sigt,siga) - 2*curRight
  endif
  end subroutine diffsourceset


  function momsl(a,xl,xr,sigt,siga)
  real(8) :: a,xl,xr,momsl,c2,c1,c0,a2,a1,a0,Aden,sigt,siga
  Aden =  9.0 *sigt*sigt*a*a-32.0
  a0   = -4.0 *(9.0*sigt*sigt*a - 6.0*sigt*a + 8.0)/Aden
  a1   =  12.0*sigt*(3.0*sigt*a - 4.0)/Aden
  a2   = -36.0*sigt*sigt/Aden
  c2   = a2*siga
  c1   = a1*siga
  c0   = a0*siga-2.0*a2/(3.0*sigt)
  momsl = (xr-xl)*((xr*xr+2.0*xl*xr+3.0*xl*xl)*c2/12.0 + (2.0*xl+xr)*c1/6.0 + c0/2.0)
  end function momsl
 

  function momsr(a,xl,xr,sigt,siga)
  real(8) :: a,xl,xr,momsr,c2,c1,c0,a2,a1,a0,Aden,sigt,siga
  Aden =  9.0 *sigt*sigt*a*a-32.0
  a0   = -4.0 *(9.0*sigt*sigt*a - 6.0*sigt*a + 8.0)/Aden
  a1   =  12.0*sigt*(3.0*sigt*a - 4.0)/Aden
  a2   = -36.0*sigt*sigt/Aden
  c2   = a2*siga
  c1   = a1*siga
  c0   = a0*siga-2.0*a2/(3.0*sigt)
  momsr = (xr-xl)*((3.0*xr*xr+2.0*xl*xr+xl*xl)*c2/12.0 + (xl+2.0*xr)*c1/6.0 + c0/2.0)
  end function momsr


  subroutine FEMomsPosited( numcells,sigt,siga,a,x,phipos,qtype )
  !Creates posited flux profile that method of manufactured solutions
  !will create a source that also yields.  A form of benchmark

  !I think the coeffs are derived wrongly. Try phipos(0)=A0/Aden, it is neg.
  !I think need to rederive a_is from Marshak BCs and normalization
  integer :: numcells
  real(8) :: sigt,siga,a,x(:)
  character(5) :: qtype
  real(8),allocatable :: phipos(:)

  real(8) :: Aden,a0,a1,a2
  integer :: i

  allocate(phipos(numcells+1))
  phipos = 0.0
  
  Aden =  9.0 *sigt*sigt*a*a - 32.0
  a0   = -4.0 *(9.0*sigt*sigt*a - 6.0*sigt*a + 8.0)/Aden
  a1   =  12.0*sigt*(3.0*sigt*a - 4.0)/Aden
  a2   = -36.0*sigt*sigt/Aden

  do i=1,numcells+1
    phipos(i) = a2*x(i)*x(i) + a1*x(i) + a0
    print *,x(i),phipos(i)
  enddo
  end subroutine FEMomsPosited


!----------------- Used in Diffusion and Sn -----------------!

  subroutine meshgen( numcells,a,x,dx )
  !Generates length mesh
  integer :: numcells
  real(8) :: a,dx
  real(8),allocatable :: x(:)

  integer :: i

  allocate(x(numcells+1))
  x = 0.0
  
  dx = a/numcells

  x(1) = 0.0
  do i=2,numcells+1
    x(i) = x(i-1) + dx
  enddo
  end subroutine meshgen


  function quadql(a,xl,xr)
  real(8) :: a,xl,xr,quadql,c2,c1,c0
  c2 = -4.0 / (a*a)
  c1 =  4.0 / a
  c0 =  0.0
  quadql = (xr-xl)*((xr*xr+2.0*xl*xr+3.0*xl*xl)*c2/12.0 + (2.0*xl+xr)*c1/6.0 + c0/2.0)
  end function quadql
 

  function quadqr(a,xl,xr)
  real(8) :: a,xl,xr,quadqr,c2,c1,c0
  c2 = -4.0 / (a*a)
  c1 =  4.0 / a
  c0 =  0.0
  quadqr = (xr-xl)*((3.0*xr*xr+2.0*xl*xr+xl*xl)*c2/12.0 + (xl+2.0*xr)*c1/6.0 + c0/2.0)
  end function quadqr


!----------------- Used Solely in Sn -----------------!

  subroutine psiinit( numcells,psi,numangs )
  !Initializes psi to psi=0.0
  integer :: numcells,numangs
  real(8),allocatable :: psi(:,:)

  allocate(psi(numangs,numcells))
  psi = 0.0
  end subroutine psiinit



  subroutine Snsourceset( numcells,a,x,dx,constq,qtype,qr,ql,sigt,siga,SnBCs,&
                          psiBCl,psiBCr,numangs )
  !Set q values and BCs (part of psi)
  integer :: numcells,numangs
  real(8) :: a,x(:),dx,constq,sigt,siga
  character(5) :: qtype
  character(7) :: SnBCs
  real(8),allocatable :: qr(:),ql(:),psiBCl(:),psiBCr(:)

  integer :: i,m

  allocate(qr(numcells))
  allocate(ql(numcells))
  qr = 0.0
  ql = 0.0

  if(qtype=='evend') then !constant source
    do i=1,numcells
      qr(i) = constq*dx/2.0
      ql(i) = constq*dx/2.0
    enddo
  elseif(qtype=='func1') then !quadratic source
    do i=1,numcells
      qr(i) = quadqr(a,x(i),x(i+1))
      ql(i) = quadql(a,x(i),x(i+1))
    enddo
  endif

  allocate(psiBCl(numangs/2))
  allocate(psiBCr(numangs/2))
  psiBCl = 0.0
  psiBCr = 0.0

  if(SnBCs=='vacuum') then
    do m=1,numangs/2
      psiBCl(m) = 0.0
      psiBCr(m) = 0.0
    enddo
  endif
  end subroutine



!----------------- Used In Sn and DSA -----------------!

  subroutine fluxinitdisc(  numcells,phi,phistart )
  !Initializes flux for diffusion and Sn
  integer :: numcells
  real(8) :: phistart
  real(8),allocatable :: phi(:)

  integer :: i

  allocate(phi(numcells))
  phi=0.0

  do i=1,numcells
    phi(i) = phistart
  enddo
  end subroutine fluxinitdisc


!----------------- Used Solely in DSA -----------------!

  subroutine DSAdiffsource( numcells,a,x,dx,q,sigt,siga,phiHsl,phiHsr,&
                               phiDSAlold,phiDSArold,iter )
  !Use flux from Sn solve to define source for diff solve as part of DSA
  !This is what couples the discontinuous flux values of Sn to the cont flux values of diff
  integer :: numcells,iter
  real(8) :: a,x(:),dx,sigt,siga
  real(8),allocatable :: q(:),phiHsl(:),phiHsr(:),phiDSAlold(:),phiDSArold(:)

  integer :: i
  real(8) :: sigs

  sigs = sigt - (siga/sigt)
  if(iter==1) allocate(q(numcells+1))
  q = 0.0

  !DSA source
  q(1)          = sigs*dx*( (phiHsl(1)       -phiDSAlold(1))       /3.0 + &
                            (phiHsr(1)       -phiDSArold(1))       /6.0 )
  do i=2,numcells
    q(i)        = sigs*dx*( (phiHsl(i-1)     -phiDSAlold(i-1))     /6.0 + &
                            (phiHsr(i-1)     -phiDSArold(i-1))     /3.0 + &
                            (phiHsl(i)       -phiDSAlold(i)  )     /3.0 + &
                            (phiHsr(i)       -phiDSArold(i)  )     /6.0 )
  enddo
  q(numcells+1) = sigs*dx*( (phiHsl(numcells)-phiDSAlold(numcells))/6.0 + &
                            (phiHsr(numcells)-phiDSArold(numcells))/3.0 )
  end subroutine DSAdiffsource



  subroutine alpha_beta( numangs,mu,wgts,alpha,beta )
  !Solves for angle constants alpha and beta used in SWLA part of DSA
  integer :: numangs
  real(8) :: mu(:),wgts(:),alpha,beta

  integer :: m
  alpha = 0.0
  beta  = 0.0
  do m=1,numangs
    if(mu(m)>0) then
      alpha = alpha + mu(m)      *wgts(m)
      beta  = beta  + mu(m)*mu(m)*wgts(m)
    endif
  enddo
  alpha = alpha/2.0
  beta  = beta /2.0 *3.0
  end subroutine alpha_beta


!----------------- Used In Diff and DSA -----------------!

  subroutine fluxinitcont(  numcells,phi,phistart )
  !Initializes flux for diffusion and Sn
  integer :: numcells
  real(8) :: phistart
  real(8),allocatable :: phi(:)

  integer :: i

  allocate(phi(numcells+1))
  phi=0.0

  do i=1,numcells+1
    phi(i) = phistart
  enddo
  end subroutine fluxinitcont


  !---------------------------------------------------------------------------------
  !
  !  Here begins what was before in the FEDiff1D module
  !
  !---------------------------------------------------------------------------------


  subroutine FEDiff(   numcells,sigt,c,a,siga,dx,phi,q,alpha,beta )
  !Loads and inverts B upon q to solve for phi in 1D Diffusion,
  !currently with Marshak BCs, using the Thomas algorithm for tri-diagonal solves
  integer :: numcells
  real(8) :: sigt,siga,c,a,dx,alpha,beta,phi(:),q(:)

  real(8),allocatable :: BB(:,:)

  integer :: i,j
  allocate(BB(numcells+1,numcells+1))
  BB = 0.0
  BB(1,1)                   = siga*dx/3.0     + 1.0/(3.0*sigt*dx) + alpha/beta !Marshak vac BCs
  BB(1,2)                   = siga*dx/6.0     - 1.0/(3.0*sigt*dx)              !a/b=0.5 for diff
  do i=2,numcells
    BB(i,i-1)               = siga*dx/6.0     - 1.0/(3.0*sigt*dx)
    BB(i,i)                 = siga*dx/3.0*2.0 + 2.0/(3.0*sigt*dx)
    BB(i,i+1)               = siga*dx/6.0     - 1.0/(3.0*sigt*dx)
  enddo
  BB(numcells+1,numcells)   = siga*dx/6.0     - 1.0/(3.0*sigt*dx)              !a/b~0.5 for DSA
  BB(numcells+1,numcells+1) = siga*dx/3.0     + 1.0/(3.0*sigt*dx) + alpha/beta !Marshak vac BCs


  !Time to invert same as phi = B \ q
  call Thomas_alg_interface( numcells,phi,BB,q )

  end subroutine FEDiff





  subroutine Thomas_alg_interface( numcells,phi,BB,q )
  !This subroutine allows me to use my previously verified Thomas algorithm,
  !but also keep input as general as possible for later use
  !with other methods of inversion
  use utilities, only: Thomas_alg
  integer :: numcells
  real(8) :: phi(:),BB(:,:),q(:)

  integer :: i,n
  real(8),allocatable :: a(:),b(:),c(:),d(:),x(:)

  n=numcells+1
  allocate(a(n))
  allocate(b(n))
  allocate(c(n))
  allocate(d(n))
  allocate(x(n))
  a = 0.0
  b = 0.0
  c = 0.0
  d = 0.0
  x = 0.0

  b(1) = BB(1,1)
  c(1) = BB(1,2)
  d(1) = q(1)
  x(1) = phi(1)
  do i=2,n-1
    a(i) = BB(i,i-1)
    b(i) = BB(i,i)
    c(i) = BB(i,i+1)
    d(i) = q(i)
    x(i) = phi(i)
  enddo
  a(n) = BB(n,n-1)
  b(n) = BB(n,n)
  d(n) = q(n)
  x(n) = phi(n)

  call Thomas_alg( n,a,b,c,d,x )

  do i=1,numcells+1
    phi(i) = x(i)
  enddo
  end subroutine Thomas_alg_interface


  !---------------------------------------------------------------------------------
  !
  !  Here begins what was before in the FESn1D module
  !
  !---------------------------------------------------------------------------------

  subroutine FESn( numcells,phiSnl,phiSnr,psil,psir,x,dx,qr,ql,sigt,siga,mu,&
                   psiBCl,psiBCr,numangs )
  integer :: numcells,numangs
  real(8) :: phiSnl(:),phiSnr(:),psil(:,:),psir(:,:),x(:),dx,qr(:),ql(:),sigt,siga,mu(:)
  real(8) :: psiBCl(:),psiBCr(:)

  integer :: i,m,j,mm
  real(8) :: A,B,C,D,Sr,Sl,sigs

  sigs = sigt-siga
  do m=1,numangs
    if(mu(m)>0.0) then
      A =  mu(m)/2.0 + sigt*dx/3.0
      B =  mu(m)/2.0 + sigt*dx/6.0
      C = -mu(m)/2.0 + sigt*dx/6.0
      D =  mu(m)/2.0 + sigt*dx/3.0
      do i=1,numcells
        if(i==1) mm=ceiling(real(m,8)/2.0)
        if(i==1) Sl=sigs*dx/2.0*(1.0/3.0*phiSnl(i)+1.0/6.0*phiSnr(i)) + 0.5*ql(i) + mu(m)*psiBCl(mm)
        if(i/=1) Sl=sigs*dx/2.0*(1.0/3.0*phiSnl(i)+1.0/6.0*phiSnr(i)) + 0.5*ql(i) + mu(m)*psir(m,i-1)
        Sr = sigs*dx/2.0*(1.0/6.0*phiSnl(i)+1.0/3.0*phiSnr(i)) + 0.5*qr(i)
        psir(m,i) = (Sr - C*Sl/A) / (-C*B/A + D)
        psil(m,i) = (Sl - B*psir(m,i)) / A
      enddo
    elseif(mu(m)<0.0) then
      A = -mu(m)/2.0 + sigt*dx/3.0
      B =  mu(m)/2.0 + sigt*dx/6.0
      C = -mu(m)/2.0 + sigt*dx/6.0
      D = -mu(m)/2.0 + sigt*dx/3.0
      do i=1,numcells
        j=numcells-i+1
        Sl = sigs*dx/2.0*(1.0/3.0*phiSnl(j)+1.0/6.0*phiSnr(j)) + 0.5*ql(j)
        if(j==numcells) mm=ceiling(real(m,8)/2.0)
        if(j==numcells) Sr=sigs*dx/2.0*(1.0/6.0*phiSnl(j)+1.0/3.0*phiSnr(j))&
                           + 0.5*qr(j) - mu(m)*psiBCr(mm)
        if(j/=numcells) Sr=sigs*dx/2.0*(1.0/6.0*phiSnl(j)+1.0/3.0*phiSnr(j))&
                           + 0.5*qr(j) - mu(m)*psil(m,j+1)
        psir(m,j) = (Sr - C*Sl/A) / (-C*B/A + D)
        psil(m,j) = (Sl - B*psir(m,j)) / A
      enddo
    endif
  enddo

  end subroutine FESn



  subroutine FESncalcphi( numcells,phiSnl,phiSnr,psil,psir,wgts,numangs )
  integer :: numcells,numangs
  real(8) :: phiSnl(:),phiSnr(:),psil(:,:),psir(:,:),wgts(:)

  integer :: i,m

  phiSnl = 0.0
  phiSnr = 0.0

  do i=1,numcells
    do m=1,numangs
      phiSnl(i) = phiSnl(i) + psil(m,i) * wgts(m)
      phiSnr(i) = phiSnr(i) + psir(m,i) * wgts(m)
    enddo
  enddo
  end subroutine FESncalcphi


  subroutine FESnerror( numcells,phiSnl,phiSnr,phiSnlold,phiSnrold,error )
  !Checks error in norm of flux at current iteration compared to last
  integer :: numcells
  real(8) :: error
  real(8) :: phiSnl(:),phiSnr(:),phiSnlold(:),phiSnrold(:)

  integer :: i
  real(8) :: rightsum,leftsum

  rightsum = 0.0
  leftsum  = 0.0
  do i=1,numcells
    rightsum = rightsum + (phiSnr(i) - phiSnrold(i))*(phiSnr(i) - phiSnrold(i))
    leftsum  = leftsum  + (phiSnl(i) - phiSnlold(i))*(phiSnl(i) - phiSnlold(i))
  enddo
  error = sqrt( (rightsum + leftsum) / numcells )
  end subroutine FESnerror



  subroutine FESnNewToOld( numcells,phiSnl,phiSnr,phiSnlold,phiSnrold )
  !Moves new flux values to "old" flux values, keeping for next error calc
  integer :: numcells
  real(8) :: phiSnl(:),phiSnr(:),phiSnlold(:),phiSnrold(:)

  integer :: i

  do i=1,numcells !save current as old
    phiSnlold(i) = phiSnl(i)
    phiSnrold(i) = phiSnr(i)
  enddo
  end subroutine FESnNewToOld


  !---------------------------------------------------------------------------------
  !
  !  Here begins what was before in the FEDSA1D module
  !
  !---------------------------------------------------------------------------------

  subroutine SWLA( numcells,phiDSAdiff,phiDSAlold,phiDSArold,phiHsl,phiHsr,&
                   phiresl,phiresr,sigt,siga,dx,alpha,beta )
  !Uses different definitions of current to combine half-step discontinuous values
  !and diffusion continuous values resulting the residual (Simplified 
  integer :: numcells
  real(8) :: sigt,siga,dx,alpha,beta
  real(8) :: phiDSAdiff(:),phiDSAlold(:),phiDSArold(:),phiHsl(:),phiHsr(:),phiresl(:),phiresr(:)

  integer :: i
  real(8) :: sigs,F,G,H,K,El,Er

  sigs = sigt*(sigt-siga)
  F = siga*dx/3.0 + alpha/beta
  G = siga*dx/6.0
  H = siga*dx/6.0
  K = siga*dx/3.0 + alpha/beta

  do i=1,numcells
    El = sigs*dx*( (phiHsl(i)-phiDSAlold(i))/3.0 + (phiHsr(i)-phiDSArold(i))/6.0 ) + &
                   alpha/beta*phiDSAdiff(i)
    Er = sigs*dx*( (phiHsl(i)-phiDSAlold(i))/6.0 + (phiHsr(i)-phiDSArold(i))/3.0 ) + &
                   alpha/beta*phiDSAdiff(i+1)
    phiresr(i) = (Er - H/F*El)/(K-G*H/F)
    phiresl(i) = (El - G*phiresr(i))/F
  enddo



  end subroutine SWLA





  subroutine HalfResSum( numcells,phiDSAl,phiDSAr,phiresl,phiresr,phiHsl,phiHsr ) !hs + res
  !Adds half-step flux values and residual for new flux values
  integer :: numcells
  real(8) :: phiDSAl(:),phiDSAr(:),phiresl(:),phiresr(:),phiHsl(:),phiHsr(:)

  integer :: i

  do i=1,numcells
    phiDSAl(i) = phiresl(i) + phiHsl(i)
    phiDSAr(i) = phiresr(i) + phiHsr(i)
  enddo
  end subroutine HalfResSum




  subroutine FEDiffSn_internaldeallocate( phipos,q,x,phiSnlold, &
                                 phiSnrold,psil,psir,psiBCl,psiBCr,qr,ql,mu,wgts, &
                                 phiHsl,phiHsr,phiresl,phiresr, &
                                 phiDSAlold,phiDSArold,phiDSAdiff )

  real(8),allocatable :: phipos(:),q(:),                x(:)
  real(8),allocatable :: phiSnlold(:),phiSnrold(:)
  real(8),allocatable :: psil(:,:),psir(:,:),psiBCl(:),psiBCr(:),qr(:),ql(:),mu(:),wgts(:)
  real(8),allocatable :: phiHsl(:),phiHsr(:),phiresl(:),phiresr(:)
  real(8),allocatable :: phiDSAlold(:),phiDSArold(:),phiDSAdiff(:)

  if(allocated(phipos)) deallocate(phipos)
  if(allocated(q)) deallocate(q)
  if(allocated(x)) deallocate(x)
  if(allocated(phiSnlold)) deallocate(phiSnlold)
  if(allocated(phiSnrold)) deallocate(phiSnrold)
  if(allocated(psil)) deallocate(psil)
  if(allocated(psir)) deallocate(psir)
  if(allocated(psiBCl)) deallocate(psiBCl)
  if(allocated(psiBCr)) deallocate(psiBCr)
  if(allocated(qr)) deallocate(qr)
  if(allocated(ql)) deallocate(ql)
  if(allocated(mu)) deallocate(mu)
  if(allocated(wgts)) deallocate(wgts)
  if(allocated(phiHsl)) deallocate(phiHsl)
  if(allocated(phiHsr)) deallocate(phiHsr)
  if(allocated(phiresl)) deallocate(phiresl)
  if(allocated(phiresl)) deallocate(phiresl)
  if(allocated(phiDSAlold)) deallocate(phiDSAlold)
  if(allocated(phiDSArold)) deallocate(phiDSArold)
  if(allocated(phiDSAdiff)) deallocate(phiDSAdiff)

  end subroutine FEDiffSn_internaldeallocate



  subroutine FEDiffSn_externaldeallocate
  !Deallocate module (not pass-by-reference) variables
  if(allocated(solve)) deallocate(solve)
  if(allocated(phidiff)) deallocate(phidiff)
  if(allocated(phiSnl)) deallocate(phiSnl)
  if(allocated(phiSnr)) deallocate(phiSnr)
  if(allocated(phiDSAl)) deallocate(phiDSAl)
  if(allocated(phiDSAr)) deallocate(phiDSAr)
  end subroutine FEDiffSn_externaldeallocate

  subroutine setflvarspassedtrue
  !For use when module is used by an external solver
  flvarspassed = .true.
  end subroutine setflvarspassedtrue



end module FEDiffSn
