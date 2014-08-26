module radtransMC
  use mcnp_random
  use utilities
  use timeman
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine radtrans_MCsim( j,o )
  use timevars, only: time
  use genRealzvars, only: sig, scatrat, numRealz, nummatSegs, matType, matLength, s
  use MCvars, only: numParts, radtrans_int, pfnumcells, rodOrplanar, sourceType, &
                    plotflux, pltflux, reflect, transmit, absorb, initcur, fluxfaces, &
                    fflux, bflux, flux
  integer  :: j,o
  real(8)  :: tt1,tt2

  integer  :: i,z
  real(8)  :: mu,db,dc,position,sc_ab,sigma

  call cpu_time(tt1)
  if( j==1 ) then
    allocate(transmit(numRealz))
    allocate(reflect(numRealz))
    allocate(absorb(numRealz))
    allocate(initcur(numRealz))
    transmit     =0.0d0
    reflect      =0.0d0
    absorb       =0.0d0
    initcur      =0.0d0
    radtrans_int =0
  endif
  do o=1,numParts !o-loop, through per particle

    if( sourceType=='left' ) then
      position = 0.0d0
      i        = 1
      mu       = 1.0d0
      if(rodOrplanar=='planar') mu = isoboundmu()
    elseif( sourceType=='intern') then
      position = s * rang()
      i        = internal_init_i(position,nummatSegs)
      mu       = merge(1.0d0,-1.0d0,rang()>0.5d0)
      if(rodOrplanar=='planar') mu = newmu()
    endif

                                initcur(j) = initcur(j) + 1d0    !initial current 
    do            ! through per part interaction
      radtrans_int=radtrans_int+1 !used to calc num of interactions

      db = merge(matLength(i+1)-position,position-matLength(i),mu>=0)/abs(mu)  !calculate db

      dc=-log(1-rang())/sig(matType(i))      !calculate dc
      sc_ab=rang()
      if( dc>db ) then                                  !escape
        if( mu>=0 ) then
          i=i+1
!print *,"particle: ",o,"  transmit"
          if(plotflux(2)=='tot') call adv_pos_col_flux(position,matLength(i),flux,j,mu)
          if(plotflux(2)=='fb') call col_fbflux(position,matLength(i),fflux,bflux,j,mu)
          if(i==nummatSegs+1) transmit(j)=transmit(j)+1      !transmit
          if(i==nummatSegs+1) exit
        else
!print *,"particle: ",o,"  reflect"
          if(plotflux(2)=='tot') call adv_pos_col_flux(position,matLength(i),flux,j,mu)
          if(plotflux(2)=='fb') call col_fbflux(position,matLength(i),fflux,bflux,j,mu)
          if(i==1)              reflect(j)=reflect(j)+1        !reflect
          if(i==1)              exit
          i=i-1
        endif
      elseif( sc_ab<scatrat(matType(i)) ) then          !scatter
!print *,"particle: ",o,"  scatter"
        if(plotflux(2)=='tot') call adv_pos_col_flux(position,position+dc*mu,flux,j,mu)
        if(plotflux(2)=='fb') call col_fbflux(position,position+dc*mu,fflux,bflux,j,mu)
        if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0) !dir of scatter
        if(rodOrplanar=='planar') mu = newmu()
      else
!print *,"particle: ",o,"  absorb"
        if(plotflux(2)=='tot') call adv_pos_col_flux(position,position+dc*mu,flux,j,mu)
                                absorb(j)=absorb(j)+1.0d0 !absorb
        if(plotflux(2)=='fb') call col_fbflux(position,position+dc*mu,fflux,bflux,j,mu)
                                exit
      endif

    enddo

  enddo

  call cpu_time(tt2)
  time(2) = time(2) + (tt2-tt1)

  end subroutine radtrans_MCsim








  subroutine radtrans_MCoutstats
  use genRealzvars, only: Adamscase, sig, scatrat, lam, s, numRealz, P
  use MCvars, only: numParts, radtrans_int, pfnumcells, rodOrplanar, plotflux, &
                    results, pltflux, reflect, transmit, absorb, initcur, fluxfaces, &
                    fflux, bflux, flux

  integer :: j,i
  real(8) :: reflection,transmission,absorption
  real(8) :: reflectionvar,transmissionvar,absorptionvar
  real(8),allocatable :: fluxave(:),fluxvar(:),fluxinput(:)
  real(8),allocatable :: fluxavemat1(:),fluxvarmat1(:),fluxavemat2(:),fluxvarmat2(:)

  333 format("reflection        :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  334 format("transmission      :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  335 format("absorption        :   ",f8.5," +-    ",f8.5," +-  ",f10.7)
  336 format(f4.1,"     ",i9,"     ",i9,"     ",f10.7,"     ",f10.7)
  337 format("     ",f10.7,"     ",f10.7,"     ",f10.7,"     ",f10.7)
  338 format("     ",f10.7,"     ",f10.7,"     ",f10.7)
  339 format("absorption        :   ",f8.5)

  do j=1,numRealz
    reflect(j) = reflect(j) / initcur(j)
    transmit(j)= transmit(j)/ initcur(j)
  enddo

  call mean_and_var_s( reflect,numRealz,reflection,reflectionvar )
  call mean_and_var_s( transmit,numRealz,transmission,transmissionvar )
  call mean_and_var_s( absorb/numParts,numRealz,absorption,absorptionvar )
  if(rodOrplanar=='planar') absorption = 1.0d0 - reflection - transmission

  write(*,*)
  write(*,*) "--- Transport: radtrans ---"
  write(*,*) "                        mean           stdev         relerr   "
  write(*,333) reflection,sqrt(reflectionvar),sqrt(reflectionvar)/reflection
  write(*,334) transmission,sqrt(transmissionvar),sqrt(transmissionvar)/transmission
  if(rodOrplanar=='rod') write(*,335)    absorption,sqrt(absorptionvar),&
                                                    sqrt(absorptionvar)/absorption
  if(rodOrplanar=='planar') write(*,339) absorption

  340 format(" interactions     :  ",f5.1," % accept  ",f5.2," ac/p    ",f5.2," rj/p")

  write(*,340) 1.0d0*100,real(radtrans_int,8)/numParts/numRealz,0.0d0


  open(unit=16, file="radtrans.txt")
  341 format(f5.2,"   ",f10.4,f10.4,f10.4,f10.4,"   ",f5.2)

  write(16,*) "#Adamscase refl      refldev   trans     transdev  plot"
  write(16,341) Adamscase+0.01d0,reflection,sqrt(reflectionvar),&
                transmission,sqrt(transmissionvar),Adamscase+0.07d0
  close(unit=16)
  call system("mv radtrans.txt plots")



  !print to file in order to plot
  open(unit=2, file="transout.txt")
  write(2,336,advance="no") Adamscase,numRealz,numParts,reflection,reflectionvar
  write(2,337,advance="no") transmission,transmissionvar,absorption,absorptionvar
  write(2,337,advance="no") sig(1),sig(2),scatrat(1),scatrat(2)
  write(2,338,advance="no") lam(1),lam(2),s
  write(2,*)
  CALL system("mv transout.txt texts")
  if( results=='remove' ) then
    CALL system("cat texts/transout_key.txt texts/transout.txt &
                 > texts/transout_results.txt")
    write(*,*) "Radtrans results written over (texts/transout_results.txt)"
  elseif( results=='add' ) then
    CALL system("cat texts/transout_results.txt texts/transout.txt &
                 > texts/transout_results2.txt")
    CALL system("mv texts/transout_results2.txt texts/transout_results.txt")
    write(*,*) "Radtrans results added to old (texts/transout_results.txt)"
  endif


  !stats for flux
  if(pltflux(1)/='noplot') then

    if(plotflux(1)=='cell') then  !for cell flux calculations
      allocate(fluxvar(pfnumcells))
      allocate(fluxave(pfnumcells))
      allocate(fluxinput(numRealz))

      if(plotflux(2)=='tot') then
        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = flux(j,i)/numParts
          enddo
          call mean_and_var_s( fluxinput,numRealz,fluxave(i),fluxvar(i) )
        enddo
      endif

      if(plotflux(2)=='fb') then
        allocate(fluxvarmat1(pfnumcells))
        allocate(fluxvarmat2(pfnumcells))
        allocate(fluxavemat1(pfnumcells))
        allocate(fluxavemat2(pfnumcells))

        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = fflux(j,i)/numParts/P(1)
          enddo
          call mean_and_var_s( fluxinput,numRealz,fluxavemat1(i),fluxvarmat1(i) )
        enddo

        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = bflux(j,i)/numParts/P(2)
          enddo
          call mean_and_var_s( fluxinput,numRealz,fluxavemat2(i),fluxvarmat2(i) )
        enddo
do i=1,4
  do j=1,numRealz
    10000 format("fflux(",i3,",",i6,"): ",f18.2,"  bflux(",i3,",",i6,"): ",f18.2)
    write(*,10000) j,i,fflux(j,i),j,i,bflux(j,i)
  enddo
enddo


!print *,"fluxavemat1: ",fluxavemat1
!print *,"fluxavemat2: ",fluxavemat2

        do i=1,pfnumcells
          do j=1,numRealz
            fluxinput(j) = fflux(j,i)/numParts + bflux(j,i)/numParts
          enddo
          call mean_and_var_s( fluxinput,numRealz,fluxave(i),fluxvar(i) )
        enddo
      endif

    
      350 format("#cell center,       ave flux,        flux dev")
      351 format(f15.7,f15.7,f15.7)
      open(unit=22, file="radMC_cellflux.txt")
      write(22,350)
      do i=1,pfnumcells
        write(22,351) (fluxfaces(i+1)+fluxfaces(i))/2.0d0,fluxave(i),sqrt(fluxvar(i))
      enddo    
      close(unit=22)

      if(plotflux(2)=='fb') then
        360 format("#cell center,    ave flux mat1,   flux dev mat1,    ave flux mat2,   flux dev mat2")
        361 format(f15.7,f15.7,f15.7,f15.7,f15.7)
        open(unit=23, file="radMC_fbcellflux.txt")
        write(23,360)
        do i=1,pfnumcells
          write(23,361) (fluxfaces(i+1)+fluxfaces(i))/2.0d0,fluxavemat1(i),sqrt(fluxvarmat1(i)),&
                                                          fluxavemat2(i),sqrt(fluxvarmat2(i))
        enddo    
        close(unit=23)
      endif

      deallocate(fluxvar)
      deallocate(fluxave)
      deallocate(fluxinput)
      if(plotflux(2)=='fb') then
        deallocate(fluxvarmat1)
        deallocate(fluxvarmat2)
        deallocate(fluxavemat1)
        deallocate(fluxavemat2)
      endif

      call system("mv radMC_cellflux.txt plots")
      if(plotflux(2)=='fb') call system("mv radMC_fbcellflux.txt plots")
    endif



  endif

  end subroutine radtrans_MCoutstats






  subroutine initialize_fluxplot
  use genRealzvars, only: s, numRealz
  use KLvars, only: KLrnumRealz
  use MCvars, only: pfnumcells, plotflux, fluxfaces, fflux, bflux, fradWoodf, &
                    bradWoodf, fKLWoodf, bKLWoodf, radWoodf, KLWoodf, Woodf, &
                    flux, radMC, radWood, KLWood

  integer :: i
  allocate(fluxfaces(pfnumcells+1))
  fluxfaces = 0.0d0
  do i=1,pfnumcells+1
    fluxfaces(i) = (s/pfnumcells) * (i-1)
  enddo
  if(radMC=='yes') then
    !i.e. different length for surface flux tally
    if(plotflux(1)=='cell' .AND. plotflux(2)=='tot') allocate(flux(numRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='tot') flux  = 0.0d0
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb')  allocate(fflux(numRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb')  fflux = 0.0d0
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb')  allocate(bflux(numRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb')  bflux = 0.0d0
  endif

  if(radWood=='yes') then
    if(plotflux(1)=='cell' .AND. plotflux(2)=='tot') allocate(radWoodf(numRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='tot') radWoodf  = 0.0d0
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') allocate(fradWoodf(numRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') fradWoodf = 0.0d0
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') allocate(bradWoodf(numRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') bradWoodf = 0.0d0
  endif

  if(KLWood=='yes') then
    if(plotflux(1)=='cell' .AND. plotflux(2)=='tot') allocate(KLWoodf(KLrnumRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='tot') KLWoodf  = 0.0d0
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') allocate(fKLWoodf(KLrnumRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') fKLWoodf = 0.0d0
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') allocate(bKLWoodf(KLrnumRealz,pfnumcells))
    if(plotflux(1)=='cell' .AND. plotflux(2)=='fb') bKLWoodf = 0.0d0
  endif

  end subroutine initialize_fluxplot





  subroutine plot_flux
  use MCvars, only: plotflux, pltflux, radMC, radWood, KLWood

  if(plotflux(1)=='cell') then

    if(pltflux(1)=='preview') then
      if(radMC=='yes' .AND. radWood=='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM_rW_KW.gnu")
      if(radMC=='yes' .AND. radWood=='yes' .AND. KLWood/='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM_rW.gnu")
      if(radMC=='yes' .AND. radWood/='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM_KW.gnu")
      if(radMC/='yes' .AND. radWood=='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rW_KW.gnu")
      if(radMC=='yes' .AND. radWood/='yes' .AND. KLWood/='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM.gnu")
      if(radMC/='yes' .AND. radWood=='yes' .AND. KLWood/='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rW.gnu")
      if(radMC/='yes' .AND. radWood/='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_KW.gnu")
    elseif(pltflux(1)=='plot') then
      if(radMC=='yes' .AND. radWood=='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM_rW_KW.p.gnu")
      if(radMC=='yes' .AND. radWood=='yes' .AND. KLWood/='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM_rW.p.gnu")
      if(radMC=='yes' .AND. radWood/='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM_KW.p.gnu")
      if(radMC/='yes' .AND. radWood=='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rW_KW.p.gnu")
      if(radMC=='yes' .AND. radWood/='yes' .AND. KLWood/='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rM.p.gnu")
      if(radMC/='yes' .AND. radWood=='yes' .AND. KLWood/='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_rW.p.gnu")
      if(radMC/='yes' .AND. radWood/='yes' .AND. KLWood=='yes') &
        call system("gnuplot plots/fluxgnus/cellflux_KW.p.gnu")
    endif

    if(pltflux(1)=='preview' .AND. plotflux(2)=='fb') then
      if(radMC=='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat1fluxgnus/cellflux_rM_rW.gnu")
      if(radMC=='yes' .AND. radWood/='yes') &
        call system("gnuplot plots/fluxgnus/mat1fluxgnus/cellflux_rM.gnu")
      if(radMC/='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat1fluxgnus/cellflux_rW.gnu")
    elseif(pltflux(1)=='plot' .AND. plotflux(2)=='fb') then
      if(radMC=='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat1fluxgnus/cellflux_rM_rW.p.gnu")
      if(radMC=='yes' .AND. radWood/='yes') &
        call system("gnuplot plots/fluxgnus/mat1fluxgnus/cellflux_rM.p.gnu")
      if(radMC/='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat1fluxgnus/cellflux_rW.p.gnu")
    endif

    if(pltflux(1)=='preview' .AND. plotflux(2)=='fb') then
      if(radMC=='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat2fluxgnus/cellflux_rM_rW.gnu")
      if(radMC=='yes' .AND. radWood/='yes') &
        call system("gnuplot plots/fluxgnus/mat2fluxgnus/cellflux_rM.gnu")
      if(radMC/='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat2fluxgnus/cellflux_rW.gnu")
    elseif(pltflux(1)=='plot' .AND. plotflux(2)=='fb') then
      if(radMC=='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat2fluxgnus/cellflux_rM_rW.p.gnu")
      if(radMC=='yes' .AND. radWood/='yes') &
        call system("gnuplot plots/fluxgnus/mat2fluxgnus/cellflux_rM.p.gnu")
      if(radMC/='yes' .AND. radWood=='yes') &
        call system("gnuplot plots/fluxgnus/mat2fluxgnus/cellflux_rW.p.gnu")
    endif

    call system("ps2pdf cellflux.ps cellflux.pdf")
    call system("ps2eps cellflux.ps")

    call system("mv cellflux.ps  plots/cellflux.ps")
    call system("mv cellflux.pdf plots/cellflux.pdf")
    call system("rm plots/cellflux.eps")
    call system("mv cellflux.eps plots/cellflux.eps")

    if(plotflux(2)=='fb') then
      call system("ps2pdf mat1cellflux.ps mat1cellflux.pdf")
      call system("ps2eps mat1cellflux.ps")

      call system("mv mat1cellflux.ps  plots/mat1cellflux.ps")
      call system("mv mat1cellflux.pdf plots/mat1cellflux.pdf")
      call system("rm plots/mat1cellflux.eps")
      call system("mv mat1cellflux.eps plots/mat1cellflux.eps")

      call system("ps2pdf mat2cellflux.ps mat2cellflux.pdf")
      call system("ps2eps mat2cellflux.ps")

      call system("mv mat2cellflux.ps  plots/mat2cellflux.ps")
      call system("mv mat2cellflux.pdf plots/mat2cellflux.pdf")
      call system("rm plots/mat2cellflux.eps")
      call system("mv mat2cellflux.eps plots/mat2cellflux.eps")
    endif

  endif

  end subroutine plot_flux









  function newmu()
  !use this for sampling a new mu within a slab
  real(8) :: newmu
 
  newmu = 2*rang() - 1 

  end function newmu



  function isoboundmu()
  !use this for an isotropic source incident on left boundary
  real(8) :: isoboundmu

  isoboundmu = sqrt(rang())

  end function isoboundmu



  function internal_init_i( position,numArrSz )
  use genRealzvars, only: matLength
  integer :: numArrSz
  real(8) :: position

  integer :: i,internal_init_i

  do i=1,numArrSz
    if( matLength(i)<=position .AND. matLength(i+1)>position ) then
      internal_init_i = i
      exit
    endif
  enddo


  end function internal_init_i



!! Radtrans and Woodcock funcs and subs


  subroutine adv_pos_col_flux( oldpos,newpos,flux,j,mu )
  use MCvars, only: pfnumcells, plotflux, pltflux, fluxfaces
  integer :: j
  real(8) :: oldpos,newpos,mu,absmu
  real(8) :: flux(:,:)

  integer :: i
  real(8) :: smallpos,larpos,dx
  character(6) :: status !'before','middle'

  if(pltflux(1)/='noplot') then
    status   = 'before'
    smallpos = merge(oldpos,newpos,oldpos<newpos)
    larpos   = merge(oldpos,newpos,oldpos>newpos)

!print *, "oldpos  : ",oldpos,"     newpos: ",newpos
!print *,"mu          : ",mu
!print *,"fluxfaces(1): ",fluxfaces(1),"   fluxfaces(2): ",fluxfaces(2)
!print *,"flux(j,1)   : ",flux(j,1),"    flux(j,2)   : ",flux(j,2)
    absmu    = abs(mu)
    dx       = fluxfaces(2)-fluxfaces(1)
    do i=1,pfnumcells

      if(plotflux(1)=='cell') then  !cell flux calculation, so far only one
        if(status=='before') then
          if(fluxfaces(i)<=smallpos .AND. smallpos<fluxfaces(i+1)) then
            if(larpos<=fluxfaces(i+1)) then
              flux(j,i)=flux(j,i)+(larpos-smallpos)*1.0d0/absmu/dx
              exit
            else
              flux(j,i)=flux(j,i)+(fluxfaces(i+1)-smallpos)*1.0d0/absmu/dx
              status = 'middle'
            endif
          endif
        elseif(status=='middle') then
          if(fluxfaces(i+1)<=larpos) then
              flux(j,i)=flux(j,i)+1.0d0/absmu
          else
              flux(j,i)=flux(j,i)+(larpos-fluxfaces(i))*1.0d0/absmu/dx
              exit
          endif
        endif
      endif 

    enddo
  endif
!print *,"flux(j,1)   : ",flux(j,1),"    flux(j,2)   : ",flux(j,2)
!print *
  oldpos = newpos !set 'position' in the outer loops to the new position

  end subroutine adv_pos_col_flux



  subroutine col_fbflux( oldpos,newpos,fflux,bflux,j,mu )
  !tallies flux contribution in each material withing fluxface bins 
  use genRealzvars, only: matType, matLength
  use MCvars, only: pfnumcells, plotflux, pltflux, fluxfaces
  integer :: j
  real(8) :: fflux(:,:),bflux(:,:)
  real(8) :: oldpos,newpos,mu,absmu

  integer :: i,k
  real(8) :: smallpos,larpos,dx,fhit,mhit,lasthit
  logical :: endflag
  character(6) :: status !'before','middle'
  if(pltflux(1)/='noplot') then
    endflag  = .false.
    status   = 'before'
    smallpos = merge(oldpos,newpos,oldpos<newpos)
    larpos   = merge(oldpos,newpos,oldpos>=newpos)

!print *, "oldpos  : ",oldpos,"     newpos: ",newpos
!print *,"mu          : ",mu
!print *,"fluxfaces(1): ",fluxfaces(1),"   fluxfaces(2): ",fluxfaces(2)
!print *,"flux(j,1)   : ",flux(j,1),"    flux(j,2)   : ",flux(j,2)
    absmu    = abs(mu)
    dx       = fluxfaces(2)-fluxfaces(1)
!    do i=1,pfnumcells
      if(plotflux(1)=='cell') then  !cell flux calculation, so far only one
        k=0 !find initial fluxface
        do
          k=k+1
          if(fluxfaces(k)>smallpos) then
            fhit = fluxfaces(k)
            exit
          endif
        enddo
        i=0 !find initial material face
        do
          i=i+1
          if(matLength(i)>smallpos) then
            mhit = matLength(i)
            exit
          endif
        enddo
!print *,"matLength()s: ",matLength(1)," ",matLength(2)," ",matLength(3)," ",matLength(4)
!print *,"fluxfaces   : ",fluxfaces(1)," ",fluxfaces(2)," ",fluxfaces(3)
!print *,"oldpos      : ",oldpos," newpos: ",newpos
!print *,"fhit        : ",fhit,"   mhit: ",mhit,"   smallpos: ",smallpos,"  larpos: ",larpos
        lasthit = smallpos
!print *
!print *
!print *,"smallpos:",smallpos,"  larpos:",larpos
!print *
        do
!print *,"    fhit:",fhit,"    mhit:",mhit
          if(fhit<mhit) then
            if(matType(i-1)==1) fflux(j,k-1) = fflux(j,k-1) + (fhit-lasthit)/absmu/dx
            if(matType(i-1)==2) bflux(j,k-1) = bflux(j,k-1) + (fhit-lasthit)/absmu/dx

!print *,"fhit: ",fhit,"    fflux(j,k-1): ",fflux(j,k-1)
            lasthit=fhit

            if(fhit<larpos) then
              k=k+1
              fhit=fluxfaces(k)
            elseif(fhit>=larpos) then
              fhit=larpos
            endif

          elseif(mhit<=fhit) then
            if(matType(i-1)==1) fflux(j,k-1) = fflux(j,k-1) + (mhit-lasthit)/absmu/dx
            if(matType(i-1)==2) bflux(j,k-1) = bflux(j,k-1) + (mhit-lasthit)/absmu/dx
            lasthit=mhit

            if(mhit<larpos) then
              i=i+1
              mhit=matLength(i)
            elseif(mhit>=larpos) then
              mhit=larpos
            endif

          endif
!print *,"fhit: ",fhit,"   larpos: ",larpos,"    mhit: ",mhit
        if(endflag) exit
        endflag = fhit .ge. larpos .AND. mhit .ge. larpos
!        if(fhit>=larpos .AND. mhit>=larpos) endflag = .true.
        enddo
      endif

    !enddo
  endif
!print *,"fflux(1,1)   : ",fflux(1,1),"    bflux(1,1)   : ",bflux(1,1)
!print *
  oldpos = newpos !set 'position' in the outer loops to the new position
!if(larpos==1.0d0) STOP
  end subroutine col_fbflux






  subroutine radtrans_resultplot
  !Plots pdfs of transmission and reflection for chosen methods
  use genRealzvars, only: numRealz
  use MCvars,       only: radMCbinplot, radWoodbinplot, KLWoodbinplot, reflect, radWoodr, &
                          KLWoodr, radWoodt, KLWoodt, transmit

  real(8) :: smrefl,lgrefl,smtran,lgtran,boundbuff

  !find reflection binning/plotting bounds
  smrefl = 1d0
  if(radMCbinplot  =='plot' .or. radMCbinplot  =='preview') smrefl = min(smrefl,minval(reflect))
  if(radWoodbinplot=='plot' .or. radWoodbinplot=='preview') smrefl = min(smrefl,minval(radWoodr))
  if(KLWoodbinplot =='plot' .or. KLWoodbinplot =='preview') smrefl = min(smrefl,minval(KLWoodr))
  lgrefl = 0d0
  if(radMCbinplot  =='plot' .or. radMCbinplot  =='preview') lgrefl = max(lgrefl,maxval(reflect))
  if(radWoodbinplot=='plot' .or. radWoodbinplot=='preview') lgrefl = max(lgrefl,maxval(radWoodr))
  if(KLWoodbinplot =='plot' .or. KLWoodbinplot =='preview') lgrefl = max(lgrefl,maxval(KLWoodr))
  boundbuff = (lgrefl-smrefl)/8d0
  smrefl = merge(smrefl-boundbuff,0d0,smrefl-boundbuff>0d0)
  lgrefl = merge(lgrefl+boundbuff,1d0,lgrefl+boundbuff<1d0)
  smrefl = smrefl - 0.0000001d0 !these to ensure binning works in case of opaque or transparent
  lgrefl = lgrefl + 0.0000001d0

  !find tranmission binning/plotting bounds
  smtran = 1d0
  if(radMCbinplot  =='plot' .or. radMCbinplot  =='preview') smtran = min(smtran,minval(transmit))
  if(radWoodbinplot=='plot' .or. radWoodbinplot=='preview') smtran = min(smtran,minval(radWoodt))
  if(KLWoodbinplot =='plot' .or. KLWoodbinplot =='preview') smtran = min(smtran,minval(KLWoodt))
  lgtran = 0d0
  if(radMCbinplot  =='plot' .or. radMCbinplot  =='preview') lgtran = max(lgtran,maxval(transmit))
  if(radWoodbinplot=='plot' .or. radWoodbinplot=='preview') lgtran = max(lgtran,maxval(radWoodt))
  if(KLWoodbinplot =='plot' .or. KLWoodbinplot =='preview') lgtran = max(lgtran,maxval(KLWoodt))  
  boundbuff = (lgtran-smtran)/8d0
  smtran = merge(smtran-boundbuff,0d0,smtran-boundbuff>0d0)
  lgtran = merge(lgtran+boundbuff,1d0,lgtran+boundbuff<1d0)
  smtran = smtran - 0.0000001d0 !these to ensure binning works in case of opaque or transparent
  lgtran = lgtran + 0.0000001d0

  !radMC binning and printing
  call system("rm plots/tranreflprofile/radMCtranreflprofile.txt")
  if(radMCbinplot  =='plot' .or. radMCbinplot  =='preview') then
    call radtrans_bin( smrefl,lgrefl,smtran,lgtran )
    call system("mv tranreflprofile.txt radMCtranreflprofile.txt")
    call system("mv radMCtranreflprofile.txt plots/tranreflprofile")
  endif

  !radWood binning and printing
  call system("rm plots/tranreflprofile/radWoodtranreflprofile.txt")
  if(radWoodbinplot=='plot' .or. radWoodbinplot=='preview') then
    call radtrans_bin( smrefl,lgrefl,smtran,lgtran )
    call system("mv tranreflprofile.txt radWoodtranreflprofile.txt")
    call system("mv radWoodtranreflprofile.txt plots/tranreflprofile")
  endif

  !KLWood binning and printing
  call system("rm plots/tranreflprofile/KLWoodtranreflprofile.txt")
  if(KLWoodbinplot =='plot' .or. KLWoodbinplot =='preview') then
    call radtrans_bin( smrefl,lgrefl,smtran,lgtran )
    call system("mv tranreflprofile.txt KLWoodtranreflprofile.txt")
    call system("mv KLWoodtranreflprofile.txt plots/tranreflprofile")
  endif

  !plot, convert, and store
  if(radMCbinplot=='preview' .or. radWoodbinplot=='preview' .or. KLWoodbinplot=='preview') then
    call system("gnuplot plots/tranreflprofile/tranprofile.p.gnu")
    call system("gnuplot plots/tranreflprofile/reflprofile.p.gnu")
  else
    call system("gnuplot plots/tranreflprofile/tranprofile.gnu")
    call system("gnuplot plots/tranreflprofile/reflprofile.gnu")
  endif
  call system("ps2pdf tranprofile.ps")
  call system("ps2pdf reflprofile.ps")
  call system("mv tranprofile.ps tranprofile.pdf plots/tranreflprofile")
  call system("mv reflprofile.ps reflprofile.pdf plots/tranreflprofile")

  end subroutine radtrans_resultplot



  subroutine radtrans_bin( smrefl,lgrefl,smtran,lgtran )
  !Heart of radtrans_resultplot, loads values to bin, and prints to generic text file
  !for each method
  use genRealzvars,         only: numRealz
  use MCvars,               only: trprofile_binnum, reflect, transmit
  real(8) :: smrefl,lgrefl,smtran,lgtran

  !local vars
  integer :: ibin
  integer,allocatable,dimension(:) :: reflcounts,trancounts  
  real(8),allocatable,dimension(:) :: reflprob,  tranprob
  real(8),allocatable,dimension(:) :: reflbounds,tranbounds

  !prepare variables
  if(allocated(reflcounts)) deallocate(reflcounts)
  if(allocated(trancounts)) deallocate(trancounts)
  if(allocated(reflprob))   deallocate(reflprob)
  if(allocated(tranprob))   deallocate(tranprob)
  if(allocated(reflbounds)) deallocate(reflbounds)
  if(allocated(tranbounds)) deallocate(tranbounds)

  if(.not. allocated(reflcounts)) allocate(reflcounts(trprofile_binnum))
  if(.not. allocated(trancounts)) allocate(trancounts(trprofile_binnum))
  if(.not. allocated(reflprob))   allocate(reflprob(trprofile_binnum))
  if(.not. allocated(tranprob))   allocate(tranprob(trprofile_binnum))
  if(.not. allocated(reflbounds)) allocate(reflbounds(trprofile_binnum+1))
  if(.not. allocated(tranbounds)) allocate(tranbounds(trprofile_binnum+1))
  reflcounts = 0
  trancounts = 0
  reflprob   = 0d0
  tranprob   = 0d0
  reflbounds = 0d0
  tranbounds = 0d0
  
  !actually store in binned fashion
  call store_in_bins( smrefl,lgrefl,trprofile_binnum,reflcounts,reflbounds,reflect,numRealz )
  call store_in_bins( smtran,lgtran,trprofile_binnum,trancounts,tranbounds,transmit,numRealz )

  tranprob = real(trancounts,8)/numRealz/((lgtran-smtran)/(trprofile_binnum-1))
  reflprob = real(reflcounts,8)/numRealz/((lgrefl-smrefl)/(trprofile_binnum-1))

  !print to generic file
  open(unit=100,file="tranreflprofile.txt")
  310 format("#         tranbounds      trancounts        reflbounds        relfcounts")
  311 format(f20.10,f20.10,f20.10,f20.10)
  write(100,310)
  do ibin=1,trprofile_binnum
    write(100,311) (tranbounds(ibin)+tranbounds(ibin+1))/2d0,tranprob(ibin),&
                   (reflbounds(ibin)+reflbounds(ibin+1))/2d0,reflprob(ibin)
  enddo
  close(unit=100)

  end subroutine radtrans_bin



end module radtransMC
