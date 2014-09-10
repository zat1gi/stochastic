module radtransMC
  use mcnp_random
  use utilities
  use timeman
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine UQ_MC( icase )
  !This subroutine perfoms Monte Carlo in the uncertain space, currently for binary mixtures.
  !'MCtransport' handles the spatial MC, but this subroutine collects data in UQ space
  use timevars, only: time
  use genRealzvars, only: numRealz
  use MCvars, only: MCcases, numParts, trannprt
  use genRealz, only: genReal
  integer :: icase

  integer :: j !which realization
  real(8) :: tt1

  call cpu_time(tt1)
  write(*,*) "Starting method: ",MCcases(icase)  
  do j=1,numRealz
!print *,"starting realization: ",j
    call genReal( j )                      !gen geometry
    if(MCcases(icase)=='radWood' .or. MCcases(icase)=='KLWood') call MCWood_setceils( j,icase )
                                           !for WMC, create ceilings

    call MCtransport( j,icase,numParts )   !transport over a realization

    if(mod( j,trannprt )==0) call radtrans_timeupdate( j,icase,tt1 )    !print time updates
  enddo !loops over realizations

  call stocMC_stats( icase )               !add routine to calc stats in stochastic space here
                                           !later make the above two loops,
                                           !batch spatial stats in between, final out here

  !add plotter here.


  end subroutine UQ_MC





  subroutine MCtransport( j,icase,tnumParts )
  !This subroutine performs MC transport over a specified geometry.
  !It can be used by a UQ_MC wrapper, and likely in the future by UQ_MLMC and/or UQ_SC.
  !'j' denotes which realization over which the MC transport is performed.
  !'icase' denotes which predefined setup over which the MC transport is performed.
  !'tnumParts' is the number of particles over which the MC transport is performed.
  !It may in the future make more sense to have flags which do or do not use certain functionality,
  !but for the time being each operation is simply chosen as a function of which case is selected.
  use genRealzvars, only: sig, scatrat, nummatSegs, matType, matLength, s
  use MCvars, only: radtrans_int, rodOrplanar, sourceType, reflect, transmit, &
                    absorb, position, oldposition, mu, areapnSamp, numpnSamp, &
                    nceilbin, Wood_rej, allowneg, distneg, MCcases, fbinmax, &
                    bbinmax, binmaxind, binmaxes

  use Woodcock, only: radWood_actsig, radWood_actscatrat, KLrxi_point
                      !above only here because of politics, rid when you finish merging!
  integer :: j,icase,tnumParts !realz number/which mode of transport/num of particles

  !local variables
  integer :: i,o
  real(8) :: db,dc,dist,  sigma,  ceilsig,woodrat,disthold
  character(9) :: fldist, flIntType, flEscapeDir, flExit


  do o=1,tnumParts                     !loop over particles
!print *," starting particle: ",o
    call genSourcePart( i,icase )      !gen source part pos, dir, and binnum (i)

    do ! simulate one pathlength of a particle
!print *,"   starting pathlength"
      flExit='clean'

      !tally number of interactions
      if(MCcases(icase)=='radMC') radtrans_int=radtrans_int+1


      !calculate distance to boundary
      select case (MCcases(icase))
        case ("radMC")
          db = merge(matLength(i+1)-position,position-matLength(i),mu>=0)/abs(mu)
        case ("radWood")
          db = merge(s-position,position,mu>=0)/abs(mu)
        case ("KLWood")
          db = merge(s-position,position,mu>=0)/abs(mu)
      end select


      !calculate distance to collision
      select case (MCcases(icase))
        case ("radMC")
          dc = -log(1-rang())/sig(matType(i))
        case ("radWood")
          ceilsig=merge(ceilsigfunc2(position,fbinmax),& !sel max sig
                        ceilsigfunc2(position,bbinmax),mu>=0)
          dc = -log(1-rang())/ceilsig                   !calc dc
        case ("KLWood")
          ceilsig=merge(ceilsigfunc2(position,fbinmax),& !sel max sig
                        ceilsigfunc2(position,bbinmax),mu>=0)
          dc = -log(1-rang())/ceilsig                   !calc dc
      end select



      !select distance limiter (add another later with LP)
      fldist = 'clean'
      dist   = min(db,dc)
      fldist = merge('boundary ','collision',db<dc)



      !if boundary chosen
      if(fldist=='boundary') then
        flEscapeDir = 'clean'
        !set direction flag
        flEscapeDir = merge('transmit','reflect ',mu>0.0d0)
        if(MCcases(icase)=='radWood' .or. MCcases(icase)=='KLWood') Wood_rej(1)=Wood_rej(1)+1 !accept path tal

        !evaluate for local or global transmission or reflection
        if(flEscapeDir=='transmit') then   !transmit
          select case (MCcases(icase))
            case ("radMC")
              i=i+1
              call MCinc_pos( matLength(i) )
              if(i==nummatSegs+1) transmit(j)=transmit(j)+1.0d0
              if(i==nummatSegs+1) flExit='exit'
            case ("radWood")
              call MCinc_pos( s )
              transmit(j) = transmit(j)
              flExit='exit'
            case ("KLWood")
              call MCinc_pos( s )
              transmit(j) = transmit(j)
              flExit='exit'
          end select
        endif

        if(flEscapeDir=='reflect') then    !reflect
          select case (MCcases(icase))
            case ("radMC")
              call MCinc_pos( matLength(i) )
              if(i==1)            reflect(j)=reflect(j)+1.0d0
              if(i==1)            flExit='exit'
              i=i-1
            case ("radWood")
              call MCinc_pos( 0.0d0 )
              reflect(j)  = reflect(j)+1
              flExit='exit'
            case ("KLWood")
              call MCinc_pos( 0.0d0 )
              reflect(j)  = reflect(j)+1
              flExit='exit'
          end select
        endif

      endif


      !if collision chosen
      if(fldist=='collision') then
        flIntType = 'clean'
        !Advance position for all outcomes
        call MCinc_pos( position + dc*mu )

        !Choose scatter, absorb, or reject interaction
        select case (MCcases(icase))
          case ("radMC")
            flIntType = merge('scatter','absorb ',rang()<scatrat(1))
          case ("radWood")
            woodrat = radWood_actsig(position,sig)/ceilsig
            if(woodrat>1.0d0) then
              !print *,"j: ",j,"  woodrat: ",woodrat
              stop 'Higher sig samples in radWood than ceiling, exiting program'
            endif
            if(woodrat<rang()) flIntType = 'reject'    !reject interaction
            if(flIntType=='clean') then                !accept interaction
              if(radWood_actscatrat(position,scatrat)>rang()) then
                flIntType = 'scatter'
              else
                flIntType = 'absorb'
              endif
            endif
          case ("KLWood")
            woodrat = KLrxi_point(j,position)/ceilsig
            if(woodrat>1.0d0) then                      !assert woodrat
              !print *,"j: ",j,"  woodrat: ",woodrat
              stop 'Higher sig samples in KLWood than ceiling, exiting program'
            endif
            if(woodrat<0.0d0 .and. allowneg=='no') stop 'Neg number sampled in KLWood, exiting program'
            if(distneg=='yes' .and. woodrat>0.0d0) then !redistribute neg option
              if(abs(disthold)>woodrat*ceilsig) then
                woodrat = 0.0d0
              else
                woodrat = (woodrat*ceilsig + disthold) / ceilsig
              endif
              disthold = 0.0d0
            endif
            if(allowneg=='yes') then                    !tallies for neg/pos if allowing neg
              if(woodrat<0.0d0) then
                numpnSamp(2)  =  numpnSamp(2)+1
                areapnSamp(2) = areapnSamp(2)+woodrat*ceilsig          
                if(woodrat*ceilsig<areapnSamp(4)) areapnSamp(4)=woodrat*ceilsig
                if(distneg=='yes') disthold = woodrat*ceilsig
                woodrat=0.0d0
              else
                numpnSamp(1)  =  numpnSamp(1)+1
                areapnSamp(1) = areapnSamp(1)+woodrat*ceilsig          
                if(woodrat*ceilsig>areapnSamp(3)) areapnSamp(3)=woodrat*ceilsig
              endif
            endif
            if(woodrat<rang()) flIntType = 'reject'    !reject interaction
            if(flIntType=='clean') then                !accept interaction
              if(scatrat(1)>rang()) then !this will need amending when scatrat is changed
                flIntType = 'scatter'
              else
                flIntType = 'absorb'
              endif
            endif
        end select

        !Tally for neg stats
        select case (MCcases(icase))
          case ("radMC")
          case ("radWood")
            if(flIntType=='reject') then
              Wood_rej(2)=Wood_rej(2)+1
            else
              Wood_rej(1)=Wood_rej(1)+1
            endif
          case ("KLWood")
            if(flIntType=='reject') then
              Wood_rej(2)=Wood_rej(2)+1
            else
              Wood_rej(1)=Wood_rej(1)+1
            endif
        end select


        !Evaluate scatter or absorb
        if(flIntType=='scatter') then     !scatter
          select case (MCcases(icase))
            case ("radMC")
              if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
              if(rodOrplanar=='planar') mu = newmu()
            case ("radWood")
              if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
              if(rodOrplanar=='planar') mu = newmu()
            case ("KLWood")
              if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
              if(rodOrplanar=='planar') mu = newmu()
          end select
        endif

        if(flIntType=='absorb') then      !absorb
          select case (MCcases(icase))
            case ("radMC")
              absorb(j)=absorb(j)+1.0d0
              flExit='exit'
            case ("radWood")
              absorb(j)=absorb(j)+1.0d0
              flExit='exit'
            case ("KLWood")
              absorb(j)=absorb(j)+1.0d0
              flExit='exit'
          end select
        endif

      endif

      if(flExit=='exit') exit
!print *,"flExit      : ",flExit
!print *,"fldist      : ",fldist
!print *,"flIntType   : ",flIntType
!print *,"flEscapeDir : ",flEscapeDir
    enddo !simulate one pathlength
  enddo !loop over particles

  end subroutine MCtransport


























  subroutine genSourcePart( i,icase )
  !This subroutine generates a position and direction, and bin if needed
  !to specify a source particle.
  use genRealzvars, only: s, nummatSegs
  use MCvars, only: position, mu, rodOrplanar, sourceType, MCcases
  integer :: i,icase

  if( sourceType=='left' ) then  !generate source particles
    position = 0.0d0
    mu       = 1.0d0
    if(rodOrplanar=='planar') mu = isoboundmu()
  elseif( sourceType=='intern' ) then
    position = s * rang()
    mu       = merge(1.0d0,-1.0d0,rang()>=0.5d0)
    if(rodOrplanar=='planar') mu = newmu()
  endif

  i = 0

  if(MCcases(icase)=='radMC') then !if bin need be set
    if(sourceType=='left') i = 1
    if(sourceType=='intern') i = internal_init_i(nummatSegs)
  endif

  end subroutine genSourcePart







  subroutine MCWood_setceils( j,icase )
  !This subroutine sets up ceiling values for Woodcock Monte Carlo.
  !These ceilings of course need to be recalculated for each new realization
  use genRealzvars, only: s, lamc, nummatSegs, sig
  use KLvars, only: numEigs
  use MCvars, only: MCcases, binmaxind, binmaxes, fbinmax, bbinmax, nceilbin

  use Woodcock, only: KLWood_binmaxes
  integer :: j,icase

  integer :: i
  real(8) :: binlength
!print *,"in MCWood_setceils"
    !select local bin maxes
    if(MCcases(icase)=='radWood') nceilbin = ceiling(s/lamc)
    if(MCcases(icase)=='KLWood' ) nceilbin = numEigs

    if(allocated(binmaxind)) deallocate(binmaxind)
    if(allocated(binmaxes)) deallocate(binmaxes)
    if(allocated(fbinmax)) deallocate(fbinmax)
    if(allocated(bbinmax)) deallocate(bbinmax)
    allocate(binmaxind(nceilbin+1))
    allocate(binmaxes(nceilbin))
    allocate(fbinmax(nceilbin))
    allocate(bbinmax(nceilbin))
    binmaxind=0.0d0
    binmaxes=0.0d0
    fbinmax=0.0d0
    bbinmax=0.0d0

    binlength=s/nceilbin
    binmaxind(1)=0.0d0
    do i=2,nceilbin+1
      binmaxind(i)=binmaxind(i-1)+binlength
    enddo

    if(MCcases(icase)=='radWood') call radWood_binmaxes2( nummatSegs )
    if(MCcases(icase)=='KLWood' ) call KLWood_binmaxes2( j )

    !create forward/backward motion max vectors
    bbinmax(1)=binmaxes(1)
    fbinmax(nceilbin)=binmaxes(nceilbin)
    do i=2,nceilbin
      bbinmax(i) = merge(bbinmax(i-1),binmaxes(i),bbinmax(i-1)>binmaxes(i))
    enddo
    do i=nceilbin-1,1,-1
      fbinmax(i) = merge(fbinmax(i+1),binmaxes(i),fbinmax(i+1)>binmaxes(i))
    enddo
!print *,"binmaxind: ",binmaxind
!print *,"bbinmax: ",bbinmax
!print *,"fbinmax: ",fbinmax
  end subroutine MCWood_setceils



  subroutine KLWood_binmaxes2( j ) ! make '1' when done with Woodcock version
  use KLvars, only: alpha, Ak, Eig, numEigs, sigave
  use MCvars, only: binmaxind, binmaxes, nceilbin

  use Woodcock, only: KLrxi_point !rid of when done with Woodcock version
  integer :: j

  integer :: i,k
  integer :: numinnersteps = 7
  integer :: numrefine     = 5
  real(8) :: safetyfactor  = 1.1d0
  real(8) :: innerstep,maxsig,maxpos,xpos,xsig,xpos1,xsig1,xpos2,xsig2

  do i=1,nceilbin
    innerstep = (binmaxind(2)-binmaxind(1)) / (numinnersteps-1)
    maxsig=0.0d0             !find initial max val
    maxpos=0.0d0
    do k=1,numinnersteps
      xpos=binmaxind(i)+(k-1)*innerstep
      xsig= KLrxi_point(j,xpos)
      if(xsig>maxsig) then
        maxsig=xsig
        maxpos=xpos
      endif
    enddo
    do k=1,numrefine     !refine
      innerstep=innerstep/2
      xpos1=maxpos-innerstep
      xsig1= KLrxi_point(j,xpos1)
      xpos2=maxpos+innerstep
      xsig2= KLrxi_point(j,xpos2)
      if(xsig1>maxsig .AND. xsig1>xsig2) then
        maxsig=xsig1
        maxpos=xpos1
      elseif(xsig2>maxsig .AND. xsig2>xsig1) then
        maxsig=xsig2
        maxpos=xpos2
      else
        maxpos=xpos
      endif
    enddo
    binmaxes(i)=maxsig*safetyfactor   !assign maxsig

  enddo

  end subroutine KLWood_binmaxes2





  subroutine radWood_binmaxes2( numArrSz ) !make simply '1' when other gone from Woodcock
  !subroutine starts to set up ceiling for WoodcockMC by mapping highest point in each bin
  use genRealzvars, only: matType, matLength, sig
  use MCvars, only: nceilbin, binmaxind, binmaxes
  integer :: numArrSz

  integer :: i,k
  real(8) :: smallersig,largersig

  smallersig=minval(sig)
  largersig =maxval(sig)

  do i=1,nceilbin
    do k=1,numArrSz

      !contains both
      if(matLength(k)>binmaxind(i) .AND. matLength(k)<binmaxind(i+1)) then
        binmaxes(i)=largersig
        exit
      endif

      !contains only one
      if(matLength(k)<=binmaxind(i) .AND. matLength(k+1)>=binmaxind(i+1)) then
        binmaxes(i)=sig(matType(k))
        exit
      endif

    enddo
  enddo
!print *,"binmaxes: ",binmaxes
  end subroutine radWood_binmaxes2




  function ceilsigfunc2(position,binmax) ! make '1' when done with Woodcock version
  use MCvars, only: nceilbin, binmaxind
  real(8) :: position,ceilsigfunc2,binmax(:)

  integer :: i
!print *,"ceilsigfunc2 spot 1"
  do i=1,nceilbin
    if( binmaxind(i)<=position .AND. binmaxind(i+1)>position ) then
      ceilsigfunc2 = binmax(i)
      exit
    endif
  enddo
!print *,"ceilsigfunc2 spot 2"
  end function ceilsigfunc2





  subroutine MCinc_pos( newposition )
  !This subroutine stores position in new position and reference passed value as position
  use MCvars, only: oldposition, position
  real(8) :: newposition
  oldposition = position
  position    = newposition
  end subroutine MCinc_pos




  subroutine stocMC_stats( icase )
  !This subroutine calculates mean and standard deviation for leakage values over stochastic
  !space for each MC transport solver.
  use genRealzvars, only: numRealz
  use MCvars, only: reflect, transmit, absorb, stocMC_reflection, &
                    stocMC_transmission, stocMC_absorption, numParts
  integer :: icase

print *,"reflect: ",reflect
print *,"transmit: ",transmit
print *,"absorb: ",absorb

  reflect  = reflect  / numParts
  transmit = transmit / numparts
  absorb   = absorb   / numParts

  call mean_and_var_s( reflect,numRealz,stocMC_reflection(icase,1),stocMC_reflection(icase,2) )
  call mean_and_var_s( transmit,numRealz,stocMC_transmission(icase,1),stocMC_transmission(icase,2) )
  call mean_and_var_s( absorb,numRealz,stocMC_absorption(icase,1),stocMC_absorption(icase,2) )

print *,"stocMC_reflection: ",stocMC_reflection
print *,"stocMC_transmission: ",stocMC_transmission
print *,"stocMC_absorption: ",stocMC_absorption
print *
print *

  end subroutine stocMC_stats




  subroutine MCallocate( icase )
  !This subroutine allocates and initializes variables that will be passed
  !through generic MCtransport subroutine.  These values will later be
  !stored in different arrays so that the variables can be re-used in
  !MCtransport if multiple cases were selected.
  use genRealzvars, only: numRealz
  use MCvars, only: transmit, reflect, absorb, radtrans_int, MCcases, &
                    numpnSamp, areapnSamp, disthold, Wood_rej
  integer :: icase

  if(.not.allocated(transmit)) allocate(transmit(numRealz))
  if(.not.allocated(reflect)) allocate(reflect(numRealz))
  if(.not.allocated(absorb)) allocate(absorb(numRealz))
  transmit     =0.0d0
  reflect      =0.0d0
  absorb       =0.0d0
  radtrans_int =0

  if(MCcases(icase)=='radWood' .or. MCcases(icase)=='KLWood') Wood_rej = 0

  if(MCcases(icase)=='KLWood') then
    numpnSamp =0
    areapnSamp=0.0d0
    disthold  =0.0d0
  endif

  end subroutine MCallocate








  subroutine radtrans_MCsim( j )
!this is the old subroutine for radMC, doomed to be deleted when UQ_MC and MCtransport overshadow it
  use timevars, only: time
  use genRealzvars, only: sig, scatrat, numRealz, nummatSegs, matType, matLength, s
  use MCvars, only: numParts, radtrans_int, pfnumcells, rodOrplanar, sourceType, &
                    plotflux, pltflux, reflect, transmit, absorb, initcur, fluxfaces, &
                    fflux, bflux, flux
  integer  :: j
  real(8)  :: tt1,tt2

  integer  :: i,z,o
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
      i        = internal_init_i(nummatSegs)
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



  function internal_init_i( numArrSz )
  use genRealzvars, only: matLength
  use MCvars, only: position
  integer :: numArrSz

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






  subroutine MCLeakage_pdfplot
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

  end subroutine MCLeakage_pdfplot



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






  subroutine MCprintstats
  !This subroutine prints reflection, transmission, and absorption stats to a '.out' file,
  !then prints that file to the screen for user friendliness.
  !Stats are from Adams89, Brantley11, and those generated here!
  use genRealzvars, only: Adamscase
  use MCvars, only: ABreflection, ABtransmission, rodOrplanar, stocMC_reflection, &
                    stocMC_transmission, stocMC_absorption, MCcases, MCcaseson
  integer :: icase

  320 format(" |AdamsMC:  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  321 format(" |BrantMC:  |",f8.5,"                 |",f8.5,"             |")
  322 format(" |AdamsLP:  |",f8.5,"                 |",f8.5,"             |")
  323 format(" |BrantLP:  |",f8.5,"                 |",f8.5,"             |")
  324 format(" |BrAtMix:  |",f8.5,"                 |",f8.5,"             |")
  325 format(" |-case:",f3.1,"-|---- Reflection and Transmission Results ------|")

  326 format(" |radMC  :  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  327 format(" |radWood:  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  328 format(" |KLWood :  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")

  !print to file
  open(unit=100,file="MCleakage.out")
  !print headings/benchmark solutions for full problem
  write(100,*)
  if(Adamscase/=0) write(100,325) Adamscase
  if(Adamscase==0) write(100,*) "|--------------- Reflection and Transmission Results ------|"
  write(100,*) "|Method    | reflave      refldev    | tranave      trandev|"
  write(100,*) "|----------|-------------------------|---------------------|"
  if(Adamscase/=0) then
    write(100,320) ABreflection(1,1),ABreflection(2,1),ABtransmission(1,1),ABtransmission(2,1)
    if(rodOrplanar=='planar') write(100,321) ABreflection(1,3),ABtransmission(1,3)
  endif


print *,"stocMC_reflection:",stocMC_reflection
print *,"stocMC_transmission: ",stocMC_transmission
  !print my solutions for radMC, radWood, KLWood
  do icase = 1,3
    if(MCcaseson(icase)==1) then
      if(MCcases(icase)=='radMC')   write(100,326) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))
      if(MCcases(icase)=='radWood') write(100,327) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))
      if(MCcases(icase)=='KLWood')  write(100,328) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))
    endif
  enddo

  !print benchmark LP solutions
  if(Adamscase/=0) then
    !write(100,*) "|          |                         |                     |" 
    write(100,*) "|----------|-------------------------|---------------------|"
    write(100,322) ABreflection(1,2),ABtransmission(1,2)
    if(rodOrplanar=='planar') write(100,323) ABreflection(1,4),ABtransmission(1,4)
  endif

  !print benchmark atomic mix solutions
  if(Adamscase/=0) then
    !write(100,*) "|          |                         |                     |" 
    write(100,*) "|----------|-------------------------|---------------------|"
    if(rodOrplanar=='planar') write(100,324) ABreflection(1,5),ABtransmission(1,5)
  endif
  write(100,*) "|----------------------------------------------------------|"
  write(100,*)
  close(unit=100)

  !move file and print to screen
  call system("mv MCleakage.out texts")

  end subroutine MCprintstats




end module radtransMC
