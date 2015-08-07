module radtransMC
  use mcnp_random
  use utilities
  implicit none

CONTAINS
  ! print statements in this module use # 300-399

  subroutine UQ_MC( icase )
  !This subroutine perfoms Monte Carlo in the uncertain space, currently for binary mixtures.
  !'MCtransport' handles the spatial MC, but this subroutine collects data and performs stats
  !in UQ space.
  use timevars, only: time
  use genRealzvars, only: numRealz, flGBgeom, posRealz
  use MCvars, only: MCcases, MCcaseson, numParts, trannprt, flfluxplotmat, refsigMode, &
                    flnegxs
  use genRealz, only: genReal
  use KLresearch, only: KL_eigenvalue
  use KLreconstruct, only: KLreconstructions
  use timeman, only: radtrans_timeupdate
  integer :: icase

  integer :: j,tnumParts,tnumRealz !'j' is which realization
  real(8) :: tt1

  call cpu_time(tt1)
  write(*,*) "Starting method: ",MCcases(icase)  
  call MCallocate( icase,tnumParts,tnumRealz )!allocate/initialize tallies

  if(  MCcases(icase)=='GaussKL' .and. flGBgeom) &
    call KL_eigenvalue

  if(  MCcases(icase)=='KLWood'   .or. &
      (MCcases(icase)=='WAMC' .and. MCcaseson(2)==1 .and. MCcases(icase)=='KLWood' ) .or. &
       MCcases(icase)=='GaussKL'     ) &
    call KLreconstructions(icase)       !create KL realz for cases that need them

  do j=1,tnumRealz

    if(MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood' .or. &
       MCcases(icase)=='LPMC'  .or. MCcases(icase)=='atmixMC' .or. flnegxs .or. &
       (.not.flnegxs .and. posRealz(j)==1)) then !if not WMC, neg ok, or realz positive

        if(MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood') &
      call genReal( j,'binary ',icase )         !gen binary geometry
        if(MCcases(icase)=='atmixMC') &
      call genReal( j,'atmixMC',icase )         !gen atomic mix geometry

        if(flfluxplotmat .and. (MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood')) &
      call MCprecalc_fluxmatnorm( j )           !collect normalization for flux in cells

        if(MCcases(icase)=='radWood' .or. MCcases(icase)=='KLWood' .or. &
          (MCcases(icase)=='WAMC' .and. (refsigMode==2 .or. refsigMode==3)) .or. &
           MCcases(icase)=='GaussKL') &
      call MCWood_setceils( j,icase )           !for WMC, create ceilings

      call MCtransport( j,icase,tnumParts )     !transport over a realization
    endif

      if(mod( j,trannprt )==0 .or. MCcases(icase)=='LPMC' .or. MCcases(icase)=='atmixMC') &
    call radtrans_timeupdate( j,icase,tt1 )   !print time updates

  enddo !loops over realizations

  call stocMC_stats( icase,tnumRealz )        !calc stats in stochastic space here
                                              !later make the above two loops,
                                              !batch spatial stats in between, final out here
  if(MCcases(icase)=='KLWood' .or. MCcases(icase)=='WAMC' .or. MCcases(icase)=='GaussKL' ) &
    call Woodnegstats(icase)

  end subroutine UQ_MC





  subroutine MCtransport( j,icase,tnumParts )
  !This subroutine performs MC transport over a specified geometry and method.
  !It can be used by a UQ_MC wrapper, and likely in the future by UQ_MLMC and/or UQ_SC.
  !'j' denotes which realization over which the MC transport is performed.
  !'icase' denotes which predefined transport method over which the MC transport is performed.
  !'tnumParts' is the number of particles over which the MC transport is performed.
  use rngvars, only: rngappnum, rngstride, setrngappnum
  use genRealzvars, only: sig, scatrat, nummatSegs, matType, matLength, s, lam, &
                          atmixsig, atmixscatrat
  use MCvars, only: radtrans_int, rodOrplanar, sourceType, reflect, transmit, &
                    absorb, position, oldposition, mu, areapnSamp, numpnSamp, &
                    nceilbin, Wood_rej, flnegxs, fldistneg, MCcases, fbinmax, &
                    bbinmax, binmaxind, binmaxes, LPamMCsums, flfluxplot, weight, &
                    refsig, wgtmax, wgtmin, wgtmaxmin, refsigMode, negwgtsigs, flCorrMC
  use genRealz, only: genReal
  use KLreconstruct, only: KLrxi_point
  use mcnp_random, only: RN_init_particle

  integer :: j,icase,tnumParts !realz number/which mode of transport/num of particles

  !local variables
  integer :: i,o,curbin
  real(8) :: db,dc,di,dist,  newpos,  sigma,  ceilsig,woodrat,disthold, tempscatrat
  real(8) :: cursigt, cursigs
  character(9) :: fldist, flIntType, flEscapeDir, flExit

  do o=1,tnumParts                     !loop over particles
    if(flCorrMC) then
      call setrngappnum(MCcases(1))
    else
      call setrngappnum(MCcases(icase))
    endif
    call RN_init_particle( int(rngappnum*rngstride+j*tnumParts+o,8) )

    call genSourcePart( i,icase )      !gen source part pos, dir, and binnum (i), init weight
    if(MCcases(icase)=='LPMC') call genReal( j,'LPMC   ',icase ) !for LP, choose starting material
    do ! simulate one pathlength of a particle
      fldist      = 'clean'
      flIntType   = 'clean'
      flEscapeDir = 'clean'
      flExit      = 'clean'
      newpos      = 101010.0d0 !later throw error if still equal this
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
        case ("LPMC")
          db = merge(s-position,position,mu>=0)/abs(mu)
        case ("atmixMC")
          db = merge(s-position,position,mu>=0)/abs(mu)
        case ("WAMC")
          curbin =solvecurbin(position)
          if(refsigMode==2) then
            db = merge(binmaxind(curbin+1)-position,position-binmaxind(curbin),mu>=0)/abs(mu)
          elseif(refsigMode==1 .or. refsigMode==3) then
            db = merge(s-position,position,mu>=0)/abs(mu)
          endif
        case ("GaussKL")
          db = merge(s-position,position,mu>=0)/abs(mu)
      end select

      !calculate distance to collision
      select case (MCcases(icase))
        case ("radMC")
          dc = -log(rang())/sig(matType(i))
        case ("radWood")
!          ceilsig=merge(ceilsigfunc(position,fbinmax),& !sel max sig
!                        ceilsigfunc(position,bbinmax),mu>=0)
          curbin =solvecurbin(position)
          ceilsig=merge(fbinmax(curbin),bbinmax(curbin),mu>=0)
          dc = -log(rang())/ceilsig                   !calc dc
        case ("KLWood")
!          ceilsig=merge(ceilsigfunc(position,fbinmax),& !sel max sig
!                        ceilsigfunc(position,bbinmax),mu>=0)
          curbin =solvecurbin(position)
          ceilsig=merge(fbinmax(curbin),bbinmax(curbin),mu>=0)
          dc = -log(rang())/ceilsig                   !calc dc
        case ("LPMC")
          dc = -log(rang())/sig(matType(1))
        case ("atmixMC")
          dc = -log(rang())/atmixsig
        case ("WAMC")
          if(refsigMode==1) then
            dc = -log(rang())/refsig 
          elseif(refsigMode==2) then
            curbin =solvecurbin(position)
            refsig =binmaxes(curbin)
            dc = -log(rang())/refsig
          elseif(refsigMode==3) then
            ceilsig=merge(fbinmax(curbin),bbinmax(curbin),mu>=0)
            dc = -log(rang())/ceilsig
          endif
        case ("GaussKL")
          curbin =solvecurbin(position)
          ceilsig=merge(fbinmax(curbin),bbinmax(curbin),mu>=0)
          dc = -log(rang())/ceilsig                   !calc dc
      end select

      !calculate distance to interface
      if(MCcases(icase)=='LPMC') di = -log(rang())*lam(matType(1))/abs(mu)

      !select distance limiter
      dist   = min(db,dc)
      fldist = merge('boundary ','collision',db<dc)
      if(MCcases(icase)=='LPMC') then
        fldist = merge('interface',fldist,di<dist)
        dist   = min(di,dist)
      endif



      !If 'boundary' event, and evaluate/flag transmission and reflection effects
      if(fldist=='boundary') then
        !set direction flag
        flEscapeDir = merge('transmit','reflect ',mu>0.0d0)
        if(MCcases(icase)=='radWood' .or. MCcases(icase)=='KLWood' .or. &
           MCcases(icase)=='GaussKL') Wood_rej(1)=Wood_rej(1)+1           !accept path tal

        !evaluate for local or global transmission or reflection
        if(flEscapeDir=='transmit') then   !transmit
          select case (MCcases(icase))
            case ("radMC")
              newpos = matLength(i+1)
              if(i==nummatSegs) transmit(j)=transmit(j) + 1.0d0
              if(i==nummatSegs) flExit='exit'
            case ("radWood")
              newpos = s
              transmit(j) = transmit(j) + 1.0d0
              flExit='exit'
            case ("KLWood")
              newpos = s
              transmit(j) = transmit(j) + 1.0d0
              flExit='exit'
            case ("LPMC")
              newpos = s
              LPamMCsums(2) = LPamMCsums(2) + 1.0d0
              flExit='exit'
            case ("atmixMC")
              newpos = s
              LPamMCsums(2) = LPamMCsums(2) + 1.0d0
              flExit='exit'
            case ("WAMC")
              if(refsigMode==2) then
                newpos = binmaxind(curbin+1) + 0.00000000001d0
              elseif(refsigMode==1 .or. refsigMode==3) then  
                newpos = s
              endif
              if(refsigMode/=2 .or. curbin+1==size(binmaxind)) then
                if(wgtmaxmin=='yes') then
                  if(    weight>wgtmax) then !option to truncate large weight values
                    weight=wgtmax
                  elseif(weight<wgtmin) then
                    weight=wgtmin
                  endif
                endif
                transmit(j) = transmit(j) + weight
                flExit='exit'
              endif
            case ("GaussKL")
              newpos = s
              transmit(j) = transmit(j) + 1.0d0
              flExit='exit'
          end select
        endif

        if(flEscapeDir=='reflect') then    !reflect
          select case (MCcases(icase))
            case ("radMC")
              newpos = matLength(i)
              if(i==1)            reflect(j)=reflect(j) + 1.0d0
              if(i==1)            flExit='exit'
            case ("radWood")
              newpos = 0.0d0
              reflect(j)  = reflect(j) + 1.0d0
              flExit='exit'
            case ("KLWood")
              newpos = 0.0d0
              reflect(j)  = reflect(j) + 1.0d0
              flExit='exit'
            case ("LPMC")
              newpos = 0.0d0
              LPamMCsums(1)  = LPamMCsums(1) + 1.0d0
              flExit='exit'
            case ("atmixMC")
              newpos = 0.0d0
              LPamMCsums(1)  = LPamMCsums(1) + 1.0d0
              flExit='exit'
            case ("WAMC")
              if(refsigMode==2) then
                newpos = binmaxind(curbin) - 0.00000000001d0
              elseif(refsigMode==1 .or. refsigMode==3) then
                newpos = 0.0d0
              endif
              if(refsigMode/=2 .or. curbin==1) then
                if(wgtmaxmin=='yes') then
                  if(   weight>wgtmax) then !option to truncate large weight values
                    weight=wgtmax
                  elseif(weight<wgtmin) then
                    weight=wgtmin
                  endif
                endif
                reflect(j) = reflect(j) + weight
                flExit='exit'
              endif
            case ("GaussKL")
              newpos = 0.0d0
              reflect(j)  = reflect(j) + 1.0d0
              flExit='exit'
          end select
        endif

      endif !endif fldist=='boundary'


      !If 'collision' event, set position and flag 'absorb', 'scatter', or 'reject'
      if(fldist=='collision') then
        !Set new position at which to advance
        newpos = position + dc*mu

        !Choose scatter, absorb, or reject (Woodcock) interaction
        select case (MCcases(icase))
          case ("radMC")
            flIntType = merge('scatter','absorb ',rang()<scatrat(matType(i)))
          case ("radWood")
            woodrat = radWood_actsig(newpos,sig)/ceilsig
            if(woodrat>1.0d0) then
              stop 'Higher sig samples in radWood than ceiling, exiting program'
            endif
            if(woodrat<rang()) flIntType = 'reject'    !reject interaction
            if(flIntType=='clean') then                !accept interaction
              if(radWood_actscatrat(newpos,scatrat)>rang()) then
                flIntType = 'scatter'
              else
                flIntType = 'absorb'
              endif
            endif
          case ("KLWood")
            !load woodcock ratio for this position and ceiling
            woodrat = KLrxi_point(j,newpos,chxstype='total')/ceilsig
            !assert within bounds, tally negstats
            if(woodrat>1.0d0) then                      !assert woodrat
              stop 'Higher sig samples in KLWood than ceiling, exiting program'
            endif
            if(woodrat<0.0d0 .and. .not.flnegxs) stop 'Neg number sampled in KLWood, exiting program'
            if(fldistneg .and. woodrat>0.0d0) then !redistribute neg option
              if(abs(disthold)>woodrat*ceilsig) then
                woodrat = 0.0d0
              else
                woodrat = (woodrat*ceilsig + disthold) / ceilsig
              endif
              disthold = 0.0d0
            endif
            if(flnegxs) then                    !tallies for neg/pos if allowing neg
              if(woodrat<0.0d0) then
                numpnSamp(2)  =  numpnSamp(2)+1
                areapnSamp(2) = areapnSamp(2)+woodrat*ceilsig          
                if(woodrat*ceilsig<areapnSamp(4)) areapnSamp(4)=woodrat*ceilsig
                if(fldistneg) disthold = woodrat*ceilsig
                woodrat=0.0d0
              else
                numpnSamp(1)  =  numpnSamp(1)+1
                areapnSamp(1) = areapnSamp(1)+woodrat*ceilsig          
                if(woodrat*ceilsig>areapnSamp(3)) areapnSamp(3)=woodrat*ceilsig
              endif
            endif
            !decide fate of particle
            if(woodrat<rang()) flIntType = 'reject'    !reject interaction
            if(flIntType=='clean') then                !accept interaction
              tempscatrat = KLWood_actscatrat(j,newpos)
              flIntType = merge('scatter','absorb ',tempscatrat>rang())
            endif
          case ("LPMC")
              flIntType = merge('scatter','absorb ',rang()<scatrat(matType(1)))
          case ("atmixMC")
              flIntType = merge('scatter','absorb ',rang()<atmixscatrat)
          case ("WAMC")
            if(refsigMode==3) then !for adaptive, decide woodcock possible or rejected interaction
              !use the next two lines for refsig from the bin
              curbin = solvecurbin(newpos)
              refsig = binmaxes(curbin)
              !use the next 7 lines for refsig locally
              refsig = localrefsig(KLrxi_point(j,newpos,chxstype='scatter'),&
                                   KLrxi_point(j,newpos,chxstype='total'))
              if(refsig>ceilsig) then
                !write(*,'(A,es9.2,A,es9.2,A,es9.2,A)') "refsig:",refsig,&
                !      " ceilsig:",ceilsig," truncation %:",(ceilsig-refsig)/refsig*100," %"
                refsig = ceilsig
              endif

              woodrat = refsig/ceilsig
              if(woodrat<rang()) flIntType='reject'
            endif
            if(flIntType=='clean') then !if refsigMode==1,2 or (==3 and not rejected yet)
              !adjust weights, choose type of interaction, flip weight sign if needed
              cursigt = KLrxi_point(j,newpos,chxstype='total')
              cursigs = KLrxi_point(j,newpos,chxstype='scatter')
              weight  = (abs(cursigs)+abs(-cursigt+refsig)) / refsig  *  weight

              if(rang()<abs(cursigs)/(abs(cursigs)+abs(-cursigt+refsig))) then
                flIntType = 'scatter'
                if(cursigs<0d0)         weight = -1d0 * weight
              else
                flIntType = 'reject '
                if(-cursigt+refsig<0d0) weight = -1d0 * weight
              endif
            endif

            if(cursigt<0.0d0) then  !negativity stats, currently overrides KLWood negstats
              numpnSamp(2)  =  numpnSamp(2)+1
              areapnSamp(2) = areapnSamp(2)+cursigt          
              if(cursigt<areapnSamp(4)) areapnSamp(4)=cursigt
            else
              numpnSamp(1)  =  numpnSamp(1)+1
              areapnSamp(1) = areapnSamp(1)+cursigt
              if(cursigt>areapnSamp(3)) areapnSamp(3)=cursigt
          endif
          case ("GaussKL")
            !load woodcock ratio for this position and ceiling
            woodrat = KLrxi_point(j,newpos,chxstype='total')/ceilsig
            !assert within bounds, tally negstats
            if(woodrat>1.0d0) then                      !assert woodrat
              stop 'Higher sig samples in KLWood than ceiling, exiting program'
            endif
            if(woodrat<0.0d0 .and. .not.flnegxs) stop 'Neg number sampled in GaussKL, exiting program'
            if(fldistneg .and. woodrat>0.0d0) then !redistribute neg option
              if(abs(disthold)>woodrat*ceilsig) then
                woodrat = 0.0d0
              else
                woodrat = (woodrat*ceilsig + disthold) / ceilsig
              endif
              disthold = 0.0d0
            endif
            if(flnegxs) then                    !tallies for neg/pos if allowing neg
              if(woodrat<0.0d0) then
                numpnSamp(2)  =  numpnSamp(2)+1
                areapnSamp(2) = areapnSamp(2)+woodrat*ceilsig          
                if(woodrat*ceilsig<areapnSamp(4)) areapnSamp(4)=woodrat*ceilsig
                if(fldistneg) disthold = woodrat*ceilsig
                woodrat=0.0d0
              else
                numpnSamp(1)  =  numpnSamp(1)+1
                areapnSamp(1) = areapnSamp(1)+woodrat*ceilsig          
                if(woodrat*ceilsig>areapnSamp(3)) areapnSamp(3)=woodrat*ceilsig
              endif
            endif
            !decide fate of particle
            if(woodrat<rang()) flIntType = 'reject'    !reject interaction
            if(flIntType=='clean') then                !accept interaction
              tempscatrat = KLWood_actscatrat(j,newpos)
              flIntType = merge('scatter','absorb ',tempscatrat>rang())
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
          case ("LPMC")
          case ("atmixMC")
          case ("WAMC")  !!!WAMCchange
          case ("GaussKL")
            if(flIntType=='reject') then
              Wood_rej(2)=Wood_rej(2)+1
            else
              Wood_rej(1)=Wood_rej(1)+1
            endif
        end select

      endif !endif fldist=='collision'

      !If 'interface' interaction, set new position (LP)
      if(fldist=='interface') then    !!!WAMCchange careful not to do this at interfaces
        newpos = position + di*mu
      endif


      !increment position
      if(newpos==101010.0d0) stop 'newpos was not set in MCtransport'
      call MCinc_pos( newpos )

      !tally flux
      if(flfluxplot) call MCfluxtallywrapper( j,icase )

      !Evaluate scatter
      if(flIntType=='scatter') then     !scatter
        if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0) !all currently do same
        if(rodOrplanar=='planar') mu = newmu()
!        select case (MCcases(icase))
!          case ("radMC")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!          case ("radWood")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!          case ("KLWood")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!          case ("LPMC")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!          case ("atmixMC")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!          case ("WAMC")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!          case ("GaussKL")
!            if(rodOrplanar=='rod')    mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
!            if(rodOrplanar=='planar') mu = newmu()
!        end select
      endif

      !Evaluate absorption
      if(flIntType=='absorb') then      !absorb
        select case (MCcases(icase))
          case ("radMC")
            absorb(j)     = absorb(j)     + 1.0d0
            flExit='exit'
          case ("radWood")
            absorb(j)     = absorb(j)     + 1.0d0
            flExit='exit'
          case ("KLWood")
            absorb(j)     = absorb(j)     + 1.0d0
            flExit='exit'
          case ("LPMC")
            LPamMCsums(3) = LPamMCsums(3) + 1.0d0
            flExit='exit'
          case ("atmixMC")
            LPamMCsums(3) = LPamMCsums(3) + 1.0d0
            flExit='exit'
          case ("WAMC")
            absorb(j)     = absorb(j)     + weight
            flExit='exit'
          case ("GaussKL")
            absorb(j)     = absorb(j)     + 1.0d0
            flExit='exit'
        end select
      endif

      !increment material type indices
      if(MCcases(icase)=='radMC') then
        if(flEscapeDir=='transmit') i = i+1
        if(flEscapeDir=='reflect')  i = i-1
      endif
      if(fldist=='interface') then
        matType(1) = merge(1,2,matType(1)==2)
      endif

      !Exit if particle history ended
      if(flExit=='exit') exit

    enddo !simulate one pathlength
  enddo !loop over particles

  end subroutine MCtransport


























  subroutine genSourcePart( i,icase )
  !This subroutine generates a position and direction, and bin index if needed (radMC)
  !to specify a source particle.
  use genRealzvars, only: s
  use MCvars, only: position, mu, rodOrplanar, sourceType, MCcases, weight
  integer :: i,icase

  if( sourceType=='left' ) then  !generate source particles
    position = 0.0d0
    mu       = isoboundmu()
    if(rodOrplanar=='rod') mu = 1.0d0
  elseif( sourceType=='intern' ) then
    position = s * rang()
    mu       = newmu()
    if(rodOrplanar=='rod') mu = merge(1.0d0,-1.0d0,rang()>=0.5d0)
  endif

  i = 0  !why is this here?  test if I can get rid of it...

  if(MCcases(icase)=='radMC') then !if bin need be set
    if(sourceType=='left')   i = 1
    if(sourceType=='intern') i = internal_init_i(position)
  elseif(MCcases(icase)=='WAMC') then !initialize weight
    weight = 1d0
  endif

  end subroutine genSourcePart







  subroutine MCWood_setceils( j,icase )
  !This subroutine sets up ceiling values for Woodcock Monte Carlo.
  !These ceilings of course need to be recalculated for each new realization
  use genRealzvars, only: s, lamc, nummatSegs, sig, posRealz
  use KLvars, only: numEigs
  use MCvars, only: MCcases, binmaxind, binmaxes, fbinmax, bbinmax, nceilbin, &
                    refsigMode, negwgtbinnum
  integer :: j,icase

  integer :: i
  real(8) :: binlength
    !select local bin maxes
    select case (MCcases(icase))
      case ("radWood")
        nceilbin = ceiling(s/lamc)
        if(nceilbin>6) nceilbin = 6
      case ("KLWood")
        nceilbin = numEigs
      case ("WAMC")
        if(refsigMode==2 .or. refsigMode==3) nceilbin = negwgtbinnum
      case ("GaussKL")
        nceilbin = numEigs
    end select

    if(j==1 .or. sum(posRealz(1:j))==1) then
      if(allocated(binmaxind)) deallocate(binmaxind)
      if(allocated(binmaxes)) deallocate(binmaxes)
      if(allocated(fbinmax)) deallocate(fbinmax)
      if(allocated(bbinmax)) deallocate(bbinmax)
    endif

    if(.not. allocated(binmaxind)) allocate(binmaxind(nceilbin+1))
    if(.not. allocated(binmaxes))  allocate(binmaxes(nceilbin))
    if(.not. allocated(fbinmax))   allocate(fbinmax(nceilbin))
    if(.not. allocated(bbinmax))   allocate(bbinmax(nceilbin))
    binmaxind=0.0d0
    binmaxes=0.0d0
    fbinmax=0.0d0
    bbinmax=0.0d0

    !set bin indices
    binlength=s/nceilbin
    binmaxind(1)=0.0d0
    do i=2,nceilbin+1
      binmaxind(i)=binmaxind(i-1)+binlength
    enddo

    !set individual bin maxes
    select case (MCcases(icase))
      case ("radWood")
        call radWood_binmaxes( nummatSegs )
      case ("KLWood")
        call KLWood_binmaxes( j )
      case ("WAMC")
        if(refsigMode==2 .or. refsigMode==3) call WAMC_binmaxes( j )
      case ("GaussKL")
        call KLWood_binmaxes( j )
    end select

    !create forward/backward motion max vectors
    bbinmax(1)=binmaxes(1)
    fbinmax(nceilbin)=binmaxes(nceilbin)
    do i=2,nceilbin
      bbinmax(i) = merge(bbinmax(i-1),binmaxes(i),bbinmax(i-1)>binmaxes(i))
    enddo
    do i=nceilbin-1,1,-1
      fbinmax(i) = merge(fbinmax(i+1),binmaxes(i),fbinmax(i+1)>binmaxes(i))
    enddo

  end subroutine MCWood_setceils



  subroutine WAMC_binmaxes( j )
  use MCvars, only: binmaxind, binmaxes, nceilbin, negwgtsigs, nwvalsperbin
  use KLreconstruct, only: KLrxi_point
  integer, intent(in) :: j

  integer :: i, anc	
  real(8) :: pos

  !load sigs, sigt, and solve refsig at each location
  do i=1,nceilbin*nwvalsperbin+1
    if( mod(i-1,nwvalsperbin)==0 ) then
      pos = binmaxind(i/nwvalsperbin+1)
    else
      anc = (i-1)/nwvalsperbin+1
      pos = binmaxind(anc) + real(mod(i-1,nwvalsperbin),8)/real(nwvalsperbin,8) &
                           * (binmaxind(anc+1)-binmaxind(anc))
    endif
    negwgtsigs(i,1) = KLrxi_point(j,pos,chxstype='scatter')
    negwgtsigs(i,2) = KLrxi_point(j,pos,chxstype='total')
    negwgtsigs(i,3) = localrefsig(negwgtsigs(i,1),negwgtsigs(i,2))
  enddo

  !take max of values in/on bin in each bin as bin max refsig value
  do i=1,nceilbin
    binmaxes(i) = maxval(negwgtsigs((i-1)*nwvalsperbin+1:i*nwvalsperbin+1,3))
  enddo
  end subroutine WAMC_binmaxes



  function localrefsig(sigs,sigt)
  !based on maxind, find optimal refsig choice.  See notes on 1-23-15.
  use MCvars, only: maxratio
  real(8), intent(in) :: sigs,sigt
  real(8) :: localrefsig

  if(sigt<0.0d0) then
    localrefsig = (abs(sigs)+abs(sigt))/(maxratio-1d0)
  elseif(abs(sigs)<sigt) then
    localrefsig = sigt
  elseif(maxratio<abs(sigs)/sigt) then
    localrefsig = (abs(sigs)+sigt)/(maxratio+1d0)
  else
    localrefsig = (abs(sigs)-sigt)/(maxratio-1d0)
  endif
  end function localrefsig




  subroutine KLWood_binmaxes( j )
  use KLvars, only: alpha, Ak, Eig, numEigs
  use MCvars, only: binmaxind, binmaxes, nceilbin
  use KLreconstruct, only: KLrxi_point

  integer, intent(in) :: j

  integer :: i,k
  integer :: numinnersteps = 11
  integer :: numrefine     = 7
  real(8) :: safetyfactor  = 1.4d0
  real(8) :: innerstep,maxsig,maxpos,xpos,xsig,xpos1,xsig1,xpos2,xsig2

  do i=1,nceilbin
    innerstep = (binmaxind(2)-binmaxind(1)) / (numinnersteps-1)
    maxsig=0.0d0             !find initial max val
    maxpos=0.0d0
    do k=1,numinnersteps
      xpos=binmaxind(i)+(k-1)*innerstep
      xsig= KLrxi_point(j,xpos,chxstype='total')
      if(xsig>maxsig) then
        maxsig=xsig
        maxpos=xpos
      endif
    enddo
    do k=1,numrefine     !refine
      innerstep=innerstep/2
      xpos1=maxpos-innerstep
      xsig1= KLrxi_point(j,xpos1,chxstype='total')
      xpos2=maxpos+innerstep
      xsig2= KLrxi_point(j,xpos2,chxstype='total')
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

  end subroutine KLWood_binmaxes





  subroutine radWood_binmaxes( numArrSz ) !make simply '1' when other gone from Woodcock
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
  end subroutine radWood_binmaxes



  function radWood_actsig(position,sig)
  use genRealzvars, only: matType, matLength
  real(8) :: position,sig(2),radWood_actsig
  integer :: i
  i = internal_init_i( position )
  radWood_actsig=sig(matType(i))
  end function radWood_actsig




  function ceilsigfunc(position,binmax) ! make '1' when done with Woodcock version
  !This function being held, if solvecurbin proves reliable, then out with this one
  use MCvars, only: nceilbin, binmaxind
  real(8) :: position,ceilsigfunc,binmax(:)
  integer :: i
  do i=1,nceilbin
    if( binmaxind(i)<=position .AND. binmaxind(i+1)>position ) then
      ceilsigfunc = binmax(i)
      exit
    endif
  enddo
  end function ceilsigfunc



  function solvecurbin(position)
  use genRealzvars, only: s
  use MCvars, only: nceilbin
  integer :: solvecurbin
  real(8) :: position

  solvecurbin = floor(position/s*nceilbin)+1
  end function solvecurbin



  function radWood_actscatrat(position,scatrat)
  use genRealzvars, only: matType, matLength
  real(8) :: position,scatrat(2),radWood_actscatrat
  integer :: i
  i = internal_init_i( position )
  radWood_actscatrat=scatrat(matType(i))
  end function radWood_actscatrat


  function KLWood_actscatrat(j,xpos)
  !For 'totxs' type xi generation, this function returns the shared scattering ratio of the materials.
  !For 'material' type xi generation, it generates tot, scat, and abs xs values,
  !uses them to assert the same total xs and produce the local scattering ratio.  
  use genRealzvars, only: scatrat
  use KLvars, only: flmatbasedxs
  use KLreconstruct, only: KLrxi_point
  use MCvars, only: numcSamp
  integer :: j
  real(8) :: eps = 0.1d0
  real(8) :: KLWood_actscatrat, xpos, scatxs, absxs, totxs

  if(.not.flmatbasedxs) then
    KLWood_actscatrat = scatrat(1)
  elseif(flmatbasedxs) then
    numcSamp(2) = numcSamp(2) + 1 !tally all scatrat samples
    scatxs = KLrxi_point(j,xpos,chxstype='scatter')
    absxs  = KLrxi_point(j,xpos,chxstype='absorb')
    !totxs  = KLrxi_point(j,xpos,chxstype='total')
    if( scatxs*absxs<0d0 ) numcSamp(1) = numcSamp(1) + 1 !tally neg scatrat samples
                           !print *,"in KLWood_actscatrat, scatxs & absxs not same sign"
    if( scatxs<0d0 ) scatxs = 0d0
    if( absxs <0d0 ) absxs  = 0d0

    KLWood_actscatrat = scatxs/(scatxs+absxs)
  endif

  end function KLWood_actscatrat



  subroutine MCinc_pos( newposition )
  !This subroutine stores position in old position and reference passed value as position
  use MCvars, only: oldposition, position
  real(8) :: newposition
  oldposition = position
  position    = newposition
  end subroutine MCinc_pos



  subroutine MCprecalc_fluxmatnorm( j )
  !This subroutine calculates how much of each material type fills each bin in a realization.
  !This is determined once geometry is created, so for 'radMC' and 'radWood', this can be
  !done right after each initial geometry calculation.  It marches through the geometry of a 
  !realization, deciding if a fluxface is next, or a material face, calculates that contribution
  !and continues, terminating when the last sub-cell has been accounted for.
  !For 'LPMC' it must be collected while the algorithm runs, since the geometry is generated
  !as part of the algorithm ('MCLPcalc_fluxmatnorm').
  use genRealzvars, only: matLength, matType
  use MCvars, only: fluxmatnorm, fluxfaces, fluxnumcells
  integer :: i, ibin, j, matnumsegs
  real(8) :: nextboundary,lastboundary
  character(8) :: flnextboundary,fllastboundary

  fllastboundary = 'matface '
  matnumsegs = size(matType)

  i    = 1
  ibin = 1
  do
    !which comes next?
    if(matLength(i+1)<=fluxfaces(ibin+1)) then
      flnextboundary = 'matface '
    else
      flnextboundary = 'fluxface'
    endif
    nextboundary = merge(matLength(i+1),fluxfaces(ibin+1),flnextboundary=='matface')
    lastboundary = merge(matLength(i)  ,fluxfaces(ibin)  ,fllastboundary=='matface')

    !add contribution of 'sub-cell', which is the distance between faces of either kind
    fluxmatnorm(ibin,j,matType(i)) = fluxmatnorm(ibin,j,matType(i)) + nextboundary - lastboundary

    !test if over
    if(ibin==fluxnumcells .and. i==matnumsegs) exit

    !update parameters
    if(flnextboundary=='fluxface') then
      ibin = ibin + 1
    elseif(flnextboundary=='matface') then
      i    = i    + 1
    endif
    fllastboundary = flnextboundary

  enddo
  end subroutine MCprecalc_fluxmatnorm





  subroutine MCfluxtallywrapper( j,icase )
  !This subroutine is the outer wrapper for flux tallies.
  !First it calculates quantities related to this tally.
  !Then for most methods it passes that info the the tallier, but for 'radWood'
  !it steps through the tracklength passing on segments at a time to be tallied.
  use genRealzvars, only: matLength, matType
  use MCvars, only: position, oldposition, flfluxplotall, flfluxplotmat, MCcases, &
                    fluxfaces
  integer :: j, i, ibin, icase
  real(8) :: minpos, maxpos, dx

  minpos = min(oldposition,position)
  maxpos = max(oldposition,position)
  dx     = fluxfaces(2)-fluxfaces(1)
  !set value for 'i'
  i=1
  if( MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood' ) &
    i = internal_init_i( minpos )

  !fluxtally irrespective of mat
  if(flfluxplotall) then
    call MCfluxtally( j,i,minpos,maxpos,'irrespective' )
  endif

  !fluxtally respective of mat
  if(MCcases(icase)=='radWood' .and. flfluxplotmat) then
    !tally length through cells
    if(matLength(i)<maxpos .and. maxpos<=matLength(i+1)) then       !first cell = last cell
      call MCfluxtally( j,i,minpos,maxpos,'respective  ' )                  !tally first/last cell
    else                                                            !first cell /= last cell
      call MCfluxtally( j,i,minpos,matlength(i+1),'respective  ' )          !tally first cell
      do
        i = i+1
        if(matLength(i)<maxpos .and. maxpos<=matLength(i+1)) then   !last cell
          call MCfluxtally( j,i,matLength(i),maxpos,'respective  ' )        !tally last cell
          exit
        else                                                        !not first or last cell
          call MCfluxtally( j,i,matLength(i),matLength(i+1),'respective  ' )!tally middle cell
        endif
      enddo
    endif
  elseif(flfluxplotmat) then
    call MCfluxtally( j,i,minpos,maxpos,'respective  ' )
    if(MCcases(icase)=='LPMC') call MCLPcalc_fluxmatnorm( minpos,maxpos,dx )
  endif

  end subroutine MCfluxtallywrapper



  subroutine MCfluxtally( j,i,minpos,maxpos,flfluxtallytype )
  !This subroutine accepts arguments from its wrapper, then tallies flux contributions accordingly.
  !It can keep full 'track' length tallies, or choose one bin (at a 'point' on the tracklength) in
  !which to place the entire contribution.
  !Flux can be calculated for the tracklength w/ or w/o respect to which material the tallies are in.
  use MCvars, only: mu, fluxall, fluxmat1, fluxmat2, fluxfaces, pltfluxtype, fluxnumcells
  integer :: j, i, ibin
  real(8) :: minpos,maxpos,absmu,dx, length,point
  character(7) :: flcontribtype
  character(12) :: flfluxtallytype
  absmu  = abs(mu)
  dx     = fluxfaces(2)-fluxfaces(1)

  !material irrespective
  if( flfluxtallytype=='irrespective' ) then
    if( pltfluxtype=='point' ) then       !point selection
      point  = rang()*(maxpos-minpos) + minpos
      length = (maxpos-minpos) / absmu
      ibin   = ceiling(point/dx)
      if(ibin==0) ibin=1  !adjust if at ends
      fluxall(ibin,j) = fluxall(ibin,j) + length
    elseif( pltfluxtype=='track' ) then   !whole tracklength
      ibin   = ceiling(minpos/dx)
      if(ibin==0) ibin=1  !adjust if at ends
      do
        call MCfluxtallysetflag( flcontribtype, ibin, minpos, maxpos )
        select case (flcontribtype)  !!!WAMCchange add seperate one for WAMC present with weight factored in
          case ("neither")
            fluxall(ibin,j) = fluxall(ibin,j) +  dx                        / absmu !niether in bin
          case ("first")
            fluxall(ibin,j) = fluxall(ibin,j) + (fluxfaces(ibin+1)-minpos) / absmu !first in bin
          case ("last")
            fluxall(ibin,j) = fluxall(ibin,j) + (maxpos-fluxfaces(ibin))   / absmu !last in bin
          case ("both")
            fluxall(ibin,j) = fluxall(ibin,j) + (maxpos-minpos)            / absmu !both in bin
        end select
        if( flcontribtype=='last' .or. flcontribtype=='both' .or. ibin==fluxnumcells ) exit
        ibin = ibin + 1
      enddo
    endif !point or track
  endif !mat irrespective

  !material respective
  if( flfluxtallytype=='respective') then
    if( pltfluxtype=='point' ) then       !point selection
      point  = rang()*(maxpos-minpos) + minpos
      length = (maxpos-minpos) / absmu
      ibin   = ceiling(point/dx)
      if(ibin==0) ibin=1  !adjust if at ends
      call MCfluxtallycontribute( j,ibin,i,length )
    elseif( pltfluxtype=='track' ) then   !whole tracklength
      ibin   = ceiling(minpos/dx)
      if(ibin==0) ibin=1  !adjust if at ends

      do
        call MCfluxtallysetflag( flcontribtype, ibin, minpos, maxpos )
        select case (flcontribtype)
          case ("neither")
            call MCfluxtallycontribute( j,ibin,i, dx                        / absmu )
          case ("first")
            call MCfluxtallycontribute( j,ibin,i,(fluxfaces(ibin+1)-minpos) / absmu )
          case ("last")
            call MCfluxtallycontribute( j,ibin,i,(maxpos-fluxfaces(ibin))   / absmu )
          case ("both")
            call MCfluxtallycontribute( j,ibin,i,(maxpos-minpos)            / absmu )
        end select
        if( flcontribtype=='last' .or. flcontribtype=='both' .or. ibin==fluxnumcells ) exit
        ibin = ibin + 1
      enddo

    endif !point or track
  endif !mat respective

  end subroutine MCfluxtally




  subroutine MCfluxtallycontribute( j,ibin,i,contribution )
  !This is used specifically in 'MCfluxtally' to add to the appropriate material bin
  use genRealzvars, only: matType
  use MCvars, only: fluxmat1, fluxmat2
  integer :: j,ibin,i
  real(8) :: contribution

  if( matType(i)==1 ) then
    fluxmat1(ibin,j) = fluxmat1(ibin,j) + contribution
  elseif( matType(i)==2 ) then
    fluxmat2(ibin,j) = fluxmat2(ibin,j) + contribution
  endif
  end subroutine MCfluxtallycontribute




  subroutine MCLPcalc_fluxmatnorm( minpos,maxpos,dx )
  !This subroutine tallies material abundances in each flux cell for LPMC so that 
  !the flux tallies can be normalized on a per material basis.  It is on the fly
  !since LP calculates geometry on the fly.
  !This routine is called from 'MCfluxtallywrapper' and highly resembles 'MCfluxtally',
  !with the only notable difference being that distances are not divided by absu(mu).
  use genRealzvars, only: matType
  use MCvars, only: fluxmatnorm, fluxfaces, fluxnumcells, position, &
                    pltfluxtype
  integer :: i, ibin
  real(8) :: minpos,maxpos,dx, point,length
  character(7) :: flcontribtype
  character(8) :: flnextboundary,fllastboundary

  if( pltfluxtype=='point' ) then       !point selection
    point  = rang()*(maxpos-minpos) + minpos
    length = (maxpos-minpos)
    ibin   = ceiling(point/dx)
    if(ibin==0) ibin=1  !adjust if at ends
    fluxmatnorm(ibin,1,matType(1)) = fluxmatnorm(ibin,1,matType(1)) + length
  elseif( pltfluxtype=='track' ) then   !whole tracklength
    ibin = ceiling(minpos/dx)
    if(ibin==0) ibin=1  !adjust if at ends
    do
      call MCfluxtallysetflag( flcontribtype, ibin, minpos, maxpos )
      select case (flcontribtype)
        case ("neither")
          fluxmatnorm(ibin,1,matType(1)) = fluxmatnorm(ibin,1,matType(1)) + &
                                            dx                                !mid bins
        case ("first")
          fluxmatnorm(ibin,1,matType(1)) = fluxmatnorm(ibin,1,matType(1)) + &
                                           (fluxfaces(ibin+1)-minpos)         !first bin
        case ("last")
          fluxmatnorm(ibin,1,matType(1)) = fluxmatnorm(ibin,1,matType(1)) + &
                                           (maxpos-fluxfaces(ibin))           !last bin
        case ("both")
          fluxmatnorm(ibin,1,matType(1)) = fluxmatnorm(ibin,1,matType(1)) + &
                                           (maxpos-minpos)                    !last bin
      end select
      if( flcontribtype=='last' .or. flcontribtype=='both' .or. ibin==fluxnumcells ) exit
      ibin = ibin + 1
    enddo

  endif !point or track

  end subroutine MCLPcalc_fluxmatnorm






  subroutine MCfluxtallysetflag( flcontribtype, ibin, minpos, maxpos )
  !This subroutine sets what type of flux contribution should be tallied.
  !The word 'first' means the first end length of the track considered is in this cell.
  !'last' means the second, or last, end of the track length is in the cell.
  !'neither' means both ends of the tracklength are not in this cell.
  !'both' means both ends of the tracklength are in the cell.
  use MCvars, only: fluxfaces
  integer :: ibin
  real(8) :: minpos, maxpos
  character(7) :: flcontribtype
  if( minpos<fluxfaces(ibin) ) then
    if( maxpos>fluxfaces(ibin+1) ) then
      flcontribtype = 'neither'
    else
      flcontribtype = 'last'
    endif
  else
    if( maxpos>fluxfaces(ibin+1) ) then
      flcontribtype = 'first'
    else
      flcontribtype = 'both'
    endif
  endif
  end subroutine MCfluxtallysetflag




  subroutine stocMC_stats( icase,tnumRealz )
  !This subroutine:
  !1) calculates mean and standard deviation for leakage values over stochastic
  !   space for each MC transport solver.
  !2) calculates ensemble averaged flux values in each cell for material
  !   irrespective and material respective flux tallies.
  !3) calls 'MCLeakage_pdfbinprint', which bins and prints leakage values for
  !   later plotting.
  use genRealzvars, only: numRealz, numPosRealz
  use MCvars, only: reflect, transmit, absorb, stocMC_reflection, LPamnumParts, &
                    stocMC_transmission, stocMC_absorption, numParts, LPamMCsums, &
                    MCcases, fluxnumcells, fluxall, flnegxs, &
                    fluxmat1, fluxmat2, stocMC_fluxall, stocMC_fluxmat1, stocMC_fluxmat2, &
                    fluxmatnorm, fluxfaces, flfluxplot, flfluxplotall, flfluxplotmat
  use timeman, only: FOM_calculation
  integer :: icase,ibin,tnumRealz,j
  real(8) :: dx,p1,p2
  real(8), allocatable :: posreflect(:),postransmit(:),posabsorb(:)
  real(8), allocatable :: posfluxall(:),posfluxmat1(:),posfluxmat2(:)

  !leakage/absorption stats
  if(MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood' .or. &
     MCcases(icase)=='KLWood'.or. MCcases(icase)=='WAMC'    .or. &
     MCcases(icase)=='GaussKL'                                     ) then
    reflect  = reflect  / numParts
    transmit = transmit / numparts
    absorb   = absorb   / numParts
    if(flnegxs .or. (MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood')) then !stats based on all realz
      call mean_and_var_s( reflect,numRealz,stocMC_reflection(icase,1),stocMC_reflection(icase,2) )
      call mean_and_var_s( transmit,numRealz,stocMC_transmission(icase,1),stocMC_transmission(icase,2) )
      call mean_and_var_s( absorb,numRealz,stocMC_absorption(icase,1),stocMC_absorption(icase,2) )
    elseif(.not.flnegxs .and. (MCcases(icase)=='KLWood' .or. & !stats based on fully-pos realz
                    MCcases(icase)=='WAMC' .or. MCcases(icase)=='GaussKL')) then
      call pos_all_pop( reflect,posreflect )
      call pos_all_pop( transmit,postransmit )
      call pos_all_pop( absorb,posabsorb )
      call mean_and_var_s( posreflect,numPosRealz,stocMC_reflection(icase,1),stocMC_reflection(icase,2) )
      call mean_and_var_s( postransmit,numPosRealz,stocMC_transmission(icase,1),stocMC_transmission(icase,2) )
      call mean_and_var_s( posabsorb,numPosRealz,stocMC_absorption(icase,1),stocMC_absorption(icase,2) )
    endif
  elseif(MCcases(icase)=='LPMC' .or. MCcases(icase)=='atmixMC') then
    stocMC_reflection(icase,1)   = LPamMCsums(1) / LPamnumParts
    stocMC_transmission(icase,1) = LPamMCsums(2) / LPamnumParts
    stocMC_absorption(icase,1)   = LPamMCsums(3) / LPamnumParts
  endif


  !flux stats
  if(flfluxplot)    dx = fluxfaces(2) - fluxfaces(1)

  if(MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood' .or. &
     MCcases(icase)=='KLWood'.or. MCcases(icase)=='GaussKL'        ) then !!!WAMCchange add here
    if(flfluxplotall) fluxall = fluxall / dx / numParts !normalize part 1
    if(flfluxplotmat) then
                      fluxmat1= fluxmat1/ dx / numParts !normalize part 1
                      fluxmat2= fluxmat2/ dx / numParts !normalize part 1
    endif    
  elseif(MCcases(icase)=='LPMC' .or. MCcases(icase)=='atmixMC') then
    if(flfluxplotall) fluxall = fluxall / dx / LPamnumParts !normalize part 1
    if(flfluxplotmat) then
                      fluxmat1= fluxmat1/ dx / LPamnumParts !normalize part 1
                      fluxmat2= fluxmat2/ dx / LPamnumParts !normalize part 1
    endif    
  endif

  if(MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood' .or. &
     MCcases(icase)=='KLWood'.or. MCcases(icase)=='GaussKL'       ) then
    if( flfluxplotall ) then  !!!WAMCchange
      do ibin=1,fluxnumcells

        if(flnegxs .or. (MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood')) then !stats from all realz
          call mean_and_var_s( fluxall(ibin,:),numRealz, &
                   stocMC_fluxall(ibin,icase,1),stocMC_fluxall(ibin,icase,2) )
        elseif(.not.flnegxs .and. (MCcases(icase)=='KLWood' .or. & !stats based on fully-pos realz
                    MCcases(icase)=='WAMC' .or. MCcases(icase)=='GaussKL')) then
          call pos_all_pop( fluxall(ibin,:),posfluxall )
          call mean_and_var_s( posfluxall(:),numPosRealz, &
                   stocMC_fluxall(ibin,icase,1),stocMC_fluxall(ibin,icase,2) )
        endif

      enddo
    endif
    if( flfluxplotmat ) then
      do ibin=1,fluxnumcells
        p1 = sum(fluxmatnorm(ibin,:,1)) / sum(fluxmatnorm(ibin,:,:))
        p2 = 1.0d0 - p1

        if(flnegxs .or. (MCcases(icase)=='radMC' .or. MCcases(icase)=='radWood')) then !stats from all realz
          fluxmat1(ibin,:) = fluxmat1(ibin,:) / p1
          fluxmat2(ibin,:) = fluxmat2(ibin,:) / p2
          call mean_and_var_s( fluxmat1(ibin,:),numRealz, &
                    stocMC_fluxmat1(ibin,icase,1),stocMC_fluxmat1(ibin,icase,2) )
          call mean_and_var_s( fluxmat2(ibin,:),numRealz, &
                    stocMC_fluxmat2(ibin,icase,1),stocMC_fluxmat2(ibin,icase,2) )
        elseif(.not.flnegxs .and. (MCcases(icase)=='KLWood' .or. & !stats based on fully-pos realz
                    MCcases(icase)=='WAMC' .or. MCcases(icase)=='GaussKL')) then
          call pos_all_pop( fluxmat1(ibin,:),posfluxmat1 )
          call pos_all_pop( fluxmat2(ibin,:),posfluxmat2 )
          posfluxmat1(:) = posfluxmat1(:) / p1
          posfluxmat2(:) = posfluxmat2(:) / p2
          call mean_and_var_s( posfluxmat1(:),numPosRealz, &
                    stocMC_fluxmat1(ibin,icase,1),stocMC_fluxmat1(ibin,icase,2) )
          call mean_and_var_s( posfluxmat2(:),numPosRealz, &
                    stocMC_fluxmat2(ibin,icase,1),stocMC_fluxmat2(ibin,icase,2) )
        endif 

      enddo
    endif

  elseif(MCcases(icase)=='LPMC') then
    if(flfluxplotall) then
      do ibin=1,fluxnumcells  !mean vals, no var, no mat specific
        stocMC_fluxall(ibin,icase,1) = fluxall(ibin,1)   !store mean
      enddo
    endif
    if(flfluxplotmat) then
      do ibin=1,fluxnumcells
        p1 = fluxmatnorm(ibin,1,1) / sum(fluxmatnorm(ibin,1,:))
        p2 = 1.0d0 - p1
        fluxmat1(ibin,1) = fluxmat1(ibin,1) / p1    !normalize part 2 for mat specific
        fluxmat2(ibin,1) = fluxmat2(ibin,1) / p2    !normalize part 2 for mat specific
        stocMC_fluxmat1(ibin,icase,1) = fluxmat1(ibin,1) !store mean 
        stocMC_fluxmat2(ibin,icase,1) = fluxmat2(ibin,1) !store mean
      enddo
    endif

  elseif(MCcases(icase)=='atmixMC') then
    if(flfluxplotall) then
      do ibin=1,fluxnumcells  !mean vals, no var, no mat specific
        stocMC_fluxall(ibin,icase,1) = fluxall(ibin,1)   !store mean
      enddo
    endif

  endif

  !leakage pdfs
  call MCLeakage_pdfbinprint( icase )         !bin and print pdf of leakage values

  !calculate FOMs
  call FOM_calculation( icase )

  end subroutine stocMC_stats



  subroutine pos_all_pop( realzarray,posrealzarray )
  !This subroutine takes an array of data and collapses it to the same array 
  !without data corresponding to realizations which at part of the domain are negative.
  use genRealzvars, only: numRealz, posRealz, numPosRealz

  real(8) :: realzarray(:)
  real(8), allocatable :: posrealzarray(:)
  integer :: j,posj

  !allocate
  if(allocated(posrealzarray)) deallocate(posrealzarray)
  allocate(posrealzarray(numPosRealz))
  posrealzarray = 0.0d0

  !populate
  posj = 0
  do j=1,numRealz
    if(posRealz(j)==1) then
      posj = posj + 1
      posrealzarray(posj) = realzarray(j)
    endif
  enddo

  end subroutine pos_all_pop




  subroutine MCfluxPrint
  !This subroutine prints the results of MC transport methods to '.out' files
  !in order to be printed to the screen.  It also clears out previous plot files.
  use MCvars, only: stocMC_fluxall, stocMC_fluxmat1, stocMC_fluxmat2, numPosMCmeths, &
                    fluxfaces, fluxnumcells, MCcases, MCcaseson, pltflux, pltmatflux, flfluxplot
  integer :: icase, ibin

  call system("test -e plots/fluxplots/radMC_fluxall.out   && rm plots/fluxplots/radMC_fluxall.out")
  call system("test -e plots/fluxplots/radMC_fluxmat.out   && rm plots/fluxplots/radMC_fluxmat.out")
  call system("test -e plots/fluxplots/radWood_fluxall.out && rm plots/fluxplots/radWood_fluxall.out")
  call system("test -e plots/fluxplots/radWood_fluxmat.out && rm plots/fluxplots/radWood_fluxmat.out")
  call system("test -e plots/fluxplots/KLWood_fluxall.out  && rm plots/fluxplots/KLWood_fluxall.out")
  call system("test -e plots/fluxplots/KLWood_fluxmat.out  && rm plots/fluxplots/KLWood_fluxmat.out")
  call system("test -e plots/fluxplots/LPMC_fluxall.out    && rm plots/fluxplots/LPMC_fluxall.out")
  call system("test -e plots/fluxplots/LPMC_fluxmat.out    && rm plots/fluxplots/LPMC_fluxmat.out")
  call system("test -e plots/fluxplots/atmixMC_fluxall.out && rm plots/fluxplots/atmixMC_fluxall.out")
  call system("test -e plots/fluxplots/atmixMC_fluxmat.out && rm plots/fluxplots/atmixMC_fluxmat.out")

  370 format("#cell center,       ave flux,        flux dev")
  371 format(f15.7,f15.7,f15.7)

  372 format("#cell center,       ave mat1 flux,   flux mat1 dev,   ave mat2 flux,   flux mat2 dev")
  373 format(f15.7,f15.7,f15.7,f15.7,f15.7)

  374 format("#cell center,       ave flux")
  375 format(f15.7,f15.7)

  376 format("#cell center,       ave mat1 flux,   ave mat2 flux")

  do icase=1,numPosMCmeths
    if(MCcaseson(icase)==1 .and. flfluxplot) then
      select case (MCcases(icase))
        case ("radMC")
          if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
            open(unit=24, file="radMC_fluxall.out")
            write(24,370)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1),sqrt(stocMC_fluxall(ibin,icase,2))
            enddo
            close(unit=24)
          endif
          if(pltmatflux=='plot' .or. pltmatflux=='preview') then
            open(unit=24, file="radMC_fluxmat.out")
            write(24,372)
            do ibin=1,fluxnumcells
              write(24,373) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxmat1(ibin,icase,1),sqrt(stocMC_fluxmat1(ibin,icase,2)),&
                            stocMC_fluxmat2(ibin,icase,1),sqrt(stocMC_fluxmat2(ibin,icase,2))
            enddo
            close(unit=24)
          endif

        case ("radWood")
          if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
            open(unit=24, file="radWood_fluxall.out")
            write(24,370)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1),sqrt(stocMC_fluxall(ibin,icase,2))
            enddo
            close(unit=24)
          endif
          if(pltmatflux=='plot' .or. pltmatflux=='preview') then
            open(unit=24, file="radWood_fluxmat.out")
            write(24,372)
            do ibin=1,fluxnumcells
              write(24,373) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxmat1(ibin,icase,1),sqrt(stocMC_fluxmat1(ibin,icase,2)),&
                            stocMC_fluxmat2(ibin,icase,1),sqrt(stocMC_fluxmat2(ibin,icase,2))
            enddo
            close(unit=24)
          endif

        case ("KLWood")
          if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
            open(unit=24, file="KLWood_fluxall.out")
            write(24,370)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1),sqrt(stocMC_fluxall(ibin,icase,2))
            enddo
            close(unit=24)
          endif
          if(pltmatflux=='plot' .or. pltmatflux=='preview') then
            open(unit=24, file="KLWood_fluxmat.out")
            write(24,370)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1),sqrt(stocMC_fluxall(ibin,icase,2))
            enddo
            close(unit=24)
          endif

        case ("LPMC")
          if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
            open(unit=24, file="LPMC_fluxall.out")
            write(24,371)
            do ibin=1,fluxnumcells
              write(24,375) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1)
            enddo
            close(unit=24)
          endif
          if(pltmatflux=='plot' .or. pltmatflux=='preview') then
            open(unit=24, file="LPMC_fluxmat.out")
            write(24,376)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxmat1(ibin,icase,1),&
                            stocMC_fluxmat2(ibin,icase,1)
            enddo
            close(unit=24)
          endif

        case ("atmixMC")
          if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
            open(unit=24, file="atmixMC_fluxall.out")
            write(24,374)
            do ibin=1,fluxnumcells
              write(24,375) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1)
            enddo
            close(unit=24)
          endif
          if(pltmatflux=='plot' .or. pltmatflux=='preview') then
            open(unit=24, file="atmixMC_fluxmat.out")
            write(24,374)
            do ibin=1,fluxnumcells
              write(24,375) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1)
            enddo
            close(unit=24)
          endif

        case ("GaussKL")
          if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
            open(unit=24, file="GaussKL_fluxall.out")
            write(24,370)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1),sqrt(stocMC_fluxall(ibin,icase,2))
            enddo
            close(unit=24)
          endif
          if(pltmatflux=='plot' .or. pltmatflux=='preview') then
            open(unit=24, file="GaussKL_fluxmat.out")
            write(24,370)
            do ibin=1,fluxnumcells
              write(24,371) (fluxfaces(ibin+1)+fluxfaces(ibin))/2.0d0,&
                            stocMC_fluxall(ibin,icase,1),sqrt(stocMC_fluxall(ibin,icase,2))
            enddo
            close(unit=24)
          endif

      end select
    endif
  enddo

  call system("test -e radMC_fluxall.out   && mv radMC_fluxall.out plots/fluxplots")
  call system("test -e radMC_fluxmat.out   && mv radMC_fluxmat.out plots/fluxplots")
  call system("test -e radWood_fluxall.out && mv radWood_fluxall.out plots/fluxplots")
  call system("test -e radWood_fluxmat.out && mv radWood_fluxmat.out plots/fluxplots")
  call system("test -e KLWood_fluxall.out  && mv KLWood_fluxall.out plots/fluxplots")
  call system("test -e KLWood_fluxmat.out  && mv KLWood_fluxmat.out plots/fluxplots")
  call system("test -e LPMC_fluxall.out    && mv LPMC_fluxall.out plots/fluxplots")
  call system("test -e LPMC_fluxmat.out    && mv LPMC_fluxmat.out plots/fluxplots")
  call system("test -e atmixMC_fluxall.out && mv atmixMC_fluxall.out plots/fluxplots")
  call system("test -e atmixMC_fluxmat.out && mv atmixMC_fluxmat.out plots/fluxplots")
  call system("test -e GaussKL_fluxall.out  && mv GaussKL_fluxall.out plots/fluxplots")
  call system("test -e GaussKL_fluxmat.out  && mv GaussKL_fluxmat.out plots/fluxplots")


  end subroutine MCfluxPrint




  subroutine MCfluxPlot
  !Builds gnus and plots with them.  Copies gnubuilding tools to local temporary 'gnu' directory,
  !which is deleted at the end.  Builds based on options like title type, plotting lines, and
  !pause or preview.  Can perform three builds, one for material irrespective flux tallying, and 
  !one for each material using material respective flux plotting.
  use MCvars, only: pltflux, radMC, radWood, KLWood, LPMC, atmixMC, pltfluxtype, pltmatflux, &
                    GaussKL

  !Clean from previous runs
  call system("test -e plots/fluxplots/fluxall.ps && rm plots/fluxplots/fluxall.ps")
  call system("test -e plots/fluxplots/fluxmat1.ps && rm plots/fluxplots/fluxmat1.ps")
  call system("test -e plots/fluxplots/fluxmat2.ps && rm plots/fluxplots/fluxmat2.ps")

  call system("test -e plots/fluxplots/fluxall.pdf && rm plots/fluxplots/fluxall.pdf")
  call system("test -e plots/fluxplots/fluxmat1.pdf && rm plots/fluxplots/fluxmat1.pdf")
  call system("test -e plots/fluxplots/fluxmat2.pdf && rm plots/fluxplots/fluxmat2.pdf")

  call system("test -e plots/fluxplots/fluxall.gnu && rm plots/fluxplots/fluxall.gnu")
  call system("test -e plots/fluxplots/fluxmat1.gnu && rm plots/fluxplots/fluxmat1.gnu")
  call system("test -e plots/fluxplots/fluxmat2.gnu && rm plots/fluxplots/fluxmat2.gnu")

  !Bring building blocks near (shorter lines of code)
  call system("cp -r plots/fluxplots/gnubuilding gnu")

  if(pltflux(1)=='plot' .or. pltflux(1)=='preview') then
    !First part of file and title
    if(pltfluxtype=='track') then
      call system("cat gnu/gen1.txt gnu/titlealltrack.txt > gnu/tempold.txt")
    elseif(pltfluxtype=='point') then
      call system("cat gnu/gen1.txt gnu/titleallpoint.txt > gnu/tempold.txt")
    endif

    !Add second general part
    call system("cat gnu/tempold.txt gnu/gen2.txt > gnu/tempnew.txt")
    call system("mv gnu/tempnew.txt gnu/tempold.txt")

    !Add everything that is to be plotted (options)
    if(radMC=='yes') then
      call system("cat gnu/tempold.txt gnu/radMCall.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(radWood=='yes') then
      call system("cat gnu/tempold.txt gnu/radWoodall.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(KLWood=='yes') then
      call system("cat gnu/tempold.txt gnu/KLWoodall.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(LPMC=='yes') then
      call system("cat gnu/tempold.txt gnu/LPMCall.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(atmixMC=='yes') then
      call system("cat gnu/tempold.txt gnu/atmixMCall.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(GaussKL=='yes') then
      call system("cat gnu/tempold.txt gnu/GaussKLall.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif

    !Add third general part
    call system("cat gnu/tempold.txt gnu/gen3.txt > gnu/tempnew.txt")
    call system("mv gnu/tempnew.txt gnu/tempold.txt")

    !Add either pause or print final line (options)
    if( pltflux(1)=='preview' ) then
      call system("cat gnu/tempold.txt gnu/preview.txt > gnu/tempnew.txt")
    elseif( pltflux(1)=='plot' ) then
      call system("cat gnu/tempold.txt gnu/plot.txt > gnu/tempnew.txt")
    endif

    !Update gnu to final name
    call system("mv gnu/tempnew.txt gnu/fluxall.gnu")
    !plot
    call system("gnuplot gnu/fluxall.gnu")

    !Store and clean up
    call system("mv flux.ps fluxall.ps")
    call system("ps2pdf fluxall.ps")
    call system("mv fluxall.ps fluxall.pdf gnu/fluxall.gnu plots/fluxplots/")

  endif !end pltflux(1)

  if(pltmatflux=='plot' .or. pltmatflux=='preview') then
    !First part of file and title
    if(pltfluxtype=='track') then
      call system("cat gnu/gen1.txt gnu/titlemat1track.txt > gnu/tempold.txt")
    elseif(pltfluxtype=='point') then
      call system("cat gnu/gen1.txt gnu/titlemat1point.txt > gnu/tempold.txt")
    endif

    !Add second general part
    call system("cat gnu/tempold.txt gnu/gen2.txt > gnu/tempnew.txt")
    call system("mv gnu/tempnew.txt gnu/tempold.txt")

    !Add everything that is to be plotted (options)
    if(radMC=='yes') then
      call system("cat gnu/tempold.txt gnu/radMCmat1.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(radWood=='yes') then
      call system("cat gnu/tempold.txt gnu/radWoodmat1.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(KLWood=='yes') then
      call system("cat gnu/tempold.txt gnu/KLWoodmat1.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(LPMC=='yes') then
      call system("cat gnu/tempold.txt gnu/LPMCmat1.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(atmixMC=='yes') then
      call system("cat gnu/tempold.txt gnu/atmixMCmat1.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(GaussKL=='yes') then
      call system("cat gnu/tempold.txt gnu/GaussKLmat1.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif

    !Add third general part
    call system("cat gnu/tempold.txt gnu/gen3.txt > gnu/tempnew.txt")
    call system("mv gnu/tempnew.txt gnu/tempold.txt")

    !Add either pause or print final line (options)
    if( pltflux(1)=='preview' ) then
      call system("cat gnu/tempold.txt gnu/preview.txt > gnu/tempnew.txt")
    elseif( pltflux(1)=='plot' ) then
      call system("cat gnu/tempold.txt gnu/plot.txt > gnu/tempnew.txt")
    endif

    !Update gnu to final name
    call system("mv gnu/tempnew.txt gnu/fluxmat1.gnu")
    !plot
    call system("gnuplot gnu/fluxmat1.gnu")

    !Store and clean up
    call system("mv flux.ps fluxmat1.ps")
    call system("ps2pdf fluxmat1.ps")
    call system("mv fluxmat1.ps fluxmat1.pdf gnu/fluxmat1.gnu plots/fluxplots/")
    
  endif !end pltmatflux for material 1

  if(pltmatflux=='plot' .or. pltmatflux=='preview') then
    !First part of file and title
    if(pltfluxtype=='track') then
      call system("cat gnu/gen1.txt gnu/titlemat2track.txt > gnu/tempold.txt")
    elseif(pltfluxtype=='point') then
      call system("cat gnu/gen1.txt gnu/titlemat2point.txt > gnu/tempold.txt")
    endif

    !Add second general part
    call system("cat gnu/tempold.txt gnu/gen2.txt > gnu/tempnew.txt")
    call system("mv gnu/tempnew.txt gnu/tempold.txt")

    !Add everything that is to be plotted (options)
    if(radMC=='yes') then
      call system("cat gnu/tempold.txt gnu/radMCmat2.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(radWood=='yes') then
      call system("cat gnu/tempold.txt gnu/radWoodmat2.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(KLWood=='yes') then
      call system("cat gnu/tempold.txt gnu/KLWoodmat2.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(LPMC=='yes') then
      call system("cat gnu/tempold.txt gnu/LPMCmat2.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(atmixMC=='yes') then
      call system("cat gnu/tempold.txt gnu/atmixMCmat2.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif
    if(GaussKL=='yes') then
      call system("cat gnu/tempold.txt gnu/GaussKLmat2.txt > gnu/tempnew.txt")
      call system("mv gnu/tempnew.txt gnu/tempold.txt")
    endif

    !Add third general part
    call system("cat gnu/tempold.txt gnu/gen3.txt > gnu/tempnew.txt")
    call system("mv gnu/tempnew.txt gnu/tempold.txt")

    !Add either pause or print final line (options)
    if( pltflux(1)=='preview' ) then
      call system("cat gnu/tempold.txt gnu/preview.txt > gnu/tempnew.txt")
    elseif( pltflux(1)=='plot' ) then
      call system("cat gnu/tempold.txt gnu/plot.txt > gnu/tempnew.txt")
    endif

    !Update gnu to final name
    call system("mv gnu/tempnew.txt gnu/fluxmat2.gnu")
    !plot
    call system("gnuplot gnu/fluxmat2.gnu")

    !Store and clean up
    call system("mv flux.ps fluxmat2.ps")
    call system("ps2pdf fluxmat2.ps")
    call system("mv fluxmat2.ps fluxmat2.pdf gnu/fluxmat2.gnu plots/fluxplots/")
    
  endif !end pltmatflux for material 2

  !Final clean up
  call system("rm -r gnu")

  end subroutine MCfluxPlot




  subroutine MCallocate( icase,tnumparts,tnumRealz )
  !This subroutine allocates and initializes variables that will be passed
  !through generic MCtransport subroutine.  These values will later be
  !stored in different arrays so that the variables can be re-used in
  !MCtransport if multiple cases were selected.
  use genRealzvars, only: numRealz, flprint, flGBgeom, GBsigave, GBsigvar, &
                          GBscatrat, GBlamc, GBs, s, sigave, lamc, scatrat, CoExp, &
                          posRealz, numPosRealz
  use MCvars, only: transmit, reflect, absorb, radtrans_int, MCcases, &
                    numpnSamp, areapnSamp, disthold, Wood_rej, LPamMCsums, &
                    numParts, LPamnumParts, fluxnumcells, fluxall, fluxmat1, &
                    fluxmat2, pltflux, pltmatflux, flfluxplotall, flfluxplotmat, &
                    fluxmatnorm, refsig, refsigMode, negwgtsigs, negwgtbinnum, &
                    nwvalsperbin, flfluxplot, fluxfaces
  use KLvars, only: flmatbasedxs
  integer :: icase,tnumParts,tnumRealz,i

  flprint = .false.

  !number of realizations allocation
  if( MCcases(icase)=='LPMC' .or. MCcases(icase)=='atmixMC' ) then
    tnumRealz = 1
  else
    tnumRealz = numRealz
  endif

  !Gauss-based input set (or not)
  if(MCcases(icase)=='GaussKL' .and. flGBgeom) then
    sigave       = GBsigave
    CoExp        = GBsigvar
    scatrat(1)   = GBscatrat
    lamc         = GBlamc
    s            = GBs
    flmatbasedxs =.false.
    if(flfluxplot) then
      if(allocated(fluxfaces)) deallocate(fluxfaces)
      allocate(fluxfaces(fluxnumcells+1))
      fluxfaces = 0.0d0
      do i=1,fluxnumcells+1
        fluxfaces(i) = (s/fluxnumcells) * (i-1)
      enddo
    endif
  endif

  !current tally allocations
  if(MCcases(icase)=='radMC'  .or. MCcases(icase)=='radWood'.or. &
     MCcases(icase)=='KLWood' .or. MCcases(icase)=='WAMC'   .or. &
     MCcases(icase)=='GaussKL'                                     ) then
    if(.not.allocated(transmit)) allocate(transmit(numRealz))
    if(.not.allocated(reflect))  allocate(reflect(numRealz))
    if(.not.allocated(absorb))   allocate(absorb(numRealz))
    transmit     = 0.0d0
    reflect      = 0.0d0
    absorb       = 0.0d0
    radtrans_int = 0
    tnumParts    = numParts
  elseif(MCcases(icase)=='LPMC' .or. MCcases(icase)=='atmixMC') then
    if(.not.allocated(LPamMCsums)) allocate(LPamMCsums(3))
    LPamMCsums   = 0.0d0
    tnumParts    = LPamnumParts
  endif

  !rejection tally allocations
  if(MCcases(icase)=='radWood' .or. MCcases(icase)=='KLWood' .or. &
     MCcases(icase)=='GaussKL'                                     ) Wood_rej = 0

  !negative xs transport tally allocations
  if(MCcases(icase)=='KLWood' .or. MCcases(icase)=='WAMC' .or. &
     MCcases(icase)=='GaussKL'                                 ) then
    if(allocated(posRealz)) deallocate(posRealz)
    allocate(posRealz(numRealz))
    posRealz   =0
    numPosRealz=0
    numpnSamp  =0
    areapnSamp =0.0d0
    disthold   =0.0d0
  else  !include this so that 'posRealz' is defined for other methods
    if(allocated(posRealz)) deallocate(posRealz)
    allocate(posRealz(numRealz))
  endif

  !set initial WAMC reference sigma value
  if(MCcases(icase)=='WAMC') then
    call setrefsig()
    if(refsigMode==2 .or. refsigMode==3) then !sig samples from which to create sigrefs
      allocate(negwgtsigs(negwgtbinnum*nwvalsperbin+1,3))
      negwgtsigs = 0d0
    endif
  endif

  !flux tally allocations
  flfluxplotall = .false.
  flfluxplotmat = .false.
  select case(MCcases(icase))
    case ("radMC")
      if(pltflux(1)=='plot' .or. pltflux(1)=='preview') flfluxplotall = .true.
      if(pltmatflux=='plot' .or. pltmatflux=='preview') flfluxplotmat = .true.
    case ("radWood")
      if(pltflux(1)=='plot' .or. pltflux(1)=='preview') flfluxplotall = .true.
      if(pltmatflux=='plot' .or. pltmatflux=='preview') flfluxplotmat = .true.
    case ("KLWood")
      if(pltflux(1)=='plot' .or. pltflux(1)=='preview' .or.&
         pltmatflux=='plot' .or. pltmatflux=='preview') flfluxplotall = .true.
    case ("LPMC")
      if(pltflux(1)=='plot' .or. pltflux(1)=='preview') flfluxplotall = .true.
      if(pltmatflux=='plot' .or. pltmatflux=='preview') flfluxplotmat = .true.
    case ("atmixMC")
      if(pltflux(1)=='plot' .or. pltflux(1)=='preview' .or.&
         pltmatflux=='plot' .or. pltmatflux=='preview') flfluxplotall = .true.
    case ("GaussKL")
      if(pltflux(1)=='plot' .or. pltflux(1)=='preview' .or.&
         pltmatflux=='plot' .or. pltmatflux=='preview') flfluxplotall = .true.
  end select
  if(flfluxplotall) then
    if(allocated(fluxall)) deallocate(fluxall)
    allocate(fluxall(fluxnumcells,tnumRealz))
    fluxall = 0.0d0
  endif
  if(flfluxplotmat) then
    if(allocated(fluxmat1)) deallocate(fluxmat1)
    if(allocated(fluxmat2)) deallocate(fluxmat2)
    allocate(fluxmat1(fluxnumcells,tnumRealz))
    allocate(fluxmat2(fluxnumcells,tnumRealz))
    fluxmat1 = 0.0d0
    fluxmat2 = 0.0d0
    if(allocated(fluxmatnorm)) deallocate(fluxmatnorm)
    allocate(fluxmatnorm(fluxnumcells,numRealz,2))
    fluxmatnorm = 0.0d0
  endif    

  end subroutine MCallocate


  subroutine setrefsig()
  !This subruotine sets the value of refsig according to the user's choice of different methods.
  !It can be used to set the value only once, or during transport simulations.
  use genRealzvars, only: sig, scatrat, P
  use MCvars, only: refsigMode,refsig,userrefsig, maxratio
  real(8) :: cursiga, cursigs

  refsig = 0.0d0
  select case(refsigMode)
    case (1)                           !set to user defined value
      refsig = userrefsig
    case (2)                           !largest sigma value
      maxratio = userrefsig
    case (3)
      maxratio = userrefsig
  end select

  end subroutine setrefsig


  subroutine MCLeakage_pdfbinprint( icase )
  !This subroutine bins leakage data, and prints this data to files to later be plotted
  !for methods which contain realizations (radMC, radWood, KLWood).
  use genRealzvars, only: numRealz
  use MCvars,       only: radMCbinplot, radWoodbinplot, KLWoodbinplot, GaussKLbinplot, &
                          reflect, transmit, MCcases, MCcaseson

  integer :: icase
  real(8) :: smrefl,lgrefl,smtran,lgtran,boundbuff

  !this nasty if statement decides whether to proceed at all.
  !from here only MCcases is needed.
  if( (MCcaseson(icase) == 1      .and. MCcases(icase) =='radMC'      .and. &  !!!WAMCchange
      (radMCbinplot     =='plot'  .or.  radMCbinplot   =='preview'))  .or.  &
      (MCcaseson(icase) == 1      .and. MCcases(icase) =='radWood'    .and. &
      (radWoodbinplot   =='plot'  .or.  radWoodbinplot =='preview'))  .or.  &
      (MCcaseson(icase) == 1      .and. MCcases(icase) =='KLWood'     .and. &
      (KLWoodbinplot    =='plot'  .or.  KLWoodbinplot  =='preview'))  .or.  &
      (MCcaseson(icase) == 1      .and. MCcases(icase) =='GaussKL'    .and. &
      (GaussKLbinplot    =='plot' .or.  GaussKLbinplot  =='preview'))          ) then

    !find reflection binning/plotting bounds
    smrefl = 1d0
    lgrefl = 0d0
    smrefl = min(smrefl,minval(reflect))
    lgrefl = max(lgrefl,maxval(reflect))
    boundbuff = (lgrefl-smrefl)/8d0
    smrefl = merge(smrefl-boundbuff,0d0,smrefl-boundbuff>0d0)
    lgrefl = merge(lgrefl+boundbuff,1d0,lgrefl+boundbuff<1d0)
    smrefl = smrefl - 0.0000001d0 !these to ensure binning works in case of opaque or transparent
    lgrefl = lgrefl + 0.0000001d0

    !find tranmission binning/plotting bounds
    smtran = 1d0
    lgtran = 0d0
    smtran = min(smtran,minval(transmit))
    lgtran = max(lgtran,maxval(transmit))
    boundbuff = (lgtran-smtran)/8d0
    smtran = merge(smtran-boundbuff,0d0,smtran-boundbuff>0d0)
    lgtran = merge(lgtran+boundbuff,1d0,lgtran+boundbuff<1d0)
    smtran = smtran - 0.0000001d0 !these to ensure binning works in case of opaque or transparent
    lgtran = lgtran + 0.0000001d0

    !remove old data files (do only first time)
    if(sum(MCcaseson(1:icase))==1) then
      call system("test -e plots/tranreflprofile/radMCtranreflprofile.txt && rm plots/tranreflprofile/radMCtranreflprofile.txt")
      call system("test -e plots/tranreflprofile/radWoodtranreflprofile.txt && rm plots/tranreflprofile/radWoodtranreflprofile.txt")
      call system("test -e plots/tranreflprofile/KLWoodtranreflprofile.txt && rm plots/tranreflprofile/KLWoodtranreflprofile.txt")
      call system("test -e plots/tranreflprofile/GaussKLtranreflprofile.txt && rm plots/tranreflprofile/GaussKLtranreflprofile.txt")
    endif

    !bin and print data
    call radtrans_bin( smrefl,lgrefl,smtran,lgtran,icase ) 

    !bin and print data, give printed data files appropriate name
    select case (MCcases(icase))
      case("radMC")
        call system("mv tranreflprofile.txt radMCtranreflprofile.txt")
        call system("mv radMCtranreflprofile.txt plots/tranreflprofile")
      case("radWood")
        call system("mv tranreflprofile.txt radWoodtranreflprofile.txt")
        call system("mv radWoodtranreflprofile.txt plots/tranreflprofile")
      case("KLWood")
        call system("mv tranreflprofile.txt KLWoodtranreflprofile.txt")
        call system("mv KLWoodtranreflprofile.txt plots/tranreflprofile")
      case("LPMC")
      case("atmixMC")
      case("WAMC")  !!!WAMCchange
        call system("mv tranreflprofile.txt KLWoodtranreflprofile.txt")
        call system("mv KLWoodtranreflprofile.txt plots/tranreflprofile")
      case("GaussKL")
        call system("mv tranreflprofile.txt GaussKLtranreflprofile.txt")
        call system("mv GaussKLtranreflprofile.txt plots/tranreflprofile")
    end select

  endif

  end subroutine MCLeakage_pdfbinprint




  subroutine radtrans_bin( smrefl,lgrefl,smtran,lgtran,icase )
  !Heart of radtrans_resultplot, for methods with realizations,
  !loads leakage values to bins (pdf), and prints to generic text file
  use genRealzvars, only: numRealz, numPosRealz, posRealz
  use MCvars, only: trprofile_binnum, reflect, transmit, flnegxs, MCcases
  integer :: icase
  real(8) :: smrefl,lgrefl,smtran,lgtran

  !local vars
  integer :: ibin
  integer,allocatable,dimension(:) :: reflcounts,trancounts  
  real(8),allocatable,dimension(:) :: reflprob,  tranprob
  real(8),allocatable,dimension(:) :: reflbounds,tranbounds
  real(8),allocatable,dimension(:) :: posreflect,postransmit

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

  if(.not.flnegxs .and. (MCcases(icase)=='KLWood' .or. &
     MCcases(icase)=='WAMC' .or. MCcases(icase)=='GaussKL') ) then
    !only use leakage data for pos realizations (if negative rejected due to method or input)
    if(allocated(posreflect)) deallocate(posreflect)
    if(allocated(postransmit)) deallocate(postransmit)
    allocate(posreflect(numPosRealz))
    allocate(postransmit(numPosRealz))
    call pos_all_pop( reflect,posreflect )
    call pos_all_pop( transmit,postransmit )

    !actually store in binned fashion
    call store_in_bins( smrefl,lgrefl,trprofile_binnum,reflcounts,reflbounds,posreflect,numPosRealz )
    call store_in_bins( smtran,lgtran,trprofile_binnum,trancounts,tranbounds,postransmit,numPosRealz )

    tranprob = real(trancounts,8)/numPosRealz/((lgtran-smtran)/(trprofile_binnum-1))
    reflprob = real(reflcounts,8)/numPosRealz/((lgrefl-smrefl)/(trprofile_binnum-1))
  else
    !actually store in binned fashion
    call store_in_bins( smrefl,lgrefl,trprofile_binnum,reflcounts,reflbounds,reflect,numRealz )
    call store_in_bins( smtran,lgtran,trprofile_binnum,trancounts,tranbounds,transmit,numRealz )

    tranprob = real(trancounts,8)/numRealz/((lgtran-smtran)/(trprofile_binnum-1))
    reflprob = real(reflcounts,8)/numRealz/((lgrefl-smrefl)/(trprofile_binnum-1))
  endif

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



  subroutine MCLeakage_pdfplot
  !This subroutine plots MC Leakage value pdfs from data files generated
  !in MCLeakage_pdfbinprint.
  use MCvars,       only: radMCbinplot, radWoodbinplot, KLWoodbinplot, GaussKLbinplot
  logical :: flplot = .false.

  !preview if at least one chose this, otherwise simply plot
  if(    radMCbinplot  =='preview' .or. radWoodbinplot=='preview' .or. &
         KLWoodbinplot =='preview' .or. GaussKLbinplot=='preview'      ) then
    call system("gnuplot plots/tranreflprofile/tranprofile.p.gnu")
    call system("gnuplot plots/tranreflprofile/reflprofile.p.gnu")
    flplot = .true.
  elseif(radMCbinplot  =='plot' .or. radWoodbinplot=='plot' .or. &
         KLWoodbinplot =='plot' .or. GaussKLbinplot=='plot'         ) then
    call system("gnuplot plots/tranreflprofile/tranprofile.gnu")
    call system("gnuplot plots/tranreflprofile/reflprofile.gnu")
    flplot = .true.
  endif
  !convert and store
  if(flplot) then
    call system("ps2pdf tranprofile.ps")
    call system("ps2pdf reflprofile.ps")
    call system("mv tranprofile.ps tranprofile.pdf plots/tranreflprofile")
    call system("mv reflprofile.ps reflprofile.pdf plots/tranreflprofile")
  endif

  end subroutine MCLeakage_pdfplot



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



  function internal_init_i( position )
  !Finds cell for 'radMC' that 'position' is located in.
  use genRealzvars, only: matLength
  real(8) :: position
  integer :: i,internal_init_i

  do i=1,size(matLength)-1
    if( matLength(i)<=position .AND. matLength(i+1)>position ) then
      internal_init_i = i
      exit
    endif
  enddo
  end function internal_init_i



  subroutine Woodnegstats(icase)
  use genRealzvars, only: numRealz, numPosRealz
  use MCvars, only: numpnSamp, areapnSamp, fldistneg, flnegxs, numcSamp, MCcases

  integer :: icase
  real(8) :: pos,neg

  open(unit=100,file="Woodnegstats.out")

  if((MCcases(icase)=='KLWood' .and. flnegxs) .or. &
     (MCcases(icase)=='WAMC'   .and. flnegxs) .or. &
     (MCcases(icase)=='GaussKL'.and. flnegxs)        ) then
    606 format("--Negative xs stats (",A8," ), keep neg xs: ",L,", neg smoothing: ",L," --")

    600 format("  Neg realz   : ",f8.5,"%, ",i21," /",i21)
    601 format("  Neg samples : ",f8.5,"%, ",i21," /",i21)
    602 format("  Neg area    : ",f8.5,"%")
    605 format("  Neg scatrat : ",f8.5,"%, ",i21," /",i21)
    603 format("  Ave neg samp: ",f11.4,"   Ave pos samp: ",f11.4)
    604 format("  Max neg samp: ",f11.4,"   Max pos samp: ",f11.4)

    write(100,606) MCcases(icase),flnegxs,fldistneg
    write(100,600) real(numRealz-numPosRealz,8)/real(numRealz,8)*100d0,numRealz-numPosRealz,numRealz
    pos = real(numpnSamp(1),8)
    neg = real(numpnSamp(2),8)
    write(100,601) neg/(pos+neg)*100d0,numpnSamp(2),numpnSamp(1)+numpnSamp(2)
    write(100,602) -areapnSamp(2)/(areapnSamp(1)-areapnSamp(2))*100d0
    write(100,605) real(numcSamp(1))/real(numcSamp(2))*100,numcSamp(1),numcSamp(2)
    write(100,603) areapnSamp(2)/numpnSamp(2),areapnSamp(1)/numpnSamp(1)
    write(100,604) areapnSamp(4),areapnSamp(3)
    write(100,*)
  elseif((MCcases(icase)=='KLWood' .and. .not.flnegxs) .or. &
         (MCcases(icase)=='WAMC'   .and. .not.flnegxs) .or. &
         (MCcases(icase)=='GaussKL'.and. .not.flnegxs)        ) then
    610 format("--Negative xs stats (",A8," ), keep neg xs: ",L," --")

    write(100,610) MCcases(icase),flnegxs
    write(100,600) real(numRealz-numPosRealz,8)/real(numRealz,8)*100d0,numRealz-numPosRealz,numRealz
    write(100,*)

  else
    write(100,*)
  endif

  close(unit=100)
  !if file already exists, cat new and old, else create, and remove old
  call system("test -e texts/Woodnegstats.out && cat texts/Woodnegstats.out Woodnegstats.out > temp.out")
  call system("test -e texts/Woodnegstats.out && mv temp.out texts/Woodnegstats.out")
  call system("test -e texts/Woodnegstats.out || cp Woodnegstats.out texts")
  call system("rm Woodnegstats.out")

  end subroutine Woodnegstats




  subroutine MCprintstats
  !This subroutine prints reflection, transmission, and absorption stats to a '.out' file,
  !then prints that file to the screen for user friendliness.
  !Stats are from Adams89, Brantley11, and those generated here.
  use genRealzvars, only: Adamscase, flGBgeom
  use KLvars, only: flMarkov, flGauss
  use MCvars, only: ABreflection, ABtransmission, rodOrplanar, stocMC_reflection, &
                    stocMC_transmission, stocMC_absorption, MCcases, MCcaseson, &
                    numPosMCmeths, LPMC, atmixMC, GaussKL
  integer :: icase

  320 format(" |AdamsMC:  |",f7.4,"   +-",f8.4,"     |",f7.4,"   +-",f8.4," |")
  321 format(" |BrantMC:  |",f8.5,"                 |",f8.5,"             |")
  322 format(" |AdamsLP:  |",f7.4,"                  |",f7.4,"              |")
  323 format(" |BrantLP:  |",f8.5,"                 |",f8.5,"             |")
  324 format(" |BrAtMix:  |",f8.5,"                 |",f8.5,"             |")
  325 format(" |-case:",f3.1,"-|---- Reflection and Transmission Results ------|")

  326 format(" |radMC  :  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  327 format(" |radWood:  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  328 format(" |KLWood :  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  331 format(" |WAMC   :  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  332 format(" |GaussKL:  |",f8.5,"  +-",f9.5,"    |",f8.5,"  +-",f9.5,"|")
  329 format(" |LPMC   :  |",f8.5,"                 |",f8.5,"             |")
  330 format(" |atmixMC:  |",f8.5,"                 |",f8.5,"             |")

  !print to file
  open(unit=100,file="MCleakage.out")

  if(flMarkov) then

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

  !print my solutions for radMC, radWood, KLWood, WAMC, GaussKL (if Markov-based geom)
  do icase = 1,numPosMCmeths
    if(MCcaseson(icase)==1) then
      if(MCcases(icase)=='radMC')   write(100,326) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))

      if(MCcases(icase)=='radWood') write(100,327) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))

      if(MCcases(icase)=='KLWood')  write(100,328) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))

      if(MCcases(icase)=='WAMC')  write(100,331) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))

      if(MCcases(icase)=='GaussKL' .and. .not.flGBgeom)  write(100,332) stocMC_reflection(icase,1),&
      sqrt(stocMC_reflection(icase,2)),stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))

    endif
  enddo

  !print for formatting if any LP solutions printed
  if(Adamscase/=0 .or. LPMC=='yes') &
    write(100,*) "|----------|-------------------------|---------------------|"

  !print benchmark LP solutions
  if(Adamscase/=0 .and. flMarkov) then
    write(100,322) ABreflection(1,2),ABtransmission(1,2)
    if(rodOrplanar=='planar') write(100,323) ABreflection(1,4),ABtransmission(1,4)
  endif

  !print my solution for LPMC
  do icase = 1,numPosMCmeths
    if(MCcaseson(icase)==1 .and. MCcases(icase)=='LPMC') write(100,329) stocMC_reflection(icase,1),&
                                                                      stocMC_transmission(icase,1)
  enddo

  !print for formatting if any atomic mix solutions printed
  if(Adamscase/=0 .or. atmixMC=='yes') &
    write(100,*) "|----------|-------------------------|---------------------|"

  !print benchmark atomic mix solutions
  if(Adamscase/=0 .and. flMarkov) then
    if(rodOrplanar=='planar') write(100,324) ABreflection(1,5),ABtransmission(1,5)
  endif

  !print my solution for atmixMC
  do icase = 1,numPosMCmeths
    if(MCcaseson(icase)==1 .and. MCcases(icase)=='atmixMC') write(100,330) stocMC_reflection(icase,1),&
                                                                      stocMC_transmission(icase,1)
  enddo

  write(100,*) "|----------------------------------------------------------|"
  write(100,*)

  endif !(flMarkov)

  if(flGauss) then

  if(GaussKL=='yes' .and. flGBgeom) then
    icase=7
    write(100,*) "|--GB-geom-|---- Reflection and Transmission Results ------|"
    write(100,*) "|Method    | reflave      refldev    | tranave      trandev|"
    write(100,*) "|----------|-------------------------|---------------------|"
    write(100,332) stocMC_reflection(icase,1),sqrt(stocMC_reflection(icase,2)),&
                   stocMC_transmission(icase,1),sqrt(stocMC_transmission(icase,2))
    write(100,*) "|----------|-------------------------|---------------------|"
    write(100,*)
  endif

  endif !(flGauss)

  close(unit=100)

  !move file and print to screen
  call system("mv MCleakage.out texts")

  end subroutine MCprintstats


end module radtransMC
