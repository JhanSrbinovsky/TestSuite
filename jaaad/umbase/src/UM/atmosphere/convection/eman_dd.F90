#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Emanuel Downdraught code
!
      SUBROUTINE EMAN_DD(n_dp,kmax_term,nlev,trlev,ntra,kterm,l_tracer  &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p, ph, timestep                            &
     &,                      th, q, qse, tracer, precip                 &
     &,                      dthbydt, dqbydt, dtrabydt                  &
     &,                      rain, snow ,Down_flux                      &
     &    )

!
! Purpose:
!   Calculation of unsaturated Down draughts using Emanuel code.
!   Note at present the 4A convectin code deals in specific humidity
!   etc but the Emanuel code says it requires mixing ratios.
!   William Ingram appears to have taken the 4.3 version of the Emanuel
!   code and imported this into the old dynamics.
!   William Ingram chose to set various constants in a way consistent
!   with the UM. These choices mean that Emanuel's original choice
!   of a varying latent heat lv=lc+(cpv-cl)*(T-273.15)is reduced
!   to lv=lc.
!
!   Called by DEEP_CONV.
!
! Current owners of code: Convection code owner
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
      integer, intent(in) ::                                            &
     &  n_dp                                                            &
                 ! No. of deep convection points
     &, nlev                                                            &
                 ! No. of model layers
     &, trlev                                                           &
                 ! No. of tracer levels
     &, ntra                                                            &
                 ! No. of tracer fields
     &, kmax_term                                                       &
                     ! highest level reached by convection
     &, kterm(n_dp)    ! level reached by convection


      logical, intent(in) ::                                            &
    & l_tracer             ! true in tracers present

      real, intent(in)    ::                                            &
     &  exner_layer_centres(n_dp,0:nlev)                                &
                                            ! Exner
     &, exner_layer_boundaries(n_dp,0:nlev)                             &
                                            ! Exner at half level above
     &, p(n_dp,0:nlev)                                                  &
                           ! Pressure  (Pa)
     &, ph(n_dp,0:nlev)                                                 &
                           ! Pressure at half level (Pa)
     &, timestep                                                        &
                           ! timestep for model physics in seconds

     &, qse(n_dp,nlev)                                                  &
                       ! Saturation specific humidity of cloud
                       ! environment (kg/kg)

     &, q (n_dp,nlev)                                                   &
                       ! specific humidity of cloud environment(kg/kg)

     &, th (n_dp,nlev)                                                  &
                           ! theta of cloud environment(K)

     &, precip(n_dp,nlev)                                               &
                           ! Precip from updraught G-R scheme(kg/m2/s)
                           ! not units for wdtrain
     &, tracer(n_dp,trlev,ntra)  !  Tracer on model levels  (kg/kg)

!
! Arguments with intent INOUT:
!

      real, intent(inout)    ::                                         &
     &  rain(n_dp)                                                      &
                               ! rainfall at surface (kg/m**2/s)
     &, snow(n_dp)             ! snowfall at surface (kg/m**2/s)

! increments
      real, intent(inout)    ::                                         &
     &  dqbydt(n_dp,nlev)                                               &
                                 ! increments to q (kg/kg/s)
     &, dthbydt(n_dp,nlev)                                              &
                             ! increments to potential temperature (K/s)
     &, dtrabydt(n_dp,nlev,ntra) ! increments to model tracers (kg/kg/s)


!
! Arguments with intent OUT:
!
      real, intent(out)    ::                                           &
     & Down_flux(n_dp,nlev) ! Downdraught mass flux (Pa/s) diagnostic



!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

      real ::                                                           &
     & wdtrain(n_dp)                                                    &
                             ! detrained precip in level k (kg/m/s3)
     &,wt(n_dp,nlev)                                                    &
                             ! fall speed for precipitation in (kg/m/s3)
     &,evap(n_dp,nlev)                                                  &
                             ! evaporation from Downdraught (/s)
     &,water(n_dp,nlev)                                                 &
                             ! rainwater in Downdraught (lp in paper)
     &,t(n_dp,nlev)                                                     &
                             ! temperature on model levels (K)
     &,h(n_dp,nlev)                                                     &
                             ! static energy (J ?)
     &,dpinv(n_dp)                                                      &
                             ! 1/dp
     &,gz(n_dp,nlev)                                                    &
                             ! gz a form of height measured from level 1
     &,dqbydt_dd1(n_dp,nlev)                                            &
                              ! increments to q (kg/kg/s)
     &,dqbydt_dd2(n_dp,nlev)                                            &
                              ! increments to q (kg/kg/s)
     &,dqbydt_dd(n_dp,nlev)                                             &
                              ! increments to q (kg/kg/s)
     &,md(n_dp,nlev)                                                    &
                             ! Downdraught mass flux Emanuel scheme
                             ! (kg/m2/s not Pa/s)
     &,qp_down (n_dp,nlev)                                              &
                             ! specific humidity of downdraught (kg/kg)
     &,trap_down(n_dp,nlev,ntra)                                        &
                                 ! tracer values of downdraught (kg/kg)
     &, lv (n_dp,nlev)       ! latent heat release on evaporation

      real ::                                                           &
     & qsm                                                              &
                             ! mean q
     &,afac                                                             &
     &,B6,C6                                                            &
     &,revap                                                            &
     &,dhdp                                                             &
                     ! dz/dp
     &,rat                                                              &
     &,ginv, lvcp                                                       &
     &,tvx,tvy, dp                                                      &
     &,coeff,fac,qstm

      Real   ::                                                         &
     & qsum(n_dp)                                                       &
                            ! qsum  conservation check
     &,qsum2(n_dp)                                                      &
                            ! qsum
     &,qsump(n_dp)                                                      &
                            ! qsum positive increments
     &,qsumn(n_dp)                                                      &
                            ! qsum negative increments
     &,qloss(n_dp)                                                      &
                            ! q lossed from atmosphere
     &,cfactor(n_dp)        ! correction factor

      integer ::                                                        &
     & i,j,k                                                            &
                     ! loop counters
     &, jtt(n_dp)                                                       &
                     ! location of 0.949*p(1)
     &, kphase(n_dp) ! level of phase change from snow to rain

! parameters for Emanuel scheme

      real, parameter ::                                                &
     &  coeffr = 1.0                                                    &
                      ! coefficient govering rate of evaporation of rain
     &, coeffs = 0.8                                                    &
                      ! coefficient govering rate of evaporation of snow
     &, omtrain = 50.0                                                  &
                       ! the assumed fall speed of rain (kg/m/s3 = Pa/s)
     &, omtsnow =  5.5                                                  &
                       ! the assumed fall speed of snow (kg/m/s3 = Pa/s)
     &, sigd = 0.05                                                     &
                      ! fractional area covered by unsaturated DD
     &, sigs = 0.12   ! fraction of precipitation falling outside cloud
!
! Model constants
!
#include "c_epslon.h"
#include "c_r_cp.h"
#include "c_lheat.h"
#include "c_0_dg_c.h"
#include "c_g.h"

!
! External routines called:
!
!   NONE

!-----------------------------------------------------------------------
! 1.0 Initialisation of arrays
!-----------------------------------------------------------------------

      lvcp = lc/cp
      ginv = 1./g


      Do i= 1,n_dp
          qsum(i) = 0.0
          qsum2(i) = 0.0
          qsump(i) = 0.0
          qsumn(i) = 0.0
      EndDo

      Do k=1,nlev
        Do i= 1,n_dp

! initialise arrays

          md(i,k) = 0.0
          wt(i,k) = 0.0
          evap(i,k) = 0.0
          water(i,k) = 0.0

          dqbydt_dd(i,k)  = 0.0
          dqbydt_dd1(i,k)  = 0.0
          dqbydt_dd2(i,k)  = 0.0

          jtt(i)      = 2
          kphase(i)   = 0

! convert from theta to T

          t(i,k)  = th(i,k) *exner_layer_centres(i,k)

          If (T(I,k) >  273.15) then  ! inconsistent T freeze
              lv(i,k) = lc
          Else
              lv(i,k) = lc+lf
          End if
        EndDo
      EndDo

! Appears to be deriving some form of height using
!
!      p=rhoRT and dp/dz=-rhog    =>   gdz = -dp RT/p
!
!  equations here are using Tv  not T
!  Taking lowest model T level as zero height.

      k = 1
        Do i= 1,n_dp
          gz(i,k)  = 0.0
          h(i,k)  = cp*T(i,k)+gz(i,k)
          qp_down(i,k) = q(i,k)
        EndDo
      Do k=2,nlev
        Do i= 1,n_dp
          tvx=T(i,k)  *(1.+c_virtual*q(i,k))
          tvy=T(i,k-1)*(1.+c_virtual*q(i,k-1))

          gz(I,k)=gz(I,k-1)+0.5*R*(TVX+TVY)*(P(I,k-1)-P(I,k))/PH(I,k-1)

          h(i,k)  = cp*T(i,k)+gz(i,k)
          qp_down(i,k) = q(i,k-1)
        EndDo
      EndDo


! initialise tracer in down draught

      If (l_tracer) then
        Do j=1,ntra
          k=1
            Do i=1,n_dp
              trap_down(i,k,j) = tracer(i,k,j)
            End do
          Do k=2,nlev
            Do i=1,n_dp
              trap_down(i,k,j) = tracer(i,k-1,j)
            End do
          End do
        End do
      End if

!-----------------------------------------------------------------------
! 2.0 Main Down draught calculation
!-----------------------------------------------------------------------
! Assuming the array precip on model levels from G-R scheme corresponds
! to the Emanuel variable wdtrain - detrained precipitation.
!  Note In Emanuel scheme
!  wdtrain(i,k)=G*EP(I,k)*M(I,k)*CLW(I,k)
!  ep(i,k) - precip efficiency
!  M(i,k)  - up draught mass flux
!  clw(i,k) - condensed cloud water
!-----------------------------------------------------------------------

      Do k=kmax_term+1,1,-1          ! level loop working Downwards

        Do i=1,n_dp           ! grid point loop
! start at level above termination
          If (k <= kterm(i)+1.and.kterm(i) /= 0) then
!-----------------------------------------------------------------------
!  Does Down draughts if precipitating efficiency >= 0.0001 at
!  top if cloud (Emanuel) => precip>0.0
!-----------------------------------------------------------------------
            If (precip(i,kterm(i)) >  0.0) then  ! test for Downdraught
!
! detrained precipitation for level k
!

               wdtrain(i) = precip(I,k)*g

! Find rain water and evaporation using provisional estimates of qp
!
! Solution of equation (9) In Emanuel's paper to find lp(k) rain water
! at level k. Note as working Down from highest level know lp(k+1)
!

!
! Value of terminal velocity and coeffecient of evaporation for rain
!
              If (T(I,k) >  273.15) then
                COEFF=COEFFR
                wt(I,k)=OMTRAIN

              Else
!
! Value of terminal velocity and coeffecient of evaporation for snow
!
                COEFF=COEFFS
                wt(I,k)=OMTSNOW

              End if

              If (T(I,k) >  273.15.and.t(i,k+1) <= 273.15) then
                kphase(i) = k       ! melting level
              Endif

! qsm - temporary estimate for water vapour of parcel at level k

              qsm=0.5*(Q(I,k)+qp_down(I,k+1))

!
! Looks like expression for evaporation without sqrt (lp(k))
!
              afac=COEFF*PH(I,k)*0.01*(qse(I,k)-qsm)                    &
     &                    /(1.0E4+2.0E3*0.01*PH(I,k)*qse(I,k))
              afac=MAX(afac,0.0)
              B6=(PH(I,k)-PH(I,k+1))*sigs*afac/wt(I,k)

              C6=(water(I,k+1)*wt(I,k+1)+wdtrain(i)/SIGD)/wt(I,k)

              If (c6 == 0.0) then
                 revap =0.0     ! set because of numeric problems
              Else
                 Revap=0.5*(-B6+SQRT(B6*B6+4.*C6))
              Endif
              If (revap <  0.0) then
                revap =0.0        ! reset
              Endif

              evap(I,k)=sigs*afac*Revap

! water - rain water (lp in Emanuel paper)

              water(I,k)=Revap*Revap

!
! Calculate precipitating Downdraught mass flux under hydrostatic approx
!
            IF(K >  1) then
              dhdp=(H(I,k)-H(I,k-1))/(P(I,k-1)-P(I,k))
              dhdp=MAX(dhdp,0.1)      ! correct units
              md(I,k)=ginv*lv(i,k)*SIGD*evap(I,k)/dhdp
              md(I,k)=MAX(md(I,k),0.0)

!
! Add a small amount of inertia to Downdraught
! (not referred to in paper?)
!
              FAC=20.0*100./(PH(I,k-2)-PH(I,k-1))

              md(I,k)=(FAC*md(I,k+1)+md(I,k))/(1.+FAC)

!
! Force Downdraught mass flux to decrease linearly to zero between about
! 950mb and the surface.
! Actually coded as decreasing linearly to zero between
! 0.949*lowest level pressure and model T level (Is this what we really
! want it to Do ?) (Did the original Emanuel code assume level 1 was
! the surface ?)
!

             IF(P(I,k) >  (0.949*P(i,1)))then
               JTT(i)=MAX(JTT(i),k)
               md(i,k)=md(i,JTT(i))*(P(i,1)-P(I,k))                     &
     &                                        /(P(i,1)-P(i,JTT(i)))

             End if

            Endif
!
! Find mixing ratio of precipitation Downdraught
!

             IF(k == 1)then
               qstm=qse(i,k)
             ELSE
               qstm=qse(i,k-1)
             End if
             IF(md(I,k) >  md(I,k+1))then
               rat=md(I,k+1)/md(I,k)
               qp_down(I,k)=qp_down(I,k+1)*rat+Q(I,k)*(1.0-rat)+ginv*   &
     &                SIGD*(PH(I,k-1)-PH(I,k))*(evap(I,k)/md(I,k))

             Else
               IF(md(I,k+1) >  0.0)then
                 qp_down(I,k)=(gz(I,k+1)-gz(I,k)+qp_down(I,k+1)*lv(i,k) &
     &                             +CP*(T(I,k+1)-T(I,k)))/lv(i,k)
               End if
             End if
!  k+1 md(i,k)= 0.0 therefore qp(i,1) remains unchanged
             qp_down(I,k)=MIN(qp_down(I,k),qstm)
             qp_down(I,k)=MAX(qp_down(I,k),0.0)

           End if
          End if
        End Do   ! loop over points
      End Do     ! loop over levels

!
! calculate tracers in downdraught plume
!
      if (l_tracer) then
        Do J=1,NTRA
          Do k=kmax_term+1,1,-1   ! level loop working Downwards
            Do i=1,n_dp           ! grid point loop

              If (k <= kterm(i)+1.and.kterm(i) /= 0) then

                IF(md(I,k) >  md(I,k+1))then
                  rat=md(I,k+1)/md(I,k)

                  TRAP_down(I,k,J)=TRAP_down(I,k+1,J)*rat               &
     &                                          +Tracer(I,k,J)*(1.-rat)
!eman orig I think this is wrong
!    &                              +TRAP_down(I,k,J)*(1.-rat)
                Else
                  IF(md(I,k+1) >  0.0)then
                     TRAP_down(I,k,J)=TRAP_down(I,k+1,J)
                  End if
                End if

              End if

            End Do   ! loop over points
          End Do     ! loop over levels
        End Do     ! loop over tracers

      Endif      ! test on tracers

! water corresponds to lp in paper
! Emanuel scheme assumes all surface precip rain. Note wt depends on
! whether rain or snow should we use this information?

!  now in mm/s = kg/m2/s  (wt in Pa/s)

      Do i=1,n_dp
        if (wt(i,1) == omtrain) then      ! falling as rain
          rain(i) = rain(i) + wt(i,1)*sigd*water(i,1)/g
        Else     ! falling as snow
          snow(i) = snow(i) + wt(i,1)*sigd*water(i,1)/g
        End if
      EndDo


! ----------------------------------------------------------------------
! Calculation of increments
! ----------------------------------------------------------------------

!  level 1

      k=1
      Do i=1,n_dp
       if (kterm(i) /= 0) then
        dpinv(i) = 1.0/(PH(i,0)-PH(i,1))

!
! increments to temperature    dtheta = dt/exner
!
! evaporation of precipitation
!

        dthbydt(i,k) = dthbydt(i,k) - (lv(i,k)/cp)*sigd*evap(i,1)       &
     &                      /exner_layer_centres(i,1)
!
! increments to q
!
! change in q due to transport by Down draught
!
        dqbydt_dd2(i,k)=dqbydt_dd2(i,k)+G*md(i,2)*(qp_down(i,2)-Q(i,1))*&
     &                dpinv(i)

!
! increase in q due to evaporation of precip
!
        dqbydt_dd1(i,k) = dqbydt_dd1(i,k)+SIGD*evap(i,1)
!
! total
!
        dqbydt_dd(i,k) = dqbydt_dd1(i,k)+dqbydt_dd2(i,k)

       End if
      End Do      ! loop over gridpoints

!
!  level where phase changes add lf*water falling
!

      Do i=1,n_dp

        if (kphase(i) /= 0) then
!
! evaporation of precipitation
!
          dthbydt(i,kphase(i)) = dthbydt(i,kphase(i))                   &
     &                       - lf*water(i,kphase(i))/cp/timestep        &
     &                               /exner_layer_centres(i,kphase(i))
        endif
      End Do      ! loop over gridpoints




!
!  Levels above k=1
!
      Do k=2,kmax_term+1        ! level loop
        Do i=1,n_dp            ! grid point loop

          If (k <= kterm(i)+1.and.kterm(i) /= 0) then

            dpinv(i)=1.0/(PH(i,k-1)-PH(I,k))

!
! Temperature increments
!
!
! evaporation of precipitation
!

            dthbydt(i,k) = dthbydt(i,k) - (lv(i,k)/cp)*sigd*evap(i,k)   &
     &                      /exner_layer_centres(i,k)

!
! change in q due to transport by Down draught
!
            dqbydt_dd2(i,k) = dqbydt_dd2(i,k) +                         &
     &                   G*dpinv(i)*(md(i,k+1)*(qp_down(i,k+1)-Q(i,k))  &
     &                              -md(I,k)*(qp_down(I,k)-Q(I,k-1)))

!
! increase in q due to evaporation of precip
!
            dqbydt_dd1(i,k) = dqbydt_dd1(i,k) + SIGD*evap(i,k)
            dqbydt_dd(i,k) = dqbydt_dd1(i,k) +  dqbydt_dd2(i,k)


          Endif           ! test on whether level in convecting column
        End Do             ! gridpoint loop
      End Do               ! level loop

      If (l_tracer) then
        Do j = 1, ntra
          k = 1
            Do i=1,n_dp
              If (kterm(i) /= 0) then
               dpinv(i) = 1.0/(PH(i,1)-PH(i,2))
               dtrabydt(i,k,j) = dtrabydt(i,k,j)                        &
     &       +  g*dpinv(i)*md(i,k+1)*(trap_down(i,k+1,j)-tracer(i,k,j))
              End if
            End Do
        End Do

        Do j = 1, ntra
          Do k= 2, kmax_term+1
            Do i=1,n_dp
              If (k >= kterm(i).and.kterm(i) /= 0) then

                dpinv(i)=1.0/(PH(i,k)-PH(I,k+1))
                dtrabydt(i,k,j) = dtrabydt(i,k,j)                       &
     &       +G*dpinv(i)*(md(I,k+1)*(TRAP_down(I,k+1,J)-tracer(I,k,J))- &
     &                    md(I,k)*(TRAP_down(I,k,J)-TRAcer(I,k-1,J)))
! eman orig I think this is wrong
!    &                    md(I,k)*(TRAP_down(I,k,J)-TRAP_down(I,k-1,J)))

              End if       ! test on whether level in convecting column
            End Do        ! gridpoint loop
          End Do           ! level loop
        End Do             ! tracer loop

      End if
!
! Conservation checks
!
      Do k= 1, nlev
        Do i=1,n_dp
           DP=(PH(i,k-1)-PH(I,k))       ! correct for layer k
           qsum(i)=qsum(i) + dqbydt_dd(i,k)*dp/g
           qsum2(i)=qsum2(i) + precip(i,k)
           if (dqbydt_dd(i,k) >  0.0) then
             qsump(i) = qsump(i) + dqbydt_dd(i,k)*dp/g
           else
             qsumn(i) = qsumn(i) + dqbydt_dd(i,k)*dp/g
           endif

! change units of Emanuel Downdraught mass flux for diagnostic output

           Down_flux(i,k) = g*md(i,k)

        EndDo           ! gridpoint loop
      EndDo             ! level
!
! Need to correct dq increments as moisture not conserved and non
! conservation is systematic leading to loss of q from atmosphere.
!
      Do i = 1, n_dp
        qloss(i) = qsum2(i)-rain(i)-qsum(i)-snow(i)
        if (qloss(i)  == 0.0) then
          cfactor(i) = 1.0
        else
          If (qsump(i) /= 0.0) then
            cfactor(i) = 1.0 + qloss(i)/qsump(i)
          Else
            If (qsumn(i) == 0.0) then
! No increments to q but no moisture balance ie qsum2 not = rain/snow
! Cannot correct increments therefore correct rain
              cfactor(i) = 1.0
              if (rain(i) >  0.0) then
                rain(i) = qsum2(i)
              else  ! assume precip snow
                snow(i) = qsum2(i)
              endif
            Else
              cfactor(i) = 1.0 - qloss(i)/qsumn(i)
            Endif
          endif
        endif
      End do
!
! Scale positive increments to increase them
!
      Do k= 1, nlev
        Do i=1,n_dp

          If (qsump(i) /= 0.0) then
           if (dqbydt_dd(i,k) >  0.0) then
             dqbydt_dd(i,k) = dqbydt_dd(i,k)*cfactor(i)
           endif
          Else       ! reduce negative increments
           if (dqbydt_dd(i,k) <  0.0) then
             dqbydt_dd(i,k) = dqbydt_dd(i,k)*cfactor(i)
           endif
          Endif

          dqbydt(i,k) = dqbydt(i,k) + dqbydt_dd(i,k)

        EndDo           ! gridpoint loop
      EndDo             ! level

! ----------------------------------------------------------------------
      Return
      END SUBROUTINE EMAN_DD

#endif
