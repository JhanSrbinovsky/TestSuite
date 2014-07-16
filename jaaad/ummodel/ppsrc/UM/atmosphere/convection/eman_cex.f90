
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Compress expand around Emanuel Downdraught code
!
      SUBROUTINE EMAN_CEX(npnts,nmid,kmax_term,nlev,trlev,ntra          &
     &,                      kterm,l_mid,l_tracer                       &
     &,                      exner_layer_centres,exner_layer_boundaries &
     &,                      p,ph,timestep,th,q,qse,tracer,precip       &
     &,                      dthbydt, dqbydt, dtrabydt                  &
     &,                      rain, snow ,Dwn_flux                       &
     &                   )

!
! Purpose:
!   First compress to points with mid-level convection.
!   Then call  Emanuel Down draughts.
!   Expand resulst back to full grid.

!
!   Called by MID_CONV.
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
     &  npnts                                                           &
                 ! No. of deep convection points
     &, nmid                                                            &
                 ! No. of mid level convection points
     &, nlev                                                            &
                 ! No. of model layers
     &, trlev                                                           &
                 ! No. of tracer levels
     &, ntra                                                            &
                 ! No. of tracer fields
     &, kmax_term                                                       &
                     ! highest level reached by convection
     &, kterm(npnts)    ! level reached by convection


      logical, intent(in) ::                                            &
     & l_tracer                                                         &
                             ! true in tracers present
     &, l_mid(npnts)         ! true  if mid level convection

      real, intent(in)    ::                                            &
     &  exner_layer_centres(npnts,0:nlev)                               &
                                             ! Exner
     &, exner_layer_boundaries(npnts,0:nlev)                            &
                                             ! Exner at half level above
     &, p(npnts,0:nlev)                                                 &
                            ! Pressure  (Pa)
     &, ph(npnts,0:nlev)                                                &
                            ! Pressure at half level (Pa)
     &, timestep                                                        &
                           ! timestep for model physics in seconds

     &, qse(npnts,nlev)                                                 &
                        ! Saturation specific humidity of cloud
                       ! environment (kg/kg)

     &, q (npnts,nlev)                                                  &
                        ! specific humidity of cloud environment(kg/kg)

     &, th (npnts,nlev)                                                 &
                            ! theta of cloud environment(K)

     &, precip(npnts,nlev)                                              &
                            ! Precip from updraught G-R scheme(kg/m2/s)
                           ! not units for wdtrain
     &, tracer(npnts,trlev,ntra)  !  Tracer on model levels  (kg/kg)

!
! Arguments with intent INOUT:
!

      real, intent(inout)    ::                                         &
     &  rain(npnts)                                                     &
                                ! rainfall at surface (kg/m**2/s)
     &, snow(npnts)             ! snowfall at surface (kg/m**2/s)

! increments
      real, intent(inout)    ::                                         &
     &  dqbydt(npnts,nlev)                                              &
                                  ! increments to q (kg/kg/s)
     &, dthbydt(npnts,nlev)                                             &
                              ! increments to potential temperature(K/s)
     &, dtrabydt(npnts,nlev,ntra) ! increments to model tracers(kg/kg/s)

!
! Arguments with intent OUT:
!
      real, intent(out)    ::                                           &
     &  Dwn_flux(npnts,nlev) ! Downdraught mass flux (Pa/s) diagnostic

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
! compressed version of most input arrays

      integer ::                                                        &
     & kterm_c(nmid)      ! compressed termination level

      real ::                                                           &
     &  rain_c(nmid)                                                    &
                           ! compressed rain
     &, snow_c(nmid)                                                    &
     &, Dwn_flux_c(nmid,nlev)                                           &
     &, dqbydt_c(nmid,nlev)                                             &
     &, dthbydt_c(nmid,nlev)                                            &
     &, dtrabydt_c(nmid,nlev,ntra)                                      &
     &, exner_layer_centres_c(nmid,0:nlev)                              &
     &, exner_layer_boundaries_c(nmid,0:nlev)                           &
     &, p_c(nmid,0:nlev)                                                &
     &, ph_c(nmid,0:nlev)                                               &
     &, qse_c(nmid,nlev)                                                &
     &, q_c(nmid,nlev)                                                  &
     &, th_c(nmid,nlev)                                                 &
     &, precip_c(nmid,nlev)                                             &
     &, tracer_c(nmid,trlev,ntra)

      Integer ::                                                        &
     & dmi(nmid)      ! index of convecting points

      Integer ::                                                        &
     & i ,j ,k      ! loop counters
!
! External routines called:
!
!   NONE

! ----------------------------------------------------------------------
! Compress input arrays
! ----------------------------------------------------------------------
      j=0
      Do i = 1,npnts
        If (l_mid(i)) Then
          j = j+1
          dmi(j) = i
        End if
      End do

      Do j = 1, nmid
        kterm_c(j) = kterm(dmi(j))
        rain_c(j) = rain(dmi(j))
        snow_c(j) = snow(dmi(j))

      End do

      Do k = 0, nlev
        Do j = 1, nmid
          exner_layer_centres_c(j,k) = exner_layer_centres(dmi(j),k)
          exner_layer_boundaries_c(j,k)=exner_layer_boundaries(dmi(j),k)
          p_c(j,k) = p(dmi(j),k)
          ph_c(j,k) = ph(dmi(j),k)
        End do
      End do

      Do k = 1, nlev
        Do j = 1, nmid
          q_c(j,k) = q(dmi(j),k)
          th_c(j,k) = th(dmi(j),k)
          qse_c(j,k) = qse(dmi(j),k)

          precip_c(j,k) = precip(dmi(j),k)

          dqbydt_c(j,k) = dqbydt(dmi(j),k)
          dthbydt_c(j,k) = dthbydt(dmi(j),k)

        End do
      End do
      If (l_tracer) then
        Do i = 1, ntra
          Do k = 1, trlev
            Do j = 1, nmid
              tracer_c(j,k,i) = tracer(dmi(j),k,i)
            End do
          End do
          Do k = 1, nlev
            Do j = 1, nmid
              dtrabydt_c(j,k,i) = dtrabydt(dmi(j),k,i)
            End do
          End do
        End do
      End if

! ----------------------------------------------------------------------
! Call Emanuel scheme for just mid level convecting points
! ----------------------------------------------------------------------
! DEPENDS ON: eman_dd
        call eman_dd (nmid,kmax_term,nlev,trlev,ntra,kterm_c,l_tracer   &
      ,               exner_layer_centres_c,exner_layer_boundaries_c    &
      ,               p_c,ph_c,timestep,th_c,q_c,qse_c,tracer_c         &
      ,               precip_c,dthbydt_c,dqbydt_c,dtrabydt_c            &
      ,               rain_c, snow_c ,dwn_flux_c)

! ----------------------------------------------------------------------
! Expand results
! ----------------------------------------------------------------------
      Do i = 1,nmid
        rain(dmi(i)) = rain_c(i)
        snow(dmi(i)) = snow_c(i)
      End do
      Do k = 1, nlev
        Do i = 1, nmid
          dwn_flux(dmi(i),k) = dwn_flux_c(i,k)
          dqbydt(dmi(i),k)   = dqbydt_c(i,k)
          dthbydt(dmi(i),k)  = dthbydt_c(i,k)
        End do
      End do

      If (l_tracer) then
        Do j = 1, ntra
          Do k = 1, nlev
            Do i = 1, nmid
              dtrabydt(dmi(i),k,j) = dtrabydt_c(i,k,j)
            End do
          End do
        End do
      End if
! ----------------------------------------------------------------------
      Return
      END SUBROUTINE EMAN_CEX

