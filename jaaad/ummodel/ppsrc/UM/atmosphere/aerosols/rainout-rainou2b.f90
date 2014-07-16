
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
      SUBROUTINE RAINOUT(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           rho_r2, Q,                                             &
     &           QCF_REMAIN, QCL_REMAIN,                                &
     &           QCF_PREVIOUS, QCL_PREVIOUS,                            &
     &           LS_RAIN3D, LS_SNOW3D,                                  &
     &           timestep,                                              &
     &           AERO_INCLOUD,                                          &
     &           AERO_ACCUM,                                            &
     &           RNOUT_AERO)
!
!---------------------------------------------------------------------
! Purpose: Makes the rain-out of dissolved aerosols during large-scale
!          precipitations. The rained-out mixing ratio depends on the
!          conversion rate of condensed water to precipitated water.
!          Re-evaporation is taken into account by transfering some
!          of the dissolved aerosol mass to the accumulation mode mass.
!          This is computed as being proportional to the amount of
!          precipitated water that re-evaporates. The routine also
!          accounts for those droplets that shrink without re-evap
!          completely by reevaporating only half of the bulk sulphate
!          concentration (unless rain is completely re-evaporated).
!
!          Called by microphys_ctl
!
! Current code owners: N Bellouin
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   6.2   26/10/05   New Deck (based on vn6.1 RAINOU2A)     N Bellouin
!                                                           & A Jones
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  UMDP 20
!
!---------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, off_x                                                           &
                                        !EW size of std. halo
     &, off_y                                                           &
                                        !NS size of std. halo
     &, halo_i                                                          &
                                        !EW extended halo
     &, halo_j                                                          &
                                        !NS extended halo
     &, model_levels                                                    &
     &, wet_model_levels
!
      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,     model_levels)            &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j, 0:model_levels)            &
     &, rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                      model_levels)  !density*r*r
!
      REAL                                                              &
     &  QCL_REMAIN(row_length,rows,wet_model_levels)                    &
     &, QCF_REMAIN(row_length,rows,wet_model_levels)                    &
     &, QCF_PREVIOUS(row_length,rows,wet_model_levels)                  &
     &, QCL_PREVIOUS(row_length,rows,wet_model_levels)                  &
     &, Q(row_length,rows,wet_model_levels)                             &
! the LS_xxxx3D arrays are the large-scale precipitation fluxes.
! The flux is defined at the bottom of the layer. It is the flux that
! falls OUT of the layer on which it is defined. Equivalently, it is the
! flux that falls INTO the layer below the one in which it is defined.
     &, LS_RAIN3D(row_length,rows,wet_model_levels)                     &
     &, LS_SNOW3D(row_length,rows,wet_model_levels)
!
      REAL timestep ! timestep in seconds
!
! Arguments with intent IN/OUT:
!   mass mixing ratio of in-cloud and accumulation/aged aerosol
      REAL                                                              &
     & AERO_INCLOUD(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &                                         model_levels),           &
     & AERO_ACCUM(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &                                         model_levels)

!
! Arguments with intent OUT (diagnostics):
      REAL                                                              &
     & RNOUT_AERO(row_length,rows)      !tracer removed kg/m2/ts
!
!  Local variables
!
      INTEGER i, j, k              !loop variables
!
      REAL FRACTION                ! fraction of droplets that
      PARAMETER(FRACTION = 0.5)    ! shrinks without re-evaporating

      REAL FRAC_AERO_SCAV            ! Dissolved aerosol is
      PARAMETER(FRAC_AERO_SCAV = 1.0)! completely dissolved

      REAL                                                              &
     &   DELTA_AERO                                                     &
                            ! amount of aerosol removed from grid box
     &,  DM_DRY                                                         &
                            ! mass p.u.area of dry air in layer
     &,  RHODR                                                          &
                            ! density multiplied by layer thick.
     &,  TOT_PRECIP3D                                                   &
                            ! total ppn in grid box
     &,  rho1                                                           &
                            ! densities at
     &,  rho2                                                           &
                            !  adjacent model levels
     &,  BETA                                                           &
                            ! see comments below
     &,  SCAV_EXP                                                       &
                            ! Exponential term in scavanged MMR
     &,  SMALLP             ! small +ve number, negligible compared to 1
!
!
! Initialisation
! ... epsilon() is defined as almost negligible, so eps/100 is negligible

      SMALLP = epsilon(1.0) / 100.0

      DO j = 1, rows
        DO i = 1, row_length
          RNOUT_AERO(i, j) = 0.0
        ENDDO
      ENDDO

      DO j = 1, rows
        DO i = 1, row_length
!
! Loop on vertical levels, from the top of the atmosphere
! to the surface. We assume no scavenging at the topmost
! level. Note that wet_model_levels <= model_levels (= in
! most cases)
!
          DO k = wet_model_levels - 1, 1, -1
!
! compute the mass of air per unit area for conversion
! of mass mixing ratio to mass
!
! first, get the density from the array rho * r^2
            rho1 = rho_r2(i,j,k)/ ( r_rho_levels(i,j,k) *               &
     &                              r_rho_levels(i,j,k) )
            rho2 = rho_r2(i,j,k+1)  / ( r_rho_levels(i,j,k+1)   *       &
     &                                  r_rho_levels(i,j,k+1)   )
!
! RHODR is the density (interpolated onto theta levels) multiplied
! by delta_r
            RHODR = rho2 * ( r_theta_levels(i,j,k) -                    &
     &                       r_rho_levels(i,j,k)   )                    &
     &              +                                                   &
     &              rho1 * ( r_rho_levels(i,j,k+1) -                    &
     &                       r_theta_levels(i,j,k) )
!
! special case for the lowest layer
            IF(k  ==  1) THEN
              RHODR = RHODR *                                           &
     &          (r_rho_levels(i,j,k+1)-r_theta_levels(i,j,k-1))         &
     &         /(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
            ENDIF
!
! conversion to dry air density
            DM_DRY = RHODR * (1.0 - Q(i,j,k) -                          &
     &        qcl_previous(i,j,k) - qcf_previous(i,j,k))

!
! compute the conversion rate of cloud water into precipitated water
            BETA = LS_RAIN3D(i,j,k)   + LS_SNOW3D(i,j,k)   -            &
     &             LS_RAIN3D(i,j,k+1) - LS_SNOW3D(i,j,k+1)

            IF(QCL_PREVIOUS(i,j,k)+QCF_PREVIOUS(i,j,k)  >   0.0)THEN
             BETA = BETA / RHODR / (QCL_PREVIOUS(i,j,k)+                &
     &                             QCF_PREVIOUS(i,j,k))
             BETA = MAX( 0.0, BETA)
!
! get the scavenged mass mixing ratio
! here, delta_aero is <= 0
! This exponential can underflow if FRAC_AERO_SCAV is too large
! Only calculate if result is significant relative to 1.0
             IF ( FRAC_AERO_SCAV * BETA * TIMESTEP < - LOG( SMALLP ) ) THEN
                SCAV_EXP = EXP( -FRAC_AERO_SCAV * BETA * TIMESTEP)
             ELSE
                SCAV_EXP = 0.0
             END IF
             DELTA_AERO = AERO_INCLOUD(i,j,k) * ( SCAV_EXP - 1.0 )

             AERO_INCLOUD(i,j,k) = AERO_INCLOUD(i,j,k) + DELTA_AERO
!
! the cumulated rain-out for the current level
! (by convention, it is positive)
             RNOUT_AERO(i,j) = RNOUT_AERO(i,j) - DELTA_AERO * DM_DRY
            ENDIF
!
! re-evaporation now.
! compute the fraction of precipitated water that re-evaporates

            BETA = LS_RAIN3D(i,j,k)   + LS_SNOW3D(i,j,k) -              &
     &             LS_RAIN3D(i,j,k+1) - LS_SNOW3D(i,j,k+1)
            IF(BETA  <   0.0) THEN
              BETA = BETA / ( LS_RAIN3D(i,j,k+1) + LS_SNOW3D(i,j,k+1) )
            ENDIF

            ! in the following, if beta > 0, it is set to 0

            IF( LS_RAIN3D(i,j,k)   + LS_SNOW3D(i,j,k)  ==  0.0) THEN
              ! total re-evaporation
              BETA = MIN( MAX(0.0, -BETA), 1.0)
            ELSE
              ! non total re-evaporation
              ! fraction is here to take into account those
              ! raindrops that shrinks without evaporating totally
              ! (the value of fraction is somehow arbitrary)
              BETA = MIN( MAX(0.0, -BETA)*FRACTION, 1.0)
            ENDIF

            ! fraction of mass mixing ratio that re-evaporates
            ! into the accumulation/aged mode
            ! DELTA_AERO is >= 0
            DELTA_AERO = BETA * RNOUT_AERO(i,j) / DM_DRY

            AERO_ACCUM(i,j,k) = AERO_ACCUM(i,j,k) + DELTA_AERO

            RNOUT_AERO(i,j) = (1.0-BETA) * RNOUT_AERO(i,j)

          ENDDO ! k

        ENDDO ! i
      ENDDO ! j

      RETURN

      END SUBROUTINE RAINOUT

