#if defined(A05_4A) || defined(A05_5A) 
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
      SUBROUTINE SCNSCV2(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           Q, QCL, QCF, TRACER,                                   &
     &           CCLDBASE, CCLDTOP,                                     &
     &           RAINRATE, SNOWRATE,                                    &
     &           L_SCAV_BELOW_CLOUD,                                    &
     &           K_RAIN, K_SNOW,                                        &
     &           ACCU_SCAV_TR                                           &
     &           )
!
! ---------------------------------------------------------------------
!  Purpose: Scavenge Sulphur Cycle tracers by convective precipitation
!
!           Called by Conv_Ctl if Sulphur Cycle is on
!
!  Current code owner:  M. Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   5.3     12/10/01  New Deck (based on SCONSC1A at vn4.5)   M Woodage
!   6.2     03/02/05  Added section 5A. R A Stratton
! 6.2      15/12/05  Correct dust diagnostics when substepping phys2.
!                                                        M. Diamantakis
!   6.4    12/12/06  Removed def 3C. R A Stratton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!  System component covered:
!
!  System task:
!
! Documentation:  UMDP 20
!
!----------------------------------------------------------------------
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
     &, wet_model_levels                                                &
     &, substep_number                                                  &
     &, CCLDBASE(row_length,rows)                                       &
                                       !convective cloud base
     &, CCLDTOP(row_length,rows)       !convective cloud top
!
      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,     model_levels)            &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j, 0:model_levels)            &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                                           !density*r*r
     &                                         model_levels)            &
     &, Q(row_length,rows,wet_model_levels)                             &
                                                   !water vapour(kg/kg)
     &, Qcl(row_length,rows,wet_model_levels)                           &
     &, Qcf(row_length,rows,wet_model_levels)                           &
     &, RAINRATE(row_length,rows)                                       &
                                       !conv rain rate at surface kg/m2
     &, SNOWRATE(row_length,rows)      !conv snow rate at surface kg/m2

!
      REAL                                                              &
     & timestep                                                         &
                                   !timestep in secs
     &, K_RAIN                                                          &
                                   !scavenging rate coeff for rain
     &, K_SNOW                     !scavenging rate coeff for snow
!
      LOGICAL                                                           &
     & L_SCAV_BELOW_CLOUD          !control for scavenging levels
!
! Arguments with intent IN/OUT:
      REAL                                                              &
     & TRACER(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &                                         model_levels)
!
! Arguments with intent OUT (diagnostics):
      REAL                                                              &
     & ACCU_SCAV_TR(row_length, rows)     !column total of scvnged trcr
!
! Local variables
!
      INTEGER                                                           &
     & START_LEVEL                                                      &
                                     ! lowest level for scavenging
     &, i, j, k                      ! Loop variables
!
      REAL TERMR,                                                       &
                                     ! to assist calcn of scav rate
     &     TERMS,                                                       &
                                     !
     &     DM,                                                          &
                                     ! mass p.u.area of air in layer
     &     DELTA_TR,                                                    &
                                     ! tracer increment due to scvnging
     &     rho1, rho2                ! air densities
!
      REAL                                                              &
     &     TOTRATE                   ! total scav rate
!
!
! Initialise ACCU_SCAV_TR array to zero before adding accumulations
          IF ( Substep_Number == 1 ) THEN
            DO J=1,ROWS
               DO I=1,ROW_LENGTH
                  ACCU_SCAV_TR(I,J)=0.0
               END DO
            END DO
          END IF
!
! Calculate total scavenging rate
!
      Do j=1,rows
        Do i=1,row_length
!
          IF (CCLDTOP(i,j) >  0) THEN
!
! Set up START_LEVEL for scavenging
            IF (L_SCAV_BELOW_CLOUD) THEN
              START_LEVEL = 1
            ELSE
              START_LEVEL = CCLDBASE(i,j)
            END IF
!
            IF (RAINRATE(i,j) <= 0.0) THEN    !check for negative ppn
              TERMR=0.0
            ELSE
              TERMR=K_RAIN*RAINRATE(i,j)
            ENDIF
!
            IF (SNOWRATE(i,j) <= 0.0) THEN
              TERMS=0.0
            ELSE
              TERMS=K_SNOW*SNOWRATE(i,j)
            ENDIF
!
! Calculate TOTRATE, *3600.0 because K_RAIN and K_SNOW are derived for
!  ppn rates in mm/hr, but model values are kg/m2/s (cf CON_SCAV)
!
            TOTRATE=(TERMR+TERMS)*3600.0*timestep
!
! Increase TOTRATE to obtain rate in cloudy part of grid box
! Assume CCA=0.05
!
            TOTRATE=TOTRATE / 0.05

! Calculate amount of tracer scavenged and add to column total
!
            DO k = START_LEVEL, CCLDTOP(i,j)
!
! Calculate proportion of tracer mixing ratio scavenged out
              DELTA_TR=TRACER(i,j,k)*(1.0-EXP(-TOTRATE))
!
! Reduce DELTA_TR to allow for non_cloudy part of grid box
              DELTA_TR = DELTA_TR * 0.05
!
! Calculate mass of air per unit area in layer for conversion of tracer
!  mixing ratio increment to mass p.u.a. for STASH.
! Avoid calculations at top model level and assume no scavenging there.
!
                IF (k  <   model_levels) THEN
! Remove the r squared factor from rho before interpolation
                  rho1=rho(i,j,k)/( r_rho_levels(i,j,k) *               &
     &                              r_rho_levels(i,j,k) )
                  rho2=rho(i,j,k+1)/( r_rho_levels(i,j,k+1) *           &
     &                                r_rho_levels(i,j,k+1) )
!
! DM = density (interpolated on to theta levels) * delta r
                  DM= rho2 * ( r_theta_levels(i,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) +                  &
     &                rho1 * ( r_rho_levels(i,j,k+1) -                  &
     &                         r_theta_levels(i,j,k) )
!
! Special case for lowest layer to get correct mass
                    If (k  ==  1) Then
                      DM= DM *                                          &
     &                  (r_rho_levels(i,j,2) - r_theta_levels(i,j,0))   &
     &                / (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
                    End If
!
! Convert DM to DRY density if level is wet
                    If (k  <=  wet_model_levels) Then
                      DM = DM * (1.0 - Q(i,j,k) - qcl(i,j,k) -          &
     &                                            qcf(i,j,k))
                    End If
!
                ELSE
                  DELTA_TR = 0.0
                  DM = 0.0
                END IF
!
! Increment column total mass p.u.a. of scavenged tracer
!
              ACCU_SCAV_TR(i,j)=ACCU_SCAV_TR(i,j)+DELTA_TR*DM
!
! Decrement tracer mixing ratio
!
              TRACER(i,j,k)=TRACER(i,j,k)-DELTA_TR
!
            END DO                     !end k loop
!
          END IF                       !END CCLDTOP > 0 TEST
!
        End Do                         !End i loop
      End Do                           !End j loop
!
      RETURN
      END SUBROUTINE SCNSCV2
!
#endif
