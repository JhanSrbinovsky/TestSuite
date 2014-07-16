
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
      SUBROUTINE SCNWSH2(                                               &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep, substep_number,                              &
     &           rho,                                                   &
     &           Q, qcl, qcf, S_MMR,                                    &
     &           CCLDBASE, CCLDTOP,                                     &
     &           RAINRATE, ACCU_SCAV_TR                                 &
     &           )
!
! ---------------------------------------------------------------------
!  Purpose: Scavenge SO2 by convective precipitation
!
!           Called by Conv_Ctl if Sulphur Cycle is on.
!
!  Current code owner:  M. Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   5.3     12/10/01  New Deck                               M Woodage
!   6.2     03/02/05  Added section 5A. R A Stratton
! 6.2      15/12/05  Correct scavenging when substepping phys2.
!   6.4     12/12/06 Removed def 3C. R A Stratton
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
     &, CCLDTOP(row_length,rows)        !convective cloud top
!
      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,     model_levels)            &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j, 0:model_levels)            &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                                           !density*r*r
     &                                        model_levels)             &
     &, Q(row_length,rows,wet_model_levels)                             &
                                                   !water vapour(kg/kg)
     &, Qcl(row_length,rows,wet_model_levels)                           &
     &, Qcf(row_length,rows,wet_model_levels)                           &
     &, RAINRATE(row_length,rows)      !conv rain rate at surface kg/m2

!
      REAL                                                              &
     & timestep                    !timestep in secs
!
! Arguments with intent IN/OUT:
      REAL                                                              &
     & S_MMR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &                               model_levels)        !mmr S in SO2
!
! Arguments with intent OUT (diagnostics):
      REAL                                                              &
     & ACCU_SCAV_TR(row_length, rows)     !column total of scvnged trcr
!
!
!  Local variables
!
      INTEGER i, j, k                !loop counters
      INTEGER START_LEVEL            !lowest level for scavenging
!
      REAL TERMR,                                                       &
                                     !scavenging rate s-1
     &     DM,                                                          &
                                     !mass p.u.area of air in layer
     &     DELTA_TR,                                                    &
                                      !tracer increment due to scvn
     &     rho1, rho2                !air densities
!
      REAL S_PPBV                    !S in SO2 in ppbv
!
      REAL RAIN_MMPH                 !RAIN rate in mm/hour
      REAL MMR_TO_PPBV        !Factor to convert S mmr to ppbv
!                                 (mol wt air / mol wt S) x 10**9
      REAL S_THOLD            !Threshold value determining form of
!                               scavenging function
      REAL L0                 !Scavenging parameter for S <= S_THOLD
      REAL L1                 !Scavenging parameter for S >  S_THOLD
!
      PARAMETER (MMR_TO_PPBV = 0.9033E9,                                &
                                            ! 28.966 E9/32.066
     &               S_THOLD = 0.3065,                                  &
                                            ! ppbv
     &                    L0 = 6.5E-5,                                  &
                                            ! s-1
     &                    L1 = 2.955E-5                                 &
                                            ! s-1
     &          )
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
! Set level for scavenging to begin
      START_LEVEL=1
!
! Calculate scavenging rate if convective cloud present
!
      Do j=1,rows
        Do i=1,row_length
!
          IF (CCLDTOP(i,j)  >   0) THEN
!
! Calculate scavenging rate using function dependent on amount of S
! present:  d(ln(S))/dt = - L1 * (R/S_THOLD)**2/3  if S  <=  S_THOLD
!           d(ln(S))/dt = - L1 * (R/S      )**2/3  if S  >   S_THOLD
!
! Calculation of scavenging rate requires S in PPBV
! and rainrate in mm/hr
!
! Convert RAINRATE to RAIN_MMPH and increase rate to obtain rate in
! cloudy part of grid box assuming CCA=0.05
!
            RAIN_MMPH = RAINRATE(i,j) * 3600.0 / 0.05
!
            DO k=START_LEVEL, CCLDTOP(i,j)
!
! Convert S_MMR to S_PPBV
              S_PPBV = S_MMR(i,j,k) * MMR_TO_PPBV
!
                IF (RAIN_MMPH <= 0.0) THEN    !check for negative ppn
                    TERMR=0.0
                ELSE
!
                  IF (S_PPBV <= S_THOLD) THEN        !Use fn for low S
                    TERMR = L0 * ( RAIN_MMPH )**0.666667
                  ELSE                               !Use fn for high S
                    TERMR = L1 * ( RAIN_MMPH/S_PPBV )**0.666667
                  ENDIF
!
                END IF
!
! Calculate amount of tracer mixing ratio scavenged out per tstep
!
              DELTA_TR = S_MMR(i,j,k)*(1.0-EXP(-TERMR*timestep))
!
! Reduce DELTA_TR to allow for non_cloudy part of grid box
! assuming CCA=0.05
!
              DELTA_TR = DELTA_TR * 0.05
!
! Calculate mass of air per unit area in layer for conversion of tracer
! mixing ratio increment to mass p.u.a. for STASH
! Avoid calculations at top model level and assume no scavenging there.
!
                IF (k  <   model_levels) THEN
! Remove the r squared factor from rho before interpolation
                  rho1= rho(i,j,k)/( r_rho_levels(i,j,k) *              &
     &                               r_rho_levels(i,j,k) )
                  rho2= rho(i,j,k+1)/( r_rho_levels(i,j,k+1) *          &
     &                                 r_rho_levels(i,j,k+1) )
!
! DM = density (interpolated on to theta levels) * delta r
                  DM = rho2 * ( r_theta_levels(i,j,k) -                 &
     &                          r_rho_levels(i,j,k) ) +                 &
     &                 rho1 * ( r_rho_levels(i,j,k+1) -                 &
     &                          r_theta_levels(i,j,k) )
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
              S_MMR(i,j,k)=S_MMR(i,j,k)-DELTA_TR
!
            END DO                    !END k loop
!
          END IF                      !END conv cloud present condn
!
        End Do                        !End i loop
      End Do                          !End j loop
!
      RETURN
      END SUBROUTINE SCNWSH2
!
