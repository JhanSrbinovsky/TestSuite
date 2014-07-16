
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
      SUBROUTINE SL3DWASH(                                              &
     &           row_length, rows,                                      &
     &           off_x, off_y, halo_i, halo_j,                          &
     &           model_levels, wet_model_levels,                        &
     &           r_rho_levels, r_theta_levels,                          &
     &           timestep,                                              &
     &           rho_r2,                                                &
     &           Q, QCL, QCF, S_MMR,                                    &
     &           LS_RAIN3D,                                             &
     &           LSCAV_TR                                               &
     &           )
!
!----------------------------------------------------------------------
! Purpose: Scavenging of SO2 by large-scale ppn.
!
!          Called by microphys_ctl if Sulphur Cycle is on.
!
! Current code owner: M. Woodage
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   06/10/01   New Deck                                M Woodage
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
     &, wet_model_levels
!
      REAL                                                              &
     &  r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j,     model_levels)            &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j, 0:model_levels)            &
     &, rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &                                       model_levels)              &
                                                           !density*r*r
!
     &, LS_RAIN3D(row_length,rows,wet_model_levels)                     &
                                                    !rain rate (kg/m2/s)
     &, Q(row_length,rows,wet_model_levels)                             &
                                                    !water vapour(kg/kg)
     &, Qcl(row_length,rows,wet_model_levels)                           &
     &, Qcf(row_length,rows,wet_model_levels)
!
      REAL                                                              &
     & timestep                                     !timestep in secs
!
! Arguments with intent IN/OUT:
      REAL                                                              &
     & S_MMR(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &                                   model_levels)    !mmr S in SO2
!
! Arguments with intent OUT (diagnostics):
      REAL LSCAV_TR(row_length,rows)          !accumulated scavenged tr
!
! Local variables
      INTEGER i, j, k       !LOOP COUNTERS
!
      REAL                                                              &
     &   TERMR                                                          &
                            !scavenging rate s-1
     &,  DELTA_TR                                                       &
                            !tracer increment due to scavnging
     &,  DM                                                             &
                            !mass p.u.area of air in layer
     &,  S_PPBV                                                         &
                            !S in SO2 in ppbv
     &,  RAIN_MMPH                                                      &
                            !RAIN rate in mm/hour
     &,  MMR_TO_PPBV                                                    &
                            !Factor to convert S mmr to ppbv
!                             (mol wt air / mol wt S) x 10**9
     &,  S_THOLD                                                        &
                            !Threshold value determining form of
!                           !   scavenging function
     &,  L0                                                             &
                            !Scavenging parameter for S <= S_THOLD
     &,  L1                                                             &
                            !Scavenging parameter for S >  S_THOLD
     &,  rho1                                                           &
                            !Densities at
     &,  rho2               !   adjacent model levels
!
      LOGICAL                                                           &
     &   L_DO_SCAV          ! for controlling scavenging
!
      PARAMETER (MMR_TO_PPBV = 0.9033E9,                                &
                                            ! 28.966E9/32.066
     &               S_THOLD = 0.3065,                                  &
                                            ! ppbv
     &                    L0 = 6.5E-5,                                  &
                                            ! s-1
     &                    L1 = 2.955E-5                                 &
                                            ! s-1
     &          )
!
!
!   Initialise LSCAV_TR to zero before doing rainout
      Do j=1,rows
        Do i=1,row_length
          LSCAV_TR(i,j) = 0.0
        End Do
      End Do
!
!
      Do j=1,rows
        Do i=1,row_length
!
          L_DO_SCAV=.TRUE.
!
! Loop over wet levels from surface up, and only scavenge at points
! from which ppn reaches the ground
!
          DO k=1,wet_model_levels
!
            IF (L_DO_SCAV) THEN
!
              IF (LS_RAIN3D(i,j,k)  >   0) THEN
!
! Calculation of scavenging rate requires S in PPBV
! and rainrate in mm/hr
! Convert S_MMR to S_PPBV
                S_PPBV = S_MMR(i,j,k)*MMR_TO_PPBV
!
! Convert LS_RAIN3D to RAIN_MMPH
                RAIN_MMPH = LS_RAIN3D(i,j,k)*3600.0

! Calculate scavenging rate using function dependent on amount of S
! present:  d(ln(S))/dt = - L1 * (R/S_THOLD)**2/3  if S  <=  S_THOLD
!           d(ln(S))/dt = - L1 * (R/S      )**2/3  if S  >   S_THOLD
!
                IF (RAIN_MMPH <= 0.0) THEN       !Check for neg ppn
                  TERMR = 0.0
                ELSE
!
                  If (S_PPBV <= S_THOLD) Then    !Use fn for low S
                    TERMR = L0 * ( RAIN_MMPH )**0.666667
                  Else                           !Use fn for high S
                    TERMR = L1 * ( RAIN_MMPH/S_PPBV )**0.666667
                  End If
!
                END IF
!
! Calculate amount of tracer scavenged out per tstep
!
                DELTA_TR=S_MMR(i,j,k)*(1.0-EXP(-TERMR*timestep))
!
! Calculate mass of air per unit area in layer for conversion of tracer
! mixing ratio increment to mass p.u.a. for STASH
! Avoid calculations at top model level and assume no scavenging there.
!
                IF (k  <   model_levels) THEN
! Remove the r squared factor from rho_r2 before interpolation
                  rho1= rho_r2(i,j,k)/( r_rho_levels(i,j,k) *           &
     &                                  r_rho_levels(i,j,k) )
                  rho2= rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1) *       &
     &                                    r_rho_levels(i,j,k+1) )
!
! DM = density (interpolated on to theta levels) * delta r
                  DM= rho2 * ( r_theta_levels(i,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) +                  &
     &                rho1 * ( r_rho_levels(i,j,k+1) -                  &
     &                         r_theta_levels(i,j,k) )
!
! Special case for lowest layer to get correct mass
                  If (k  ==  1) Then
                    DM = DM *                                           &
     &                  (r_rho_levels(i,j,2)-r_theta_levels(i,j,0))     &
     &                 /(r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
                  End If
!
! Convert DM to DRY density if level is wet
                  If (k  <=  wet_model_levels) Then
                    DM = DM * (1.0 - Q(i,j,k) - qcl(i,j,k) - qcf(i,j,k))
                  End If
!
                ELSE
                  DELTA_TR = 0.0
                  DM = 0.0
                END IF
!
! Increment accumulated scavenged tracer in column, multiplying by DM
! to convert mmr to mass per unit area.
!
                LSCAV_TR(i,j)=LSCAV_TR(i,j)+DELTA_TR*DM
!
! Decrement tracer mixing ratio
!
                S_MMR(i,j,k)=S_MMR(i,j,k)-DELTA_TR
!
              ELSE                    !level with zero ppn found
!
                L_DO_SCAV=.FALSE.     !stop scavenging in this column
!
              END IF                  !END LS_RAIN3D > 0 condition
!
            END IF                    !END L_DO_SCAV condition
!
          END DO                      !END K loop
!
        End Do                        !End i loop
      End Do                          !End j loop
!
      RETURN
      END SUBROUTINE SL3DWASH
!
