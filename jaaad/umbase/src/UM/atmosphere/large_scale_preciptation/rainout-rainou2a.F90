#if defined(A17_2A)
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
     &           TRACER,                                                &
     &           RNOUT_TRACER)
!
!---------------------------------------------------------------------
! Purpose: This subroutine removes dissolved tracer aerosol assuming
!          the amount in grid box is reduced in the same proportion as
!          the reduction in the total condensed water (liquid + ice)
!          due to ppn. (It is assumed that the concn of tracer
!          is the same in every droplet and ice particle)
!
!           Called by LSPP_CTL
!
! Current code owners: D L Roberts, M Woodage
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   06/10/01   New Deck (based on vn4.5 RAINOU1A)     M Woodage
!   6.2   11/11/05   Added pre-processor test for
!                    scientific version A17_2A              N Bellouin
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
     &, LS_RAIN3D(row_length,rows,wet_model_levels)                     &
     &, LS_SNOW3D(row_length,rows,wet_model_levels)
!
! Arguments with intent IN/OUT:
      REAL                                                              &
     & TRACER(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &                                         model_levels)
!
! Arguments with intent OUT (diagnostics):
      REAL                                                              &
     & RNOUT_TRACER(row_length,rows)      !tracer removed kg/m2/ts
!
!  Local variables
!
      INTEGER i, j, k              !loop variables
!
      REAL                                                              &
     &   FRACTION                                                       &
                                   !fraction of water remaining
     &,  DELTA_TRACER                                                   &
                                   !amount tracer removed from grid box
     &,  DM                                                             &
                                   !mass p.u.area of air in layer
     &,  TOT_PRECIP3D                                                   &
                                   !total ppn in grid box
     &,  rho1                                                           &
                                   !densities at
     &,  rho2                      !  adjacent model levels
!
      LOGICAL                                                           &
     &   L_DO_SCAV                 !for controlling scavenging
!
! General Method:
! Restrict calculations to points where there is some condensed water,
! and from which some ppn reaches the surface.
! Cases to consider:
! a)If condensed water has increased, leave dissolved tracer unchanged
! b)If condensed water has decreased but is non_neg, reduce dissolved
!   tracer in same ratio
! c)If condensed water has decreased to  <  0, remove all dissolved
!   tracer
!
!   Initialise RNOUT_TRACER to zero before doing rainout
      Do j=1,rows
        Do i=1,row_length
          RNOUT_TRACER(i,j) = 0.0
        End Do
      End Do
!
      Do j=1,rows
        Do i=1,row_length
!
          L_DO_SCAV=.TRUE.
!
          DO k=1,wet_model_levels    !loop over wet levels from surf up
!
            IF (L_DO_SCAV) THEN
!
              TOT_PRECIP3D=LS_RAIN3D(i,j,k)+LS_SNOW3D(i,j,k)
!
              IF (TOT_PRECIP3D  >   0.0)  THEN
!
                IF (QCL_PREVIOUS(i,j,k)+QCF_PREVIOUS(i,j,k)             &
     &                                                    >   0.0) THEN
!
                  IF (QCL_REMAIN(i,j,k)+QCF_REMAIN(i,j,k)               &
     &                >   QCL_PREVIOUS(i,j,k)+QCF_PREVIOUS(i,j,k)) THEN
                    FRACTION=1.0
                  ELSE IF (QCL_REMAIN(i,j,k)+QCF_REMAIN(i,j,k)          &
     &                                                    >=  0.0) THEN
                    FRACTION=(QCL_REMAIN(i,j,k)+QCF_REMAIN(i,j,k))      &
     &                       /(QCL_PREVIOUS(i,j,k)+QCF_PREVIOUS(i,j,k))
                  ELSE
                    FRACTION=0.0
                  ENDIF
!
! Calculate amount TRACER removed per grid box
                  DELTA_TRACER = TRACER(i,j,k) * (1.0-FRACTION)
!
! Calculate mass of air per unit area in layer for conversion of tracer
!  mixing ratio increment to mass p.u.a. for STASH.
! Avoid calculations at top model level and assume no scavenging there.
!
                  IF (k  <   model_levels) THEN
! Remove the r squared factor from rho_r2 before interpolation
                    rho1= rho_r2(i,j,k)/( r_rho_levels(i,j,k) *         &
     &                                    r_rho_levels(i,j,k) )
                    rho2= rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1) *     &
     &                                      r_rho_levels(i,j,k+1) )
!
! DM = density (interpolated on to theta levels) * delta r
                    DM = rho2 * ( r_theta_levels(i,j,k) -               &
     &                            r_rho_levels(i,j,k) ) +               &
     &                   rho1 * ( r_rho_levels(i,j,k+1) -               &
     &                            r_theta_levels(i,j,k) )
!
! Special case for lowest layer to get correct mass
                    If (k  ==  1) Then
                      DM = DM *                                         &
     &                  (r_rho_levels(i,j,2)-r_theta_levels(i,j,0))     &
     &                 /(r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
                    End If
!
! Convert DM to DRY density if level is wet
                    If (k  <=  wet_model_levels) Then
                      DM = DM * (1.0 - Q(i,j,k) -                       &
     &                qcl_previous(i,j,k) - qcf_previous(i,j,k))
                    End If
!
                  ELSE
                    DELTA_TRACER = 0.0
                    DM = 0.0
                  END IF
!
! Accumulate removed TRACER for each level
!
                  RNOUT_TRACER(i,j)=RNOUT_TRACER(i,j)+DELTA_TRACER*DM
!
! Decrement TRACER
!
                  TRACER(i,j,k)=TRACER(i,j,k)-DELTA_TRACER
!
                ENDIF                 !END QCLF_PREVIOUS > 0 condition
!
              ELSE                    !level with zero ppn found
!
                L_DO_SCAV=.FALSE.     !stop scavenging in this column
!
              ENDIF                   !END TOT_PRECIP3D condition
!
            END IF                    !END L_DO_SCAV condition
!
          END DO                      !END k loop
!
        End Do                        !End i loop
      End Do                          !End j loop
!
      RETURN
      END SUBROUTINE RAINOUT
!
#endif
