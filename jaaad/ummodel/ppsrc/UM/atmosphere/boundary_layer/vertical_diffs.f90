
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
 ! SUBROUTINE VERTICAL_DIFFS

      SUBROUTINE VERTICAL_DIFFS (                                       &
! IN Field dimensioning/pointers.
     &  rows, row_length, n_rows, halo_i, halo_j, off_x, off_y          &
     &, BL_LEVELS                                                       &
     &, TIMESTEP, LQ_MIX_BL                                             &
! IN Vertical coordinate information.
     &, R_RHO_LEVELS                                                    &
     &, R_THETA_LEVELS                                                  &
! IN Fields.
     &, RHO_P_RSQ, rho_wet, rho_dry                                     &
! OUT Vertical differences required by physics.
     &, RHO_uv, RHO_tq, RHO_DRY_tq                                      &
     &, DZL_charney, RDZ                                                &
     &, Z1_UV, Z1_TQ                                                    &
     &, RDZ_CHARNEY_GRID                                                &
     &, DTRDZ_CHARNEY_GRID                                              &
     &, dtrdz_u, dtrdz_v, rdz_u, rdz_v                                  &
     &  )
! PURPOSE:
!
! METHOD:
!
! HISTORY:
! DATE   VERSION   COMMENT
! ----   -------   -------
! 10/5/02  UM5.4   Change DTRDZ arrays to reflect spherical geometry
!                  and change rdz_u,v to remove assumption that
!                  rho-levels are midway between theta-levels.
!                                             Adrian Lock
!
! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

! LOGICAL SWITCHES

      LOGICAL LQ_MIX_BL         ! IN switch for using mixing ratios

      INTEGER                                                           &
     &  rows, row_length, n_rows, halo_i, halo_j                        &
     &, BL_LEVELS, off_x, off_y

      REAL                                                              &
     &  R_RHO_LEVELS (1-halo_i:row_length+halo_i,                       &
     &                  1-halo_j:rows+halo_j, BL_LEVELS+1)              &
     &, R_THETA_LEVELS (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:BL_LEVELS+1)            &
     &, RHO_P_RSQ(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &            BL_LEVELS+1)                                          &
!                       ! wet density times r^2 on rho levels (kg/m3)
     &, rho_wet(row_length, rows, bl_levels+1)                          &
!                       ! wet density on rho levels (kg/m3)
     &, rho_dry(row_length, rows, bl_levels+1)                          &
!                       ! dry density on rho levels (kg/m3)
     &, TIMESTEP

! ARGUMENTS WITH INTENT IN/OUT. IE: INPUT VARIABLES CHANGED ON OUTPUT.


! ARGUMENTS WITH INTENT OUT. IE: OUTPUT VARIABLES.

      REAL                                                              &
     &  RHO_uv(row_length, rows, BL_LEVELS+1)                           &
!                               ! OUT density on UV (ie. rho) levels;
!                               !    used in RHOKH so dry density if
!                               !    Lq_mix_bl is true
     &, RHO_tq(row_length, rows, BL_LEVELS)                             &
!                               ! OUT density on TQ (ie. theta) levels;
!                               !    used in RHOKM so wet density
     &, RHO_DRY_tq(row_length, rows, BL_LEVELS)                         &
!                               ! OUT density on TQ (ie. theta) levels;
!                               !    used in non-turb flux integration
!                               !    so dry density if Lq_mix_bl is true
     &, DTRDZ_CHARNEY_GRID (row_length,rows, BL_LEVELS)                 &
     &, DZL_charney(row_length,rows,BL_LEVELS)                          &
     &, RDZ(row_length,rows,BL_LEVELS)                                  &
                                ! OUT  Reciprocal of distance between
     &, RDZ_CHARNEY_GRID(row_length,rows,BL_LEVELS)                     &
     &, Z1_UV(row_length,rows)                                          &
     &, Z1_TQ(row_length,rows)

      Real                                                              &
              ! quantities interpolated to u and v grids.
     &  rdz_u (row_length, rows, 2:bl_levels)                           &
     &, rdz_v (row_length, n_rows, 2:bl_levels)                         &
     &, dtrdz_u (row_length, rows, bl_levels)                           &
     &, dtrdz_v (row_length, n_rows, bl_levels)

! LOCAL VARIABLES.

      INTEGER                                                           &
     &  K                                                               &
     &, I,J

      REAL                                                              &
     &  DZ_P                                                            &
     &, DZ_T

      Real                                                              &
     &  work1(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)    &
     &, work2(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS+1)  &
     &, work3(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)

      Real                                                              &
     &  dzl_u(row_length, rows, bl_levels)                              &
     &, dzl_v(row_length, n_rows, bl_levels)

! ----------------------------------------------------------------------
! CALCULATE rho on layer centres and layer boundaries
! ----------------------------------------------------------------------
      IF ( lq_mix_bl ) THEN
!       ! conservation of mixing ratios requires use of dry density
        DO K=1,BL_LEVELS+1
          Do j= 1,rows
          DO i= 1,row_length
            RHO_uv(i,j,k) = rho_dry(i,j,k)
          Enddo
          Enddo
        Enddo
      ELSE
        DO K=1,BL_LEVELS+1
          Do j= 1,rows
          DO i= 1,row_length
            RHO_uv(i,j,k) = rho_wet(i,j,k)
          Enddo
          Enddo
        Enddo
      ENDIF
!
! Interpolate RHO_uv to temperature levels for RHO_DRY_tq
!  - can't think of a better name but this will be wet or dry
!    depending on lq_mix_bl
!
! DEPENDS ON: p_to_t
      CALL P_TO_T (                                                     &
! IN Field dimensions and pointers.
     &  row_length, rows, halo_i, halo_j, 0, 0, BL_LEVELS               &
! IN Vertical coordinate levels.
     &, r_theta_levels, r_rho_levels                                    &
! IN field to be interpolated.
     &, RHO_uv                                                          &
! OUT Interpolated field.
     &, RHO_DRY_tq                                                      &
     & )
!
! Interpolate RHO_wet to temperature levels - RHO_tq
!
! DEPENDS ON: p_to_t
      CALL P_TO_T (                                                     &
! IN Field dimensions and pointers.
     &  row_length, rows, halo_i, halo_j, 0, 0, BL_LEVELS               &
! IN Vertical coordinate levels.
     &, r_theta_levels, r_rho_levels                                    &
! IN field to be interpolated.
     &, RHO_WET                                                         &
! OUT Interpolated field.
     &, RHO_tq                                                          &
     & )


! ----------------------------------------------------------------------
! CALCULATE DTRDZ, DZL, RDZ
! ----------------------------------------------------------------------
      k = 1
      Do J=1,rows
        DO I = 1,row_length
          DZ_P = R_theta_LEVELS(I,j,K)-R_THETA_LEVELS(I,j,k-1)
          DZ_T = R_RHO_LEVELS(I,j,K+1)-R_THETA_LEVELS(I,j,k-1)
          RDZ(I,j,K) = 1./                                              &
     &                   (R_RHO_LEVELS(I,j,K) - R_THETA_LEVELS(I,j,K-1))
          Z1_UV(I,j) = R_RHO_LEVELS(I,j,K) - R_THETA_LEVELS (I,j,k-1)
          Z1_TQ(I,j) = R_THETA_LEVELS(I,j,K) -                          &
     &                     R_THETA_LEVELS (I,j,k-1)
          RDZ_CHARNEY_GRID(I,j,K) = 1./ DZ_P
!         ! Following is dt/(rho * r^2 * Dz), used for flux-divergence
!         ! equation for scalars - hence dry density if lq_mix_bl
!         ! is true, wet otherwise
          DTRDZ_CHARNEY_GRID(I,j,K) = TIMESTEP/                         &
     &              ( R_THETA_LEVELS (I,j,K)*R_THETA_LEVELS (I,j,K)*    &
     &                                      RHO_DRY_TQ(I,j,K)*DZ_T )
! Dzl_charney( ,,1) is such that
! 0.5 * DZL_Charney = height of first temperature level
          DZL_charney(I,j,K) = 2. *                                     &
     &               (R_theta_LEVELS(I,j,1) - R_theta_LEVELS(I,j,0))
        End Do
      End Do

      DO K = 2, BL_LEVELS
        Do J=1,rows
          DO I = 1,row_length
            DZ_P = R_THETA_LEVELS (I,j,K) - R_THETA_LEVELS(I,j,K-1)
            DZ_T = R_RHO_LEVELS (I,j,K+1) - R_RHO_LEVELS(I,j,K)
            RDZ(I,j,K) = 1./                                            &
     &                (R_RHO_LEVELS(I,j,K) - R_RHO_LEVELS(I,j,K-1))
            DZL_charney(I,j,K) = DZ_T
            RDZ_CHARNEY_GRID(I,j,K) = 1./ DZ_P
!         ! Following is dt/(rho * r^2 * Dz), used for flux-divergence
!         ! equation for scalars - hence dry density if lq_mix_bl
!         ! is true, wet otherwise
            DTRDZ_CHARNEY_GRID(I,j,K) = TIMESTEP/                       &
     &              ( R_THETA_LEVELS (I,j,K)*R_THETA_LEVELS (I,j,K)*    &
     &                                      RHO_DRY_TQ(I,j,K)*DZ_T )
          END DO
        END DO
      END DO

      Do k = 1, BL_levels
        Do J=1,n_rows+off_y
          DO I = 1,row_length+off_x
            work1(I,j,K) = R_theta_LEVELS(I,j,K)-R_THETA_LEVELS(I,j,k-1)
          End Do
        End Do
      End Do

      Do k = 1, BL_levels
        Do J=1,n_rows+off_y
          DO I = 1,row_length+off_x
            work3(I,j,k) = TIMESTEP/( RHO_P_RSQ(I,J,K) * work1(I,j,K) )
          End Do
        End Do
      End Do

      Do k = 2, BL_levels
        Do J=1,n_rows+off_y
          DO I = 1,row_length+off_x
!           ! Note work1(K=1) will not be correct so dzl_u,v(K=1) below
!           ! will also be wrong but the K=1 values are not used.
            work1(I,j,K) = R_RHO_LEVELS(I,j,K)-R_RHO_LEVELS(I,j,k-1)
          End Do
        End Do
      End Do

! ----------------------------------------------------------------------
! horizontal interpolations to momentum points.
! ----------------------------------------------------------------------

! DEPENDS ON: p_to_u
      Call P_TO_U(work3,row_length,rows,bl_levels,                      &
     &            off_x, off_y, dtrdz_u)

! DEPENDS ON: p_to_u
      Call P_TO_U(work1,row_length,rows,bl_levels,                      &
     &            off_x, off_y, dzl_u)

! DEPENDS ON: p_to_v
      Call P_TO_V(work3,row_length,rows,n_rows,                         &
     &            bl_levels, off_x, off_y, dtrdz_v)

! DEPENDS ON: p_to_v
      Call P_TO_V(work1,row_length,rows,n_rows,                         &
     &            bl_levels, off_x, off_y, dzl_v)

! calculate rdz_u, rdz_v from local arrays dzl_u, dzl_v.

      do k=2,bl_levels
        do j= 1,rows
          do i= 1,row_length
            rdz_u(I,j,k) = 1.0/dzl_u(I,j,k)
          end do
        end do
      end do

      do k= 2,bl_levels
        do j= 1,n_rows
          do i= 1,row_length
            rdz_v(I,j,k) = 1.0/dzl_v(I,j,k)
          end do
        end do
      end do

      RETURN
      END SUBROUTINE VERTICAL_DIFFS
