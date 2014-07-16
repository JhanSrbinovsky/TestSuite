#if defined(A04_3B) || defined(A04_3C) || defined(A04_3D)
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine MASS_CALC
!
! Purpose:
!   To calculate mass of air in a model layer
!
!   Called by mcr_ctl
!
! Current owners of code: S.Woodward
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   5.5      12/02/03  Based on code in RAINOUT.  S Wooward
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! Documentation: UMDP20
!
!----------------------------------------------------------------------

      SUBROUTINE MASS_CALC(                                             &
     &           ROW_LENGTH, ROWS,                                      &
     &           MODEL_LEVELS, WET_MODEL_LEVELS,                        &
     &           R_RHO_LEVELS, R_THETA_LEVELS,                          &
     &           TIMESTEP,                                              &
     &           RHO_R2,                                                &
     &           Q, QCL, QCF,                                           &
     &           DM                                                     &
     &           )
!
! calculates mass of (dry) air per square metre
!

!
      IMPLICIT NONE
!
! Arguments with intent IN:
!
      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
     &, ROWS                                                            &
     &, MODEL_LEVELS                                                    &
     &, WET_MODEL_LEVELS                                                &
     &, CCLDBASE(ROW_LENGTH,ROWS)                                       &
                                       !convective cloud base
     &, CCLDTOP(ROW_LENGTH,ROWS)       !convective cloud top
!
      REAL                                                              &
     &  R_RHO_LEVELS(ROW_LENGTH,ROWS,MODEL_LEVELS)                      &
     &, R_THETA_LEVELS(ROW_LENGTH,ROWS,0:MODEL_LEVELS)                  &
     &, RHO_R2(ROW_LENGTH,ROWS,MODEL_LEVELS)                            &
                                              !density*r*r
     &, Q(ROW_LENGTH,ROWS,WET_MODEL_LEVELS)                             &
                                                   !water vapour(kg/kg)
     &, QCL(ROW_LENGTH,ROWS,WET_MODEL_LEVELS)                           &
     &, QCF(ROW_LENGTH,ROWS,WET_MODEL_LEVELS)

!
      REAL                                                              &
     & TIMESTEP                    !timestep in secs
!
!
! Arguments with intent OUT :
      REAL                                                              &
     & DM(ROW_LENGTH, ROWS, MODEL_LEVELS)     !mass of air
!
! Local variables
!
      INTEGER                                                           &
     & I, J, K                                                          &
                                    ! Loop variables
     &,TOPLEV                       ! level loop limit
!
      REAL                                                              &
     & RHO1, RHO2                ! air densities
!
!
!
!
      DO K = 1,MODEL_LEVELS-1
        DO J = 1,ROWS
          DO I = 1,ROW_LENGTH
! Remove the r squared factor from rho before interpolation
            RHO1 = RHO_R2(I,J,K) /                                      &
     &          ( R_RHO_LEVELS(I,J,K) * R_RHO_LEVELS(I,J,K) )
            RHO2 = RHO_R2(I,J,K+1)/                                     &
     &          ( R_RHO_LEVELS(I,J,K+1) *  R_RHO_LEVELS(I,J,K+1) )
! DM = density (interpolated on to theta levels) * delta r
            DM(I,J,K) = RHO2 * ( R_THETA_LEVELS(I,J,K) -                &
     &                         R_RHO_LEVELS(I,J,K) ) +                  &
     &               RHO1 * ( R_RHO_LEVELS(I,J,K+1) -                   &
     &                         R_THETA_LEVELS(I,J,K) )
!
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !MODEL_LEVELS

! Special case for lowest layer to get correct mass
       DO J = 1,ROWS
         DO I = 1,ROW_LENGTH
           DM(I,J,1) =                                                  &
     &     DM(I,J,1) * (R_RHO_LEVELS(I,J,2) - R_THETA_LEVELS(I,J,0)) /  &
     &           (R_RHO_LEVELS(I,J,2) - R_RHO_LEVELS(I,J,1))
         ENDDO !ROW_LENGTH
       ENDDO !ROWS

 ! Convert DM to DRY density if level is wet
      TOPLEV=MODEL_LEVELS - 1
      IF (WET_MODEL_LEVELS  <   TOPLEV) TOPLEV = WET_MODEL_LEVELS
      DO K = 1,TOPLEV
        DO J = 1,ROWS
          DO I = 1,ROW_LENGTH
            DM(I,J,K) = DM (I,J,K) /                                    &
     &              (1.0 + Q(I,J,K) + QCL(I,J,K) + QCF(I,J,K))
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !WET_MODEL_LEVELS

! Top level
       DO J = 1,ROWS
         DO I = 1,ROW_LENGTH
           DM(I,J,MODEL_LEVELS) = 0.0
         ENDDO !ROW_LENGTH
       ENDDO !ROWS
!
!
!
      RETURN
      END SUBROUTINE MASS_CALC
!
#endif
