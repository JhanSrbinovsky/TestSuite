
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ PC2 Cloud Scheme: Inhomogeneous forcing
! Subroutine Interface:
      SUBROUTINE PC2_INHOMOGENEOUS(                                     &
!      Array dimensions
     & levels, row_length,rows                                          &
!      Prognostic Fields
     &,CF, CFL, CFF, QCL, QCF                                           &
!      Forcing quantities for driving the inhomogeneous forcing
     &,Q4_L, Q4_F, LS, G_L, G_F)
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates the change in liquid cloud fraction,
!   ice cloud fraction and total cloud fraction as a result of
!   inhomogeneous forcing of the gridbox with Q4 increments.
!
! Method:
!   Uses the method proposed in Annex B of the PC2 cloud scheme project
!   report, partitioned between liquid and ice phases by the method
!   given in Annex D. Condensate forcing is assumed to take place in
!   conjunction with cloudy air. This subroutine does NOT update vapour,
!   temperature and condensate values, which are assumed to have been
!   calculated by the physics scheme which is responsible for the call
!   to this subroutine. NOTE: This subroutine does NOT check that
!   increments calculated are sensible.
!
! Current Owner of Code: D. R. Wilson
!
! History:
! Version   Date     Comment
!  5.4    22-07-02   Original Code (Damian Wilson)
!
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation
!
!  Global Variables:----------------------------------------------------
!      None
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER                                                           &
                        !, INTENT(IN)
     & LEVELS                                                           &
!       No. of levels being processed.
     &,row_length,rows
!       Row length and number of rows being processed.
!
      REAL                                                              &
                        !, INTENT(IN)
     & QCL(row_length,rows,LEVELS)                                      &
!       Liquid content (kg water per kg air)
     &,QCF(row_length,rows,LEVELS)                                      &
!       Ice content (kg water per kg air)
     &,Q4_L(row_length,rows,LEVELS)                                     &
!       Rate of change of liquid content with time from the forcing
!       (kg kg-1 s-1)
     &,Q4_F(row_length,rows,LEVELS)                                     &
!       Rate of change of ice content with time from forcing mechanism
!       (kg kg-1 s-1)
     &,LS(row_length,rows,LEVELS)                                       &
!       In cloud condensate content of the cloud which is forced into
!       the gridbox (liquid plus ice) (kg kg-1)
     &,G_L(row_length,rows,LEVELS)                                      &
!       Volume fraction of injected cloud which contains liquid
     &,G_F(row_length,rows,LEVELS)
!       Volume fraction of injected cloud which contains ice. Note that
!       G_L + G_F need not equal one if there is overlap between liquid
!       and ice
!
      REAL                                                              &
                        !, INTENT(INOUT)
     & CF(row_length,rows,LEVELS)                                       &
!       Total cloud fraction (no units)
     &,CFL(row_length,rows,LEVELS)                                      &
!       Liquid cloud fraction (no units)
     &,CFF(row_length,rows,LEVELS)
!       Liquid cloud fraction (no units)
!
!  External functions:
!
!  Local parameters and other physical constants------------------------
      REAL                                                              &
     & TOLERANCE
!       Tolerance to avoid a divide by zero
!
      PARAMETER(                                                        &
     &          TOLERANCE=1.0E-10                                       &
     &          )
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL                                                              &
     & DENOM                                                            &
                  ! Denominator in calculation
     &,Q4         ! Sum of Q4_L and Q4_F
!
!  (b) Others.
      INTEGER K,I,J       ! Loop counters: K - vertical level index
!                           I,J - horizontal position index
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
!
! Loop round levels to be processed
! Levels_do1:
      DO k=LEVELS,1,-1
! Rows_do2:
        DO j=1,rows
! Row_length_do2:
          DO I=1,row_length
!
! Calculate total Q4. Only perform the calculation if the total Q4 is
! non-zero.
!
            Q4 = Q4_L(i,j,k) + Q4_F(i,j,k)
            IF (Q4  /=  0.0) THEN
!
! Calculate the change in total cloud fraction
!
              DENOM = (LS(i,j,k) - QCL(i,j,k) - QCF(i,j,k) )
!
              IF ( ABS(DENOM)  >   TOLERANCE ) THEN
!
                DENOM = 1.0 / DENOM
                CF(i,j,k)  = CF(i,j,k)  +                               &
     &                       ( 1.0        - CF(i,j,k)  ) * Q4 * DENOM
                CFL(i,j,k) = CFL(i,j,k) +                               &
     &                       ( G_L(i,j,k) - CFL(i,j,k) ) * Q4 * DENOM
                CFF(i,j,k) = CFF(i,j,k) +                               &
     &                       ( G_F(i,j,k) - CFF(i,j,k) ) * Q4 * DENOM
!
! Otherwise cloud fraction will go to one. In theory, cloud fraction
! can never be reduced by this process.
!
              ELSE
                CF(i,j,k)  = 1.0
                CFL(i,j,k) = 1.0
                CFF(i,j,k) = 1.0
              END IF
!
            END IF
! Row_length_do1:
          END DO
! Rows_do1:
        END DO
! Levels_do1:
      END DO
!
! End of the subroutine
!
      RETURN
      END SUBROUTINE PC2_INHOMOGENEOUS
