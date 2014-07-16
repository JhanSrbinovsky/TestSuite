
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates constants used in large-scale precipitation scheme.
!
!  SUBROUTINE GAMMAF--------------------------------------------------
!   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
!   A POLYNOMIAL APPROXIMATION
! --------------------------------------------------------------------
         SUBROUTINE GAMMAF(Y,GAM)
           IMPLICIT NONE
           REAL                                                         &
                              !, INTENT(IN)
     &       Y
           REAL                                                         &
                              !, INTENT(OUT)
     &       GAM
! Gamma function of Y
!
! LOCAL VARIABLE
           INTEGER I,M
           REAL GG,G,PARE,X
! --------------------------------------------------------------------
           GG=1.
           M=Y
           X=Y-M
           IF (M /= 1) THEN
             DO I=1,M-1
               G=Y-I
               GG=GG*G
             END DO
           END IF
           PARE=-0.5748646*X+0.9512363*X*X-0.6998588*X*X*X              &
     &     +0.4245549*X*X*X*X-0.1010678*X*X*X*X*X+1.
           GAM=PARE*GG
           RETURN
         END SUBROUTINE GAMMAF
!
