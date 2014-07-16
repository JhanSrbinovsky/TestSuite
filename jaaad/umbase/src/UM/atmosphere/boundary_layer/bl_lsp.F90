#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE BL_LSP( ROW_LENGTH,ROWS,BL_LEVELS,                     &
     &                   QCF,Q,T )
!
!  Purpose: Convert temperature from liquid ice to liquid, and convert
!           the vapour+liquid+ice variable (Q) to vapour+liquid. This
!           subroutine is used if the mixed phase precipitation scheme
!           is selected AND a full boundary layer treatment is not
!           performed.
!
! D Wilson    <- programmer
!
!  Model            Modification history from model version 4.4:
! version  Date
! 4.4      Sept 97  Originally Coded
!                                                 Damian Wilson
! 5.5      17/04/03  Remove references to obsolete sections
!                    A03_3A,3B,5B,7A. T.White
!
        IMPLICIT NONE
!
        INTEGER                                                         &
     &    ROW_LENGTH                                                    &
                                ! IN   Length of each row
     &,   ROWS                                                          &
                                ! IN   Number of rows
     &,   BL_LEVELS             ! IN   Number of boundary layer levels
!
        REAL                                                            &
     &    QCF(ROW_LENGTH,ROWS,BL_LEVELS)                                &
                                         ! INOUT Ice water content
     &,   Q(ROW_LENGTH,ROWS,BL_LEVELS)                                  &
                                         ! INOUT
!                                  IN    Vapour+liquid+ice content
!                                  OUT   Vapour+liquid content
     &,   T(ROW_LENGTH,ROWS,BL_LEVELS)   ! INOUT
!                                  IN    Liquid ice temperature
!                                  OUT   Liquid temperature
! Temporary Space
        INTEGER                                                         &
     &          I                                                       &
                                 ! Counter over points
     &,         J                                                       &
                                 ! Counter over points
     &,         K                ! Counter over boundary layer levels
        REAL NEWQCF              ! Temporary variable for QCF
#include "c_r_cp.h"
#include "c_lheat.h"
      REAL LSRCP                 ! IN Latent heat of sublimation / Cp
      PARAMETER( LSRCP=((LC+LF)/CP) )
!
!
      DO K=1,BL_LEVELS
        DO J=1,ROWS
          DO I=1, ROW_LENGTH
! Convert Q (vapour+liquid+ice) to (vapour+liquid)
            Q(I,J,K)=Q(I,J,K)-QCF(I,J,K)
! Check that Q is not negative
            IF (Q(I,J,K)  <   0.0) THEN
! Evaporate ice to keep Q positive, but don't let ice go negative
! itself
              NEWQCF=MAX(QCF(I,J,K)+Q(I,J,K),0.0)
              Q(I,J,K)=Q(I,J,K)+(QCF(I,J,K)-NEWQCF)
              QCF(I,J,K)=NEWQCF
            ENDIF
! Adjust T from T liquid ice to T liquid
            T(I,J,K)=T(I,J,K)+LSRCP*QCF(I,J,K)
          END DO
        END DO
      END DO
! End the subroutine
      RETURN
      END SUBROUTINE BL_LSP
#endif
