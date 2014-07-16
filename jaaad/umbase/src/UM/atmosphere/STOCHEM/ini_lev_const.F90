#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE INI_LEV_CONST(z_top_of_model,                          &
     & first_constant_r_rho_level)
!-----------------------------------------------------------------------
! Purpose:
! To initialise level dependant constants
!
! Method:
! Only for the smooth generation method
!
! Original Programmer: Colin Johnson
!
! Current code owner: Colin Johnson
!
! History:
! Date        Version     Comment
! -------     -------     -------
! 17/9/01     1.0         Original              Colin Johnson
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
!-----------------------------------------------------------------------
#include "typsize.h"
      INTEGER, INTENT(IN) :: first_constant_r_rho_level
      REAL,    INTENT(IN) :: z_top_of_model

      INTEGER :: I
      REAL :: Eta
      REAL, DIMENSION(0:NMETLEV) :: Zsea_Theta
      REAL, DIMENSION(0:NMETLEV) :: Ck_Theta
      REAL, DIMENSION(NMETLEV)   :: Zsea_Rho
      REAL, DIMENSION(NMETLEV)   :: Ck_Rho

! See RCF_SETUP_LEVDEPC_MOD for the following code:

! Theta
      Zsea_Theta=eta_theta*z_top_of_model

      DO I=0,first_constant_r_rho_level-1
        Ck_theta(I)=(1.0 - eta_theta( I )/                              &
     &     eta_rho(first_constant_r_rho_level))** 2
      ENDDO

      DO I=first_constant_r_rho_level,Model_Levels
        Ck_Theta(I) = 0.0
      ENDDO

! Rho
      Zsea_Rho=eta_rho*z_top_of_model

      DO I=1,first_constant_r_rho_level-1
        Ck_Rho(I)=(1.0 - eta_rho( I )/                                  &
     &      eta_rho(first_constant_r_rho_level))** 2
      ENDDO

      DO I=first_constant_r_rho_level,Model_Levels
        Ck_Rho(I) = 0.0
      ENDDO

! STOCHEM - full levels

      Zsea_Stochem=Eta_Stochem*z_top_of_model

      WHERE (Eta_stochem < eta_rho(first_constant_r_rho_level))
        Ck_stochem=(1.0-Eta_stochem/                                    &
     &             eta_rho(first_constant_r_rho_level))**2
      ELSEWHERE
        Ck_Stochem=0.0
      ENDWHERE

! STOCHEM - half levels

      DO I=1,NLEV
        Eta=(Eta_stochem(I)+Eta_stochem(I-1))/2.0
        Zsea_Stochem_half(I)=Eta*z_top_of_model
        IF (Eta < eta_rho(first_constant_r_rho_level)) THEN
          Ck_Stochem_half(I)=(1.0-Eta/                                  &
     &             eta_rho(first_constant_r_rho_level))**2
        ELSE
          Ck_Stochem_half(I)=0.0
        ENDIF
      ENDDO

      END SUBROUTINE INI_LEV_CONST
#endif
