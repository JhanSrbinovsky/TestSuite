#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      FUNCTION ETAtoR(Eta,Orog,z_top_of_model,                          &
     &                first_constant_r_rho_level)
! ----------------------------------------------------------------------
! Purpose:
! To convert from eta to r
!
! Method:
! Only for the smooth generation method
!
! Original Programmer: Colin Johnson
!
! Current code owner: Colin Johnson
!
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.3    24/05/01  Created.  C.E. Johnson
!  5.3    05/09/01  Converted to a function.  C.E. Johnson
!  5.5    22/01/04  Commented out error check. K. Ketelsen
!  6.1    22/08/04  Minor tidying of code. M.G. Sanderson
!
!VVV  V1.1  ETAtoR  5/9/01 - Converted to Function
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: first_constant_r_rho_level

      REAL, INTENT(IN)  :: Eta
      REAL, INTENT(IN)  :: Orog
      REAL, INTENT(IN)  :: z_top_of_model

      REAL :: ETAtoR
!kk      CHARACTER(LEN=72):: cmessage

      REAL :: Ck_eta

!kk      IF (Eta < 0.0 .OR. Eta > 1.0 ) THEN
!kk        cmessage='Eta out of range '
!kk        WRITE(6,*) cmessage,' Eta: ',Eta
!kk        CALL EREPORT('ETAtoR',1,cmessage)
!kk      ENDIF

      IF (Eta < Eta_Rho(first_constant_r_rho_level)) THEN
        Ck_eta=(1.0-Eta/eta_rho(first_constant_r_rho_level))**2
      ELSE
        Ck_eta=0.0
      ENDIF

      ETAtoR=Eta*z_top_of_model+Ck_eta*Orog+Earth_radius

      END FUNCTION ETAtoR
#endif
