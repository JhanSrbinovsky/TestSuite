#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      FUNCTION RTOETA(r,orog,z_top_of_model,first_constant_r_rho_level)
! ---------------------------------------------------------------------
! Purpose:
! To convert from r to eta
!
! Method:
!
! Original Programmer: Colin Johnson
!
! Current code owner: Colin Johnson
!
! Version   Date                    Comment
!  5.3    24/05/01  Created.  C.E. Johnson
!  5.5    22/01/04  Commented out error check. K. Ketelsen
!  6.1    22/08/04  Minor tidying of code. M.G. Sanderson
!  6.2    28/03/06  Minor changes for vn6.2  M.G. Sanderson
! ---------------------------------------------------------------------
!
!VVV  V1.0  RtoETA 24/5/01 - Original version
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      IMPLICIT NONE
! ----------------------------------------------------------------------
      INTEGER, INTENT(IN)  :: first_constant_r_rho_level

      REAL, INTENT(IN)  :: z_top_of_model
      REAL, INTENT(IN)  :: Orog
      REAL, INTENT(IN)  :: R

      REAL :: RtoEta
      REAL :: a
      REAL :: b
      REAL :: c
      REAL :: r_a
      CHARACTER(LEN=72):: cmessage

      r_a=r-Earth_radius
!     IF (r_a < orog .or. r_a-z_top_of_model > 1.0) THEN
!       cmessage='R is out of range'
!       WRITE(6,*) cmessage,'R: ',R
!       CALL EREPORT('RtoETA',1,cmessage)
!     ENDIF

      RtoEta=r_a/z_top_of_model

      IF (Orog>1.0E-05 .AND. RtoEta<Eta_Rho(first_constant_r_rho_level))&
     & THEN
        a=Orog/(Eta_Rho(first_constant_r_rho_level))**2
        b=z_top_of_model-(2.0*orog/                                     &
     &      Eta_Rho(first_constant_r_rho_level))
        c=orog-r_a
        RtoEta=(-b+SQRT(b**2-4.0*a*c))/(2.0*a)
      ENDIF

      END FUNCTION RTOETA
#endif
