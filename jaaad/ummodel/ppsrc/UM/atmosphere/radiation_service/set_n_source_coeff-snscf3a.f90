


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to set number of source coefficients.
!
! Method:
!       The two-stream approximation is examined and the number
!       of coefficients is set accordingly.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION SET_N_SOURCE_COEFF(ISOLIR, L_IR_SOURCE_QUAD              &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     ISOLIR
!             SPECTRAL REGION
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             FLAG FOR QUADRATIC INFRA-RED SOURCE
!
      INTEGER                                                           &
                !, INTENT(OUT)
     &     SET_N_SOURCE_COEFF
!             RETURNED NUMBER OF SOURCE COEFFICIENTS
!
!
!
      IF (ISOLIR == IP_SOLAR) THEN
         SET_N_SOURCE_COEFF=2
      ELSE
         IF (L_IR_SOURCE_QUAD) THEN
            SET_N_SOURCE_COEFF=2
         ELSE
            SET_N_SOURCE_COEFF=1
         ENDIF
      ENDIF
!
!
!
      RETURN
      END FUNCTION SET_N_SOURCE_COEFF
