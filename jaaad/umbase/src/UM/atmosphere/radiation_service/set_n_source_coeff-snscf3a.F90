#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
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
#include "spcrg3a.h"
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
#endif
#endif
