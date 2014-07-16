#if defined(A70_1Z)
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
! Current owner of code: James Manners
!
! Description of code:
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
!     Include header files.
#include "c_kinds.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic infra-red source
!
      INTEGER ::                                                        &
     &    SET_N_SOURCE_COEFF
!           Returned number of source coefficients
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
