#if defined(A70_1Z)
#if defined(A01_3Z) || defined(A02_3Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the scattering flag in this a band.
!
! Method:
!       Straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                  Comment
!       6.2             13-02-06              Include code into UM
!                                             (J.-C. Thelen)
!
! Description of code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_SCATTERING(I_SCATTER_METHOD                        &
     &   , L_SWITCH_SCATTER                                             &
     &   , I_SCATTER_METHOD_BAND                                        &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     Include header files.
#include "c_kinds.h"
#include "scatter_method_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    I_SCATTER_METHOD
!           Method of treating scattering
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SWITCH_SCATTER
!           Scattering switch for the band
      INTEGER, INTENT(OUT) ::                                           &
     &    I_SCATTER_METHOD_BAND
!           Scattering flag in this band
!
!
!
      IF (I_SCATTER_METHOD == IP_SCATTER_FULL) THEN
!
!       Perform a full scattering calculation in this band
        I_SCATTER_METHOD_BAND=IP_SCATTER_FULL
!
      ELSE IF (I_SCATTER_METHOD == IP_NO_SCATTER_ABS) THEN
!
!       Scattering extinction is to be neglected if specified in this
!       band. Otherwise a full scattering calculation is required.
        IF (L_SWITCH_SCATTER) THEN
          I_SCATTER_METHOD_BAND=IP_NO_SCATTER_ABS
        ELSE
          I_SCATTER_METHOD_BAND=IP_SCATTER_FULL
        ENDIF
!
      ELSE IF (I_SCATTER_METHOD == IP_NO_SCATTER_EXT) THEN
!
!       Scattering extinction is to be treated as absorption if
!       specified in this band. otherwise a full scattering
!       calculation is required.
        IF (L_SWITCH_SCATTER) THEN
          I_SCATTER_METHOD_BAND=IP_NO_SCATTER_ABS
        ELSE
          I_SCATTER_METHOD_BAND=IP_SCATTER_FULL
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SET_SCATTERING
#endif
#endif
