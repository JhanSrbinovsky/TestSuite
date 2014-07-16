#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to determine whether densities are required for clouds.
!
! Method:
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-06-96                L_CLOUD_DENSITY set
!                                               as .FALSE. initially
!                                               to provide a default.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      FUNCTION L_CLOUD_DENSITY(N_CONDENSED, I_PHASE_CMP, L_CLOUD_CMP    &
     &   , I_CONDENSED_PARAM                                            &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     INCLUDE COMDECKS
#include "dimfix3a.h"
#include "iclprm3a.h"
#include "phase3a.h"
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_CONDENSED                                                  &
!             NUMBER OF TYPES OF CONDENSATE
     &   , I_PHASE_CMP(NPD_CLOUD_COMPONENT)                             &
!             PHASES OF COMPONENTS
     &   , I_CONDENSED_PARAM(NPD_CLOUD_COMPONENT)
!             PARAMETRIZATIONS OF COMPONENTS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_CLOUD_CMP(NPD_CLOUD_COMPONENT)
!             FLAGS FOR ENABLED COMPONENTS
      LOGICAL                                                           &
                !, INTENT(OUT)
     &     L_CLOUD_DENSITY
!             RETURNED FLAG FOR CALCULATING DENSITY
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     K
!             LOOP VARIABLE
!
!
      L_CLOUD_DENSITY=.FALSE.
!
!     DENSITIES MUST BE CALCULATED IF SUN & SHINE'S PARAMETRIZATIONS
!     ARE USED.
      DO K=1, N_CONDENSED
         L_CLOUD_DENSITY=L_CLOUD_DENSITY.OR.                            &
     &      (L_CLOUD_CMP(K).AND.(I_PHASE_CMP(K) == IP_PHASE_ICE).AND.   &
     &      ( (I_CONDENSED_PARAM(K) == IP_SUN_SHINE_VN2_VIS).OR.        &
     &        (I_CONDENSED_PARAM(K) == IP_SUN_SHINE_VN2_IR) ) )
      ENDDO
!
!
!
      RETURN
      END FUNCTION L_CLOUD_DENSITY
#endif
#endif
