#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find the optical depth and single scattering albedo.
!
! Method:
!       Depending on the treatment of scattering, the optical and
!       and single scattering albedo are determined from the
!       extinctions supplied.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       5.1             04-04-00                Tolerances replaced
!                                               by F90 intrinsics.
!                                               (J. M. Edwards)
!       5.3             04-10-01                Obsolete code removed.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SINGLE_SCATTERING(I_SCATTER_METHOD_BAND                &
     &   , N_PROFILE, N_LAYER, I_TOP                                    &
     &   , D_MASS                                                       &
     &   , K_GREY_TOT, K_EXT_SCAT, K_GAS_ABS                            &
     &   , TAU, OMEGA                                                   &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
#include "sctmth3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             TREATMENT OF SCATTERING IN THIS BAND
!
!                       Atmospheric Properties
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , I_TOP
!             TOP LAYER TO CONSIDER
      REAL                                                              &
                !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
!
!                       Optical Propeties
      REAL                                                              &
                !, INTENT(IN)
     &     K_GREY_TOT(NPD_PROFILE, NPD_LAYER)                           &
!             ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT(NPD_PROFILE, NPD_LAYER)                           &
!             SCATTERING EXTINCTION
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS EXTINCTION
!
!                       Single Scattering Propeties
      REAL                                                              &
                !, INTENT(OUT)
     &     TAU(NPD_PROFILE, NPD_LAYER)                                  &
!             OPTICAL DEPTH
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             SINGLE SCATTERING ALBEDO
!
!
!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     L                                                            &
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
!
!
!
!     THE MACHINE TOLERANCE IS ADDED TO THE DENOMINATOR IN THE
!     EXPRESSION FOR OMEGA TO PREVENT DIVISION BY ZERO: THIS IS
!     SIGNIFICANT ONLY IF THE TOTAL EXTINCTION IS SMALL, AND THUS
!     WILL NOT SENSIBLY AFFECT ANY PHYSICAL RESULTS.
!
      IF (I_SCATTER_METHOD_BAND == IP_SCATTER_FULL) THEN
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I))             &
     &            *D_MASS(L, I)
               OMEGA(L, I)=K_EXT_SCAT(L, I)                             &
     &            /(K_GREY_TOT(L, I)+K_GAS_ABS(L, I)+TINY(OMEGA))
            ENDDO
         ENDDO
!
      ELSE IF (I_SCATTER_METHOD_BAND == IP_NO_SCATTER_ABS) THEN
!
!        THE SCATTERING EXTINCTION IS IGNORED COMPLETELY, SO
!        ONLY THE ABSORPTIVE CONTRIBUTIONS TO THE SINGLE
!        SCATTERING PROPETIES ARE INCLUDED.
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I)              &
     &            -K_EXT_SCAT(L, I))*D_MASS(L, I)
               OMEGA(L, I)=0.0
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE SINGLE_SCATTERING
#endif
#endif
