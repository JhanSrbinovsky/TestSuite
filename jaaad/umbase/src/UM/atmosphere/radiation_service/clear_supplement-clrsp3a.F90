#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate clear-sky fluxes.
!
! Method:
!       This subroutine is called after fluxes including clouds have
!       been calculated to find the corresponding clear-sky fluxes.
!       The optical properties of the column are already known.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             10-04-96                New solver added
!                                               (J. M. Edwards)
!       4.5             18-05-98                Obsolete solvers
!                                               removed.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER              &
     &   , I_SOLVER_CLEAR                                               &
     &   , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE    &
     &   , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN                       &
     &   , S_DOWN_FREE, S_UP_FREE                                       &
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                      &
     &   , SOURCE_GROUND                                                &
     &   , L_SCALE_SOLAR, ADJUST_SOLAR_KE                               &
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                          &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_PROFILE                                                  &
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
#include "stdio3a.h"
#include "spcrg3a.h"
#include "solver3a.h"
#include "error3a.h"
#include "dimfix3a.h"
!
!     DUMMY VARIABLES.
      INTEGER                                                           &
                !, INTENT(OUT)
     &     IERR
!             ERROR FLAG
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , I_SOLVER_CLEAR
!             SOLVER FOR CLEAR FLUXES
      LOGICAL                                                           &
                !, INTENT(IN)
     &     L_SCALE_SOLAR
!             SCALING APPLIED TO SOLAR BEAM
      REAL                                                              &
            !, INTENT(IN)
     &     TRANS_FREE(NPD_PROFILE, NPD_LAYER)                           &
!             TRANSMISSION COEFFICIENTS
     &   , REFLECT_FREE(NPD_PROFILE, NPD_LAYER)                         &
!             REFLECTION COEFFICIENTS
     &   , TRANS_0_FREE(NPD_PROFILE, NPD_LAYER)                         &
!             DIRECT TRANSMISSION COEFFICIENTS
     &   , SOURCE_COEFF_FREE(NPD_PROFILE, NPD_LAYER, NPD_SOURCE_COEFF)  &
!             COEFFICIENTS IN SOURCE TERMS
     &   , S_DOWN_FREE(NPD_PROFILE, NPD_LAYER)                          &
!             DOWNWARD SOURCE
     &   , S_UP_FREE(NPD_PROFILE, NPD_LAYER)                            &
!             UPWARD SOURCE
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE ALBEDO
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)                              &
!             DIRECT ALBEDO
     &   , FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT TOTAL FLUX
     &   , FLUX_INC_DIRECT(NPD_PROFILE)                                 &
!             INCIDENT DIRECT FLUX
     &   , SOURCE_GROUND(NPD_PROFILE)                                   &
!             GROUND SOURCE FUNCTION
     &   , ADJUST_SOLAR_KE(NPD_PROFILE, NPD_LAYER)
!             SCALING OF SOLAR BEAM
!
!
      REAL                                                              &
            !, INTENT(OUT)
     &     FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)                 &
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUXES
!
!
!     DUMMY VARIABALES.
      INTEGER                                                           &
     &     N_EQUATION
!             NUMBER OF EQUATIONS
      REAL                                                              &
     &     A3(NPD_PROFILE, 3, 2*NPD_LAYER+2)                            &
!             TRIDIAGONAL MATRIX
     &   , A5(NPD_PROFILE, 5, 2*NPD_LAYER+2)                            &
!             PENTADIAGONAL MATRIX
     &   , B(NPD_PROFILE, 2*NPD_LAYER+2)                                &
!             RHS OF MATRIX EQUATION
     &   , WORK_1(NPD_PROFILE, 2*NPD_LAYER+2)                           &
!             WORKING ARRAY FOR SOLVER
     &   , WORK_2(NPD_PROFILE, 2*NPD_LAYER+2)
!             WORKING ARRAY FOR SOLVER
!
!     SUBROUTINES CALLED:
      EXTERNAL                                                          &
     &     SOLAR_SOURCE                                                 &
     &   , SET_MATRIX_PENTADIAGONAL                                     &
     &   , BAND_SOLVER, SOLVER_HOMOGEN_DIRECT
!
!
!     THE SOURCE FUNCTIONS ONLY NEED TO BE RECALCULATED IN THE VISIBLE.
      IF (ISOLIR == IP_SOLAR) THEN
! DEPENDS ON: solar_source
         CALL SOLAR_SOURCE(N_PROFILE, N_LAYER                           &
     &      , FLUX_INC_DIRECT                                           &
     &      , TRANS_0_FREE, SOURCE_COEFF_FREE                           &
     &      , L_SCALE_SOLAR, ADJUST_SOLAR_KE                            &
     &      , FLUX_DIRECT_CLEAR                                         &
     &      , S_DOWN_FREE, S_UP_FREE                                    &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
      ENDIF
!
!
!     SELECT AN APPROPRIATE SOLVER FOR THE EQUATIONS OF TRANSFER.
!
      IF (I_SOLVER_CLEAR == IP_SOLVER_PENTADIAGONAL) THEN
!
!        CALCULATE THE ELEMENTS OF THE MATRIX EQUATIONS.
! DEPENDS ON: set_matrix_pentadiagonal
         CALL SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER               &
     &      , TRANS_FREE, REFLECT_FREE                                  &
     &      , S_DOWN_FREE, S_UP_FREE                                    &
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                   &
     &      , FLUX_DIRECT_CLEAR(1, N_LAYER), FLUX_INC_DOWN              &
     &      , SOURCE_GROUND                                             &
     &      , A5, B                                                     &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
         N_EQUATION=2*N_LAYER+2
!
! DEPENDS ON: band_solver
         CALL BAND_SOLVER(N_PROFILE, N_EQUATION                         &
     &      , 2, 2                                                      &
     &      , A5, B                                                     &
     &      , FLUX_TOTAL_CLEAR                                          &
     &      , NPD_PROFILE, 2*NPD_LAYER+2                                &
     &      , WORK_1                                                    &
     &      )
!
      ELSE IF (I_SOLVER_CLEAR == IP_SOLVER_HOMOGEN_DIRECT) THEN
!
!        SOLVE FOR THE FLUXES IN THE COLUMN DIRECTLY.
! DEPENDS ON: solver_homogen_direct
         CALL SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER                  &
     &      , TRANS_FREE, REFLECT_FREE                                  &
     &      , S_DOWN_FREE, S_UP_FREE                                    &
     &      , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                   &
     &      , FLUX_DIRECT_CLEAR(1, N_LAYER), FLUX_INC_DOWN              &
     &      , SOURCE_GROUND                                             &
     &      , FLUX_TOTAL_CLEAR                                          &
     &      , NPD_PROFILE, NPD_LAYER                                    &
     &      )
!
      ELSE
!
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** ERROR: THE SOLVER SPECIFIED IS NOT VALID '             &
     &      //'FOR CLEAR FLUXES.'
         IERR=I_ERR_FATAL
         RETURN
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE CLEAR_SUPPLEMENT
#endif
#endif
