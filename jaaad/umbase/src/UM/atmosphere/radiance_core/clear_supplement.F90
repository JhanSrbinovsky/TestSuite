#if defined(A70_1Z)
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
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CLEAR_SUPPLEMENT(IERR, N_PROFILE, N_LAYER              &
     &  , I_SOLVER_CLEAR                                                &
     &  , TRANS_FREE, REFLECT_FREE, TRANS_0_FREE, SOURCE_COEFF_FREE     &
     &  , ISOLIR, FLUX_INC_DIRECT, FLUX_INC_DOWN                        &
     &  , S_DOWN_FREE, S_UP_FREE                                        &
     &  , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                       &
     &  , SOURCE_GROUND                                                 &
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
     &  , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR                           &
     &  , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                         &
     &  )
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_LAYER                                                      &
!           Size allocated for atmospheric layers
     &  , ND_SOURCE_COEFF
!           Size allocated for layers
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "spectral_region_pcf3z.h"
#include "solver_pcf3z.h"
#include "error_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , ISOLIR                                                        &
!           Spectral region
     &  , I_SOLVER_CLEAR
!           Solver for clear fluxes
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR
!           Scaling applied to solar beam
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TRANS_FREE(ND_PROFILE, ND_LAYER)                              &
!           Transmission coefficients
     &  , REFLECT_FREE(ND_PROFILE, ND_LAYER)                            &
!           Reflection coefficients
     &  , TRANS_0_FREE(ND_PROFILE, ND_LAYER)                            &
!           Direct transmission coefficients
     &  , SOURCE_COEFF_FREE(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)      &
!           Coefficients in source terms
     &  , S_DOWN_FREE(ND_PROFILE, ND_LAYER)                             &
!           Downward source
     &  , S_UP_FREE(ND_PROFILE, ND_LAYER)                               &
!           Upward source
     &  , ALBEDO_SURFACE_DIFF(ND_PROFILE)                               &
!           Diffuse albedo
     &  , ALBEDO_SURFACE_DIR(ND_PROFILE)                                &
!           Direct albedo
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident total flux
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , SOURCE_GROUND(ND_PROFILE)                                     &
!           Ground source function
     &  , ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)
!           Scaling of solar beam
!
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT_CLEAR(ND_PROFILE, 0: ND_LAYER)                    &
!           Clear direct flux
     &  , FLUX_TOTAL_CLEAR(ND_PROFILE, 2*ND_LAYER+2)
!           Clear total fluxes
!
!
!     Dummy variabales.
      INTEGER                                                           &
     &    N_EQUATION
!           Number of equations
      REAL  (Real64) ::                                                 &
     &    A5(ND_PROFILE, 5, 2*ND_LAYER+2)                               &
!           Pentadiagonal matrix
     &  , B(ND_PROFILE, 2*ND_LAYER+2)                                   &
!           Rhs of matrix equation
     &  , WORK_1(ND_PROFILE, 2*ND_LAYER+2)
!           Working array for solver
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    SOLAR_SOURCE, SET_MATRIX_PENTADIAGONAL                        &
     &  , BAND_SOLVER, SOLVER_HOMOGEN_DIRECT
!
!
!
!     The source functions only need to be recalculated in the visible.
      IF (ISOLIR == IP_SOLAR) THEN
! DEPENDS ON: solar_source
        CALL SOLAR_SOURCE(N_PROFILE, N_LAYER                            &
     &    , FLUX_INC_DIRECT                                             &
     &    , TRANS_0_FREE, SOURCE_COEFF_FREE                             &
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
     &    , FLUX_DIRECT_CLEAR                                           &
     &    , S_DOWN_FREE, S_UP_FREE                                      &
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
      ENDIF
!
!
!     Select an appropriate solver for the equations of transfer.
!
      IF (I_SOLVER_CLEAR == IP_SOLVER_PENTADIAGONAL) THEN
!
!       Calculate the elements of the matrix equations.
! DEPENDS ON: set_matrix_pentadiagonal
        CALL SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER                &
     &    , TRANS_FREE, REFLECT_FREE                                    &
     &    , S_DOWN_FREE, S_UP_FREE                                      &
     &    , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                     &
     &    , FLUX_DIRECT_CLEAR(1, N_LAYER), FLUX_INC_DOWN                &
     &    , SOURCE_GROUND                                               &
     &    , A5, B                                                       &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
        N_EQUATION=2*N_LAYER+2
!
! DEPENDS ON: band_solver
        CALL BAND_SOLVER(N_PROFILE, N_EQUATION                          &
     &    , 2, 2                                                        &
     &    , A5, B                                                       &
     &    , FLUX_TOTAL_CLEAR                                            &
     &    , WORK_1                                                      &
     &    , ND_PROFILE, 5, 2*ND_LAYER+2                                 &
     &    )
!
      ELSE IF (I_SOLVER_CLEAR == IP_SOLVER_HOMOGEN_DIRECT) THEN
!
!       Solve for the fluxes in the column directly.
! DEPENDS ON: solver_homogen_direct
        CALL SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER                   &
     &    , TRANS_FREE, REFLECT_FREE                                    &
     &    , S_DOWN_FREE, S_UP_FREE                                      &
     &    , ISOLIR, ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR             &
     &    , FLUX_DIRECT_CLEAR(1, N_LAYER), FLUX_INC_DOWN                &
     &    , SOURCE_GROUND                                               &
     &    , FLUX_TOTAL_CLEAR                                            &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
!
      ELSE
!
         WRITE(IU_ERR, '(/A)')                                          &
     &      '*** Error: The solver specified for clear-sky fluxes '     &
     &      //'is not valid.'
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
