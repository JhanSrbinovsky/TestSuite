#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve the two-stream equations in a column.
!
! Method:
!       The coefficients of the two-stream equations are calculated.
!       From these we obtain the transmission and reflection
!       coefficients and the source terms. Depending on the solver
!       selected, an appropriate set of matrix equations is formulated
!       and solved to give the fluxes.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE TWO_STREAM(IERR                                        &
!                     Atmospheric Properties
     &  , N_PROFILE, N_LAYER                                            &
!                     Two-stream Scheme
     &  , I_2STREAM                                                     &
!                     Options for Solver
     &  , I_SOLVER                                                      &
!                     Options for Equivalent Extinction
     &  , L_SCALE_SOLAR, ADJUST_SOLAR_KE                                &
!                     Spectral Region
     &  , ISOLIR                                                        &
!                     Infra-red Properties
     &  , DIFF_PLANCK                                                   &
     &  , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                               &
!                     Conditions at TOA
     &  , FLUX_INC_DOWN, FLUX_INC_DIRECT, SEC_0                         &
!                     Surface Conditions
     &  , DIFFUSE_ALBEDO, DIRECT_ALBEDO, D_PLANCK_FLUX_SURFACE          &
!                     Single Scattering Properties
     &  , TAU, OMEGA, ASYMMETRY                                         &
!                     Fluxes Calculated
     &  , FLUX_DIRECT, FLUX_TOTAL                                       &
!                     Dimensions
     &  , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                         &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of dummy arrays.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for profiles
     &  , ND_LAYER                                                      &
!           Size allocated for layers
     &  , ND_SOURCE_COEFF
!           Size allocated for source coefficients
!
!     Include header files.
#include "c_kinds.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
#include "solver_pcf3z.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy variables.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of profiles
     &  , N_LAYER                                                       &
!           Number of layers
     &  , ISOLIR                                                        &
!           Spectral region
     &  , I_SOLVER                                                      &
!           Solver employed
     &  , I_2STREAM
!           Two-stream scheme
      INTEGER, INTENT(INOUT) ::                                         &
     &    IERR
!           Error flag
      LOGICAL, INTENT(IN) ::                                            &
     &    L_SCALE_SOLAR                                                 &
!           Scaling applied to solar flux
     &  , L_IR_SOURCE_QUAD
!           Use quadratic source term
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE, ND_LAYER)                                     &
!           Optical depth
     &  , OMEGA(ND_PROFILE, ND_LAYER)                                   &
!           Albedo of single scattering
     &  , ASYMMETRY(ND_PROFILE, ND_LAYER)                               &
!           Asymmetry
     &  , SEC_0(ND_PROFILE)                                             &
!           Secants of solar zenith angles
     &  , DIFFUSE_ALBEDO(ND_PROFILE)                                    &
!           Diffuse albedo
     &  , DIRECT_ALBEDO(ND_PROFILE)                                     &
!           Direct albedo
     &  , FLUX_INC_DOWN(ND_PROFILE)                                     &
!           Incident total flux
     &  , FLUX_INC_DIRECT(ND_PROFILE)                                   &
!           Incident direct flux
     &  , DIFF_PLANCK(ND_PROFILE, ND_LAYER)                             &
!           Difference in Planckian fluxes across layers
     &  , D_PLANCK_FLUX_SURFACE(ND_PROFILE)                             &
!           Ground source function
     &  , ADJUST_SOLAR_KE(ND_PROFILE, ND_LAYER)                         &
!           Adjustment of solar beam with equivalent extinction
     &  , DIFF_PLANCK_2(ND_PROFILE, ND_LAYER)
!           2x2nd differences of Planckian
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    FLUX_DIRECT(ND_PROFILE, 0: ND_LAYER)                          &
!           Direct flux
     &  , FLUX_TOTAL(ND_PROFILE, 2*ND_LAYER+2)
!           Total fluxes
!
!
!     Local variables.
      INTEGER                                                           &
     &    N_EQUATION
!           Number of equations
      REAL  (Real64) ::                                                 &
     &    TRANS(ND_PROFILE, ND_LAYER)                                   &
!           Transmission of layer
     &  , REFLECT(ND_PROFILE, ND_LAYER)                                 &
!           Reflectance of layer
     &  , TRANS_0(ND_PROFILE, ND_LAYER)                                 &
!           Direct transmittance
     &  , SOURCE_COEFF(ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF)           &
!           Source coefficients
     &  , S_DOWN(ND_PROFILE, ND_LAYER)                                  &
!           Downward source
     &  , S_UP(ND_PROFILE, ND_LAYER)
!           Upward source
      REAL  (Real64) ::                                                 &
     &    A5(ND_PROFILE, 5, 2*ND_LAYER+2)                               &
!           Pentadigonal matrix
     &  , B(ND_PROFILE, 2*ND_LAYER+2)                                   &
!           RHS of matrix equation
     &  , WORK_1(ND_PROFILE, 2*ND_LAYER+2)
!           Working array for solver
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    TWO_COEFF, SOLAR_SOURCE, IR_SOURCE                            &
     &  , SET_MATRIX_PENTADIAGONAL                                      &
     &  , BAND_SOLVER, SOLVER_HOMOGEN_DIRECT
!
!
!
!     Calculate the two-stream coefficients.
! DEPENDS ON: two_coeff
      CALL TWO_COEFF(IERR                                               &
     &  , N_PROFILE, 1, N_LAYER                                         &
     &  , I_2STREAM, L_IR_SOURCE_QUAD                                   &
     &  , ASYMMETRY, OMEGA, TAU                                         &
     &  , ISOLIR, SEC_0                                                 &
     &  , TRANS, REFLECT, TRANS_0                                       &
     &  , SOURCE_COEFF                                                  &
     &  , ND_PROFILE, 1, ND_LAYER, 1, ND_LAYER, ND_SOURCE_COEFF         &
     &  )
      IF (IERR /= I_NORMAL) RETURN
!
!     Calculate the appropriate source terms.
      IF (ISOLIR == IP_SOLAR) THEN
! DEPENDS ON: solar_source
        CALL SOLAR_SOURCE(N_PROFILE, N_LAYER                            &
     &    , FLUX_INC_DIRECT                                             &
     &    , TRANS_0, SOURCE_COEFF                                       &
     &    , L_SCALE_SOLAR, ADJUST_SOLAR_KE                              &
     &    , FLUX_DIRECT                                                 &
     &    , S_DOWN, S_UP                                                &
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
! DEPENDS ON: ir_source
        CALL IR_SOURCE(N_PROFILE, 1, N_LAYER                            &
     &    , SOURCE_COEFF, DIFF_PLANCK                                   &
     &    , L_IR_SOURCE_QUAD, DIFF_PLANCK_2                             &
     &    , S_DOWN, S_UP                                                &
     &    , ND_PROFILE, ND_LAYER, ND_SOURCE_COEFF                       &
     &    )
      ENDIF
!
!     Select an appropriate solver for the equations of transfer.
!
      IF (I_SOLVER == IP_SOLVER_PENTADIAGONAL) THEN
! DEPENDS ON: set_matrix_pentadiagonal
        CALL SET_MATRIX_PENTADIAGONAL(N_PROFILE, N_LAYER                &
     &    , TRANS, REFLECT                                              &
     &    , S_DOWN, S_UP                                                &
     &    , DIFFUSE_ALBEDO, DIRECT_ALBEDO                               &
     &    , FLUX_DIRECT(1, N_LAYER), FLUX_INC_DOWN                      &
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , A5, B                                                       &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
        N_EQUATION=2*N_LAYER+2
!
! DEPENDS ON: band_solver
        CALL BAND_SOLVER(N_PROFILE, N_EQUATION                          &
     &    , 2, 2                                                        &
     &    , A5, B                                                       &
     &    , FLUX_TOTAL                                                  &
     &    , WORK_1                                                      &
     &    , ND_PROFILE, 5, 2*ND_LAYER+2                                 &
     &    )
!
      ELSE IF (I_SOLVER == IP_SOLVER_HOMOGEN_DIRECT) THEN
!
! DEPENDS ON: solver_homogen_direct
        CALL SOLVER_HOMOGEN_DIRECT(N_PROFILE, N_LAYER                   &
     &    , TRANS, REFLECT                                              &
     &    , S_DOWN, S_UP                                                &
     &    , ISOLIR, DIFFUSE_ALBEDO, DIRECT_ALBEDO                       &
     &    , FLUX_DIRECT(1, N_LAYER), FLUX_INC_DOWN                      &
     &    , D_PLANCK_FLUX_SURFACE                                       &
     &    , FLUX_TOTAL                                                  &
     &    , ND_PROFILE, ND_LAYER                                        &
     &    )
!
      ELSE
!
        WRITE(IU_ERR, '(/A)')                                           &
     &    '***Error: The solver and the cloud scheme are incompatiable.'
        IERR=I_ERR_FATAL
        RETURN
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE TWO_STREAM
#endif
