


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
! STDIO3A defines unit numbers for standard i/o in two-stream radiation
! code.
      INTEGER,PARAMETER:: IU_STDIN=5
      INTEGER,PARAMETER:: IU_STDOUT=6
      INTEGER,PARAMETER:: IU_ERR=6
! STDIO3A end
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
! SOLVER3A defines reference numbers for solvers for two-stream
! radiation code.

      ! pentadiagonal scheme
      INTEGER,PARAMETER:: IP_SOLVER_PENTADIAGONAL=1

      ! mixed column scheme using full endekadiagonal matrix
      INTEGER,PARAMETER:: IP_SOLVER_MIX_11=6

      ! mixed column scheme with approximate scattering
      INTEGER,PARAMETER:: IP_SOLVER_MIX_APP_SCAT=9

      ! direct mixed column scheme for full fluxes
      INTEGER,PARAMETER:: IP_SOLVER_MIX_DIRECT=11

      ! direct solver for a homogeneous column
      INTEGER,PARAMETER:: IP_SOLVER_HOMOGEN_DIRECT=13

      ! direct solver for triple column
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE=14

      ! direct solver for triple column approximating scattering
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE_APP_SCAT=15

      ! direct mixed column scheme for full fluxes (modified
      !   for correct treatment of shadowing by Robin Hogan)
      INTEGER,PARAMETER:: IP_SOLVER_MIX_DIRECT_HOGAN=16

      ! direct solver for triple column (modified for
      !   correct treatment of shadowing by Robin Hogan)
      INTEGER,PARAMETER:: IP_SOLVER_TRIPLE_HOGAN=17

! SOLVER3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
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
