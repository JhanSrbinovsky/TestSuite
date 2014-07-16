


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the properties of the surface.
!
! Method:
!       The albedo of the surface and the infra-red source
!       function are set according to the parametrization in use.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_SURFACE_PROPERTIES(N_POINT_TYPE, INDEX_SURFACE     &
     &   , I_SPEC_SURFACE                                               &
     &   , ISOLIR, IB                                                   &
     &   , SURFACE_ALBEDO, ALBEDO_FIELD_DIFF, ALBEDO_FIELD_DIR          &
     &   , N_DIR_ALBEDO_FIT, DIRECT_ALBEDO_PARM, SEC_0                  &
     &   , EMISSIVITY_GROUND, EMISSIVITY_FIELD                          &
     &   , ALBEDO_SURFACE_DIFF, ALBEDO_SURFACE_DIR                      &
     &   , THERMAL_GROUND_BAND                                          &
     &   , NPD_PROFILE, NPD_BAND, NPD_SURFACE, NPD_ALBEDO_PARM          &
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
     &   , NPD_SURFACE                                                  &
!             MAXIMUM NUMBER OF SURFACES
     &   , NPD_BAND                                                     &
!             MAXIMUM NUMBER OF BANDS
     &   , NPD_ALBEDO_PARM
!             MAXIMUM NUMBER OF ALBEDO PARAMETERS
!
!     COMDECK CALLS
! SRFSP3A defines permitted methods of specifying the  surface albedo and
! emissivity for two-stream radiation code.

      ! properties specified by surface type
      INTEGER,PARAMETER:: IP_SURFACE_SPECIFIED=1

      ! properties passed into code
      INTEGER,PARAMETER:: IP_SURFACE_INTERNAL=2

      ! direct albedo fitted as polynomial
      INTEGER,PARAMETER:: IP_SURFACE_POLYNOMIAL=3

      ! fit in the functional form used by payne
      INTEGER,PARAMETER:: IP_SURFACE_PAYNE=4
! SRFSP3A end
! SPCRG3A defines flags for different portions of the spectrum in
! two-stream radiation code.
      INTEGER,PARAMETER:: IP_SOLAR=1
      INTEGER,PARAMETER:: IP_INFRA_RED=2
! SPCRG3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_POINT_TYPE(NPD_SURFACE)                                    &
!             NUMBER OF POINTS OF EEACH TYPE
     &   , INDEX_SURFACE(NPD_PROFILE, NPD_SURFACE)                      &
!             LIST OF POINTS OF EACH TYPE
     &   , I_SPEC_SURFACE(NPD_SURFACE)                                  &
!             METHOD OF SPECIFYING SURFACES
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , IB                                                           &
!             NUMBER OF BAND
     &   , N_DIR_ALBEDO_FIT(NPD_SURFACE)
!             NUMBER OF PARAMETERS IN FIT TO DIRECT ALBEDO
      REAL                                                              &
                !, INTENT(IN)
     &     SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE)                        &
!             BAND-DEPENDENT SURFACE ALBEDOS
     &   , ALBEDO_FIELD_DIFF(NPD_PROFILE)                               &
!             SPECIFIED DIFFUSE ALBEDO FIELD
     &   , ALBEDO_FIELD_DIR(NPD_PROFILE)                                &
!             SPECIFIED DIRECT ALBEDO FIELD
     &   , DIRECT_ALBEDO_PARM(0: NPD_ALBEDO_PARM, NPD_BAND, NPD_SURFACE)&
!             COEFFICIENTS FOR DIRECT ALBEDOS
     &   , SEC_0(NPD_PROFILE)                                           &
!             SECANTS OF ZENITH ANGLES
     &   , EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)                     &
!             BAND-DEPENDENT EMISSIVITIES
     &   , EMISSIVITY_FIELD(NPD_PROFILE)
!             SPECIFIED EMISSIVITIES
      REAL                                                              &
                !, INTENT(INOUT)
     &     THERMAL_GROUND_BAND(NPD_PROFILE)
!             THERMAL SOURCE FROM GROUND
      REAL                                                              &
                !, INTENT(OUT)
     &     ALBEDO_SURFACE_DIFF(NPD_PROFILE)                             &
!             DIFFUSE SURFACE ALBEDOS
     &   , ALBEDO_SURFACE_DIR(NPD_PROFILE)
!             DIRECT SURFACE ALBEDOS

!
!     LOCAL VARIABLES.
      INTEGER                                                           &
     &     J                                                            &
!             LOOP VARIABLE
     &   , K                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , IC
!             LOOP VARIABLE
!
!
!     SET THE SURFACE ALBEDOS.
      DO K=1, NPD_SURFACE
!
         IF (I_SPEC_SURFACE(K) == IP_SURFACE_SPECIFIED) THEN
            DO J=1, N_POINT_TYPE(K)
               L=INDEX_SURFACE(J, K)
               ALBEDO_SURFACE_DIFF(L)                                   &
     &            =SURFACE_ALBEDO(IB, K)
               ALBEDO_SURFACE_DIR(L)                                    &
     &            =SURFACE_ALBEDO(IB, K)
            ENDDO
         ENDIF
!
         IF (I_SPEC_SURFACE(K) == IP_SURFACE_INTERNAL) THEN
            DO J=1, N_POINT_TYPE(K)
               L=INDEX_SURFACE(J, K)
               ALBEDO_SURFACE_DIFF(L)=ALBEDO_FIELD_DIFF(L)
               ALBEDO_SURFACE_DIR(L)=ALBEDO_FIELD_DIR(L)
            ENDDO
         ENDIF
!
         IF (I_SPEC_SURFACE(K) == IP_SURFACE_POLYNOMIAL) THEN
            DO J=1, N_POINT_TYPE(K)
               L=INDEX_SURFACE(J, K)
               ALBEDO_SURFACE_DIFF(L)                                   &
     &            =SURFACE_ALBEDO(IB, K)
               ALBEDO_SURFACE_DIR(L)                                    &
     &            =DIRECT_ALBEDO_PARM(N_DIR_ALBEDO_FIT(K), IB, K)
            ENDDO
            DO IC=N_DIR_ALBEDO_FIT(K), 1, -1
               DO J=1, N_POINT_TYPE(K)
                  L=INDEX_SURFACE(J, K)
                  ALBEDO_SURFACE_DIR(L)                                 &
     &               =ALBEDO_SURFACE_DIR(L)/SEC_0(L)                    &
     &               +DIRECT_ALBEDO_PARM(IC-1, IB, K)
               ENDDO
            ENDDO
         ENDIF
!
         IF (I_SPEC_SURFACE(K) == IP_SURFACE_PAYNE) THEN
            DO J=1, N_POINT_TYPE(K)
               L=INDEX_SURFACE(J, K)
               ALBEDO_SURFACE_DIFF(L)=0.06E+00
               ALBEDO_SURFACE_DIR(L)=DIRECT_ALBEDO_PARM(1, IB, K)       &
     &           /(DIRECT_ALBEDO_PARM(2, IB, K)                         &
     &           +DIRECT_ALBEDO_PARM(3, IB, K)                          &
     &           *EXP(-DIRECT_ALBEDO_PARM(4, IB, K)*LOG(SEC_0(L))))
            ENDDO
         ENDIF
      ENDDO
!
!     SET THE EMISSIVITY AND MULTIPLY THE SOURCE FUNCTION
!     IN THE INFRA-RED.
      IF (ISOLIR == IP_INFRA_RED) THEN
         DO K=1, NPD_SURFACE
!
            IF (I_SPEC_SURFACE(K) == IP_SURFACE_SPECIFIED) THEN
               DO J=1, N_POINT_TYPE(K)
                  L=INDEX_SURFACE(J, K)
                  THERMAL_GROUND_BAND(L)=EMISSIVITY_GROUND(IB, K)       &
     &               *THERMAL_GROUND_BAND(L)
               ENDDO
            ENDIF
!
            IF (I_SPEC_SURFACE(K) == IP_SURFACE_INTERNAL) THEN
               DO J=1, N_POINT_TYPE(K)
                  L=INDEX_SURFACE(J, K)
                  THERMAL_GROUND_BAND(L)=EMISSIVITY_FIELD(L)            &
     &               *THERMAL_GROUND_BAND(L)
               ENDDO
            ENDIF
!
            IF (I_SPEC_SURFACE(K) == IP_SURFACE_POLYNOMIAL) THEN
               DO J=1, N_POINT_TYPE(K)
                  L=INDEX_SURFACE(J, K)
                  THERMAL_GROUND_BAND(L)=EMISSIVITY_GROUND(IB, K)       &
     &               *THERMAL_GROUND_BAND(L)
               ENDDO
            ENDIF
!
            IF (I_SPEC_SURFACE(K) == IP_SURFACE_PAYNE) THEN
               DO J=1, N_POINT_TYPE(K)
                  L=INDEX_SURFACE(J, K)
!                 SINCE THE EMISSIVITY IS 1.0 CONTINUE.
                  CONTINUE
               ENDDO
            ENDIF
!
         ENDDO
      ENDIF
!
!
      RETURN
      END SUBROUTINE SET_SURFACE_PROPERTIES
