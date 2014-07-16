


! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find maximally overlapped energy transfer coefficients.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             17-05-96                Add code for new
!                                               solvers.
!       4.5             18-05-98                Reference to obsolete
!                                               solvers removed.
!                                               (J. M. Edwards)
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE OVERLAP_MIX_RANDOM(N_PROFILE, N_LAYER, N_CLOUD_TOP     &
     &   , ISOLIR, I_SOLVER                                             &
     &   , W_CLOUD, W_FREE                                              &
     &   , CLOUD_OVERLAP                                                &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
!
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
! CLCFPT3A defines pointers in CLOUD_OVERLAP for two-stream radiation
! code.
!
! note that several pointers are identical since only certain
! groups of coefficients are relevant to a particular scheme.

      INTEGER,PARAMETER:: IP_CLOVLP_GFF=1 ! GAMMA-FREE-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_GFC=2 ! GAMMA-FREE-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_GCF=3 ! GAMMA-CLOUD-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_GCC=4 ! GAMMA-CLOUD-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_BFF=5 ! BETA-FREE-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_BFC=6 ! BETA-FREE-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_BCF=7 ! BETA-CLOUD-FREE
      INTEGER,PARAMETER:: IP_CLOVLP_BCC=8 ! BETA-CLOUD-CLOUD
      INTEGER,PARAMETER:: IP_CLOVLP_GFM=5 ! GAMMA_F-
      INTEGER,PARAMETER:: IP_CLOVLP_GFP=6 ! GAMMA_F+
      INTEGER,PARAMETER:: IP_CLOVLP_BFM=7 ! BETA_F-
      INTEGER,PARAMETER:: IP_CLOVLP_BFP=8 ! BETA_F+
      INTEGER,PARAMETER:: IP_CLOVLP_GM=5  ! GAMMA_-
      INTEGER,PARAMETER:: IP_CLOVLP_GP=6  ! GAMMA_+
      INTEGER,PARAMETER:: IP_CLOVLP_BM=7  ! BETA_-
      INTEGER,PARAMETER:: IP_CLOVLP_BP=8  ! BETA_+

      ! pointers for triple overlaps:

      INTEGER,PARAMETER:: IP_CLOVLP_V11=1
      INTEGER,PARAMETER:: IP_CLOVLP_V12=2
      INTEGER,PARAMETER:: IP_CLOVLP_V13=3
      INTEGER,PARAMETER:: IP_CLOVLP_V21=4
      INTEGER,PARAMETER:: IP_CLOVLP_V22=5
      INTEGER,PARAMETER:: IP_CLOVLP_V23=6
      INTEGER,PARAMETER:: IP_CLOVLP_V31=7
      INTEGER,PARAMETER:: IP_CLOVLP_V32=8
      INTEGER,PARAMETER:: IP_CLOVLP_V33=9
      INTEGER,PARAMETER:: IP_CLOVLP_U11=10
      INTEGER,PARAMETER:: IP_CLOVLP_U12=11
      INTEGER,PARAMETER:: IP_CLOVLP_U13=12
      INTEGER,PARAMETER:: IP_CLOVLP_U21=13
      INTEGER,PARAMETER:: IP_CLOVLP_U22=14
      INTEGER,PARAMETER:: IP_CLOVLP_U23=15
      INTEGER,PARAMETER:: IP_CLOVLP_U31=16
      INTEGER,PARAMETER:: IP_CLOVLP_U32=17
      INTEGER,PARAMETER:: IP_CLOVLP_U33=18

! CLCFPT3A end
!
!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , N_LAYER                                                      &
!             NUMBER OF LAYERS
     &   , N_CLOUD_TOP                                                  &
!             TOPMOST CLOUDY LAYER
     &   , ISOLIR                                                       &
!             SPECTRAL REGION
     &   , I_SOLVER
!             SOLVER TO BE USED
      REAL                                                              &
                !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             CLOUD AMOUNTS
!
      REAL                                                              &
                !, INTENT(OUT)
     &     W_FREE(NPD_PROFILE, NPD_LAYER)                               &
!             CLOUD-FREE AMOUNTS
     &   , CLOUD_OVERLAP(NPD_PROFILE, 0: NPD_LAYER, NPD_OVERLAP_COEFF)
!             COEFFICIENTS FOR TRANSFER OF ENERGY AT INTERFACE
!
!
!     LOCAL ARGUMENTS.
      INTEGER                                                           &
     &     I                                                            &
     &   , L
!
!
!     SET THE FREE FRACTIONS IN EACH LAYER.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00
         ENDDO
      ENDDO
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            W_FREE(L, I)=1.0E+00-W_CLOUD(L, I)
         ENDDO
      ENDDO
!
!     EVALUATE THE EXTENT OF OVERLAP BETWEEN LAYERS OF CLOUD
!     AT THE INTERFACE BETWEEN THE ITH AND (I+1)ST LAYER ON THE
!     ASSUMPTION OF RANDOM OVERLAP BETWEEN ADJACENT LAYERS.
!     THE TOP AND BOTTOM BOUNDARIES ARE EXCEPTIONAL.
!
!     IN THE SOLAR REGION COEFFICIENTS FOR DOWNWARD COUPLING OF THE
!     FLUXES ARE REQUIRED. THESE COEFFICIENTS ARE ALSO NEEDED FOR
!     INFRA-RED CALCULATIONS WITH APPROXIMATE SCATTERING.
!
      IF ( (I_SOLVER == IP_SOLVER_MIX_DIRECT).OR.                       &
     &     (I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN).OR.                 &
     &     (ISOLIR == IP_SOLAR).OR.                                     &
     &     ( (ISOLIR == IP_INFRA_RED).AND.                              &
     &       (I_SOLVER == IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO I=N_CLOUD_TOP-1, N_LAYER-1
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFF)=W_FREE(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GFC)=W_CLOUD(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCF)=W_FREE(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GCC)=W_CLOUD(L, I+1)
            ENDDO
         ENDDO
!
         DO L=1, N_PROFILE
!           BOTTOM BOUNDARY:
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GFF)=1.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GFC)=0.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GCF)=1.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GCC)=0.0E+00
         ENDDO
!
      ENDIF
!
!     WITH APPROXIMATE SCATTERING IN THE LONGWAVE THE CORRESPONDING
!     UPWARD COEFFICIENTS ARE NEEDED.
!
      IF ( (I_SOLVER == IP_SOLVER_MIX_DIRECT).OR.                       &
     &     (I_SOLVER == IP_SOLVER_MIX_DIRECT_HOGAN).OR.                 &
     &     ( (ISOLIR == IP_INFRA_RED).AND.                              &
     &       (I_SOLVER == IP_SOLVER_MIX_APP_SCAT) ) ) THEN
!
         DO L=1, N_PROFILE
!           TOP CLOUDY BOUNDARY:
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BFF)=1.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BFC)=1.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BCF)=0.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BCC)=0.0E+00
         ENDDO
!
         DO I=N_CLOUD_TOP, N_LAYER
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFF)=W_FREE(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BFC)=W_FREE(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCF)=W_CLOUD(L, I)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BCC)=W_CLOUD(L, I)
            ENDDO
         ENDDO
      ENDIF
!
!
      IF (I_SOLVER == IP_SOLVER_MIX_11) THEN
!
         DO L=1, N_PROFILE
!
!           TOP BOUNDARY:
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_GM)=0.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_GP)               &
     &         =2.0E+00-4.0E+00*W_CLOUD(L, N_CLOUD_TOP)
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BM)=0.0E+00
            CLOUD_OVERLAP(L, N_CLOUD_TOP-1, IP_CLOVLP_BP)=2.0E+00
!
!           BOTTOM BOUNDARY:
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GM)=0.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_GP)=2.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BM)=0.0E+00
            CLOUD_OVERLAP(L, N_LAYER, IP_CLOVLP_BP)                     &
     &         =2.0E+00-4.0E+00*W_CLOUD(L, N_LAYER)
!
         ENDDO
!
         DO I=N_CLOUD_TOP, N_LAYER-1
            DO L=1, N_PROFILE
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GM)=0.0E+00
               CLOUD_OVERLAP(L, I, IP_CLOVLP_GP)                        &
     &            =2.0E+00-4.0E+00*W_CLOUD(L, I+1)
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BM)=0.0E+00
               CLOUD_OVERLAP(L, I, IP_CLOVLP_BP)                        &
     &            =2.0E+00-4.0E+00*W_CLOUD(L, I)
            ENDDO
         ENDDO
!
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE OVERLAP_MIX_RANDOM
