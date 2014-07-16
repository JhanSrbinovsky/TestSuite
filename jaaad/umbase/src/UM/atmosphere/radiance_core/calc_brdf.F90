#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the arrays of BRDF terms.
!
! Purpose:
!   This routine is called to calculate a set of arrays related to
!   the BRDF for later efficiency.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used.
!
! Symmetries of the BRDF and storage:
!   Since the BRDF is defined only for downwelling incident
!   radiances and upwelling reflected radiances it cannot be
!   uniquely defined as a double expension in spherical harmonics.
!   To make a unique expansion we stipulate that only harmonics of
!   even parity will be used: if odd harmonics were chosen we would
!   get into difficulties with the Gibb's phenomenon, as all odd
!   harmonics vanish on the horizontal.
!   F(j, l, l', m) will therefore be 0 unless l+m and l'+m are
!   both even, so the indices of storage for l and l' are set to
!   l/2 and l'/2. There are more efficient schemes of storage
!   that depend on the fact that F vanishes if l<m or l'<m, but
!   such a scheme has not been selected because of its extra
!   complexity.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_BRDF(ISOLIR, MS_MIN, MS_MAX                       &
     &  , IA_SPH_MM                                                     &
     &  , UPLM_SOL, UPLM_ZERO                                           &
     &  , N_BRDF_BASIS_FNC, LS_BRDF_TRUNC, F_BRDF                       &
     &  , N_PROFILE, N_DIRECTION, DIRECTION                             &
     &  , BRDF_SOL, BRDF_HEMI                                           &
     &  , ND_PROFILE, ND_RADIANCE_PROFILE, ND_DIRECTION                 &
     &  , ND_MAX_ORDER, ND_SPH_COEFF                                    &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC                              &
     &  )
!
!
!
      IMPLICIT NONE
!
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_RADIANCE_PROFILE                                           &
!           Size allocated for profiles where radiances are calculated
     &  , ND_DIRECTION                                                  &
!           Size allocated for viewing directions
     &  , ND_MAX_ORDER                                                  &
!           Size allcoated for polar orders
     &  , ND_SPH_COEFF                                                  &
!           Size allocated for spherical coefficients
     &  , ND_BRDF_BASIS_FNC                                             &
!           Size allocated for BRDF basis functions
     &  , ND_BRDF_TRUNC
!           Size allocated for truncation of BRDFs
!
!     Include header files
#include "c_kinds.h"
#include "c_pi.h"
#include "spectral_region_pcf3z.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE
!           Number of atmospheric profiles
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Spectral region
      INTEGER, INTENT(IN) ::                                            &
     &    MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , IA_SPH_MM(0: ND_MAX_ORDER)
!           Address of spherical coefficient for (m, m) for each m
      REAL  (Real64), INTENT(IN) ::                                     &
     &    UPLM_ZERO(ND_SPH_COEFF)
!           Array of Upsilon_l^m and derivatives at polar angles of pi/2
      INTEGER, INTENT(IN) ::                                            &
     &    N_BRDF_BASIS_FNC                                              &
!           Number of basis functions for BRDFs
     &  , LS_BRDF_TRUNC
!           Order of truncation applied to BRDFs
      REAL  (Real64), INTENT(IN) ::                                     &
     &    F_BRDF(ND_BRDF_BASIS_FNC, 0: ND_BRDF_TRUNC/2                  &
     &      , 0: ND_BRDF_TRUNC/2, 0: ND_BRDF_TRUNC)
!           Array of moments of BRDF basis functions
      INTEGER, INTENT(IN) ::                                            &
     &    N_DIRECTION
!           Number of viewing directions
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIRECTION(ND_RADIANCE_PROFILE, ND_DIRECTION, 2)               &
!           Cosines of polar viewing angles and actual azimuthal
!           viewing angles
     &  , UPLM_SOL(ND_PROFILE, ND_SPH_COEFF)
!           Upsilon terms for solar radiation
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    BRDF_SOL(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)         &
!           The BRDF evaluated for scattering from the solar
!           beam into the viewing direction
     &  , BRDF_HEMI(ND_PROFILE, ND_BRDF_BASIS_FNC, ND_DIRECTION)
!           The BRDF evaluated for scattering from isotropic
!           radiation into the viewing direction
!
!
!     Local variables
      INTEGER                                                           &
     &    LS                                                            &
!           Polar order of harmonic
     &  , LSR                                                           &
!           Reduced polar order of harmonic
     &  , LS_P                                                          &
!           Polar order of harmonic
     &  , LSR_P                                                         &
!           Reduced polar order of harmonic
     &  , MS                                                            &
!           Azimuthal order of spherical harmonic
     &  , J                                                             &
!           Loop variable
     &  , L                                                             &
!           Loop variable
     &  , ID
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    UP_LM(ND_PROFILE, ND_BRDF_TRUNC+1, ND_DIRECTION)              &
!           Polar parts of spherical harmonics in the viewing
!           directions
     &  , SS1(ND_PROFILE, 0: ND_BRDF_TRUNC)                             &
!           Products of the BRDF and the solar harmonics
     &  , AZIM_FACTOR(ND_PROFILE, ND_DIRECTION)                         &
!           Azimuthal factors
     &  , KAPPA(ND_BRDF_TRUNC+1)                                        &
!           Hemispherical quadratic integrals of spherical harmonics
!           (reduced storage does not seem worth the effort here)
     &  , FK(ND_BRDF_TRUNC+1)
!           Sum of products of the BRDF and KAPPA over l'
!
!     Subroutines called:
      EXTERNAL                                                          &
     &    EVAL_UPLM
!
!
!
      IF (ISOLIR == IP_SOLAR) THEN

!       Initialize the BRDF for solar radiation
        DO ID=1, N_DIRECTION
          DO J=1, N_BRDF_BASIS_FNC
            DO L=1, N_PROFILE
              BRDF_SOL(L, J, ID)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
!
!       Loop over azimuthal orders.
        DO MS=MS_MIN, MIN(MS_MAX, LS_BRDF_TRUNC)
!
!         Caclulate the azimuthal factors.
          IF (MS == 0) THEN
            DO ID=1, N_DIRECTION
              DO L=1, N_PROFILE
                AZIM_FACTOR(L, ID)=1.0E+00_Real64
              ENDDO
            ENDDO
          ELSE
            DO ID=1, N_DIRECTION
              DO L=1, N_PROFILE
                AZIM_FACTOR(L, ID)                                      &
     &          =2.0E+00_Real64*COS(REAL(MS,Real64)*DIRECTION(L,ID,2))
              ENDDO
            ENDDO
          ENDIF
!
!         Calculate spherical harmonics in the viewing directions
!         at this azimuthal order.
          DO ID=1, N_DIRECTION
! DEPENDS ON: eval_uplm
            CALL EVAL_UPLM(MS, LS_BRDF_TRUNC                            &
     &        , N_PROFILE, DIRECTION(1, ID, 1), UP_LM(1, 1, ID)         &
     &        , ND_PROFILE)
          ENDDO
!
!         Now loop over basis functions.
          DO J=1, N_BRDF_BASIS_FNC
!
!           The array SS1 pulls in the solar dependence, which is
!           independent of the viewing direction. At this stage both
!           MS and J are fixed.


            DO LS=MS, LS_BRDF_TRUNC, 2
              DO L=1, N_PROFILE
                SS1(L, LS)=F_BRDF(J, LS/2, MS/2, MS)                    &
     &           *UPLM_SOL(L, IA_SPH_MM(MS))
              ENDDO

              DO LS_P=MS+2, LS_BRDF_TRUNC, 2
                DO L=1, N_PROFILE
                  SS1(L, LS)=F_BRDF(J, LS/2, LS_P/2, MS)                &
     &             *UPLM_SOL(L, IA_SPH_MM(MS)+LS_P-MS)
                ENDDO
              ENDDO
            ENDDO
!
!           Now consider each direction, incrementing the solar
!           BRDF.
            DO ID=1, N_DIRECTION
              DO LS=MS, LS_BRDF_TRUNC, 2
                DO L=1, N_PROFILE
                  BRDF_SOL(L, J, ID)=BRDF_SOL(L, J, ID)                 &
     &              +SS1(L, LS)*UP_LM(L, LS+1-MS, ID)                   &
     &              *AZIM_FACTOR(L, ID)
                ENDDO
              ENDDO
            ENDDO
!
          ENDDO
!
        ENDDO
!
      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
!       Initialize.
        DO ID=1, N_DIRECTION
          DO J=1, N_BRDF_BASIS_FNC
            DO L=1, N_PROFILE
              BRDF_HEMI(L, J, ID)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDDO
!
!       Only azimuthally symmetric terms contribute.
        DO MS=MS_MIN, 0
!
          DO LSR_P=1, LS_BRDF_TRUNC-MS+1, 2
            KAPPA(LSR_P)=2.0E+00_Real64*PI                              &
     &        *UPLM_ZERO(IA_SPH_MM(0)+LSR_P-1)*UPLM_ZERO(2)             &
     &        /REAL((LSR_P-2)*(LSR+1+2*MS), Real64)
          ENDDO
!
!         Calculate spherical harmonics in the viewing directions
!         at this azimuthal order.
          DO ID=1, N_DIRECTION
! DEPENDS ON: eval_uplm
            CALL EVAL_UPLM(MS, LS_BRDF_TRUNC                            &
     &        , N_PROFILE, DIRECTION(1, ID, 1), UP_LM(1, 1, ID)         &
     &        , ND_PROFILE)
          ENDDO
!
!         Now loop over basis functions.
          DO J=1, N_BRDF_BASIS_FNC
!
            DO LSR=1, LS_BRDF_TRUNC-MS+1, 2
              FK(LSR)=0.0E+00_Real64
              DO LSR_P=1, LS_BRDF_TRUNC-MS+1, 2
                FK(LSR)=FK(LSR)+KAPPA(LSR_P)                            &
     &            *F_BRDF(J, (LSR-1+MS)/2, (LSR_P-1+MS)/2, MS)
              ENDDO
              DO ID=1, N_DIRECTION
                DO L=1, N_PROFILE
                  BRDF_HEMI(L, J, ID)=BRDF_HEMI(L, J, ID)               &
     &              +FK(LSR)*UP_LM(L, LSR, ID)
                ENDDO
              ENDDO
            ENDDO
!
          ENDDO
!
        ENDDO
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE CALC_BRDF
#endif
