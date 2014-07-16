#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate basis functions for the diffuse albedo.
!
! Purpose:
!   This routine takes the BRDF supplied and calculates a diffuse
!   albedo for isotropic radiation, which is used to define an
!   equivalent extinction.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used. See calc_brdf.f for a note on
!   the symmetries of the BRDF and storage.
!
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE DIFF_ALBEDO_BASIS(N_BRDF_BASIS_FNC                     &
     &  , LS_BRDF_TRUNC, F_BRDF                                         &
     &  , UPLM_ZERO                                                     &
     &  , DIFFUSE_ALB_BASIS                                             &
     &  , ND_BRDF_BASIS_FNC, ND_BRDF_TRUNC, ND_SPH_COEFF                &
     &  )
!
!
      IMPLICIT NONE
!
! Include Header files
#include "c_kinds.h"
!
!     Sizes of arrays
      INTEGER, INTENT(IN) ::                                            &
     &    ND_BRDF_BASIS_FNC                                             &
!           Size allocated for BRDF basis functions
     &  , ND_BRDF_TRUNC                                                 &
!           Size allocated for truncation of BRDFs
     &  , ND_SPH_COEFF
!           Size allocated for spherical coefficients (dimensioning
!           as ND_BRDF_TRUNC+1 might seem natural, but this can
!           lead to problems at low orders of truncation if
!           ND_BRDF_TRUNC is set too large higher in the program.
!
!     Include header files
#include "c_pi.h"
!
!
!     Dummy arguments.
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
!
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    DIFFUSE_ALB_BASIS(ND_BRDF_BASIS_FNC)
!           The diffuse albedo for isotropic radiation calculated
!           from the appropriate BRDF basis function
!
!
!     Local variables
      INTEGER                                                           &
     &    LSR                                                           &
!           Reduced polar order of harmonic
     &  , LSR_P                                                         &
!           Reduced polar order of harmonic
     &  , J
!           Loop variable
!
      REAL  (Real64) ::                                                 &
     &    FACTOR(ND_BRDF_TRUNC+1)                                       &
!           Term involving a sum over l' calculated for speed.
     &  , SUM_P(ND_BRDF_TRUNC+1)
!           Sum of products of the BRDF and factors over l'
!
!
!
!
      DO J=1, N_BRDF_BASIS_FNC
        DIFFUSE_ALB_BASIS(J)=0.0E+00_Real64
!
        DO LSR_P=1, LS_BRDF_TRUNC+1, 2
          FACTOR(LSR_P)=UPLM_ZERO(LSR_P)                                &
     &      *REAL(1-2*MOD(LSR_P-1, 2), Real64)                          &
     &      /REAL((LSR_P-2)*(LSR_P+1), Real64)
        ENDDO
!
        DO LSR=1, LS_BRDF_TRUNC+1, 2
          SUM_P(LSR)=0.0E+00_Real64
          DO LSR_P=1, LS_BRDF_TRUNC+1, 2
            SUM_P(LSR)=SUM_P(LSR)+FACTOR(LSR_P)                         &
     &        *F_BRDF(J, (LSR-1)/2, (LSR_P-1)/2, 0)
          ENDDO
          DIFFUSE_ALB_BASIS(J)=DIFFUSE_ALB_BASIS(J)                     &
     &      +SUM_P(LSR)*UPLM_ZERO(LSR)                                  &
     &      /REAL((LSR-2)*(LSR+1), Real64)
        ENDDO
!
        DIFFUSE_ALB_BASIS(J)=DIFFUSE_ALB_BASIS(J)*4.0E+00_Real64*PI
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE DIFF_ALBEDO_BASIS
#endif
