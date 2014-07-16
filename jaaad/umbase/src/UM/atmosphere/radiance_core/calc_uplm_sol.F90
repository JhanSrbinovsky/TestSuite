#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate Upsilon_l^m(0) for the solar direction.
!
! Purpose:
!   This routine is called to determine the values of spherical
!   harmonics in the solar direction.
!
! Method:
!   As this routine is called only once speed is not too critical
!   so direct calculation is used.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE CALC_UPLM_SOL(N_PROFILE, MS_MIN, MS_MAX, IA_SPH_MM     &
     &  , LS_LOCAL_TRUNC, ZEN_0, UPLM_SOL                               &
     &  , ND_PROFILE, ND_MAX_ORDER, ND_SPH_COEFF)
!
!
      IMPLICIT NONE
!
!
!     Include header files
#include "c_kinds.h"
#include "c_pi.h"
!
!
!     Dummy arguments.
      INTEGER, INTENT(IN) ::                                            &
     &    ND_PROFILE                                                    &
!           Size allocated for atmospheric profiles
     &  , ND_SPH_COEFF                                                  &
!           Number of spherical coefficients
     &  , ND_MAX_ORDER
!           Maximum order of calculation
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric profiles
     &  , MS_MIN                                                        &
!           Lowest azimuthal order calculated
     &  , MS_MAX                                                        &
!           Highest azimuthal order calculated
     &  , IA_SPH_MM(0: ND_MAX_ORDER)                                    &
!           Address of spherical oefficient for (m, m) for each m
     &  , LS_LOCAL_TRUNC(0: ND_MAX_ORDER)
!           Local truncation at this order
      REAL  (Real64), INTENT(IN) ::                                     &
     &    ZEN_0(ND_PROFILE)
!           Cosines of solar zenith angles
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    UPLM_SOL(ND_PROFILE, ND_SPH_COEFF)
!           Array of Upsilon_l^m evaluated in the solar direction.
!
!
!     Local variables
      INTEGER                                                           &
     &    LS                                                            &
!           Order of harmonic
     &  , MS                                                            &
!           Azimuthal quantum number of harmonic
     &  , J                                                             &
!           Temporary address
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    PRODUCT                                                       &
!           Factorial terms of Y_lm
     &  , LSR                                                           &
!           Real polar order of harmonic
     &  , MSR
!           Real azimuthal order of harmonic
!
!
!
!     Note here that ZEN_0 holds the cosine of the zenith angle, so
!     the cosine of the solar direction is actually -ZEN_0.
      DO MS=MS_MIN, MS_MAX
!
!       Calculate Upsilon_m^m(n_sol) to start the recurrence.
        J=IA_SPH_MM(MS)
!
        PRODUCT=1.0E+00_Real64
        MSR=REAL(MS, Real64)
        IF (MS >  0) THEN
          DO K=1, MS
            PRODUCT=(1.0E+00_Real64-5.0E-01_Real64/REAL(K, Real64))     &
     &                                                    *PRODUCT
          ENDDO
          DO L=1, N_PROFILE
            UPLM_SOL(L, J)=(-1.0E+00_Real64)**MS                        &
     &        *SQRT((1.0E+00_Real64-ZEN_0(L)*ZEN_0(L))**MS*PRODUCT      &
     &        *(2.0E+00_Real64*MSR+1.0E+00_Real64)/(4.0E+00_Real64*PI))
          ENDDO
        ELSE
          DO L=1, N_PROFILE
            UPLM_SOL(L, J)=1.0E+00_Real64/SQRT(4.0E+00_Real64*PI)
          ENDDO
        ENDIF
!
!       Calculate the next polar order to enable the recurrence to
!       start.
        IF (MS <= LS_LOCAL_TRUNC(MS)+1) THEN
          DO L=1, N_PROFILE
            UPLM_SOL(L, J+1)                                            &
     &        =-ZEN_0(L)*SQRT(2.0E+00_Real64*MSR                        &
     &        +3.0E+00_Real64)*UPLM_SOL(L, J)
          ENDDO
        ENDIF
!
!       Complete the recurrence on l.
        DO LS=MS+2, LS_LOCAL_TRUNC(MS)+1
          J=IA_SPH_MM(MS)+LS-MS
          LSR=REAL(LS, Real64)
          DO L=1, N_PROFILE
            UPLM_SOL(L, J)                                              &
     &        =SQRT(((2.0E+00_Real64*LSR-1.0E+00_Real64)                &
     &        *(2.0E+00_Real64*LSR+1.0E+00_Real64))                     &
     &        /((LSR+MSR)*(LSR-MSR)))                                   &
     &        *(-ZEN_0(L))*UPLM_SOL(L, J-1)                             &
     &        -SQRT(((2.0E+00_Real64*LSR+1.0E+00_Real64)                &
     &        *(LSR-1.0E+00_Real64-MSR)*(LSR-1.0E+00_Real64+MSR))       &
     &        /((2.0E+00_Real64*LSR-3.0E+00_Real64)                     &
     &        *(LSR-MSR)*(LSR+MSR)))*UPLM_SOL(L, J-2)
          ENDDO
        ENDDO
!
      ENDDO
!
!
!
      RETURN
      END SUBROUTINE CALC_UPLM_SOL
#endif
