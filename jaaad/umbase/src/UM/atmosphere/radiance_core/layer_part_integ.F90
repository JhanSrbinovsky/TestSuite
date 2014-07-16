#if defined(A70_1Z)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set the particular integral for one layer
!
! Purpose:
!   This routine calculates the particular integral in the
!   requested spectral region for the current layer.
!
! Method:
!
!   The solar particular integral is calculated using a recurrence,
!   while the particular integral in the infra-red is calculated
!   from the Planckian terms. A complementary function is added to
!   the naive form of the particular integral to maintain
!   numerical conditioning.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE LAYER_PART_INTEG(                                      &
!                     Basic sizes
     &    N_PROFILE, LS_TRUNC, MS, N_RED_EIGENSYSTEM                    &
!                     Numerical arrays of spherical terms
     &  , CG_COEFF, MU, EIG_VEC, THETA                                  &
!                     Solar variables
     &  , ISOLIR, I_DIRECT_TOP, I_DIRECT_BOTTOM, MU_0, UPLM_SOL         &
!                     Infra-red variables
     &  , DIFF_PLANCK, L_IR_SOURCE_QUAD, DIFF_PLANCK_2                  &
!                     Optical properies
     &  , TAU, SQS2                                                     &
!                     Output variables
     &  , SOURCE_TOP, SOURCE_BOTTOM, UPM_C, K_SOL, Z_SOL, Q_0, Q_1      &
!                     Dimensions
     &  , ND_PROFILE, ND_MAX_ORDER, ND_RED_EIGENSYSTEM                  &
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
     &  , ND_MAX_ORDER                                                  &
!           Size allocated for orders of spherical harmonics
     &  , ND_RED_EIGENSYSTEM
!           Size allocated for the reduced eigensystem
!
!     Include header files.
#include "c_kinds.h"
#include "c_pi.h"
#include "spectral_region_pcf3z.h"
!
!     Dummy arguments
      INTEGER, INTENT(IN) ::                                            &
     &    N_PROFILE                                                     &
!           Number of atmospheric layers
     &  , N_RED_EIGENSYSTEM
!           Size of the reduced eigensystem
      INTEGER, INTENT(IN) ::                                            &
     &    MS                                                            &
!           Azimuthal order
     &  , LS_TRUNC
!           The truncating order of the system of equations
      INTEGER, INTENT(IN) ::                                            &
     &    ISOLIR
!           Flag for spectral region
      REAL  (Real64), INTENT(IN) ::                                     &
     &    CG_COEFF(LS_TRUNC+1-MS)                                       &
!           Clebsch-Gordan coefficients
     &  , MU(ND_PROFILE, ND_RED_EIGENSYSTEM)                            &
!           (Positive) Eigenvalues
     &  , EIG_VEC(ND_PROFILE, 2*ND_RED_EIGENSYSTEM, ND_RED_EIGENSYSTEM) &
!           Eigenvectors of the full systems for positive eigenvalues
!           (these are scaled by the s-coefficients in the routine
!           EIG_SYS)
     &  , THETA(ND_PROFILE, ND_RED_EIGENSYSTEM)
!           Exponentials of optical depths along slant paths defined
!           by the eigenvalues
      REAL  (Real64), INTENT(IN) ::                                     &
     &    TAU(ND_PROFILE)                                               &
!           Optical depths of the layers
     &  , SQS2(ND_PROFILE, 0: ND_MAX_ORDER)
!           S-coefficients
      REAL  (Real64), INTENT(IN) ::                                     &
     &    MU_0(ND_PROFILE)                                              &
!           Cosine of solar zenith angle
     &  , I_DIRECT_TOP(ND_PROFILE)                                      &
!           The direct solar radiance at the top of the current layer
     &  , I_DIRECT_BOTTOM(ND_PROFILE)                                   &
!           The direct solar radiance at the bottom of
!           the current layer
     &  , UPLM_SOL(ND_PROFILE, LS_TRUNC+2-MS)
!           Spherical harmonics of the solar direction
      LOGICAL, INTENT(IN) ::                                            &
     &    L_IR_SOURCE_QUAD
!           Flag for quadratic source function in the IR
      REAL  (Real64), INTENT(IN) ::                                     &
     &    DIFF_PLANCK(ND_PROFILE)                                       &
!           Differences in the hemispheric Planckian FLUX (bottom-top)
!           across the layer
     &  , DIFF_PLANCK_2(ND_PROFILE)
!           Twice the second differences in the hemispheric Planckian
!           FLUX
!
      INTEGER, INTENT(OUT) ::                                           &
     &    K_SOL(ND_PROFILE)
!           Index of eigenvalue used for solar conditioning
      REAL  (Real64), INTENT(OUT) ::                                    &
     &    SOURCE_TOP(ND_PROFILE, LS_TRUNC+1-MS)                         &
!           Source function at the top of the layer
     &  , SOURCE_BOTTOM(ND_PROFILE, LS_TRUNC+1-MS)                      &
!           Source function at the bottom of the layer
     &  , Z_SOL(ND_PROFILE, LS_TRUNC+1-MS)                              &
!           Coefficient of the solar particular integral
     &  , Q_0(ND_PROFILE)                                               &
!           Term for thermal particular integral
     &  , Q_1(ND_PROFILE)                                               &
!           Term for thermal particular integral
     &  , UPM_C(ND_PROFILE, 2*ND_RED_EIGENSYSTEM)
!           Arrays for coefficients of the complementary function
!           used to condition the particular integral
!
!
!     Local variables
      INTEGER                                                           &
     &    LSR                                                           &
!           Reduced polar order
     &  , K                                                             &
!           Loop variable
     &  , L
!           Loop variable
      REAL  (Real64) ::                                                 &
     &    V_SOL(ND_PROFILE, LS_TRUNC+2-MS)                              &
!           Solar particular integral
     &  , V_DIF(ND_PROFILE, LS_TRUNC+2-MS)                              &
!           Difference between particular integral and
!           eigenvector
     &  , GAMMA(ND_PROFILE)                                             &
!           Constant used in the solar particular integral
     &  , X(ND_PROFILE)                                                 &
!           Temporary variable
     &  , Z(ND_PROFILE, LS_TRUNC+1-MS)                                  &
!           Another Temporary variable
     &  , M1LS                                                          &
!           -1 raised to the power l+m
     &  , EIG_SEP(ND_PROFILE)                                           &
!           Separation of the eigenvalue from the cosine of the
!           solar zenith angle
     &  , EIG_SEP_TMP(ND_PROFILE)                                       &
!           Temporary version of the above used in searching
     &  , EIG_DIFF                                                      &
!           Difference between eigenvalue and cosine of zenith
!           angle
     &  , CONST
!
      REAL (Real64) ::                                                  &
     &    NORM1                                                         &
     &  , NORM2                                                         &
     &  , PROD                                                          &
     &  , RES
!
!
      UPM_C=0.0

      IF (ISOLIR == IP_SOLAR) THEN
!
!       If an eigenvalue approaches the cosine of the zenith angle
!       the solution will become ill-conditioned. We reduce this
!       ill conditioning by subtracting a multiple of the
!       eigensolution with the eigenvalue closest to the cosine
!       of the zenith angle.
!
!       Fine the closest eigenvalue.
        DO L=1, N_PROFILE
          K_SOL(L)=1
          EIG_SEP(L)=ABS(MU(L, 1)-MU_0(L))
        ENDDO
        DO K=1, N_RED_EIGENSYSTEM
          DO L=1, N_PROFILE
            EIG_SEP_TMP(L)=ABS(MU(L, K)-MU_0(L))
            IF (EIG_SEP_TMP(L) <  EIG_SEP(L)) THEN
              K_SOL(L)=K
              EIG_SEP(L)=EIG_SEP_TMP(L)
            ENDIF
          ENDDO
        ENDDO
!
!       Determine the particular integral for this layer.
!       Upward recurrence is stable here.
        DO L=1, N_PROFILE

          V_DIF(L, 1)=0.0E+00_Real64

          M1LS=REAL(1-2*MOD(1, 2), Real64)

          V_DIF(L, 2)=( -MU_0(L)*SQS2(L, MS)*V_DIF(L, 1)                &
     &               + SQS2(L, MS)*M1LS*EIG_VEC(L,1,K_SOL(L)) )         &
     &      /CG_COEFF(1)

        ENDDO
!       V_SOL is required one order beyond the truncation to
!       complete the solution.
        DO LSR=3, LS_TRUNC+2-MS
          DO L=1, N_PROFILE

            M1LS=REAL(1-2*MOD(LSR-1, 2), Real64)

            V_DIF(L, LSR)                                               &
     &        =(-MU_0(L)*SQS2(L, LSR+MS-2)*V_DIF(L, LSR-1)              &
     &        -CG_COEFF(LSR-2)*V_DIF(L, LSR-2)                          &
     &        +SQS2(L, LSR+MS-2)*M1LS*EIG_VEC(L,LSR-1,K_SOL(L)))        &
     &        /CG_COEFF(LSR-1)

          ENDDO
        ENDDO
        DO L=1, N_PROFILE

          GAMMA(L)=UPLM_SOL(L, LS_TRUNC+2-MS)                           &
     &            /V_DIF(L, LS_TRUNC+2-MS)

        ENDDO
!
!       Set the solution to remove ill-conditioning. Note that the
!       first element of the eigenvector cannot be zero, since the
!       recurrence would then force all elements to be 0.
!
!        DO L=1, N_PROFILE
!           UPM_C(L, K_SOL(L))=-I_DIRECT_TOP(L)*GAMMA(L)
!     &      *V_SOL(L, 1)/EIG_VEC(L, 1, K_SOL(L))
!           UPM_C(L, K_SOL(L))=0.0

!          UPM_C(L, K_SOL(L))=-I_DIRECT_TOP(L)
!     &                       * (GAMMA(L)*V_SOL(L,1)-UPLM_SOL(L,1))
!     &                       /EIG_VEC(L, 1, K_SOL(L))
!        ENDDO

!          norm1=0.0
!          norm2=0.0
!          prod=0.0
!          do lsr=1,ls_trunc+1-ms
!             norm1=norm1+v_sol(3281,lsr)*v_sol(3281,lsr)*
!     &                  sqs2(l,lsr+ms-1)

!             norm2=norm2+EIG_VEC(3281, lsr, l) *
!     &                   EIG_VEC(3281, lsr, l) *
!     &                   sqs2(l,lsr+ms-1)
!
!             prod=prod+v_sol(3281,lsr)*
!     &                   EIG_VEC(3281, lsr, l)*
!     &                   sqs2(l,lsr+ms-1)*
!     &                   real(1-2*mod(lsr-1, 2))
!
!          enddo
!          norm1=sqrt(norm1)
!          norm2=sqrt(norm2)
!          if (norm1*norm2 == 0.0) then
!            res=10.0
!          else
!            res=prod/(norm1*norm2)
!          endif
!        enddo

!       Calculate the source function at the top and bottom
!       of this layer.
        DO LSR=1, LS_TRUNC+1-MS
          DO L=1, N_PROFILE

            Z_SOL(L, LSR)=I_DIRECT_TOP(L)                               &
     &        *(GAMMA(L)*V_DIF(L, LSR)-UPLM_SOL(L, LSR))

            SOURCE_TOP(L, LSR)=I_DIRECT_TOP(L)                          &
     &        *(GAMMA(L)*V_DIF(L, LSR)-UPLM_SOL(L, LSR))


            IF( EIG_SEP(L) <  1.0E-06) THEN


              EIG_DIFF=TAU(L)/(MU_0(L)*MU(L,K_SOL(L)))

              CONST = - GAMMA(L)*EXP(-TAU(L)/MU_0(L))                   &
     &                *( EIG_DIFF+0.5*EIG_DIFF**2                       &
     &                *(MU(L,K_SOL(L))-MU_0(L)))

            ELSE

              EIG_DIFF=TAU(L)*(1.0/MU(L,K_SOL(L))-1.0/MU_0(L))

              IF (EIG_DIFF <  0.0 ) THEN

                 CONST=GAMMA(L)*EXP(-TAU(L)/MU(L,K_SOL(L)))             &
     &                *(EXP(EIG_DIFF)-1.0)                              &
     &                /(MU(L,K_SOL(L))-MU_0(L))

              ELSE

                 CONST=GAMMA(L)*EXP(-TAU(L)/MU_0(L))                    &
     &                *(1.0-EXP(-EIG_DIFF))                             &
     &                /(MU(L,K_SOL(L))-MU_0(L))

              ENDIF
            ENDIF

            M1LS=REAL(1-2*MOD(LSR-1, 2), Real64)

            SOURCE_BOTTOM(L, LSR)                                       &
     &        = I_DIRECT_BOTTOM(L) * EXP(-TAU(L)/MU_0(L))               &
     &        *( GAMMA(L)*V_DIF(L, LSR)-UPLM_SOL(L, LSR))               &
     &        + I_DIRECT_BOTTOM(L) *EIG_VEC(L, LSR, K_SOL(L))           &
     &        *M1LS*CONST

          ENDDO
        ENDDO
!
! Second method to remove ill-conditioning. Test for the distance
! between the cosine of the zenith angle and the eigenvalue. If
! distance is less than epsilon we look for a particular integral
! of the form TAU*EXP(_TAU/MU) insetad of EXP(_TAU/MU).
!

!        DO L=1,N_PROFILE


!             V_SOL(L,2) =
!     &         (SQS2(L,MS)-1.0)*UPLM_SOL(L,1)/CG_COEFF(2)

!             V_SOL(L,LS_TRUNC-MS) =
!     &         (SQS2(L,LS_TRUNC-MS)-1.0)*UPLM_SOL(L,LS_TRUNC+1-MS)/
!     &          CG_COEFF(LS_TRUNC+1-MS)


! Recursion backwards

!             DO LSR=LS_TRUNC-2-MS,MS+1,-2
!                V_SOL(L,LSR)=
!     &           ( (SQS2(L,LSR)-1.0)*UPLM_SOL(L,LSR+1) -
!     &              CG_COEFF(LSR+2)*V_SOL(L,LSR+2) ) /
!     &              CG_COEFF(LSR+1)
!             ENDDO

! Recursion forwards

!             DO LSR=2,LS_TRUNC-MS,2
!                V_SOL(L,LSR+2)=
!     &           ( (SQS2(L,LSR)-1.0)*UPLM_SOL(L,LSR+1) -
!     &              CG_COEFF(LSR+1)*V_SOL(L,LSR) ) /
!     &              CG_COEFF(LSR+2)
!             ENDDO
!
!             DO LSR=1, LS_TRUNC+1-MS
!               SOURCE_TOP(L, LSR)=0.0
!               SOURCE_BOTTOM(L, LSR)=I_DIRECT_BOTTOM(L)*V_SOL(L,LSR)
!     &                              *TAU(L)*THETA(L, K_SOL(L))
!             ENDDO
!

!           ENDIF
!        ENDDO

      ELSE IF (ISOLIR == IP_INFRA_RED) THEN
!
!       The variation of the Planckian across the layer can be either
!       linear or quadratic in the optical depth. The particular
!       integrals tend to infinity as the optical depth tends to 0, so
!       a particular form of the complementary function must be added
!       to cancel off the singularity; otherwise ill-conditioning will
!       arise. Since the Planckian is azimuthally symmetric only terms
!       with m=0 are affected. Linear variations in the Planckian
!       produce a term in the particular integral with polar order 1.
!       More complicated variations produce terms at higher orders.
!       Note that ill-conditioning has been removed only in the case of
!       linear variations so far as the quadratic case is more
!       complicated. To deal with the case when TAU is 0, we add a
!       tolerance: it is therefore essential that Q_0 should be used
!       consistently to avoid errors in this limit.
!
        IF (MS == 0) THEN
!
          DO L=1, N_PROFILE
            Q_0(L)=SQRT(4.0E+00_Real64/(3.0E+00_Real64*PI))             &
     &        *DIFF_PLANCK(L)/(SQS2(L, 1)*TAU(L)+EPSILON(TAU))
          ENDDO
!
          IF (L_IR_SOURCE_QUAD) THEN
!
            DO L=1, N_PROFILE
              Q_1(L)=2.0E+00_Real64                                     &
     &          *SQRT(4.0E+00_Real64/(3.0E+00_Real64*PI))               &
     &          *DIFF_PLANCK_2(L)/(SQS2(L, 1)*TAU(L)**2+EPSILON(TAU))
              SOURCE_TOP(L, 1)                                          &
     &          =CG_COEFF(1)*Q_1(L)/SQS2(L, 0)
              SOURCE_BOTTOM(L, 1)=SOURCE_TOP(L, 1)
              SOURCE_TOP(L, 2)=Q_0(L)-0.5E+00_Real64*Q_1(L)
              SOURCE_BOTTOM(L, 2)=Q_0(L)+0.5E+00_Real64*Q_1(L)
            ENDDO
            IF (LS_TRUNC >  1) THEN
              DO L=1, N_PROFILE
                SOURCE_TOP(L, 3)=CG_COEFF(2)*Q_1(L)/SQS2(L, 2)
                SOURCE_BOTTOM(L, 3)=SOURCE_TOP(L, 3)
              ENDDO
            ENDIF
!
          ELSE
!
            DO L=1, N_PROFILE
              SOURCE_TOP(L, 1)=0.0E+00_Real64
              SOURCE_BOTTOM(L, 1)=0.0E+00_Real64
              SOURCE_TOP(L, 2)=Q_0(L)
              SOURCE_BOTTOM(L, 2)=SOURCE_TOP(L, 2)
            ENDDO
            IF (LS_TRUNC >  1) THEN
              DO L=1, N_PROFILE
                SOURCE_TOP(L, 3)=0.0E+00_Real64
                SOURCE_BOTTOM(L, 3)=0.0E+00_Real64
              ENDDO
            ENDIF
!
          ENDIF
!
!         Higher orders are unaffected.
          DO LSR=4, LS_TRUNC+1-MS
            DO L=1, N_PROFILE
              SOURCE_TOP(L, LSR)=0.0E+00_Real64
              SOURCE_BOTTOM(L, LSR)=0.0E+00_Real64
            ENDDO
          ENDDO
!
!         Now define the part of the complementary function to
!         restore conditioning.
          DO K=1, N_RED_EIGENSYSTEM
            DO L=1, N_PROFILE
              UPM_C(L, K+N_RED_EIGENSYSTEM)                             &
     &          =-Q_0(L)*SQS2(L, 1)*EIG_VEC(L, 2, K)
              UPM_C(L, K)=-UPM_C(L, K+N_RED_EIGENSYSTEM)
            ENDDO
          ENDDO
!
!         We take advantage of the relationship between the formats
!         of the positive and negative exponentials to reduce the
!         number of operations.
          DO LSR=1, LS_TRUNC+1-MS
            M1LS=REAL(1-2*MOD(LSR-1, 2), Real64)
            DO L=1, N_PROFILE
              X(L)=UPM_C(L, 1+N_RED_EIGENSYSTEM)*(THETA(L, 1)-M1LS)     &
     &          *EIG_VEC(L, LSR, 1)
            ENDDO
            DO K=2, N_RED_EIGENSYSTEM
              DO L=1, N_PROFILE
                X(L)=X(L)                                               &
     &            +UPM_C(L, K+N_RED_EIGENSYSTEM)*(THETA(L, K)-M1LS)     &
     &            *EIG_VEC(L, LSR, K)
              ENDDO
            ENDDO
            DO L=1, N_PROFILE
              SOURCE_TOP(L, LSR)=SOURCE_TOP(L, LSR)+X(L)
              SOURCE_BOTTOM(L, LSR)=SOURCE_BOTTOM(L, LSR)-M1LS*X(L)
            ENDDO
          ENDDO
!
        ELSE
!
!         This code should never be executed as non-zero azimuthal
!         orders are not relevant in the IR, but it is here for
!         safety.
          DO LSR=1, LS_TRUNC+1-MS
            DO L=1, N_PROFILE
              SOURCE_TOP(L, LSR)=0.0E+00_Real64
              SOURCE_BOTTOM(L, LSR)=0.0E+00_Real64
            ENDDO
          ENDDO
        ENDIF
!
      ENDIF
!
!
!
      RETURN
      END SUBROUTINE LAYER_PART_INTEG
#endif
