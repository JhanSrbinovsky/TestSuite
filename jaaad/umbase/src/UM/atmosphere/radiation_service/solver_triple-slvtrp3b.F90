#if defined(A70_1B) || defined(A70_1C)
#if defined(A01_3A) || defined(A02_3A) \
 || defined(A01_3C) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to solve for mixed fluxes scattering without a matrix.
!
! Method:
!       Gaussian elimination in an upward direction is employed to
!       determine effective albedos for lower levels of the atmosphere.
!       This allows a downward pass of back-substitution to be carried
!       out to determine the upward and downward fluxes.
!
!       Following testing, the original code, which was optimised for the
!       NEC, has been put under a VECTOR definition while revised code,
!       optimised for the IBM has been put a NOT VECTOR definition.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE SOLVER_TRIPLE(N_PROFILE, N_LAYER, N_CLOUD_TOP          &
     &   , T, R, S_DOWN, S_UP                                           &
     &   , T_STRAT, R_STRAT, S_DOWN_STRAT, S_UP_STRAT                   &
     &   , T_CONV, R_CONV, S_DOWN_CONV, S_UP_CONV                       &
     &   , V11, V12, V13, V21, V22, V23, V31, V32, V33                  &
     &   , U11, U12, U13, U21, U22, U23, U31, U32, U33                  &
     &   , L_NET                                                        &
     &   , FLUX_INC_DOWN                                                &
     &   , SOURCE_GROUND_FREE, SOURCE_GROUND_STRAT                      &
     &   , SOURCE_GROUND_CONV, ALBEDO_SURFACE_DIFF                      &
     &   , FLUX_TOTAL                                                   &
     &   , NPD_PROFILE, NPD_LAYER                                       &
     &   )
 
      IMPLICIT NONE
 
!     SIZES OF DUMMY ARRAYS.
      INTEGER, INTENT(IN) :: NPD_PROFILE    ! MAXIMUM NUMBER OF PROFILES
      INTEGER, INTENT(IN) :: NPD_LAYER      ! MAXIMUM NUMBER OF LAYERS  
 
!     DUMMY ARGUMENTS.
      INTEGER, INTENT(IN) :: N_PROFILE      ! NUMBER OF PROFILES
      INTEGER, INTENT(IN) :: N_LAYER        ! NUMBER OF LAYERS  
      INTEGER, INTENT(IN) :: N_CLOUD_TOP    ! TOPMOST CLOUDY LAYER
      LOGICAL, INTENT(IN) :: L_NET   !FLAG FOR CALCULATION OF NET FLUXES
 
      REAL                                                              &
                !, INTENT(IN)
     &     T(NPD_PROFILE, NPD_LAYER)                                    &
!             CLEAR-SKY TRANSMISSION
     &   , R(NPD_PROFILE, NPD_LAYER)                                    &
!             CLEAR-SKY REFLECTION
     &   , S_DOWN(NPD_PROFILE, NPD_LAYER)                               &
!             CLEAR-SKY DOWNWARD SOURCE FUNCTION
     &   , S_UP(NPD_PROFILE, NPD_LAYER)                                 &
!             CLEAR-SKY UPWARD SOURCE FUNCTION
     &   , T_STRAT(NPD_PROFILE, NPD_LAYER)                              &
!             STRATFIFORM TRANSMISSION
     &   , R_STRAT(NPD_PROFILE, NPD_LAYER)                              &
!             STRATFIFORM REFLECTION
     &   , S_DOWN_STRAT(NPD_PROFILE, NPD_LAYER)                         &
!             DOWNWARD STRATFIFORM SOURCE FUNCTION
     &   , S_UP_STRAT(NPD_PROFILE, NPD_LAYER)                           &
!             UPWARD STRATFIFORM SOURCE FUNCTION
     &   , T_CONV(NPD_PROFILE, NPD_LAYER)                               &
!             CONVECTIVE TRANSMISSION
     &   , R_CONV(NPD_PROFILE, NPD_LAYER)                               &
!             CONVECTIVE REFLECTION
     &   , S_DOWN_CONV(NPD_PROFILE, NPD_LAYER)                          &
!             DOWNWARD CONVECTIVE SOURCE FUNCTION
     &   , S_UP_CONV(NPD_PROFILE, NPD_LAYER)
!             UPWARD CONVECTIVE SOURCE FUNCTION

      REAL                                                              &
                !, INTENT(IN)
     &     V11(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V12(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V13(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V21(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V22(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V23(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V31(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V32(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
     &   , V33(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR DOWNWARD RADIATION
      REAL                                                              &
     &     U11(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U12(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U13(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U21(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U22(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U23(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U31(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U32(NPD_PROFILE, 0: NPD_LAYER)                               &
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
     &   , U33(NPD_PROFILE, 0: NPD_LAYER)
!             ENERGY TRANSFER COEFFICIENT FOR UPWARD RADIATION
      REAL                                                              &
                !, INTENT(IN)
     &     FLUX_INC_DOWN(NPD_PROFILE)                                   &
!             INCIDENT FLUX
     &   , SOURCE_GROUND_FREE(NPD_PROFILE)                              &
!             SOURCE FROM GROUND (CLEAR SKY)
     &   , SOURCE_GROUND_STRAT(NPD_PROFILE)                             &
!             SOURCE FROM GROUND (CLOUDY REGION)
     &   , SOURCE_GROUND_CONV(NPD_PROFILE)                              &
!             SOURCE FROM GROUND (CLOUDY REGION)
     &   , ALBEDO_SURFACE_DIFF(NPD_PROFILE)
!             DIFFUSE ALBEDO

      REAL, INTENT(OUT) :: FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
                         ! TOTAL FLUX
 
!     LOCAL VARIABLES.
      INTEGER :: I    ! LOOP VARIABLE
      INTEGER :: L    ! LOOP VARIABLE
#if defined(VECTOR)
 
!     EFFECTIVE COUPLING ALBEDOS AND SOURCE FUNCTIONS:
      REAL                                                              &
     &     ALPHA(N_PROFILE,3,3,NPD_LAYER+1),G(N_PROFILE,3,NPD_LAYER+1)

!     TERMS FOR DOWNWARD PROPAGATION:
      REAL                                                              &
     &     GAMMA(N_PROFILE,3,3,NPD_LAYER), H(N_PROFILE,3,NPD_LAYER)     &
     &   , BETA(N_PROFILE,3,3,NPD_LAYER)

!     AUXILIARY NUMERICAL VARIABLES REQUIRED ONLY IN THE CURRENT LAYER:
      REAL                                                              &
     &     THETA(N_PROFILE,3,3,0:NPD_LAYER)                             &
     &   , LAMBDA(N_PROFILE,0:3,0:NPD_LAYER)
 
!     TEMPORARY FLUXES
      REAL :: FLUX_DOWN(3,0: NPD_LAYER)
      REAL :: FLUX_UP(3,0: NPD_LAYER)
      REAL :: FLUX_TEMP(N_PROFILE,6,0:NPD_LAYER)

      INTEGER  :: n_UV_vars    ! number of variables in UV array
      INTEGER  :: n_T_R_vars   ! number of variables in T_R array

      PARAMETER                                                         &
     & (n_UV_vars  = 18,                                                &
     &  n_T_R_vars = 12)

!     New temporary arrays for the optimised version
      REAL                                                              &
     &     UV(n_UV_vars,0:N_LAYER,NPD_PROFILE)                          &
     &   , T_R(n_T_R_vars,N_LAYER,NPD_PROFILE)

!     THIS ROUTINE IS SPECIFIC TO CASES OF THREE REGIONS AND IT IS
!     ASSUMED THAT 1 REPRESENTS CLEAR SKIES, 2 REPRESENTS STRATIFORM
!     CLOUDS AND 3 REPRESENTS CONVECTIVE CLOUD.

! Write U and V into a scalar array

      DO I = 0,N_LAYER
!*DIR$ CACHE_BYPASS UV,V11,V21,V31,V12
         DO L=1,N_PROFILE
            UV(1,I,L) = V11(L,I)
            UV(2,I,L) = V21(L,I)
            UV(3,I,L) = V31(L,I)
            UV(4,I,L) = V12(L,I)
         END DO
      END DO
      DO I=0,N_LAYER

!*DIR$ CACHE_BYPASS UV,V22,V32,V13,V23
         DO L=1,N_PROFILE
            UV(5,I,L) = V22(L,I)
            UV(6,I,L) = V32(L,I)
            UV(7,I,L) = V13(L,I)
            UV(8,I,L) = V23(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,U11,U21,U31,V33
         DO L=1,N_PROFILE
            UV(9,I,L) = V33(L,I)
            UV(10,I,L) = U11(L,I)
            UV(11,I,L) = U21(L,I)
            UV(12,I,L) = U31(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,U12,U22,U32,U13
         DO L=1,N_PROFILE
            UV(13,I,L) = U12(L,I)
            UV(14,I,L) = U22(L,I)
            UV(15,I,L) = U32(L,I)
            UV(16,I,L) = U13(L,I)
         END DO
      END DO
      DO I=0,N_LAYER
!*DIR$ CACHE_BYPASS UV,U23,U33
         DO L=1,N_PROFILE
            UV(17,I,L) = U23(L,I)
            UV(18,I,L) = U33(L,I)
         END DO
      END DO
      DO I=1,N_LAYER
!*DIR$ CACHE_BYPASS T_R,T,T_STRAT,T_CONV,R
         DO L=1,N_PROFILE
            T_R(1,I,L) = T(L,I)
            T_R(2,I,L) = T_STRAT(L,I)
            T_R(3,I,L) = T_CONV(L,I)
            T_R(4,I,L) = R(L,I)
         END DO
      END DO
      DO I=1,N_LAYER
!*DIR$ CACHE_BYPASS T_R,R_STRAT,R_CONV,S_DOWN,S_UP
         DO L=1,N_PROFILE
            T_R(5,I,L) = R_STRAT(L,I)
            T_R(6,I,L) = R_CONV(L,I)
            T_R(7,I,L) = S_DOWN(L,I)
            T_R(8,I,L) = S_UP(L,I)
         END DO
      END DO
      DO I=1,N_LAYER
!*DIR$ CACHE_BYPASS T_R,S_DOWN_STRAT,S_UP_STRAT,S_UP_CONV,S_DOWN_CONV
         DO L=1,N_PROFILE
            T_R(9,I,L) = S_DOWN_STRAT(L,I)
            T_R(10,I,L) = S_UP_STRAT(L,I)
            T_R(11,I,L) = S_DOWN_CONV(L,I)
            T_R(12,I,L) = S_UP_CONV(L,I)
         END DO
      END DO
 
!     INITIALIZE AT THE BOTTOM OF THE COLUMN FOR UPWARD ELIMINATION.
      DO L=1, N_PROFILE
         alpha(L,1,1,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         alpha(L,1,2,N_LAYER+1)=0.0E+00
         alpha(L,1,3,N_LAYER+1)=0.0E+00
         alpha(L,2,1,N_LAYER+1)=0.0E+00
         alpha(L,2,2,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         alpha(L,2,3,N_LAYER+1)=0.0E+00
         alpha(L,3,1,N_LAYER+1)=0.0E+00
         alpha(L,3,2,N_LAYER+1)=0.0E+00
         alpha(L,3,3,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         g(L,1,N_LAYER+1)=SOURCE_GROUND_FREE(L)
         g(L,2,N_LAYER+1)=SOURCE_GROUND_STRAT(L)
         g(L,3,N_LAYER+1)=SOURCE_GROUND_CONV(L)
 
!     UPWARD ELIMINATION THROUGH THE CLOUDY LAYERS.

         DO I=N_LAYER, N_CLOUD_TOP, -1
            theta(L,1,1,I)=alpha(L,1,1,I+1)*UV(1,I,L)                   &
     &          +alpha(L,1,2,I+1)*UV(2,I,L)                             &
     &          +alpha(L,1,3,I+1)*UV(3,I,L)
            theta(L,1,2,I)=alpha(L,1,1,I+1)*UV(4,I,L)                   &
     &           +alpha(L,1,2,I+1)*UV(5,I,L)                            &
     &           +alpha(L,1,3,I+1)*UV(6,I,L)
            theta(L,1,3,I)=alpha(L,1,1,I+1)*UV(7,I,L)                   &
     &           +alpha(L,1,2,I+1)*UV(8,I,L)                            &
     &           +alpha(L,1,3,I+1)*UV(9,I,L)
            theta(L,2,1,I)=alpha(L,2,1,I+1)*UV(1,I,L)                   &
     &            +alpha(L,2,2,I+1)*UV(2,I,L)                           &
     &            +alpha(L,2,3,I+1)*UV(3,I,L)
            theta(L,2,2,I)=alpha(L,2,1,I+1)*UV(4,I,L)                   &
     &           +alpha(L,2,2,I+1)*UV(5,I,L)                            &
     &           +alpha(L,2,3,I+1)*UV(6,I,L)
            theta(L,2,3,I)=alpha(L,2,1,I+1)*UV(7,I,L)                   &
     &           +alpha(L,2,2,I+1)*UV(8,I,L)                            &
     &           +alpha(L,2,3,I+1)*UV(9,I,L)
            theta(L,3,1,I)=alpha(L,3,1,I+1)*UV(1,I,L)                   &
     &           +alpha(L,3,2,I+1)*UV(2,I,L)                            &
     &           +alpha(L,3,3,I+1)*UV(3,I,L)
            theta(L,3,2,I)=alpha(L,3,1,I+1)*UV(4,I,L)                   &
     &           +alpha(L,3,2,I+1)*UV(5,I,L)                            &
     &           +alpha(L,3,3,I+1)*UV(6,I,L)
            theta(L,3,3,I)=alpha(L,3,1,I+1)*UV(7,I,L)                   &
     &           +alpha(L,3,2,I+1)*UV(8,I,L)                            &
     &           +alpha(L,3,3,I+1)*UV(9,I,L)

            beta(L,3,1,I)=-theta(L,3,1,I)*T_R(4,I,L)
            beta(L,3,2,I)=-theta(L,3,2,I)*T_R(5,I,L)
            beta(L,3,3,I)=1.0E+00/(1.0E+00-theta(L,3,3,I)*T_R(6,I,L))
            gamma(L,3,1,I)=theta(L,3,1,I)*T_R(1,I,L)
            gamma(L,3,2,I)=theta(L,3,2,I)*T_R(2,I,L)
            gamma(L,3,3,I)=theta(L,3,3,I)*T_R(3,I,L)
            h(L,3,I)=g(L,3,I+1)+theta(L,3,1,I)*T_R(7,I,L)               &
     &         +theta(L,3,2,I)*T_R(9,I,L)                               &
     &         +theta(L,3,3,I)*T_R(11,I,L)
 
            lambda(L,3,I)=theta(L,2,3,I)*T_R(6,I,L)*beta(L,3,3,I)
            beta(L,2,2,I)=1.0E+00                                       &
     &         /(1.0E+00-theta(L,2,2,I)*T_R(5,I,L)                      &
     &           +lambda(L,3,I)*beta(L,3,2,I))
            beta(L,2,1,I)=-theta(L,2,1,I)*T_R(4,I,L)                    &
     &                  +lambda(L,3,I)*beta(L,3,1,I)
            gamma(L,2,1,I)=theta(L,2,1,I)*T_R(1,I,L)                    &
     &           +lambda(L,3,I)*gamma(L,3,1,I)
            gamma(L,2,2,I)=theta(L,2,2,I)*T_R(2,I,L)+                   &
     &           lambda(L,3,I)*gamma(L,3,2,I)
            gamma(L,2,3,I)=theta(L,2,3,I)*T_R(3,I,L)+lambda(L,3,I)*     &
     &           gamma(L,3,3,I)
            h(L,2,I)=g(L,2,I+1)+theta(L,2,1,I)*T_R(7,I,L)               &
     &         +theta(L,2,2,I)*T_R(9,I,L)                               &
     &         +theta(L,2,3,I)*T_R(11,I,L)                              &
     &         +lambda(L,3,I)*h(L,3,I)
 
            lambda(L,3,I)=theta(L,1,3,I)*T_R(6,I,L)*beta(L,3,3,I)
            lambda(L,2,I)=(theta(L,1,2,I)*T_R(5,I,L)                    &
     &         -lambda(L,3,I)*beta(L,3,2,I))*beta(L,2,2,I)
            beta(L,1,1,I)=1.0E+00                                       &
     &       /(1.0E+00-theta(L,1,1,I)*T_R(4,I,L)                        &
     &      +lambda(L,3,I)*beta(L,3,1,I)                                &
     &         +lambda(L,2,I)*beta(L,2,1,I))
            gamma(L,1,1,I)=theta(L,1,1,I)*T_R(1,I,L)                    &
     &           +lambda(L,3,I)*gamma(L,3,1,I)                          &
     &          +lambda(L,2,I)*gamma(L,2,1,I)
            gamma(L,1,2,I)=theta(L,1,2,I)*T_R(2,I,L)                    &
     &           +lambda(L,3,I)*gamma(L,3,2,I)                          &
     &           +lambda(L,2,I)*gamma(L,2,2,I)
            gamma(L,1,3,I)=theta(L,1,3,I)*T_R(3,I,L)                    &
     &           +lambda(L,3,I)*gamma(L,3,3,I)                          &
     &           +lambda(L,2,I)*gamma(L,2,3,I)
            h(L,1,I)=g(L,1,I+1)+theta(L,1,1,I)*T_R(7,I,L)               &
     &         +theta(L,1,2,I)*T_R(9,I,L)                               &
     &         +theta(L,1,3,I)*T_R(11,I,L)                              &
     &         +lambda(L,3,I)*h(L,3,I)+lambda(L,2,I)*h(L,2,I)
 
            lambda(L,3,I)=UV(18,I-1,L)*T_R(3,I,L)*beta(L,3,3,I)
            lambda(L,2,I)=(UV(15,I-1,L)*T_R(2,I,L)+lambda(L,3,I)        &
     &         *beta(L,3,2,I))*beta(L,2,2,I)
            lambda(L,1,I)=(UV(12,I-1,L)*T_R(1,I,L)                      &
     &           +lambda(L,3,I)*beta(L,3,1,I)                           &
     &         +lambda(L,2,I)*beta(L,2,1,I))*beta(L,1,1,I)
            alpha(L,3,1,I)=UV(12,I-1,L)*T_R(4,I,L)                      &
     &           +lambda(L,3,I)*gamma(L,3,1,I)                          &
     &           +lambda(L,2,I)*gamma(L,2,1,I)                          &
     &           +lambda(L,1,I)*gamma(L,1,1,I)
            alpha(L,3,2,I)=UV(15,I-1,L)*T_R(5,I,L)                      &
     &         +lambda(L,3,I)*gamma(L,3,2,I)                            &
     &         +lambda(L,2,I)*gamma(L,2,2,I)                            &
     &         +lambda(L,1,I)*gamma(L,1,2,I)
            alpha(L,3,3,I)=UV(18,I-1,L)*T_R(6,I,L)                      &
     &         +lambda(L,3,I)*gamma(L,3,3,I)                            &
     &         +lambda(L,2,I)*gamma(L,2,3,I)                            &
     &         +lambda(L,1,I)*gamma(L,1,3,I)

              g(L,3,I)=UV(12,I-1,L)*T_R(8,I,L)                          &
     &         + UV(15,I-1,L)*T_R(10,I,L)                               &
     &         +UV(18,I-1,L)*T_R(12,I,L)                                &
     &         +lambda(L,3,I)*h(L,3,I)+lambda(L,2,I)*h(L,2,I)           &
     &         +lambda(L,1,I)*h(L,1,I)
 
            lambda(L,3,I)=UV(17,I-1,l)*T_R(3,I,L)*beta(L,3,3,I)
            lambda(L,2,I)=(UV(14,I-1,l)*T_R(2,I,L)                      &
     &           +lambda(L,3,I)*beta(L,3,2,I))*beta(L,2,2,I)
            lambda(L,1,I)=(UV(11,I-1,l)*T_R(1,I,L)                      &
     &           +lambda(L,3,I)*beta(L,3,1,I)                           &
     &           +lambda(L,2,I)*beta(L,2,1,I))*beta(L,1,1,I)
            alpha(L,2,1,I)=UV(11,I-1,L)*T_R(4,I,L)                      &
     &           +lambda(L,3,I)*gamma(L,3,1,I)                          &
     &           +lambda(L,2,I)*gamma(L,2,1,I)                          &
     &           +lambda(L,1,I)*gamma(L,1,1,I)
            alpha(L,2,2,I)=UV(14,I-1,l)*T_R(5,I,L)                      &
     &           +lambda(L,3,I)*gamma(L,3,2,I)                          &
     &           +lambda(L,2,I)*gamma(L,2,2,I)                          &
     &           +lambda(L,1,I)*gamma(L,1,2,I)
            alpha(L,2,3,I)=UV(17,I-1,L)*T_R(6,I,L)                      &
     &         +lambda(L,3,I)*gamma(L,3,3,I)                            &
     &         +lambda(L,2,I)*gamma(L,2,3,I)                            &
     &         +lambda(L,1,I)*gamma(L,1,3,I)
            g(L,2,I)=UV(11,I-1,L)*T_R(8,I,L)                            &
     &           +UV(14,I-1,L)*T_R(10,I,L)                              &
     &           +UV(17,I-1,L)*T_R(12,I,L)                              &
     &           +lambda(L,3,I)*h(L,3,I)+lambda(L,2,I)*h(L,2,I)         &
     &           +lambda(L,1,I)*h(L,1,I)
 
            lambda(L,3,I)=UV(16,I-1,L)*T_R(3,I,L)*beta(L,3,3,I)
            lambda(L,2,I)=(UV(13,I-1,L)*T_R(2,I,L)                      &
     &         +lambda(L,3,I)*beta(L,3,2,I))                            &
     &         *beta(L,2,2,I)
            lambda(L,1,I)=(UV(10,I-1,L)*T_R(1,I,L)                      &
     &           +lambda(L,3,I)*beta(L,3,1,I)                           &
     &           +lambda(L,2,I)*beta(L,2,1,I))*beta(L,1,1,I)
            alpha(L,1,1,I)=UV(10,I-1,L)*T_R(4,I,L)                      &
     &           +lambda(L,3,I)*gamma(L,3,1,I)                          &
     &           +lambda(L,2,I)*gamma(L,2,1,I)                          &
     &           +lambda(L,1,I)*gamma(L,1,1,I)
            alpha(L,1,2,I)=UV(13,I-1,L)*T_R(5,I,L)                      &
     &         +lambda(L,3,I)*gamma(L,3,2,I)                            &
     &         +lambda(L,2,I)*gamma(L,2,2,I)                            &
     &         +lambda(L,1,I)*gamma(L,1,2,I)
            alpha(L,1,3,I)=UV(16,I-1,L)*T_R(6,I,L)                      &
     &         +lambda(L,3,I)*gamma(L,3,3,I)                            &
     &         +lambda(L,2,I)*gamma(L,2,3,I)                            &
     &         +lambda(L,1,I)*gamma(L,1,3,I)
            g(L,1,I)=UV(10,I-1,L)*T_R(8,I,L)                            &
     &           +UV(13,I-1,L)*T_R(10,I,L)                              &
     &           +UV(16,I-1,L)*T_R(12,I,L)                              &
     &           +lambda(L,3,I)*h(L,3,I)                                &
     &           +lambda(L,2,I)*h(L,2,I)+lambda(L,1,I)*h(L,1,I)
 
         ENDDO
 
         IF (N_CLOUD_TOP >  1) THEN
 
!           THE LAYER ABOVE THE CLOUD: ONLY ONE SET OF ALPHAS
!           IS NOW NEEDED.
 
            I=N_CLOUD_TOP-1
 
            IF (N_CLOUD_TOP <  N_LAYER) THEN
!              If all columns are clear down to the surface, the
!              coefficients UV will not be set, so the case when
!              N_CLOUD_TOP equals N_LAYER must be treated separately.
               theta(L,1,1,I)=alpha(L,1,1,I+1)*UV(1,I,L)                &
     &            +alpha(L,1,2,I+1)*UV(2,I,L)                           &
     &            +alpha(L,1,3,I+1)*UV(3,I,L)
            ELSE
               theta(L,1,1,I)=alpha(L,1,1,I+1)
            ENDIF
 
            beta(L,1,1,I)=1.0E+00/(1.0E+00-theta(L,1,1,I)*T_R(4,I,L))
            gamma(L,1,1,I)=theta(L,1,1,I)*T_R(1,I,L)
            h(L,1,I)=g(L,1,I+1)+theta(L,1,1,I)*T_R(7,I,L)
 
            lambda(L,0,I)=T_R(1,I,L)*beta(L,1,1,I)
            alpha(L,1,1,I)=T_R(4,I,L)+lambda(L,0,I)*gamma(L,1,1,I)
            g(L,1,I)=T_R(8,I,L)+lambda(L,0,I)*h(L,1,I)
 
         ENDIF
 
 
         DO I=N_CLOUD_TOP-2, 1, -1
 
            beta(L,1,1,I)=1.0E+00/(1.0E+00-alpha(L,1,1,I+1)*T_R(4,I,L))
            gamma(L,1,1,I)=alpha(L,1,1,I+1)*T_R(1,I,L)
            h(L,1,I)=g(L,1,I+1)+alpha(L,1,1,I+1)*T_R(7,I,L)
 
            lambda(L,1,I)=T_R(1,I,L)*beta(L,1,1,I)
            alpha(L,1,1,I)=T_R(4,I,L)+lambda(L,1,I)*gamma(L,1,1,I)
            g(L,1,I)=T_R(8,I,L)+lambda(L,1,I)*h(L,1,I)
 
         ENDDO

!     Initialize for downward back-substitution. If there is cloud
!     in the top layer the upward flux must be calculated allowing
!     for reflection from all three regions.
      FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
      IF (N_CLOUD_TOP >  1) THEN
         FLUX_TOTAL(L, 1)=g(L,1,1)+alpha(L,1,1,1)*FLUX_INC_DOWN(L)
      ELSE
         FLUX_TOTAL(L, 1)=g(L,1,1)+FLUX_INC_DOWN(L)                     &
     &      *(UV(1, 0, L)*alpha(L,1,1,1)+UV(2, 0, L)*alpha(L,1,2,1)     &
     &      +UV(3, 0, L)*alpha(L,1,3,1))
      ENDIF
 
!     SWEEP DOWNWARD THROUGH THE CLEAR-SKY REGION, FINDING THE DOWNWARD
!     FLUX AT THE TOP OF THE LAYER AND THE UPWARD FLUX AT THE BOTTOM.
      DO I=1, N_CLOUD_TOP-1
            FLUX_TOTAL(L, 2*I+1)=(gamma(L,1,1,I)*FLUX_TOTAL(L, 2*I)     &
     &         +h(L,1,I))*beta(L,1,1,I)
            FLUX_TOTAL(L, 2*I+2)=T_R(1,I,L)*FLUX_TOTAL(L, 2*I)          &
     &        +T_R(4,I,L)*FLUX_TOTAL(L, 2*I+1)+T_R(7,I,L)
      ENDDO
 
!     PASS INTO THE TOP CLOUDY LAYER. USE FLUX_DOWN_[1,2,3] TO HOLD,
!     PROVISIONALLY, THE DOWNWARD FLUXES JUST BELOW THE TOP OF THE
!     LAYER, THEN CALCULATE THE UPWARD FLUXES AT THE BOTTOM AND
!     FINALLY THE DOWNWARD FLUXES AT THE BOTTOM OF THE LAYER.
      I=N_CLOUD_TOP
         flux_temp(L,1,I)=UV(1,I-1,L)*FLUX_TOTAL(L, 2*I)
         flux_temp(L,2,I)=UV(2,I-1,L)*FLUX_TOTAL(L, 2*I)
         flux_temp(L,3,I)=UV(3,I-1,L)*FLUX_TOTAL(L, 2*I)
         flux_temp(L,4,I)=(gamma(L,1,1,I)*flux_temp(L,1,I)              &
     &      +gamma(L,1,2,I)*flux_temp(L,2,I)                            &
     &      +gamma(L,1,3,I)*flux_temp(L,3,I)                            &
     &      +h(L,1,I))*beta(L,1,1,I)
         flux_temp(L,5,I)=(gamma(L,2,1,I)*flux_temp(L,1,I)              &
     &      +gamma(L,2,2,I)*flux_temp(L,2,I)                            &
     &      +gamma(L,2,3,I)*flux_temp(L,3,I)+h(L,2,I)                   &
     &      -beta(L,2,1,I)*flux_temp(L,4,I))*beta(L,2,2,I)
         flux_temp(L,6,I)=(gamma(L,3,1,I)*flux_temp(L,1,I)              &
     &      +gamma(L,3,2,I)*flux_temp(L,2,I)                            &
     &      +gamma(L,3,3,I)*flux_temp(L,3,I)+h(L,3,I)                   &
     &      -beta(L,3,1,I)*flux_temp(L,4,I)-beta(L,3,2,I)               &
     &      *flux_temp(L,5,I))                                          &
     &      *beta(L,3,3,I)
         flux_temp(L,1,I)=T_R(1,I,L)*flux_temp(L,1,I)                   &
     &      +T_R(4,I,L)*flux_temp(L,4,I)+T_R(7,I,L)
         flux_temp(L,2,I)=T_R(2,I,L)*flux_temp(L,2,I)                   &
     &      +T_R(5,I,L)*flux_temp(L,5,I)+T_R(9,I,L)
         flux_temp(L,3,I)=T_R(3,I,L)*flux_temp(L,3,I)                   &
     &      +T_R(6,I,L)*flux_temp(L,6,I)+T_R(11,I,L)
 
!     THE MAIN LOOP OF BACK-SUBSTITUTION. THE PROVISIONAL USE OF THE
!     DOWNWARD FLUXES IS AS ABOVE.
         DO I=N_CLOUD_TOP+1, N_LAYER
            flux_temp(L,1,I)=UV(1,I-1,L)*flux_temp(L,1,I-1)             &
     &           +UV(4,I-1,L)*flux_temp(L,2,I-1)                        &
     &           +UV(7,I-1,L)*flux_temp(L,3,I-1)
            flux_temp(L,2,I)=UV(2,I-1,L)*flux_temp(L,1,I-1)             &
     &           +UV(5,I-1,L)*flux_temp(L,2,I-1)                        &
     &           +UV(8,I-1,L)*flux_temp(L,3,I-1)
            flux_temp(L,3,I)=UV(3,I-1,l)*flux_temp(L,1,I-1)             &
     &           +UV(6,I-1,L)*flux_temp(L,2,I-1)                        &
     &           +UV(9,I-1,L)*flux_temp(L,3,I-1)
            flux_temp(L,4,I)=(gamma(L,1,1,I)*flux_temp(L,1,I)           &
     &         +gamma(L,1,2,I)*flux_temp(L,2,I)                         &
     &         +gamma(L,1,3,I)*flux_temp(L,3,I)                         &
     &         +h(L,1,I))*beta(L,1,1,I)
            flux_temp(L,5,I)=(gamma(L,2,1,I)*flux_temp(L,1,I)           &
     &         +gamma(L,2,2,I)*flux_temp(L,2,I)                         &
     &         +gamma(L,2,3,I)*flux_temp(L,3,I)+h(L,2,I)                &
     &         -beta(L,2,1,I)*flux_temp(L,4,I))*beta(L,2,2,I)
            flux_temp(L,6,I)=(gamma(L,3,1,I)*flux_temp(L,1,I)           &
     &         +gamma(L,3,2,I)*flux_temp(L,2,I)                         &
     &         +gamma(L,3,3,I)*flux_temp(L,3,I)+h(L,3,I)                &
     &         -beta(L,3,1,I)*flux_temp(L,4,I)                          &
     &         -beta(L,3,2,I)*flux_temp(L,5,I))                         &
     &         *beta(L,3,3,I)
            flux_temp(L,1,I)=T_R(1,I,L)*flux_temp(L,1,I)                &
     &         +T_R(4,I,L)*flux_temp(L,4,I)+T_R(7,I,L)
            flux_temp(L,2,I)=T_R(2,I,L)*flux_temp(L,2,I)                &
     &         +T_R(5,I,L)*flux_temp(L,5,I)+T_R(9,I,L)
            flux_temp(L,3,I)=T_R(3,I,L)*flux_temp(L,3,I)                &
     &         +T_R(6,I,L)*flux_temp(L,6,I)+T_R(11,I,L)
         ENDDO

!     CALCULATE THE OVERALL FLUX.
         DO I=N_CLOUD_TOP, N_LAYER
            FLUX_TOTAL(L, 2*I+1)=flux_temp(L,4,I)+flux_temp(L,5,I)      &
     &         +flux_temp(L,6,I)
            FLUX_TOTAL(L, 2*I+2)=flux_temp(L,1,I)+flux_temp(L,2,I)      &
     &         +flux_temp(L,3,I)
         ENDDO

      END DO
#else

!     EFFECTIVE COUPLING ALBEDOS AND SOURCE FUNCTIONS:
      REAL :: ALPHA11(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA12(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA13(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA21(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA22(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA23(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA31(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA32(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: ALPHA33(NPD_PROFILE, NPD_LAYER+1) 
      REAL :: G1(NPD_PROFILE, NPD_LAYER+1)
      REAL :: G2(NPD_PROFILE, NPD_LAYER+1)
      REAL :: G3(NPD_PROFILE, NPD_LAYER+1)

!     TERMS FOR DOWNWARD PROPAGATION:
      REAL :: GAMMA11(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA12(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA13(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA21(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA22(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA23(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA31(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA32(NPD_PROFILE, NPD_LAYER)
      REAL :: GAMMA33(NPD_PROFILE, NPD_LAYER)
      REAL :: BETA11(NPD_PROFILE, NPD_LAYER)
      REAL :: BETA21(NPD_PROFILE, NPD_LAYER)
      REAL :: BETA22(NPD_PROFILE, NPD_LAYER)
      REAL :: BETA31(NPD_PROFILE, NPD_LAYER)
      REAL :: BETA32(NPD_PROFILE, NPD_LAYER)
      REAL :: BETA33(NPD_PROFILE, NPD_LAYER)
      REAL :: H1(NPD_PROFILE, NPD_LAYER)
      REAL :: H2(NPD_PROFILE, NPD_LAYER)
      REAL :: H3(NPD_PROFILE, NPD_LAYER)
 
!     AUXILAIRY NUMERICAL VARIABLES REQUIRED ONLY IN THE CURRENT LAYER:
      REAL :: THETA11
      REAL :: THETA12
      REAL :: THETA13
      REAL :: THETA21
      REAL :: THETA22
      REAL :: THETA23
      REAL :: THETA31
      REAL :: THETA32
      REAL :: THETA33
      REAL :: LAMBDA3
      REAL :: LAMBDA2
      REAL :: LAMBDA1
      REAL :: LAMBDA
 
!     TEMPORARY FLUXES
      REAL                                                              &
     &     FLUX_DOWN_1(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DOWNWARD FLUXES OUTSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_2(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_DOWN_3(NPD_PROFILE, 0: NPD_LAYER)                       &
!             DOWNWARD FLUXES INSIDE CLOUDS JUST BELOW I'TH LEVEL
     &   , FLUX_UP_1(NPD_PROFILE, 0: NPD_LAYER)                         &
!             UPWARD FLUXES OUTSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_2(NPD_PROFILE, 0: NPD_LAYER)                         &
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
     &   , FLUX_UP_3(NPD_PROFILE, 0: NPD_LAYER)
!             UPWARD FLUXES INSIDE CLOUDS JUST ABOVE I'TH LEVEL
         
        
!     THIS ROUTINE IS SPECIFIC TO CASES OF THREE REGIONS AND IT IS      
!     ASSUMED THAT 1 REPRESENTS CLEAR SKIES, 2 REPRESENTS STRATIFORM    
!     CLOUDS AND 3 REPRESENTS CONVECTIVE CLOUD.        
        
!     INITIALIZE AT THE BOTTOM OF THE COLUMN FOR UPWARD ELIMINATION.    
      DO L=1, N_PROFILE
         ALPHA11(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA12(L, N_LAYER+1)=0.0E+00
         ALPHA13(L, N_LAYER+1)=0.0E+00
         ALPHA21(L, N_LAYER+1)=0.0E+00
         ALPHA22(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         ALPHA23(L, N_LAYER+1)=0.0E+00
         ALPHA31(L, N_LAYER+1)=0.0E+00
         ALPHA32(L, N_LAYER+1)=0.0E+00
         ALPHA33(L, N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
         G1(L, N_LAYER+1)=SOURCE_GROUND_FREE(L)
         G2(L, N_LAYER+1)=SOURCE_GROUND_STRAT(L)
         G3(L, N_LAYER+1)=SOURCE_GROUND_CONV(L)
      ENDDO

!     UPWARD ELIMINATION THROUGH THE CLOUDY LAYERS.
      DO I=N_LAYER, N_CLOUD_TOP, -1        
         DO L=1, N_PROFILE
            THETA11=ALPHA11(L, I+1)*V11(L, I)+ALPHA12(L, I+1)*V21(L, I) &
     &         +ALPHA13(L, I+1)*V31(L, I)
            THETA12=ALPHA11(L, I+1)*V12(L, I)+ALPHA12(L, I+1)*V22(L, I) &
     &         +ALPHA13(L, I+1)*V32(L, I)
            THETA13=ALPHA11(L, I+1)*V13(L, I)+ALPHA12(L, I+1)*V23(L, I) &
     &         +ALPHA13(L, I+1)*V33(L, I)
            THETA21=ALPHA21(L, I+1)*V11(L, I)+ALPHA22(L, I+1)*V21(L, I) &
     &         +ALPHA23(L, I+1)*V31(L, I)
            THETA22=ALPHA21(L, I+1)*V12(L, I)+ALPHA22(L, I+1)*V22(L, I) &
     &         +ALPHA23(L, I+1)*V32(L, I)
            THETA23=ALPHA21(L, I+1)*V13(L, I)+ALPHA22(L, I+1)*V23(L, I) &
     &         +ALPHA23(L, I+1)*V33(L, I)
            THETA31=ALPHA31(L, I+1)*V11(L, I)+ALPHA32(L, I+1)*V21(L, I) &
     &         +ALPHA33(L, I+1)*V31(L, I)
            THETA32=ALPHA31(L, I+1)*V12(L, I)+ALPHA32(L, I+1)*V22(L, I) &
     &         +ALPHA33(L, I+1)*V32(L, I)
            THETA33=ALPHA31(L, I+1)*V13(L, I)+ALPHA32(L, I+1)*V23(L, I) &
     &         +ALPHA33(L, I+1)*V33(L, I)

            BETA31(L, I)=-THETA31*R(L, I)
            BETA32(L, I)=-THETA32*R_STRAT(L, I)
            BETA33(L, I)=1.0E+00/(1.0E+00-THETA33*R_CONV(L, I))
            GAMMA31(L, I)=THETA31*T(L, I)
            GAMMA32(L, I)=THETA32*T_STRAT(L, I)
            GAMMA33(L, I)=THETA33*T_CONV(L, I)
            H3(L, I)=G3(L, I+1)+THETA31*S_DOWN(L, I)                    &
     &         +THETA32*S_DOWN_STRAT(L, I)                              &
     &         +THETA33*S_DOWN_CONV(L, I)
 
            LAMBDA3=THETA23*R_CONV(L, I)*BETA33(L, I)
            BETA22(L, I)=1.0E+00                                        &
     &         /(1.0E+00-THETA22*R_STRAT(L, I)+LAMBDA3*BETA32(L, I))
            BETA21(L, I)=-THETA21*R(L, I)+LAMBDA3*BETA31(L, I)
            GAMMA21(L, I)=THETA21*T(L, I)+LAMBDA3*GAMMA31(L, I)
            GAMMA22(L, I)=THETA22*T_STRAT(L, I)+LAMBDA3*GAMMA32(L, I)
            GAMMA23(L, I)=THETA23*T_CONV(L, I)+LAMBDA3*GAMMA33(L, I)
            H2(L, I)=G2(L, I+1)+THETA21*S_DOWN(L, I)                    &
     &         +THETA22*S_DOWN_STRAT(L, I)+THETA23*S_DOWN_CONV(L, I)    &
     &         +LAMBDA3*H3(L, I)
 
            LAMBDA3=THETA13*R_CONV(L, I)*BETA33(L, I)
            LAMBDA2=(THETA12*R_STRAT(L, I)-LAMBDA3*BETA32(L, I))        &
     &         *BETA22(L, I)
            BETA11(L, I)=1.0E+00                                        &
     &         /(1.0E+00-THETA11*R(L, I)+LAMBDA3*BETA31(L, I)           &
     &         +LAMBDA2*BETA21(L, I))
            GAMMA11(L, I)=THETA11*T(L, I)+LAMBDA3*GAMMA31(L, I)         &
     &         +LAMBDA2*GAMMA21(L, I)
            GAMMA12(L, I)=THETA12*T_STRAT(L, I)+LAMBDA3*GAMMA32(L, I)   &
     &         +LAMBDA2*GAMMA22(L, I)
            GAMMA13(L, I)=THETA13*T_CONV(L, I)+LAMBDA3*GAMMA33(L, I)    &
     &         +LAMBDA2*GAMMA23(L, I)
            H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)                    &
     &         +THETA12*S_DOWN_STRAT(L, I)+THETA13*S_DOWN_CONV(L, I)    &
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)
 
            LAMBDA3=U33(L, I-1)*T_CONV(L, I)*BETA33(L, I)
            LAMBDA2=(U32(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))    &
     &         *BETA22(L, I)
            LAMBDA1=(U31(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)           &
     &         +LAMBDA2*BETA21(L, I))*BETA11(L, I)
            ALPHA31(L, I)=U31(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)     &
     &         +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
            ALPHA32(L, I)=U32(L, I-1)*R_STRAT(L, I)                     &
     &         +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)             &
     &         +LAMBDA1*GAMMA12(L, I)
            ALPHA33(L, I)=U33(L, I-1)*R_CONV(L, I)                      &
     &         +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)             &
     &         +LAMBDA1*GAMMA13(L, I)
 
            G3(L, I)=U31(L, I-1)*S_UP(L, I)+U32(L, I-1)*S_UP_STRAT(L, I)&
     &         +U33(L, I-1)*S_UP_CONV(L, I)                             &
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
 
            LAMBDA3=U23(L, I-1)*T_CONV(L, I)*BETA33(L, I)
            LAMBDA2=(U22(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))    &
     &         *BETA22(L, I)
            LAMBDA1=(U21(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)           &
     &         +LAMBDA2*BETA21(L, I))*BETA11(L, I)
            ALPHA21(L, I)=U21(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)     &
     &         +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
            ALPHA22(L, I)=U22(L, I-1)*R_STRAT(L, I)                     &
     &         +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)             &
     &         +LAMBDA1*GAMMA12(L, I)
            ALPHA23(L, I)=U23(L, I-1)*R_CONV(L, I)                      &
     &         +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)             &
     &         +LAMBDA1*GAMMA13(L, I)
            G2(L, I)=U21(L, I-1)*S_UP(L, I)+U22(L, I-1)*S_UP_STRAT(L, I)&
     &         +U23(L, I-1)*S_UP_CONV(L, I)                             &
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
 
            LAMBDA3=U13(L, I-1)*T_CONV(L, I)*BETA33(L, I)
            LAMBDA2=(U12(L, I-1)*T_STRAT(L, I)+LAMBDA3*BETA32(L, I))    &
     &         *BETA22(L, I)
            LAMBDA1=(U11(L, I-1)*T(L, I)+LAMBDA3*BETA31(L, I)           &
     &         +LAMBDA2*BETA21(L, I))*BETA11(L, I)
            ALPHA11(L, I)=U11(L, I-1)*R(L, I)+LAMBDA3*GAMMA31(L, I)     &
     &         +LAMBDA2*GAMMA21(L, I)+LAMBDA1*GAMMA11(L, I)
            ALPHA12(L, I)=U12(L, I-1)*R_STRAT(L, I)                     &
     &         +LAMBDA3*GAMMA32(L, I)+LAMBDA2*GAMMA22(L, I)             &
     &         +LAMBDA1*GAMMA12(L, I)
            ALPHA13(L, I)=U13(L, I-1)*R_CONV(L, I)                      &
     &         +LAMBDA3*GAMMA33(L, I)+LAMBDA2*GAMMA23(L, I)             &
     &         +LAMBDA1*GAMMA13(L, I)
            G1(L, I)=U11(L, I-1)*S_UP(L, I)+U12(L, I-1)*S_UP_STRAT(L, I)&
     &         +U13(L, I-1)*S_UP_CONV(L, I)                             &
     &         +LAMBDA3*H3(L, I)+LAMBDA2*H2(L, I)+LAMBDA1*H1(L, I)
 
         ENDDO
      ENDDO
 
      IF (N_CLOUD_TOP > 1) THEN
 
!     THE LAYER ABOVE THE CLOUD: ONLY ONE SET OF ALPHAS IS NOW NEEDED.
 
         I=N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            IF (N_CLOUD_TOP < N_LAYER) THEN
!              If all columns are clear down to the surface, the
!              coefficients UV will not be set, so the case when
!              N_CLOUD_TOP equals N_LAYER must be treated separately.
               THETA11=ALPHA11(L, I+1)*V11(L, I)                        &
     &                +ALPHA12(L, I+1)*V21(L, I)                        &
     &                +ALPHA13(L, I+1)*V31(L, I)
            ELSE
               THETA11=ALPHA11(L, I+1)
            ENDIF
 
            BETA11(L, I)=1.0E+00/(1.0E+00-THETA11*R(L, I))
            GAMMA11(L, I)=THETA11*T(L, I)
            H1(L, I)=G1(L, I+1)+THETA11*S_DOWN(L, I)
 
            LAMBDA=T(L, I)*BETA11(L, I)
            ALPHA11(L, I)=R(L, I)+LAMBDA*GAMMA11(L, I)
            G1(L, I)=S_UP(L, I)+LAMBDA*H1(L, I)
         ENDDO

      ENDIF
 
      DO I=N_CLOUD_TOP-2, 1, -1
         DO L=1, N_PROFILE
            BETA11(L, I)=1.0E+00/(1.0E+00-ALPHA11(L, I+1)*R(L, I))
            GAMMA11(L, I)=ALPHA11(L, I+1)*T(L, I)
            H1(L, I)=G1(L, I+1)+ALPHA11(L, I+1)*S_DOWN(L, I)
 
            LAMBDA1=T(L, I)*BETA11(L, I)
            ALPHA11(L, I)=R(L, I)+LAMBDA1*GAMMA11(L, I)
            G1(L, I)=S_UP(L, I)+LAMBDA1*H1(L, I)
         ENDDO
      ENDDO

!     Initialize for downward back-substitution. If there is cloud
!     in the top layer the upward flux must be calculated allowing
!     for reflection from all three regions.
      IF (N_CLOUD_TOP > 1) THEN
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
          FLUX_TOTAL(L, 1)=G1(L, 1)+ALPHA11(L, 1)*FLUX_INC_DOWN(L)
        ENDDO
      ELSE
        DO L=1, N_PROFILE
          FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
          FLUX_TOTAL(L, 1)=G1(L, 1)+FLUX_INC_DOWN(L)                    &
     &       *(V11(L, 0)*ALPHA11(L, 1)+V21(L, 0)*ALPHA12(L, 1)          &
     &       + V31(L, 0)*ALPHA13(L, 1))
        ENDDO
      ENDIF
 
!     SWEEP DOWNWARD THROUGH THE CLEAR-SKY REGION, FINDING THE DOWNWARD
!     FLUX AT THE TOP OF THE LAYER AND THE UPWARD FLUX AT THE BOTTOM.
      DO I=1, N_CLOUD_TOP-1
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=(GAMMA11(L, I)*FLUX_TOTAL(L, 2*I)      &
     &         +H1(L, I))*BETA11(L, I)
            FLUX_TOTAL(L, 2*I+2)=T(L, I)*FLUX_TOTAL(L, 2*I)             &
     &         +R(L, I)*FLUX_TOTAL(L, 2*I+1)+S_DOWN(L, I)
         ENDDO
      ENDDO
 
!     PASS INTO THE TOP CLOUDY LAYER. USE FLUX_DOWN_[1,2,3] TO HOLD,
!     PROVISIONALLY, THE DOWNWARD FLUXES JUST BELOW THE TOP OF THE
!     LAYER, THEN CALCULATE THE UPWARD FLUXES AT THE BOTTOM AND
!     FINALLY THE DOWNWARD FLUXES AT THE BOTTOM OF THE LAYER.
      I=N_CLOUD_TOP
      DO L=1, N_PROFILE
         FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_DOWN_3(L, I)=V31(L, I-1)*FLUX_TOTAL(L, 2*I)
         FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)               &
     &      +GAMMA12(L, I)*FLUX_DOWN_2(L, I)                            &
     &      +GAMMA13(L, I)*FLUX_DOWN_3(L, I)                            &
     &      +H1(L, I))*BETA11(L, I)
         FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)               &
     &      +GAMMA22(L, I)*FLUX_DOWN_2(L, I)                            &
     &      +GAMMA23(L, I)*FLUX_DOWN_3(L, I)+H2(L, I)                   &
     &      -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22(L, I)
         FLUX_UP_3(L, I)=(GAMMA31(L, I)*FLUX_DOWN_1(L, I)               &
     &      +GAMMA32(L, I)*FLUX_DOWN_2(L, I)                            &
     &      +GAMMA33(L, I)*FLUX_DOWN_3(L, I)+H3(L, I)                   &
     &      -BETA31(L, I)*FLUX_UP_1(L, I)-BETA32(L, I)*FLUX_UP_2(L, I)) &
     &      *BETA33(L, I)
         FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                    &
     &      +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
         FLUX_DOWN_2(L, I)=T_STRAT(L, I)*FLUX_DOWN_2(L, I)              &
     &      +R_STRAT(L, I)*FLUX_UP_2(L, I)+S_DOWN_STRAT(L, I)
         FLUX_DOWN_3(L, I)=T_CONV(L, I)*FLUX_DOWN_3(L, I)               &
     &      +R_CONV(L, I)*FLUX_UP_3(L, I)+S_DOWN_CONV(L, I)
      ENDDO
 
!     THE MAIN LOOP OF BACK-SUBSTITUTION. THE PROVISIONAL USE OF THE
!     DOWNWARD FLUXES IS AS ABOVE.
      DO I=N_CLOUD_TOP+1, N_LAYER
         DO L=1, N_PROFILE
            FLUX_DOWN_1(L, I)=V11(L, I-1)*FLUX_DOWN_1(L, I-1)           &
     &         +V12(L, I-1)*FLUX_DOWN_2(L, I-1)                         &
     &         +V13(L, I-1)*FLUX_DOWN_3(L, I-1)
            FLUX_DOWN_2(L, I)=V21(L, I-1)*FLUX_DOWN_1(L, I-1)           &
     &         +V22(L, I-1)*FLUX_DOWN_2(L, I-1)                         &
     &         +V23(L, I-1)*FLUX_DOWN_3(L, I-1)
            FLUX_DOWN_3(L, I)=V31(L, I-1)*FLUX_DOWN_1(L, I-1)           &
     &         +V32(L, I-1)*FLUX_DOWN_2(L, I-1)                         &
     &         +V33(L, I-1)*FLUX_DOWN_3(L, I-1)
            FLUX_UP_1(L, I)=(GAMMA11(L, I)*FLUX_DOWN_1(L, I)            &
     &         +GAMMA12(L, I)*FLUX_DOWN_2(L, I)                         &
     &         +GAMMA13(L, I)*FLUX_DOWN_3(L, I)                         &
     &         +H1(L, I))*BETA11(L, I)
            FLUX_UP_2(L, I)=(GAMMA21(L, I)*FLUX_DOWN_1(L, I)            &
     &         +GAMMA22(L, I)*FLUX_DOWN_2(L, I)                         &
     &         +GAMMA23(L, I)*FLUX_DOWN_3(L, I)+H2(L, I)                &
     &         -BETA21(L, I)*FLUX_UP_1(L, I))*BETA22(L, I)
            FLUX_UP_3(L, I)=(GAMMA31(L, I)*FLUX_DOWN_1(L, I)            &
     &         +GAMMA32(L, I)*FLUX_DOWN_2(L, I)                         &
     &         +GAMMA33(L, I)*FLUX_DOWN_3(L, I)+H3(L, I)                &
     &         -BETA31(L, I)*FLUX_UP_1(L, I)                            &
     &         -BETA32(L, I)*FLUX_UP_2(L, I))                           &
     &         *BETA33(L, I)
            FLUX_DOWN_1(L, I)=T(L, I)*FLUX_DOWN_1(L, I)                 &
     &         +R(L, I)*FLUX_UP_1(L, I)+S_DOWN(L, I)
            FLUX_DOWN_2(L, I)=T_STRAT(L, I)*FLUX_DOWN_2(L, I)           &
     &         +R_STRAT(L, I)*FLUX_UP_2(L, I)+S_DOWN_STRAT(L, I)
            FLUX_DOWN_3(L, I)=T_CONV(L, I)*FLUX_DOWN_3(L, I)            &
     &         +R_CONV(L, I)*FLUX_UP_3(L, I)+S_DOWN_CONV(L, I)
         ENDDO
      ENDDO
 
!     CALCULATE THE OVERALL FLUX.
      DO I=N_CLOUD_TOP, N_LAYER
         DO L=1, N_PROFILE
            FLUX_TOTAL(L, 2*I+1)=FLUX_UP_1(L, I)+FLUX_UP_2(L, I)        &
     &         +FLUX_UP_3(L, I)
            FLUX_TOTAL(L, 2*I+2)=FLUX_DOWN_1(L, I)+FLUX_DOWN_2(L, I)    &
     &         +FLUX_DOWN_3(L, I)
         ENDDO
      ENDDO
#endif
 
!     REDUCE TO NET FLUXES IF REQUIRED.
      IF (L_NET) THEN
         DO I=0, N_LAYER
            DO L=1, N_PROFILE
               FLUX_TOTAL(L, I+1)                                       &
     &            =FLUX_TOTAL(L, 2*I+2)-FLUX_TOTAL(L, 2*I+1)
            ENDDO
         ENDDO
      ENDIF
 
      RETURN
      END SUBROUTINE SOLVER_TRIPLE
#endif
#endif
