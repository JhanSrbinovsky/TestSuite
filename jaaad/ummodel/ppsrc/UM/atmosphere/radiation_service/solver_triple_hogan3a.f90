


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
!       This version has been modified by Robin Hogan to fix asymmetry 
!       between stratiform and convective clouds, and to allow 
!       shadowing, as documented in Shonk & Hogan, 2007, J. Climate.
!
! Current Owner of Code: James Manners
!
! Description of Code:
!   FORTRAN 95
!
!- ---------------------------------------------------------------------
  SUBROUTINE SOLVER_TRIPLE_HOGAN(N_PROFILE, N_LAYER, N_CLOUD_TOP    &
     , T, R, S_DOWN, S_UP                                           &
     , T_STRAT, R_STRAT, S_DOWN_STRAT, S_UP_STRAT                   &
     , T_CONV, R_CONV, S_DOWN_CONV, S_UP_CONV                       &
     , V11, V12, V13, V21, V22, V23, V31, V32, V33                  &
     , U11, U12, U13, U21, U22, U23, U31, U32, U33                  &
     , L_NET                                                        &
     , FLUX_INC_DOWN                                                &
     , SOURCE_GROUND_FREE, SOURCE_GROUND_STRAT                      &
     , SOURCE_GROUND_CONV, ALBEDO_SURFACE_DIFF                      &
     , FLUX_TOTAL                                                   &
     , NPD_PROFILE, NPD_LAYER                                       &
     )
  
  
  IMPLICIT NONE


! Sizes of dummy arrays:

! Maximum number of profiles
  INTEGER, INTENT(IN) :: NPD_PROFILE
! Maximum number of layers
  INTEGER, INTENT(IN) :: NPD_LAYER
  
  
! Dummy arguments:
  
! Number of profiles
  INTEGER, INTENT(IN) :: N_PROFILE
! Number of layers
  INTEGER, INTENT(IN) :: N_LAYER
! Topmost cloudy layer
  INTEGER, INTENT(IN) :: N_CLOUD_TOP
  
! Flag for calculation of net fluxes
  LOGICAL, INTENT(IN) :: L_NET
  
! Clear-sky transmission
  REAL, INTENT(IN) :: T(NPD_PROFILE, NPD_LAYER)
! Clear-sky reflection
  REAL, INTENT(IN) :: R(NPD_PROFILE, NPD_LAYER)
! Clear-sky downward source function
  REAL, INTENT(IN) :: S_DOWN(NPD_PROFILE, NPD_LAYER)
! Clear-sky upward source function
  REAL, INTENT(IN) :: S_UP(NPD_PROFILE, NPD_LAYER)
! Stratfiform transmission
  REAL, INTENT(IN) :: T_STRAT(NPD_PROFILE, NPD_LAYER)
! Stratfiform reflection
  REAL, INTENT(IN) :: R_STRAT(NPD_PROFILE, NPD_LAYER)
! Downward stratfiform source function
  REAL, INTENT(IN) :: S_DOWN_STRAT(NPD_PROFILE, NPD_LAYER)
! Upward stratfiform source function
  REAL, INTENT(IN) :: S_UP_STRAT(NPD_PROFILE, NPD_LAYER)
! Convective transmission
  REAL, INTENT(IN) :: T_CONV(NPD_PROFILE, NPD_LAYER)
! Convective reflection
  REAL, INTENT(IN) :: R_CONV(NPD_PROFILE, NPD_LAYER)
! Downward convective source function
  REAL, INTENT(IN) :: S_DOWN_CONV(NPD_PROFILE, NPD_LAYER)
! Upward convective source function
  REAL, INTENT(IN) :: S_UP_CONV(NPD_PROFILE, NPD_LAYER)
  
! Energy transfer coefficients for downward radiation
  REAL, INTENT(IN) :: V11(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V12(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V13(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V21(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V22(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V23(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V31(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V32(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: V33(NPD_PROFILE, 0: NPD_LAYER)
  
! Energy transfer coefficients for upward radiation
  REAL, INTENT(IN) :: U11(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U12(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U13(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U21(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U22(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U23(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U31(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U32(NPD_PROFILE, 0: NPD_LAYER)
  REAL, INTENT(IN) :: U33(NPD_PROFILE, 0: NPD_LAYER)
  
! Incident flux
  REAL, INTENT(IN) :: FLUX_INC_DOWN(NPD_PROFILE)
! Source from ground (clear sky)
  REAL, INTENT(IN) :: SOURCE_GROUND_FREE(NPD_PROFILE)
! Source from ground (cloudy region)
  REAL, INTENT(IN) :: SOURCE_GROUND_STRAT(NPD_PROFILE)
! Source from ground (convective cloud region)
  REAL, INTENT(IN) :: SOURCE_GROUND_CONV(NPD_PROFILE)
! Diffuse albedo
  REAL, INTENT(IN) :: ALBEDO_SURFACE_DIFF(NPD_PROFILE)
  
! Total flux
  REAL, INTENT(OUT) :: FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
  
! Local variables:
  
! Loop variables
  INTEGER :: I, L
  
! Effective coupling albedos and source functions:
  REAL :: alpha(N_PROFILE,3,NPD_LAYER+1)
  REAL :: g(N_PROFILE,3,NPD_LAYER+1)
  
! Terms for downward propagation:
  REAL :: gamma(N_PROFILE,3,NPD_LAYER)
  REAL :: h(N_PROFILE,3,NPD_LAYER)
  REAL :: beta(N_PROFILE,3,NPD_LAYER)
  
! Auxiliary numerical variables required only in the current layer:
  REAL :: theta(N_PROFILE,3,0:NPD_LAYER)
  REAL :: lambda(N_PROFILE,0:3,0:NPD_LAYER)
  
! Temporary fluxes
  REAL :: flux_temp(N_PROFILE,6,0:NPD_LAYER)
  
! Number of variables in UV array
  INTEGER, PARAMETER :: n_UV_vars = 18
! Number of variables in T_R array
  INTEGER, PARAMETER :: n_T_R_vars = 12
  
! New temporary arrays for the optimised version
  REAL :: UV(n_UV_vars,0:N_LAYER,NPD_PROFILE)
  REAL :: T_R(n_T_R_vars,N_LAYER,NPD_PROFILE)


! We do a local transpose of variables, such that the first dimension
! is N_LAYER and the second dimension is N_PROFILE, as this better 
! suits the way the algorithm accesses the data. (Actually a little 
! more complicated than this as all the variables are copied into a 
! single array, where the first dimension is the variable index, and 
! the second and third dimensions are N_LAYER and N_PROFILE)

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

! This routine is specific to cases of three regions and it is
! assumed that 1 represents clear skies, 2 represents stratiform
! clouds and 3 represents convective cloud.

! Initialize at the bottom of the column for upward elimination.
  DO L=1, N_PROFILE
     alpha(L,1,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
     alpha(L,2,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
     alpha(L,3,N_LAYER+1)=ALBEDO_SURFACE_DIFF(L)
     g(L,1,N_LAYER+1)=SOURCE_GROUND_FREE(L)
     g(L,2,N_LAYER+1)=SOURCE_GROUND_STRAT(L)
     g(L,3,N_LAYER+1)=SOURCE_GROUND_CONV(L)

! Upward elimination through the cloudy layers.
     DO I=N_LAYER, N_CLOUD_TOP, -1

        theta(L,1,I)=alpha(L,1,I+1)*UV(1,I,L)                       &
           +alpha(L,2,I+1)*UV(2,I,L)                                &
           +alpha(L,3,I+1)*UV(3,I,L)
        theta(L,2,I)=alpha(L,1,I+1)*UV(4,I,L)                       &
           +alpha(L,2,I+1)*UV(5,I,L)                                &
           +alpha(L,3,I+1)*UV(6,I,L)
        theta(L,3,I)=alpha(L,1,I+1)*UV(7,I,L)                       &
           +alpha(L,2,I+1)*UV(8,I,L)                                &
           +alpha(L,3,I+1)*UV(9,I,L)

        beta(L,1,I)=1.0e+00/(1.0e+00-theta(L,1,I)*T_R(4,I,L))
        beta(L,2,I)=1.0e+00/(1.0e+00-theta(L,2,I)*T_R(5,I,L))
        beta(L,3,I)=1.0e+00/(1.0e+00-theta(L,3,I)*T_R(6,I,L))
        gamma(L,1,I)=theta(L,1,I)*T_R(1,I,L)
        gamma(L,2,I)=theta(L,2,I)*T_R(2,I,L)
        gamma(L,3,I)=theta(L,3,I)*T_R(3,I,L)
        h(L,1,I)=g(L,1,I+1)+theta(L,1,I)*T_R(7,I,L)
        h(L,2,I)=g(L,2,I+1)+theta(L,2,I)*T_R(9,I,L)
        h(L,3,I)=g(L,3,I+1)+theta(L,3,I)*T_R(11,I,L)

        lambda(L,1,I)=T_R(8,I,L) +h(L,1,I)*T_R(1,I,L)*beta(L,1,I)
        lambda(L,2,I)=T_R(10,I,L)+h(L,2,I)*T_R(2,I,L)*beta(L,2,I)
        lambda(L,3,I)=T_R(12,I,L)+h(L,3,I)*T_R(3,I,L)*beta(L,3,I)

        alpha(L,1,I)=T_R(4,I,L)                                     &
           +theta(L,1,I)*T_R(1,I,L)*T_R(1,I,L)*beta(L,1,I)
        alpha(L,2,I)=T_R(5,I,L)                                     &
           +theta(L,2,I)*T_R(2,I,L)*T_R(2,I,L)*beta(L,2,I)
        alpha(L,3,I)=T_R(6,I,L)                                     &
           +theta(L,3,I)*T_R(3,I,L)*T_R(3,I,L)*beta(L,3,I)

        g(L,1,I)=UV(10,I-1,L)*lambda(L,1,I)                         &
           +UV(13,I-1,L)*lambda(L,2,I)                              &
           +UV(16,I-1,L)*lambda(L,3,I)
        g(L,2,I)=UV(11,I-1,L)*lambda(L,1,I)                         &
           +UV(14,I-1,L)*lambda(L,2,I)                              &
           +UV(17,I-1,L)*lambda(L,3,I)
        g(L,3,I)=UV(12,I-1,L)*lambda(L,1,I)                         &
           +UV(15,I-1,L)*lambda(L,2,I)                              &
           +UV(18,I-1,L)*lambda(L,3,I)

     ENDDO


     IF (N_CLOUD_TOP >  1) THEN

! The layer above the cloud: only one set of alphas is now needed.

        I=N_CLOUD_TOP-1

        IF (N_CLOUD_TOP <  N_LAYER) THEN
! If all columns are clear down to the surface, the coefficients UV
! will not be set, so the case when N_CLOUD_TOP equals N_LAYER must 
! be treated separately.
           theta(L,1,I)=alpha(L,1,I+1)*UV(1,I,L)                    &
              +alpha(L,2,I+1)*UV(2,I,L)                             &
              +alpha(L,3,I+1)*UV(3,I,L)
        ELSE
           theta(L,1,I)=alpha(L,1,I+1)
        ENDIF

        beta(L,1,I)=1.0e+00/(1.0e+00-theta(L,1,I)*T_R(4,I,L))
        gamma(L,1,I)=theta(L,1,I)*T_R(1,I,L)
        h(L,1,I)=g(L,1,I+1)+theta(L,1,I)*T_R(7,I,L)

        lambda(L,0,I)=T_R(1,I,L)*beta(L,1,I)
        alpha(L,1,I)=T_R(4,I,L)+lambda(L,0,I)*gamma(L,1,I)
        g(L,1,I)=T_R(8,I,L)+lambda(L,0,I)*h(L,1,I)

     ENDIF


     DO I=N_CLOUD_TOP-2, 1, -1

        beta(L,1,I)=1.0e+00/(1.0e+00-alpha(L,1,I+1)*T_R(4,I,L))
        gamma(L,1,I)=alpha(L,1,I+1)*T_R(1,I,L)
        h(L,1,I)=g(L,1,I+1)+alpha(L,1,I+1)*T_R(7,I,L)

        lambda(L,1,I)=T_R(1,I,L)*beta(L,1,I)
        alpha(L,1,I)=T_R(4,I,L)+lambda(L,1,I)*gamma(L,1,I)
        g(L,1,I)=T_R(8,I,L)+lambda(L,1,I)*h(L,1,I)

     ENDDO


! Initialize for downward back-substitution. If there is cloud
! in the top layer the upward flux must be calculated allowing
! for reflection from all three regions.
     FLUX_TOTAL(L, 2)=FLUX_INC_DOWN(L)
     IF (N_CLOUD_TOP >  1) THEN
        FLUX_TOTAL(L, 1)=g(L,1,1)+alpha(L,1,1)*FLUX_INC_DOWN(L)
     ELSE
        FLUX_TOTAL(L, 1)=g(L,1,1)+FLUX_INC_DOWN(L)                  &
           *(UV(1,0,L)*alpha(L,1,1)+UV(2,0,L)*alpha(L,2,1)          &
           + UV(3,0,L)*alpha(L,3,1))
     ENDIF

! Sweep downward through the clear-sky region, finding the downward
! flux at the top of the layer and the upward flux at the bottom.
     DO I=1, N_CLOUD_TOP-1
        FLUX_TOTAL(L, 2*I+1)=(gamma(L,1,I)*FLUX_TOTAL(L, 2*I)       &
           +h(L,1,I))*beta(L,1,I)
        FLUX_TOTAL(L, 2*I+2)=T_R(1,I,L)*FLUX_TOTAL(L, 2*I)          &
           +T_R(4,I,L)*FLUX_TOTAL(L, 2*I+1)+T_R(7,I,L)
     ENDDO

! Pass into the top cloudy layer. Use flux_temp(1,2,3) to hold,
! provisionally, the downward fluxes just below the top of the layer,
! then calculate the upward fluxes at the bottom, flux_temp(4,5,6),
! and finally the downward fluxes at the bottom of the layer, again
! in flux_temp(1,2,3).
     I=N_CLOUD_TOP
     flux_temp(L,1,I)=UV(1,I-1,L)*FLUX_TOTAL(L, 2*I)
     flux_temp(L,2,I)=UV(2,I-1,L)*FLUX_TOTAL(L, 2*I)
     flux_temp(L,3,I)=UV(3,I-1,L)*FLUX_TOTAL(L, 2*I)
     flux_temp(L,4,I)=(gamma(L,1,I)*flux_temp(L,1,I)                &
        +h(L,1,I))*beta(L,1,I)
     flux_temp(L,5,I)=(gamma(L,2,I)*flux_temp(L,2,I)                &
        +h(L,2,I))*beta(L,2,I)
     flux_temp(L,6,I)=(gamma(L,3,I)*flux_temp(L,3,I)                &
        +h(L,3,I))*beta(L,3,I)
     flux_temp(L,1,I)=T_R(1,I,L)*flux_temp(L,1,I)                   &
        +T_R(4,I,L)*flux_temp(L,4,I)+T_R(7,I,L)
     flux_temp(L,2,I)=T_R(2,I,L)*flux_temp(L,2,I)                   &
        +T_R(5,I,L)*flux_temp(L,5,I)+T_R(9,I,L)
     flux_temp(L,3,I)=T_R(3,I,L)*flux_temp(L,3,I)                   &
        +T_R(6,I,L)*flux_temp(L,6,I)+T_R(11,I,L)

! The main loop of back-substitution. The provisional use of
! flux_temp(1,2,3) is as above.
     DO I=N_CLOUD_TOP+1, N_LAYER
        flux_temp(L,1,I)=UV(1,I-1,L)*flux_temp(L,1,I-1)             &
           +UV(4,I-1,L)*flux_temp(L,2,I-1)                          &
           +UV(7,I-1,L)*flux_temp(L,3,I-1)
        flux_temp(L,2,I)=UV(2,I-1,L)*flux_temp(L,1,I-1)             &
           +UV(5,I-1,L)*flux_temp(L,2,I-1)                          &
           +UV(8,I-1,L)*flux_temp(L,3,I-1)
        flux_temp(L,3,I)=UV(3,I-1,l)*flux_temp(L,1,I-1)             &
           +UV(6,I-1,L)*flux_temp(L,2,I-1)                          &
           +UV(9,I-1,L)*flux_temp(L,3,I-1)
        flux_temp(L,4,I)=(gamma(L,1,I)*flux_temp(L,1,I)             &
           +h(L,1,I))*beta(L,1,I)
        flux_temp(L,5,I)=(gamma(L,2,I)*flux_temp(L,2,I)             &
           +h(L,2,I))*beta(L,2,I)
        flux_temp(L,6,I)=(gamma(L,3,I)*flux_temp(L,3,I)             &
           +h(L,3,I))*beta(L,3,I)
        flux_temp(L,1,I)=T_R(1,I,L)*flux_temp(L,1,I)                &
           +T_R(4,I,L)*flux_temp(L,4,I)+T_R(7,I,L)
        flux_temp(L,2,I)=T_R(2,I,L)*flux_temp(L,2,I)                &
           +T_R(5,I,L)*flux_temp(L,5,I)+T_R(9,I,L)
        flux_temp(L,3,I)=T_R(3,I,L)*flux_temp(L,3,I)                &
           +T_R(6,I,L)*flux_temp(L,6,I)+T_R(11,I,L)
     ENDDO


! Calculate the overall flux.
     DO I=N_CLOUD_TOP, N_LAYER
        FLUX_TOTAL(L, 2*I+1)=flux_temp(L,4,I)+flux_temp(L,5,I)      &
           +flux_temp(L,6,I)
        FLUX_TOTAL(L, 2*I+2)=flux_temp(L,1,I)+flux_temp(L,2,I)      &
           +flux_temp(L,3,I)
     ENDDO

  END DO

! Reduce to net fluxes if required.
  IF (L_NET) THEN
     DO I=0, N_LAYER
        DO L=1, N_PROFILE
           FLUX_TOTAL(L, I+1)                                       &
              =FLUX_TOTAL(L, 2*I+2)-FLUX_TOTAL(L, 2*I+1)
        ENDDO
     ENDDO
  ENDIF

  RETURN
  END SUBROUTINE SOLVER_TRIPLE_HOGAN
