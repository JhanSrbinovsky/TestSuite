
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE TR_MIX_OPT --------------------------------------------
!LL
!LL Purpose:Calculate tracer flux and pass through to IMP_MIX to solve
!LL
!LL  SDJ  <- Programmers of some or all of previous code or changes
!LL
!LL  Model           Modification history:
!LL version  Date
!
!        6.2  15/02/06 New deck. Optimised version of TRMIX2C at vn 6.1
!                            Vector optimisation by Klaus Ketelman:
!                            Split i,j loop into 3 and move them into
!                            IENT loop             Stephanie Woodward
!
!        6.4   Nov 06  Make tests on "a /= 0" safe        A P Lock
!LL
!LL
!LL  Programming standard: UM Documentation Paper No 3
!LL
!LL  Documentation: UM Documentation Paper No 24.
!LL
!*---------------------------------------------------------------------
!*L  Arguments :-
      SUBROUTINE TR_MIX_OPT (                                           &
     & halo_i, halo_j, row_length, rows, BL_LEVELS, off_x, off_y        &
     &,gamma                                                            &
     &,RHOKH_RDZ,RHOKH_1                                                &
     &,DTRDZ, r_rho_levels, r_theta_levels                              &
     &,TIMESTEP                                                         &
     &,F_FIELD,FIELD                                                    &
     &,SURF_EM,RES_FACTOR,SURF_DEP_FLUX                                 &
     &,KENT, WE_LIM, T_FRAC, ZRZI                                       &
     &,KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                       &
     &,ZH ,ZHSC, Z_HALF                                                 &
     &,ERROR,LTIMER                                                     &
     &)

      IMPLICIT NONE
      INTEGER                                                           &
     & row_length                                                       &
     &, rows                                                            &
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j                                                          &
                    ! Size of halo in j direction.
     &, off_x                                                           &
                    ! Size of small halo in i
     &, off_y                                                           &
                    ! Size of small halo in j.
     &, BL_LEVELS   ! IN No. of atmospheric levels for boundary layer

! Co-ordinate arrays
      Real                                                              &
                                   ! IN local vertical co-ordinates
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:BL_LEVELS)              &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j,BL_LEVELS)

      REAL                                                              &
     & gamma(bl_levels)                                                 &
     &,RHOKH_RDZ(row_length,rows,2:BL_LEVELS)                           &
!                                  ! IN Mixing coeff. above surface
!                                  !    = RHOKH(,K)*RDZ(K)
!                                  !    for K>=2 (from IMP_SOLVER).
     &,RHOKH_1(row_length,rows)                                         &
                                   ! IN  Surface exchange coeff.
!                                  !     from P243 (SF_EXCH)
     &,DTRDZ(row_length,rows,BL_LEVELS)                                 &
!                                  ! IN  dt/(rho*r*r*dz) for scalar
!                                  !     flux divergence
     &,TIMESTEP                                                         &
                                   ! IN Timestep in seconds.
     &,SURF_EM(row_length,rows)                                         &
                                   ! IN, Surface emissions in kg/m2/s
     &,RES_FACTOR(row_length,rows)                                      &
                                   ! IN, dry dep coeff=Ra/(Ra+Rb+Rc)
     &, WE_LIM(row_length,rows,3)                                       &
                                  ! IN rho*entrainment rate implied by
!                                  !     placing of subsidence
     &, ZRZI(row_length,rows,3)                                         &
                                   ! IN (z-z_base)/(z_i-z_base)
     &, T_FRAC(row_length,rows,3)                                       &
                                   ! IN a fraction of the timestep
     &, WE_LIM_DSC(row_length,rows,3)                                   &
!                                 ! IN rho*entrainment rate implied by
!                                  !     placing of subsidence
     &, ZRZI_DSC(row_length,rows,3)                                     &
                                   ! IN (z-z_base)/(z_i-z_base)
     &, T_FRAC_DSC(row_length,rows,3)                                   &
!                                  ! IN a fraction of the timestep
     &, Z_HALF(row_length,rows,BL_LEVELS)                               &
                                   ! IN Z_HALF(*,K) is height of half
!                                  !    level k-1/2.
     &, ZHSC(row_length,rows)                                           &
                                   ! IN Top of decoupled layer
     &, ZH(row_length,rows)        ! IN Top of surface mixed layer
!
      REAL                                                              &
     & F_FIELD(row_length,rows,BL_LEVELS)                               &
                                   ! OUT Flux of tracer in kg/m2/s.
     &,FIELD(row_length,rows,BL_LEVELS)                                 &
                                   ! INOUT Tracer amount in kg/kg.
     &,SURF_DEP_FLUX(row_length,rows)
                                  ! OUT, surface deposn flux (kg/m2/s)
!
      INTEGER                                                           &
     &  KENT(row_length,rows)                                           &
                                   ! IN grid-level of SML inversion
     &, KENT_DSC(row_length,rows)                                       &
                                   ! IN grid-level of DSC inversion
     &,ERROR                       ! OUT 1 if bad arguments, else 0.
!
      LOGICAL                                                           &
     & LTIMER                      ! IN Logical switch for TIMER
!                                  !    diagnostics
!*
!*L  External references :-
      EXTERNAL IMP_MIX,TIMER
!*
!*L  Local and other symbolic constants :-
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------

      REAL,PARAMETER:: ONE=1.0
      REAL,PARAMETER:: SMALLP=TINY(ONE)

!*
!*L Workspace :-
!
      REAL                                                              &
     & RHOK_DEP(row_length, rows)                                       &
                                    ! Surface deposition coeficient
     &,GAMMA_RHOKH_RDZ(row_length,rows,2:BL_LEVELS)                     &
!                            ! gamma*RHOKH_RDZ
     &,DZ_DISC                                                          &
                             ! Temporary in subgrid zi calculation
     &,DZLKP1                                                           &
                             ! Thickness of theta-level K+1
     &,F_FIELD_ENT                                                      &
                             ! Time-level n entrainment flux
     &,DFIELD_INV            ! inversion jump
!*
!  Local scalars :-
      INTEGER                                                           &
     & I,J                                                              &
                ! Loop counter (horizontal field index).
     &,K                                                                &
                ! Loop counter (vertical index).
     &,KM1                                                              &
                ! Max(K-1,2)
     &,IENT     ! Loop counter for entrainment

!
      real,dimension(row_length,rows)         :: DFIELD_SML,DFIELD_DSC
!
!---------------------------------------------------------------------
!L  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!---------------------------------------------------------------------
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('TRMIX   ',3)
      END IF
!

      ERROR=0
!
!
!---------------------------------------------------------------------
!L 1.  Calculate flux of tracer:
!---------------------------------------------------------------------
!L 1.1 Above the surface
!---------------------------------------------------------------------
!
      DO 1 K=2,BL_LEVELS
        DO J=1, rows
          DO I=1, row_length
            GAMMA_RHOKH_RDZ(I,J,K) = GAMMA(K) * RHOKH_RDZ(I,J,K)
            F_FIELD(I,J,K) = - RHOKH_RDZ(I,J,K) *                       &
     &                         (FIELD(I,J,K) - FIELD(I,J,K-1))
          End Do
        End Do
   1  CONTINUE
!
!---------------------------------------------------------------------
!L 1.2 At the surface: (i) set surface flux equal to input emissions
!L                  (should be passed in as ancillary file, else ZERO)
!L                     (ii) Use input resistance factors to calculate
!L                 surface deposition (if ZERO then no dry deposition)
!---------------------------------------------------------------------
!
        DO J=1, rows
          DO I=1, row_length
            F_FIELD(I,J,1) = SURF_EM(I,J)  ! Inject surface emissions
            RHOK_DEP(I,J) = RES_FACTOR(I,J) * RHOKH_1(I,J)
            SURF_DEP_FLUX(I,J) = -RHOK_DEP(I,J) * FIELD(I,J,1)
            RHOK_DEP(I,J) = GAMMA(1) * RHOK_DEP(I,J)
          End Do
        End Do

!
!----------------------------------------------------------------------
! Add on explicit entrainment fluxes
! These are calculated from the time-level n profile and parametrized
! entrainment rate to give a time-level n flux (F_FIELD_ENT).  This is
! then implemented using an equivalent RHOKH as the diagnosed inversion
! jump can change significantly across the timestep and so, therefore,
! should the parametrized flux.
!---------------------------------------------------------------------
      DO J=1, rows
      DO I=1, row_length
        K = KENT(I,j)-1     ! equal to originally diagnosed NTML
!       !----------------------------
!       ! diagnose SML inversion jump
!       !----------------------------
        IF (K  ==  BL_LEVELS-1) THEN
          DFIELD_SML(i,j) = FIELD(I,j,K+1) - FIELD(I,j,K)
        ELSE
          DZLKP1  = Z_HALF(I,j,K+2) - Z_HALF(I,j,K+1)
          DZ_DISC = Z_HALF(I,j,K+2) - ZH(I,j)
          IF (DZ_DISC/DZLKP1  >   0.1) THEN
            DFIELD_SML(i,j) = (FIELD(I,j,K+1)-FIELD(I,j,K)) *           &
     &                                                 DZLKP1 /DZ_DISC

            IF ( FIELD(I,j,K+2)  >   FIELD(I,j,K+1) .AND.               &
     &                 FIELD(I,j,K+1)  >   FIELD(I,j,K) ) THEN
              DFIELD_SML(i,j) = MIN( FIELD(I,j,K+2)-FIELD(I,j,K),       &
     &                                                DFIELD_SML(i,j))
            ELSEIF ( FIELD(I,j,K+2)  <   FIELD(I,j,K+1) .AND.           &
     &                 FIELD(I,j,K+1)  <   FIELD(I,j,K) ) THEN
              DFIELD_SML(i,j) = MAX( FIELD(I,j,K+2)-FIELD(I,j,K),       &
     &                                                DFIELD_SML(i,j))
            ELSE  ! FIELD non-monotonic
              DFIELD_SML(i,j) = FIELD(I,j,K+1)-FIELD(I,j,K)
            END IF
          ELSE
            DFIELD_SML(i,j) = FIELD(I,j,K+2) - FIELD(I,j,K)
          END IF
        END IF

        K = KENT_DSC(I,j)-1     ! equal to originally diagnosed NTDSC

!       !----------------------------
!       ! diagnose DSC inversion jump
!       !----------------------------
        IF (K  ==  BL_LEVELS-1) THEN
          DFIELD_DSC(i,j) = FIELD(I,j,K+1) - FIELD(I,j,K)
        ELSEIF ( KENT_DSC(I,j)  >   2 ) THEN

          DZLKP1  = Z_HALF(I,j,K+2) - Z_HALF(I,j,K+1)
          DZ_DISC = Z_HALF(I,j,K+2) - ZHSC(I,j)
          IF (DZ_DISC/DZLKP1  >   0.1) THEN
            DFIELD_DSC(i,j) = (FIELD(I,j,K+1)-FIELD(I,j,K)) *           &
     &                                                DZLKP1 /DZ_DISC

            IF ( FIELD(I,j,K+2)  >   FIELD(I,j,K+1) .AND.               &
     &                 FIELD(I,j,K+1)  >   FIELD(I,j,K) ) THEN
              DFIELD_DSC(i,j) = MIN( FIELD(I,j,K+2)-FIELD(I,j,K),       &
     &                                                DFIELD_DSC(i,j))
            ELSEIF ( FIELD(I,j,K+2)  <   FIELD(I,j,K+1) .AND.           &
     &                 FIELD(I,j,K+1)  <   FIELD(I,j,K) ) THEN
              DFIELD_DSC(i,j) = MAX( FIELD(I,j,K+2)-FIELD(I,j,K),       &
     &                                                DFIELD_DSC(i,j))
            ELSE  ! FIELD non-monotonic
              DFIELD_DSC(i,j) = FIELD(I,j,K+1)-FIELD(I,j,K)
            END IF
          ELSE
            DFIELD_DSC(i,j) = FIELD(I,j,K+2) - FIELD(I,j,K)
          END IF
        END IF
      END DO
      END DO                       !close I and J loops
!       !------------------------------------
!       ! calculate entrainment fluxes and KH
!       !------------------------------------
      DO IENT = 1,3
        DO J=1, rows
!CDIR nodep
        DO I=1, row_length
          K = KENT(I,j)-2+IENT
          IF ( K > 1 .AND. K <= BL_LEVELS .AND.                         &
     &                     T_FRAC(I,j,IENT) > 0.0) THEN
          IF ( ABS(FIELD(I,j,K)-FIELD(I,j,K-1)) >= SMALLP ) THEN
           DFIELD_INV = DFIELD_SML(i,j)
           IF ( DFIELD_SML(i,j)/(FIELD(I,j,K)-FIELD(I,j,K-1)) <  0.0 )  &
!             ! DFIELD_INV must have same sign as
!             ! local gradient to get right sign of flux
     &        DFIELD_INV = FIELD(I,j,K)-FIELD(I,j,K-1)
           F_FIELD_ENT = T_FRAC(I,j,IENT) * ( F_FIELD(I,j,1)            &
     &        - ( WE_LIM(I,j,IENT) * DFIELD_INV + F_FIELD(I,j,1) )      &
     &          * ZRZI(I,j,IENT) )
!           ! interpolation to surface flux must not change the sign
!           ! of the entrainment flux otherwise KH will be <0!
           IF ( F_FIELD_ENT / (FIELD(I,j,K)-FIELD(I,j,K-1))  >   0.0 )  &
     &        F_FIELD_ENT = - T_FRAC(I,j,IENT) *                        &
     &                  WE_LIM(I,j,IENT) * DFIELD_INV * ZRZI(I,j,IENT)
!           ! Restrict size of RHOKH for numerical safety
            KM1 = MAX( K-1, 2 )
            GAMMA_RHOKH_RDZ(I,j,K) = GAMMA_RHOKH_RDZ(I,j,K) +           &
     &              MIN( - GAMMA(K)*F_FIELD_ENT/                        &
     &                     ( FIELD(I,j,K)-FIELD(I,j,K-1) ),             &
     &                   10.0*GAMMA_RHOKH_RDZ(I,j,KM1) )
!           ! Recalculate explicit flux using entrainment KH
            F_FIELD(I,J,K) = - (GAMMA_RHOKH_RDZ(I,J,K) / GAMMA(K)) *    &
     &                         (FIELD(I,J,K) - FIELD(I,J,K-1))
          END IF
          END IF
        end do
        end do
      END DO

      DO IENT = 1,3
        DO J=1, rows
!CDIR nodep
        DO I=1, row_length
          K = KENT_DSC(I,j)-2+IENT
          IF ( KENT_DSC(I,j) >= 3 .AND. K <= BL_LEVELS ) THEN
            IF ( T_FRAC_DSC(I,j,IENT)  >   0.0                          &
     &           .AND. ABS(FIELD(I,j,K)-FIELD(I,j,K-1)) >= SMALLP ) THEN
             DFIELD_INV = DFIELD_DSC(i,j)
!             ! DFIELD_INV must have same sign as
!             ! local gradient to get right sign of flux
             IF (DFIELD_DSC(i,j)/(FIELD(I,j,K)-FIELD(I,j,K-1)) <  0.0)  &
     &          DFIELD_INV = FIELD(I,j,K)-FIELD(I,j,K-1)
             F_FIELD_ENT = - T_FRAC_DSC(I,j,IENT) *                     &
     &          WE_LIM_DSC(I,j,IENT) * DFIELD_INV * ZRZI_DSC(I,j,IENT)
!             ! Restrict size of RHOKH for numerical safety
             KM1 = MAX( K-1, 2 )
             GAMMA_RHOKH_RDZ(I,j,K) = GAMMA_RHOKH_RDZ(I,j,K) +          &
     &                 MIN( - GAMMA(K)*F_FIELD_ENT/                     &
     &                                ( FIELD(I,j,K)-FIELD(I,j,K-1) ),  &
     &                      10.0*GAMMA_RHOKH_RDZ(I,j,KM1) )
!             ! Recalculate explicit flux using entrainment KH
              F_FIELD(I,J,K) = - (GAMMA_RHOKH_RDZ(I,J,K) / GAMMA(K)) *  &
     &                         (FIELD(I,J,K) - FIELD(I,J,K-1))
            END IF
          END IF
        end do
        end do
      END DO


!---------------------------------------------------------------------
!L 2.  Call routine IMPL_CAL to calculate incrememnts to tracer field
!L     and suface deposition flux for output
!---------------------------------------------------------------------
!
! DEPENDS ON: imp_mix
      CALL IMP_MIX (                                                    &
     & halo_i,halo_j,row_length,rows,BL_LEVELS                          &
     &,GAMMA_RHOKH_RDZ,RHOK_DEP                                         &
     &,DTRDZ, r_rho_levels, r_theta_levels                              &
     &,TIMESTEP                                                         &
     &,F_FIELD,SURF_DEP_FLUX,FIELD                                      &
     &,ERROR,LTIMER                                                     &
     & )

!

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('TRMIX   ',4)
      END IF

      RETURN
      END SUBROUTINE TR_MIX_OPT
