
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!*LL  SUBROUTINE IMP_MIX -----------------------------------------------
!LL
!LL  Purpose: Calculate turbulent mixing increments for a passive tracer
!LL           using an implicit numerical scheme.  The tridiagonal
!LL           matices are inverted using simple Gaussian elimination.
!LL
!LL
!LL  Model           Modification history :
!LL version  Date
!LL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
!LL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
!LL                                     S J Swarbrick
!LL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!LL  5.1   20/06/00  Change references of 1.e30 to rmdi. R Rawlins
!LL  5.4   30/05/02  Remove redundant rapid mixing code, include
!LL                  spherical geometry in flux divergences and
!LL                  calculate divergences in height rather than
!LL                  pressure coordinates.   Adrian Lock
!    5.5 17/04/03 Remove references to obsolete sections
!                 A03_3A,5A,5B,7A. T.White
!LL  5.5   17/02/03  Remove duplicate declaration of DTRDZ M.Hughes
!LL
!LL  SDJ  <- Programmers of some or all of previous code or changes
!LL
!LL  Programming standard: UM Documentation Paper No. 3
!LL
!LL  Documentation: UM Documentation Paper No 24.
!LL
!*----------------------------------------------------------------------
!*L  Arguments :-
      SUBROUTINE IMP_MIX (                                              &
     & halo_i, halo_j, row_length, rows, bl_levels                      &
     &,gamma_rhokh_rdz, gamma_rhok_dep                                  &
     &,dtrdz, r_rho_levels, r_theta_levels                              &
     &,timestep                                                         &
     &,f_field,surf_dep_flux,field                                      &
     &,error,ltimer                                                     &
     & )
      IMPLICIT NONE
!
!  Inputs :-
!
      INTEGER                                                           &
     & row_length                                                       &
     &,rows                                                             &
     &,halo_i                                                           &
                                   ! Size of halo in i direction.
     &,halo_j                                                           &
                                   ! Size of halo in j direction.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
!                                  !    which boundary layer fluxes are
!                                  !    calculated.
      REAL                                                              &
     & GAMMA_RHOKH_RDZ(row_length,rows,2:BL_LEVELS)                     &
!                                  ! IN Turbulent mixing coefs. above
!                                  !    surface, =GAMMA(K)*RHOKH(,K)
!                                  !    *RDZ(K) for K>=2 (from KMKH).
     &,GAMMA_RHOK_DEP(row_length,rows)                                  &
                                       ! IN Surface exchange coefficient
!                                  !    for surface deposition*GAMMA(1)

!                                  ! IN  dt/(rho*r*r*dz) for scalar
!                                  !     flux divergence
     &,TIMESTEP                    ! IN Timestep in seconds.
!
! Co-ordinate arrays
      Real                                                              &
                                   ! IN local vertical co-ordinates
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:BL_LEVELS)              &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j,BL_LEVELS)

!  Next 2 arrays are IN as "explicit" fluxes and OUT as "implicit"
!  fluxes.
!
      REAL                                                              &
     & F_FIELD(row_length,rows,BL_LEVELS)                               &
                                           ! INOUT Flux of tracer
     &,SURF_DEP_FLUX(row_length,rows)                                   &
                                           ! INOUT surface deposition
                                           !       flux
     &,FIELD(row_length,rows,BL_LEVELS)    ! INOUT Amount of tracer

      INTEGER                                                           &
     & ERROR                       ! OUT 1 if bad arguments, else 0.

      LOGICAL                                                           &
     & LTIMER                      ! IN Logical switch for TIMER diagnos

!*
!*L  External references :-

      EXTERNAL TIMER

!*
!*L  Local and other symbolic constants :-
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!*
!*L Workspace :-
!   4*BL_LEVELS + 4 blocks of real workspace are required.
      REAL                                                              &
     & AF(row_length,rows,BL_LEVELS)                                    &
                                           ! Elements in rows in matrix
!                                  ! equation (modified during
!                                  ! Gaussian elimination calculations).
     &,D_FIELD(row_length,rows,BL_LEVELS)                               &
                                           ! Delta FIELD (tracer field)
!                                  ! elements of vector on RHS, then
!                                  ! LHS, of eqn P244.79.
     &,DTRDZ(row_length,rows,BL_LEVELS)
                                   ! -g.dt/dp for bottom BL_LEVELS
!                                  ! model layers (needed in P245).
!*
!  Local scalars :-
      REAL                                                              &
     & CF                                                               &
                ! Matrix element for local increments
     &,RBF                                                              &
                ! Reciprocal of B for local increments
     &,r_sq                                                             &
                ! r*r
     &,rr_sq    ! 1/(r*r)

      INTEGER                                                           &
     & BLM1                                                             &
                ! BL_LEVELS minus 1.
     &,I                                                                &
                ! Loop counter (horizontal field index).
     &,J                                                                &
                ! Offset version of I.
     &,K                                                                &
                ! Loop counter (vertical index).
     &,KM1                                                              &
                ! K minus 1.
     &,KP1      ! K plus 1.
!
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IMPMIX  ',3)
      ENDIF

      ERROR=0

      BLM1 = BL_LEVELS-1

!L----------------------------------------------------------------------
!L (A) Calculations on P-grid.
!-----------------------------------------------------------------------
!! 1.0 For simulations on a sphere use spherical geometry for vertical
!!     derivatives so multiply fluxes and rho*K_h by r*r
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
       do j=1,rows
       DO I = 1,row_length
         r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
         GAMMA_RHOKH_RDZ(I,J,K) = r_sq * GAMMA_RHOKH_RDZ(I,J,K)
         F_FIELD(I,j,K)         = r_sq * F_FIELD(I,j,K)
       ENDDO
       ENDDO
      ENDDO
      do j=1,rows
      DO I = 1,row_length
        r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
        GAMMA_RHOK_DEP(I,j) = r_sq * GAMMA_RHOK_DEP(I,j)
        F_FIELD(I,j,1)      = r_sq * F_FIELD(I,j,1)
        SURF_DEP_FLUX(I,j)  = r_sq * SURF_DEP_FLUX(I,j)
      ENDDO
      ENDDO
!L----------------------------------------------------------------------
!L 4.  Calculate those matrix and vector elements on the LHS of eqn
!L     which are to do with implicit solution of the tracer
!L     transport problem at the surface and between all levels.
!L     Begin with "upward sweep" through lower half of matrix).
!-----------------------------------------------------------------------
!
      DO J= 1, rows
      DO I= 1, row_length
!-----------------------------------------------------------------------
!L 4.2 Lowest atmospheric layer FIELD row of matrix.
!-----------------------------------------------------------------------
!         ! "Explicit" increment to FIELD(1)
          D_FIELD(I,j,1) = -DTRDZ(I,j,1) *                              &
     &                   ( F_FIELD(I,j,2) - F_FIELD(I,j,1)              &
     &                                  - SURF_DEP_FLUX(I,j) )

          CF = -DTRDZ(I,j,1) * GAMMA_RHOK_DEP(I,j)
          AF(I,j,1) = -DTRDZ(I,j,1) * GAMMA_RHOKH_RDZ(I,j,2)
          RBF = 1.0 / ( 1.0 - AF(I,j,1) - CF )
          D_FIELD(I,j,1) = RBF * D_FIELD(I,j,1)
          AF(I,j,1) = RBF * AF(I,j,1)
      END Do
      END Do
!-----------------------------------------------------------------------
!L 4.3 Rows of matrix applying to FIELD transport into model layers in
!L     the "middle" of the "boundary" layer, i.e. all but the bottom
!L     layer and the top "boundary" layer.
!-----------------------------------------------------------------------
      DO K=2,BLM1
        KP1 = K+1
        KM1 = K-1
        DO J= 1, rows
        DO I= 1, row_length
!
!   "Explicit" flux divergence across layer giving explicit FIELD
!   increment due to mixing
!
            D_FIELD(I,j,K) = -DTRDZ(I,j,K) *                            &
     &                       (F_FIELD(I,j,KP1) - F_FIELD(I,j,K))
            AF(I,j,K) = -DTRDZ(I,j,K) * GAMMA_RHOKH_RDZ(I,j,KP1)
            CF = -DTRDZ(I,j,K) * GAMMA_RHOKH_RDZ(I,j,K)
              RBF = 1.0 / ( 1.0 - AF(I,j,K)                             &
     &                            - CF * ( 1.0 + AF(I,j,KM1) ) )
              D_FIELD(I,j,K) = RBF * ( D_FIELD(I,j,K)                   &
     &                                    - CF*D_FIELD(I,j,KM1) )
            AF(I,j,K) = RBF * AF(I,j,K)
        END Do
        END Do
      END Do
!-----------------------------------------------------------------------
!L 4.4 Top "boundary" layer FIELD row of matrix. FIELD for this layer
!L     can then be, and is, updated.
!-----------------------------------------------------------------------
      DO J=1, rows
      DO I=1, row_length
        D_FIELD(I,j,BL_LEVELS) = DTRDZ(I,j,BL_LEVELS) *                 &
     &                           F_FIELD(I,j,BL_LEVELS)
!
        CF = -DTRDZ(I,j,BL_LEVELS) * GAMMA_RHOKH_RDZ(I,j,BL_LEVELS)
          RBF = 1.0 / ( 1.0 - CF*( 1.0 + AF(I,j,BLM1) ) )
          D_FIELD(I,j,BL_LEVELS) = RBF * ( D_FIELD(I,j,BL_LEVELS)       &
     &                                        - CF*D_FIELD(I,j,BLM1) )
        FIELD(I,j,BL_LEVELS) = FIELD(I,j,BL_LEVELS) +                   &
     &                         D_FIELD(I,j,BL_LEVELS)
      END Do
      END Do
!
!-----------------------------------------------------------------------
!L 5.  "Downward sweep" through whole matrix.  FIELD is updated when
!L     the final implicit increments have been calculated.
!-----------------------------------------------------------------------
!L 5.1 Remaining FIELD rows of matrix and add implicit increments
!-----------------------------------------------------------------------
!
      DO K=BLM1,1,-1
        DO J=1, rows
        DO I=1, row_length
            D_FIELD(I,j,K) = D_FIELD(I,j,K) -                           &
     &                       AF(I,j,K)*D_FIELD(I,j,K+1)
            FIELD(I,j,K) = FIELD(I,j,K) + D_FIELD(I,j,K)
        END Do
        END Do
          END Do
!
!-----------------------------------------------------------------------
!L 6.  Calculate final implicit flux of tracer.
!-----------------------------------------------------------------------
!! 6.0 First divide by r*r for diagnostics.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
       do j=1,rows
       DO I = 1,row_length
         rr_sq = 1.0 / ( r_rho_levels(i,j,k)*r_rho_levels(i,j,k) )
         GAMMA_RHOKH_RDZ(I,J,K) = rr_sq * GAMMA_RHOKH_RDZ(I,J,K)
         F_FIELD(I,j,K)         = rr_sq * F_FIELD(I,j,K)
       ENDDO
       ENDDO
      ENDDO
      do j=1,rows
      DO I = 1,row_length
        rr_sq = 1.0 / ( r_theta_levels(i,j,0)*r_theta_levels(i,j,0) )
        GAMMA_RHOK_DEP(I,j) = rr_sq * GAMMA_RHOK_DEP(I,j)
        F_FIELD(I,j,1)      = rr_sq * F_FIELD(I,j,1)
        SURF_DEP_FLUX(I,j)  = rr_sq * SURF_DEP_FLUX(I,j)
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
!L 6.1 Surface fluxes.
!-----------------------------------------------------------------------
!
      DO j = 1, rows
      DO i = 1, row_length
          SURF_DEP_FLUX(I,j) = SURF_DEP_FLUX(I,j)                       &
     &                        - GAMMA_RHOK_DEP(I,j) * D_FIELD(I,j,1)
        END Do
      END Do
!
!-----------------------------------------------------------------------
!L 6.2 Fluxes at layer interfaces above the surface.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        KM1 = K-1
        DO j = 1, rows
        DO i = 1, row_length
!
!  Calculate and store implicit fluxes due to local mixing.
!
            F_FIELD(I,j,K) = F_FIELD(I,j,K) - GAMMA_RHOKH_RDZ(I,j,K)    &
     &                       * ( D_FIELD(I,j,K) - D_FIELD(I,j,KM1) )
        END Do
        END Do
      END Do
!
  999 CONTINUE   ! Branch for error exit.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('IMPMIX  ',4)
      ENDIF

      RETURN
      END SUBROUTINE IMP_MIX
