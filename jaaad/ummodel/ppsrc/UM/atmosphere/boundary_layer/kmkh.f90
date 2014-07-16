
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE KMKH---------------------------------------------------
!!!
!!!  Purpose: To set the turbulent mixing coefficients KM and KH
!!!           (Note: should be used after any vertical interpolation
!!!                  but before any horizontal interpolation.)
!!!
!!!
!!!
!!!  Programming standard:
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE KMKH (                                                 &
     & row_length,rows,BL_LEVELS, off_x,off_y                           &
     &,Keep_Ri_FA,RHOKM,RHO_KM,RHOKH,RHOKMZ,RHOKHZ                      &
     &,NTML,CUMULUS,RHOKM_TOP,RHOKH_TOP,UNSTABLE                        &
     &,NTDSC,DSC                                                        &
     &,SML_DISC_INV,DSC_DISC_INV                                        &
     &,nSCMDpkgs,L_SCMDiags                                             &
     &,LTIMER                                                           &
     & )

      Use cv_run_mod, Only:                                             &
          l_conv4a

      IMPLICIT NONE

      LOGICAL LTIMER             ! IN Flag for TIMER diagnostics

      INTEGER                                                           &
     & row_length,rows,off_x,off_y                                      &
     &,BL_LEVELS                                                        &
!                  ! IN No. of atmospheric levels for which boundary 
!                       layer fluxes are calculated.
     &,Keep_Ri_FA
                   ! IN switch to keep local mixing in free atmosphere

      LOGICAL                                                           &
     & CUMULUS(row_length,rows)                                         &
                                ! IN flag for Cu in the bl
     &,UNSTABLE(row_length,rows)                                        &
                                ! IN Flag for unstable surface layer.
     &,DSC(row_length,rows)     ! IN Flag set if decoupled
!                               !    stratocumulus layer found
      INTEGER                                                           &
     & NTML(row_length,rows)                                            &
                                ! IN Number of model levels in the
!                                   turbulently mixed layer.
     &,NTDSC(row_length,rows)                                           &
                                ! IN Top level for turb mixing in
!                                   cloud layer
     &,SML_DISC_INV(row_length,rows)                                    &
!                               ! IN Flags for whether discontinuous
     &,DSC_DISC_INV(row_length,rows)
!                               ! IN inversions are diagnosed
      REAL                                                              &
     & RHOKMZ(row_length,rows,2:BL_LEVELS)                              &
!                             ! IN Non-local turbulent mixing
!                                  coefficient for momentum.
     &,RHOKHZ(row_length,rows,2:BL_LEVELS)                              &
!                             ! IN Non-local turbulent mixing
!                                  coefficient for heat and moisture
     &,RHOKM_TOP(row_length,rows,2:BL_LEVELS)                           &
!                             ! IN Non-local top-down turbulent
!                             !    mixing coefficient for momentum.
     &,RHOKH_TOP(row_length,rows,2:BL_LEVELS)
!                             ! IN Non-local top-down turbulent
!                             !    mixing coefficient for heat
!                             !    and moisture.
      REAL                                                              &
     & RHOKM(1-off_x:row_length+off_x,1-off_y:rows+off_y,BL_LEVELS)     &
!                             ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for momentum.
     &,RHOKH(row_length,rows,BL_LEVELS)
!                             ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for heat and moisture.

      REAL                                                              &
     & RHO_KM(row_length,rows,2:BL_LEVELS)
!                             ! OUT RHO * KM before horizontal
!                                   interpolation to UV-grid.
!
! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     & nSCMDpkgs              ! No of SCM diagnostics packages

      LOGICAL                                                           &
     & L_SCMDiags(nSCMDpkgs)  ! Logicals for SCM diagnostics packages
!
!!----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER

!!----------------------------------------------------------------------
! Start blopt8a

! Description:
!   Permissible settings for BL options.
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2      27/01/06 Original code.  J. M. Edwards
!
      INTEGER, PARAMETER :: Off = 0  ! Switch disabled
      INTEGER, PARAMETER :: On  = 1  ! Switch enabled
!
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!     Options for non-gradient stress following
!     Brown and Grant (1997), version 2 including a limit on its size
!
!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)
!
!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2
!
!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used a dynamic criterion in the
!       diagnosis of BL types
!
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
!
!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a

!
!  Define local storage.

      Character*(*), Parameter ::  RoutineName = 'kmkh'

      REAL                                                              &
     & RHOKH_DIAG(row_length,rows,bl_levels)                            &
     &,RHOKM_DIAG(row_length,rows,bl_levels)
!                            ! Arrays for SCM diagnostics,
!                            ! diffusivity of heat and momentum kg/(ms)
!
      INTEGER                                                           &
     & I,j                                                              &
                       ! Loop counter (horizontal field index).
     &,K             ! Loop counter (vertical level index).


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('KMKH    ',3)
      ENDIF


      IF (Keep_Ri_FA == ON) THEN 
!-----------------------------------------------------------------------
!       ! Set local K's to zero at the LCL in cumulus and at the
!       ! top of an unstable BL with a well-defined inversion
!-----------------------------------------------------------------------
        DO K=2,BL_LEVELS
          do j=1,rows
          do i=1,row_length
            IF ( ( CUMULUS(I,j) .OR. SML_DISC_INV(I,j) == 1) .AND.      &
     &           (K == NTML(I,j)+1 .OR. K == NTML(I,j)+2) ) THEN
              RHOKH(I,j,K) = 0.0
              RHOKM(I,j,K) = 0.0
            ENDIF

            IF ( DSC_DISC_INV(I,j)  ==  1 .AND.                         &
     &           (K == NTDSC(I,j)+1 .OR. K == NTDSC(I,j)+2) ) THEN
              RHOKH(I,j,K) = 0.0
              RHOKM(I,j,K) = 0.0
            ENDIF

          ENDDO ! P_POINTS,i
          ENDDO ! P_POINTS,j
        ENDDO ! BL_LEVELS

      ELSE
!-----------------------------------------------------------------------
!       ! Set local K's to zero from the LCL in cumulus and from the
!       ! top of an unstable BL with a well-defined inversion
!-----------------------------------------------------------------------
        DO K=2,BL_LEVELS
          do j=1,rows
          do i=1,row_length
            IF(CUMULUS(I,j) .AND. ( (L_conv4a .AND. K >  NTML(I,j))     &
     &         .OR. (.NOT. L_conv4a .AND. K >= NTML(I,j)) )) THEN
              RHOKH(I,j,K)=0.0
              RHOKM(I,j,K)=0.0
            ENDIF

            IF ( DSC_DISC_INV(I,j)  ==  1 .AND. K  >   NTDSC(I,j) ) THEN
              RHOKH(I,j,K) = 0.0
              RHOKM(I,j,K) = 0.0
            ENDIF

            IF ( SML_DISC_INV(I,j)  ==  1 .AND. K  >   NTML(I,j) ) THEN
!             !   This also means no local mixing within any DSC layer
              RHOKH(I,j,K) = 0.0
              RHOKM(I,j,K) = 0.0
            ENDIF

          ENDDO ! P_POINTS,i
          ENDDO ! P_POINTS,j
        ENDDO ! BL_LEVELS

      ENDIF ! test on Keep_Ri_FA

!-----------------------------------------------------------------------
!     ! Set non-local K's to zero at the LCL in cumulus layers,
!     ! including level NTML if not L_conv4a convection scheme
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        do j=1,rows
        do i=1,row_length

          IF(CUMULUS(I,j) .AND. ( (L_conv4a .AND. K == NTML(I,j)+1)     &
     &         .OR. (.NOT. L_conv4a .AND.                               &
     &                     K >= NTML(I,j).AND.K <  NTML(I,j)+2) )) THEN
           RHOKHZ(I,j,K)=0.0
           RHOKMZ(I,j,K)=0.0
           RHOKH_TOP(I,j,K)=0.0
           RHOKM_TOP(I,j,K)=0.0
          ENDIF
        ENDDO ! P_POINTS,i
        ENDDO ! P_POINTS,j
      ENDDO ! BL_LEVELS

!-----------------------------------------------------------------------
!     ! Set KM and KH to be the maximum of the local and non-local 
!     ! values andstore RHO_KM on P-grid for output.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        do j=1,rows
        do i=1,row_length

          RHOKH(I,j,K) = MAX( RHOKH(I,j,K) ,                            &
     &                                  RHOKHZ(I,j,K)+RHOKH_TOP(I,j,K) )
          RHOKM(I,j,K) = MAX( RHOKM(I,j,K) ,                            &
     &                                  RHOKMZ(I,j,K)+RHOKM_TOP(I,j,K) )

          RHO_KM(I,j,K) = RHOKM(I,j,K)
        ENDDO ! P_POINTS,i
        ENDDO ! P_POINTS,j
      ENDDO ! BL_LEVELS
!

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('KMKH    ',4)
      ENDIF

      RETURN
      END SUBROUTINE KMKH
