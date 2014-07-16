
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!
!!!  Purpose: Calculate explicit fluxes of TL and QT
!!!
!!!  Model           Modification history
!!! version  Date
!!!  4.4    09/09/97   New deck  R.N.B.Smith
!!!  5.2    14/03/00   Non-local scalar entrainment fluxes now
!!!                    parametrized explicity, rather than through
!!!                    RHOKH.     A.P.Lock
!!!  6.2    23/01/06   Improvements to Single Column Model Diagnostics
!!!                    System                          A. Kerr-Munslow
!!!  6.2    16/02/06   Optional updated parametrization for
!!!                    non-gradient scalar fluxes.     A.P.Lock
!!!
!!!  RNBS  <- Programmers of some or all of previous code or changes
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 3
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------

! SUBROUTINE EX_FLUX_TQ

      SUBROUTINE EX_FLUX_TQ (                                           &
     &  row_length,rows, BL_LEVELS, FLUX_GRAD                           &
     &, TL, QW, RDZ, FTL, FQW                                           &
     &, RHOKH, RHOKHZ                                                   &
     &, GRAD_T_ADJ, GRAD_Q_ADJ, RHOF2, RHOFSC                           &
     &, FT_NT, FQ_NT, FT_NT_DSCB, FQ_NT_DSCB                            &
     &, TOTHF_ZH, TOTHF_ZHSC, TOTQF_ZH, TOTQF_ZHSC                      &
     &, NTML, NTDSC, NBDSC                                              &
     &, nSCMDpkgs,L_SCMDiags                                            &
     &, LTIMER                                                          &
     &  )


      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL                                                           &
     &  LTIMER          ! IN Flag for TIMER diagnostics

      INTEGER                                                           &
     & row_length,rows                                                  &
     &,BL_LEVELS                                                        &
                              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed ! <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA
     &,FLUX_GRAD                                                        &
                              ! IN Switch for non-gradient formulation
     &,NTML(row_length,rows)                                            &
                              ! IN Number of model layers in turbulent
!                                  mixing layer.
     &,NTDSC(row_length,rows)                                           &
                              ! IN Top level for turb mixing in any
!                                  decoupled Sc layer
     &,NBDSC(row_length,rows) ! IN lowest flux level in DSC layer

      REAL                                                              &
     &  TL(row_length,rows, BL_LEVELS)                                  &
!                             ! IN Liquid/frozen water temperture (K)
     &, QW(row_length,rows, BL_LEVELS)                                  &
!                             ! IN Total water content (kg/kg)
     &, RHOKH(row_length,rows, BL_LEVELS)                               &
!                             ! IN Exchange coeffs for scalars
     &, RHOKHZ(row_length,rows,2:BL_LEVELS)                             &
!                             ! IN Non-local turbulent mixing
!                             !    coefficient for heat and moisture
     &, RDZ(row_length,rows, BL_LEVELS)                                 &
!                             ! IN RDZ(,1) is the reciprocal
!                             !     height of level 1, i.e. of the
!                             !     middle of layer 1.  For K > 1,
!                             !     RDZ(,K) is the reciprocal of the
!                             !     vertical distance from level
!                             !     K-1 to level K.
     &, GRAD_T_ADJ(row_length,rows)                                     &
!                             ! IN Temperature gradient adjustmenent
!                             !    for non-local mixing in unstable
!                             !    turbulent boundary layer.
     &, GRAD_Q_ADJ(row_length,rows)
!                             ! IN Humidity gradient adjustment
!                                  for non-local mixing in unstable
!                                  turbulent boundary layer.

! Additional variables for SCM diagnostics which are dummy in full UM
      INTEGER                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      LOGICAL                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!
      REAL, DIMENSION(row_length,rows, BL_LEVELS+1) ::                  &
     &  FT_NT                                                           &
                              ! IN Non-turbulent heat and moisture flux
     &, FQ_NT                 !    (on rho levels, surface flux(K=1)=0)
      REAL, DIMENSION(row_length,rows, 2:BL_LEVELS) ::                  &
     &  RHOF2                                                           &
                              ! IN f2 and fsc term shape profiles
     &, RHOFSC                !

      REAL, DIMENSION(row_length,rows) ::                               &
     &  TOTHF_ZH                                                        &
                              ! IN Total heat fluxes at inversions
     &, TOTHF_ZHSC                                                      &
                              !
     &, TOTQF_ZH                                                        &
                              ! IN Total moisture fluxes at inversions
     &, TOTQF_ZHSC                                                      &
                              !
     &, FT_NT_DSCB                                                      &
                              ! IN Non-turbulent heat and moisture flux
     &, FQ_NT_DSCB            !      at the base of the DSC layer.

! ARGUMENTS WITH INTENT INOUT.

      REAL                                                              &
     &  FTL(row_length,rows, BL_LEVELS)                                 &
!                            ! INOUT FTL(,K) contains net turb
!                                   sensible heat flux into layer K
!                                   from below; so FTL(,1) is the
!                                   surface sensible heat, H. (W/m2)
     &, FQW(row_length,rows, BL_LEVELS)
!                            ! INOUT Moisture flux between layers
!                                   (kg per square metre per sec).
!                                   FQW(,1) is total water flux
!                                   from surface, 'E'.





      INTEGER                                                           &
     &  I, J, K

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
!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------

      REAL                                                              &
     & GRCP
      PARAMETER (                                                       &
     & GRCP = G/CP                                                      &
     & )


! LOCAL VARIABLES.

       Character*(*), Parameter ::  RoutineName = 'ex_flux_tq'

       REAL                                                             &
     &  GRAD_FTL(row_length,rows,BL_LEVELS)                             &
                                            ! K*dth/dz
     &, GRAD_FQW(row_length,rows,BL_LEVELS)                             &
                                            ! K*dq/dz
     &, NON_GRAD_FTL(row_length,rows,BL_LEVELS)                         &
                                                ! Non-gradient flux
     &, NON_GRAD_FQW(row_length,rows,BL_LEVELS)                         &
                                                ! Non-gradient flux
     &, F2_FTL(row_length,rows,BL_LEVELS)                               &
                                            ! Heat flux: f2 term
     &, F2_FQW(row_length,rows,BL_LEVELS)                               &
                                            ! Moisture flux: f2 term
     &, FSC_FTL(row_length,rows,BL_LEVELS)                              &
                                            ! Heat flux: fsc term
     &, FSC_FQW(row_length,rows,BL_LEVELS)  ! Moisture flux: fsc term

      EXTERNAL TIMER

!-----------------------------------------------------------------------

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXFLUXTQ',3)
      ENDIF

      DO K=1,BL_LEVELS
        do j=1,rows
        DO I=1,row_length
          GRAD_FTL(I,j,K)=0.0
          GRAD_FQW(I,j,K)=0.0
          NON_GRAD_FTL(I,j,K)=0.0
          NON_GRAD_FQW(I,j,K)=0.0
          F2_FTL(I,j,K) =0.0
          F2_FQW(I,j,K) =0.0
          FSC_FTL(I,j,K)=0.0
          FSC_FQW(I,j,K)=0.0
        ENDDO
        ENDDO
      ENDDO

      DO K=2,BL_LEVELS
!-----------------------------------------------------------------------
!! 1. "Explicit" fluxes of TL and QW, on P-grid.
!-----------------------------------------------------------------------
        do j=1,rows
        DO I=1,row_length
          GRAD_FTL(I,j,K)= - RHOKH(I,j,K) *                             &
     &      ( ( ( TL(I,j,K) - TL(I,j,K-1) ) * RDZ(I,j,K) ) + GRCP )
          GRAD_FQW(I,j,K)= - RHOKH(I,j,K) *                             &
     &          ( QW(I,j,K) - QW(I,j,K-1) ) * RDZ(I,j,K)
!         !-------------------------------------------------------------
!         ! Entrainment fluxes were specified directly in FTL,FQW in
!         ! KMKHZ so now add on gradient-dependent fluxes (note that
!         ! RHOKH should be zero at these levels).
!         !-------------------------------------------------------------
          FTL(I,j,K) = FTL(I,j,K) + GRAD_FTL(I,j,K)
          FQW(I,j,K) = FQW(I,j,K) + GRAD_FQW(I,j,K)
!         !-------------------------------------------------------------
!         !  Add surface-drive gradient adjustment terms to fluxes
!         !  within the surface-based mixed layer.
!         !-------------------------------------------------------------
          IF (K  <=  NTML(I,j) ) THEN
            NON_GRAD_FTL(I,j,K) = RHOKHZ(I,j,K) * GRAD_T_ADJ(I,j)
            NON_GRAD_FQW(I,j,K) = RHOKHZ(I,j,K) * GRAD_Q_ADJ(I,j)
            FTL(I,j,K) = FTL(I,j,K) + NON_GRAD_FTL(I,j,K)
            FQW(I,j,K) = FQW(I,j,K) + NON_GRAD_FQW(I,j,K)
          ENDIF
        ENDDO
        ENDDO
      ENDDO
!-----------------------------------------------------------------------
!! 2. Lock and Whelan revised non-gradient formulation
!-----------------------------------------------------------------------
      IF (FLUX_GRAD  ==  LockWhelan2006) THEN
        DO K=2,BL_LEVELS
         do j=1,rows
         DO I=1,row_length

          IF ( K  <=  NTML(I,j) ) THEN
            F2_FTL(I,j,K)  = RHOF2(I,j,K)  * TOTHF_ZH(I,j)
            FSC_FTL(I,j,K) = RHOFSC(I,j,K) * TOTHF_ZH(I,j)
            FTL(I,j,K)   = FTL(I,j,K) + F2_FTL(I,j,K) + FSC_FTL(I,j,K)  &
     &                     - FT_NT(I,j,K)

            F2_FQW(I,j,K)  = RHOF2(I,j,K)  * TOTQF_ZH(I,j)
            FSC_FQW(I,j,K) = RHOFSC(I,j,K) * TOTQF_ZH(I,j)
            FQW(I,j,K)   = FQW(I,j,K) + F2_FQW(I,j,K) + FSC_FQW(I,j,K)  &
     &                     - FQ_NT(I,j,K)
          ENDIF

          IF ( K  <=  NTDSC(I,j) .AND. K >= NBDSC(I,j) ) THEN

            F2_FTL(I,j,K)  = RHOF2(I,j,K)  *                            &
     &                        ( TOTHF_ZHSC(I,j)-FT_NT_DSCB(I,j) )
            FSC_FTL(I,j,K) = RHOFSC(I,j,K) *                            &
     &                        ( TOTHF_ZHSC(I,j)-FT_NT_DSCB(I,j) )
            FTL(I,j,K)   = FTL(I,j,K) + F2_FTL(I,j,K) + FSC_FTL(I,j,K)  &
     &                     - ( FT_NT(I,j,K)-FT_NT_DSCB(I,j) )

            F2_FQW(I,j,K)  = RHOF2(I,j,K)  *                            &
     &                        ( TOTQF_ZHSC(I,j)-FQ_NT_DSCB(I,j) )
            FSC_FQW(I,j,K) = RHOFSC(I,j,K) *                            &
     &                        ( TOTQF_ZHSC(I,j)-FQ_NT_DSCB(I,j) )
            FQW(I,j,K)   = FQW(I,j,K) + F2_FQW(I,j,K) + FSC_FQW(I,j,K)  &
     &                     - ( FQ_NT(I,j,K)-FQ_NT_DSCB(I,j) )
          ENDIF
         ENDDO
         ENDDO
        ENDDO
      ENDIF   ! FLUX_GRAD
!

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXFLUXTQ',4)
      ENDIF

      RETURN
      END SUBROUTINE EX_FLUX_TQ
