#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
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

#if defined(SCMA)
! INOUT scmop is declared here
#include "s_scmop.h"
#endif
      INTEGER                                                           &
     &  I, J, K

#include "blopt8a.h"
#include "c_g.h"
#include "c_r_cp.h"

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
#if defined(SCMA)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      IF (L_SCMDiags(SCMDiag_bl)) THEN
!
        DO K=1,bl_levels
          DO j=1,rows
            DO I=1,row_length
              GRAD_FTL(I,j,K)= CP * GRAD_FTL(I,j,K)
              NON_GRAD_FTL(I,j,K) = CP * NON_GRAD_FTL(I,j,K)
            END DO ! I
          END DO ! j
        END DO ! K
!
! DEPENDS ON: scmoutput
       Call SCMoutput(GRAD_FTL,                                         &
            'Grad_ftl','Down gradient flux of TL','W/m2',               &
            t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(NON_GRAD_FTL,                                     &
            'Surf_ng_ftl','Surface driven non-gradient flux of TL',     &
            'W/m2',                                                     &
            t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(GRAD_FQW,                                         &
            'Grad_fqw','Down-gradient flux of QW','kg/m2/s',            &
            t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(NON_GRAD_FQW,                                     &
            'Surf_ng_fqw','Surface driven non-gradient flux of QW',     &
            'kg/m2/s',                                                  &
            t_avg,d_bl,default_streams,'',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(F2_FTL,                                           &
            'f2_ftl','FTL: f2 term','W/m2',                             &
            t_avg, d_bl, default_streams, '',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(F2_FQW,                                           &
            'f2_fqw','FQW: f2 term','kg/m2/s',                          &
            t_avg, d_bl, default_streams, '',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(FSC_FTL,                                          &
            'fSc_ftl','FTL: fSc term','W/m2',                           &
            t_avg, d_bl, default_streams, '',RoutineName)

! DEPENDS ON: scmoutput
       Call SCMoutput(FSC_FQW,                                          &
            'fSc_fqw','FQW: fSc term','kg/m2/s',                        &
            t_avg, d_bl, default_streams, '',RoutineName)
!
      END IF ! L_SCMDiags(SCMDiag_bl)
#endif

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EXFLUXTQ',4)
      ENDIF

      RETURN
      END SUBROUTINE EX_FLUX_TQ
#endif
