#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
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
#include "blopt8a.h"
#if defined(SCMA)
#include "s_scmop.h"
#endif

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

#if defined(SCMA)
      DO K=1,bl_levels
        DO j=1,rows
          DO I=1,row_length
            ! save local K
            RHOKH_DIAG(I,j,K) = RHOKH(I,j,K)
            RHOKM_DIAG(I,j,K) = RHOKM(I,j,K)
          END DO ! i
        END DO ! j
      END DO ! k
#endif

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
#if defined(SCMA)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
      IF (L_SCMDiags(SCMDiag_bl)) THEN
!
!     Output the diffusivities.
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokm,                                           &
             'momdif','Diffusivity of momentum','kg/(ms)',              &
             t_avg,d_bl,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokh,                                           &
             'htdiff','Diffusivity of heat','kg/(ms)',                  &
             t_avg,d_bl,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokm_diag,                                      &
             'KM_local','Diffusivity of momentum','kg/(ms)',            &
             t_avg,d_bl,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokh_diag,                                      &
             'KH_local','Diffusivity of heat','kg/(ms)',                &
             t_avg,d_bl,default_streams,'',RoutineName)
!
        DO j=1,rows
          DO I=1,row_length
            RHOKH_DIAG(I,j,1) = 0.0
            RHOKM_DIAG(I,j,1) = 0.0
            DO K=2,bl_levels
              RHOKH_DIAG(I,j,K) = RHOKHZ(I,j,K)
              RHOKM_DIAG(I,j,K) = RHOKMZ(I,j,K)
            END DO ! k
          END DO ! i
        END DO ! j
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokm_diag,                                      &
             'KM_surf','Diffusivity of momentum','kg/(ms)',             &
             t_avg,d_bl,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokh_diag,                                      &
             'KH_surf',                                                 &
             'Diffusivity of heat','kg/(ms)',                           &
             t_avg,d_bl,default_streams,'',RoutineName)
!
        DO K=2,bl_levels
          DO j=1,rows
            DO I=1,row_length
              RHOKH_DIAG(I,j,K) = RHOKH_TOP(I,j,K)
              RHOKM_DIAG(I,j,K) = RHOKM_TOP(I,j,K)
            END DO ! i
          END DO ! j
        END DO ! k
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokm_diag,                                      &
             'KM_top','Diffusivity of momentum','kg/(ms)',              &
             t_avg,d_bl,default_streams,'',RoutineName)
!
! DEPENDS ON: scmoutput
        Call SCMoutput(rhokh_diag,                                      &
             'KH_top','Diffusivity of heat','kg/(ms)',                  &
             t_avg,d_bl,default_streams,'',RoutineName)
!
      END IF ! L_SCMDiags(SCMDiag_bl)
#endif

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('KMKH    ',4)
      ENDIF

      RETURN
      END SUBROUTINE KMKH
#endif
