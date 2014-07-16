#if defined(A71_1A)
! SUBROUTINE CALC_VIS_PROB-------------------------------------------
!
! Purpose: Calculates visibility probability,
!          The visibility probability is similar to the cloud fraction
!          except it records the fraction of a grid box with RH
!          greater than that required for the critical visibility
!          (e.g 1 km), taking into account precipitation.
!
!          Suitable for single-column use.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   27/07/01   Original code.  Pete Clark.
!   5.5   17/04/03   Remove reference to obsolete section
!                    C90_1A. T.White
!   6.2   03/02/06   Move to section A71_1A. P.Selwood
!
!
! Programming standard:  Unified Model Documentation Paper No 3,
!                        Version 5, dated 08/12/92.
!
! Documentation:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!
!    NOTE: New UM Doc Paper to be produced soon (S.Cusack 5/11/01)
!
!END----------------------------------------------------------------
!
!
!Arguments:---------------------------------------------------------
      SUBROUTINE CALC_VIS_PROB(                                         &
     &  PSTAR,RHCRIT,LEVELS,POINTS,PFIELD                               &
                                                  !INPUT
     & ,T,AEROSOL,L_MURK,Q,QCL,QCF,VIS,NVIS                             &
                                                  !INPUT
     & ,LCA,CCA,PCT                                                     &
                                                  !INPUT
     & ,Beta_LS_Rain, Beta_LS_Snow                                      &
                                                  !INPUT
     & ,Beta_C_Rain, Beta_C_Snow                                        &
                                                  !INPUT
     & ,Prob_of_Vis                                                     &
                                                  !OUTPUT
     & ,ERROR                                                           &
                                                  !OUTPUT
     & )
      IMPLICIT NONE
      INTEGER                                                           &
     & LEVELS                                                           &
                           ! IN No. of levels being processed.
     &,POINTS                                                           &
                           ! IN No. of gridpoints being processed.
     &,PFIELD                                                           &
                           ! IN No. of points in global field (at one
!                          !    vertical level).
     &,NVIS                ! IN No. of visibility thresholds
      REAL                                                              &
     & PSTAR(PFIELD)                                                    &
                           ! IN Surface pressure (Pa).
     &,RHCRIT(LEVELS)                                                   &
                           ! IN Critical relative humidity.  See the
!                          !    the paragraph incorporating eqs P292.11
!                          !    to P292.14; the values need to be tuned
!                          !    for the given set of levels.
     &,Q(PFIELD,LEVELS)                                                 &
                           ! IN Specific Humidity
!                          !    (kg per kg air).
     &,QCL(PFIELD,LEVELS)                                               &
                           ! Cloud liquid water content at
!                          !     processed levels (kg per kg air).
     &,QCF(PFIELD,LEVELS)                                               &
                           ! Cloud ice content at processed levels
!                          !    (kg per kg air).
     &,T(PFIELD,LEVELS)                                                 &
                           ! IN Temperature (K).
     &,AEROSOL(PFIELD,LEVELS)                                           &
                              ! IN Aerosol mixing ratio(ug/kg)
     &,VIS(NVIS)                                                        &
                              ! Visibility thresholds
     &,LCA(PFIELD)                                                      &
                                       ! IN Total Layer Cloud.
     &,CCA(PFIELD)                                                      &
                                       ! IN Convective Cloud.
     &,Beta_LS_Rain(PFIELD,LEVELS)                                      &
                                       ! IN Scattering in LS Rain.
     &,Beta_LS_Snow(PFIELD,LEVELS)                                      &
                                       ! IN Scattering in LS Snow.
     &,Beta_C_Rain(PFIELD,LEVELS)                                       &
                                       ! IN Scattering in Conv Rain
     &,Beta_C_Snow(PFIELD,LEVELS)      ! IN Scattering in Conv Snow
      LOGICAL                                                           &
     &   L_MURK               ! IN : Aerosol present

      LOGICAL                                                           &
     & PCT                              ! IN T:Cloud amounts are in %

      REAL                                                              &
     & Prob_of_Vis(PFIELD,LEVELS,NVIS) ! OUT Vis prob at processed level
!                          !     (decimal fraction).
      INTEGER ERROR        ! OUT 0 if OK; 1 if bad arguments.
!
!---------------------------------------------------------------------
!  Workspace usage----------------------------------------------------
      REAL                                                              &
                           ! "Automatic" arrays on Cray.
     & Vis_Threshold(PFIELD,LEVELS,NVIS)                                &
     &,Prob_of_Vis_LS(PFIELD,LEVELS,NVIS)                               &
     &,Prob_of_Vis_C(PFIELD,LEVELS,NVIS)                                &
     &,P_LSP(PFIELD)                                                    &
     &,P_CP(PFIELD)                                                     &
     &,Temp ! Temporary
!  External subroutine called ----------------------------------------
      EXTERNAL FOG_FR
!
! Parameters
      REAL large_distance
!
      PARAMETER(                                                        &
     & large_distance = 1.0e6                                           &
     &         )
!
! Local, including SAVE'd, storage------------------------------------
!
      INTEGER K,I,J     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.
!--------------------------------------------------------------------
!  Local and other physical constants----------------------------------
#include "c_pi.h"
#include "c_visbty.h"
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
!-----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
!-----------------------------------------------------------------------
      ERROR=0
      IF(POINTS >  PFIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF

      IF(PCT) THEN
        DO I=1,POINTS
          P_CP(I) =MAX(CCA(I)/100.0,0.0)
          P_LSP(I)=MAX((1.0-P_CP(I))*LCA(I)/100.0,0.0)
        ENDDO
      ELSE
        DO I=1,POINTS
          P_CP(I) =MAX(CCA(I),0.0)
          P_LSP(I)=MAX((1.0-P_CP(I))*LCA(I),0.0)
        ENDDO
      ENDIF


      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            Vis_Threshold(I,K,J)=VIS(J)
          ENDDO
        ENDDO
      ENDDO

! DEPENDS ON: fog_fr
      CALL FOG_FR(                                                      &
     & PSTAR,RHCRIT,LEVELS,PFIELD,                                      &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,Vis_Threshold,Prob_of_Vis,NVIS        &
     & )

      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            IF (P_LSP(I)  >   0.0) THEN
              Temp = 1.0/VIS(J)+                                        &
     &         (Beta_LS_Rain(I,K)+Beta_LS_SNOW(I,K))/LnLiminalContrast
              IF (Temp  >   0.0) THEN
                Vis_Threshold(I,K,J)=MAX(1.0/Temp , VIS(J))
              ELSE
                Vis_Threshold(I,K,J)=large_distance
              ENDIF
            ELSE
              Vis_Threshold(I,K,J)=VIS(J)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

! DEPENDS ON: fog_fr
      CALL FOG_FR(                                                      &
     & PSTAR,RHCRIT,LEVELS,PFIELD,                                      &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,Vis_Threshold,Prob_of_Vis_LS,NVIS     &
     & )

      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            IF (P_CP(I)  >   0.0) THEN
              Temp = 1.0/VIS(J)+                                        &
     &         (Beta_C_Rain(I,K)+Beta_C_SNOW(I,K))/LnLiminalContrast
              IF (Temp  >   0.0) THEN
                Vis_Threshold(I,K,J)=MAX(1.0/Temp , VIS(J))
              ELSE
                Vis_Threshold(I,K,J)=large_distance
              ENDIF
            ELSE
              Vis_Threshold(I,K,J)=VIS(J)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

! DEPENDS ON: fog_fr
      CALL FOG_FR(                                                      &
     & PSTAR,RHCRIT,LEVELS,PFIELD,                                      &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,Vis_Threshold,Prob_of_Vis_C,NVIS      &
     & )

      DO K=1,LEVELS
        DO I=1,POINTS
          DO J=1,NVIS
            Prob_of_Vis(I,K,J)=(1.0-P_CP(I)-P_LSP(I))*                  &
     &                            Prob_of_Vis(I,K,J) +                  &
     &                          P_LSP(I)*                               &
     &                            Prob_of_Vis_LS(I,K,J) +               &
     &                          P_CP(I)*                                &
     &                            Prob_of_Vis_C(I,K,J)
          ENDDO
        ENDDO
      ENDDO
!
 9999 Continue

      RETURN
      END SUBROUTINE CALC_VIS_PROB
#endif
