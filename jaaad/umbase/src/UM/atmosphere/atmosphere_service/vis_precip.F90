#if defined(A71_1A)
!   5.5   17/04/03   Remove reference to obsolete section
!                    C90_1A. T.White
!  SUBROUTINE VIS_PRECIP ---------------------------------------------
!
!     PURPOSE:
! Process fields of precipitation intensity to give scattering coefft
! in 1/metres.
! Calculated at model level (eg bottom eta level 25m)
! or level within surface layer eg screen ht ( 1.5M )
!
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   27/07/01   Original code.  Pete Clark.
!   6.2   03/02/06   Moved to A71_1A. P.Selwood
!   6.2   22/08/05   Remove spaces from GOTOs. P.Selwood
!
!  Programming standard: U M Doc. Paper No. 4
!
!  Logical components covered :
!
!  Project task:
!
!  External documentation
!    Forecasting Research Scientific Paper NO.4
!    Diagnosis of visibility in the UK Met Office Mesoscale Model
!    and the use of a visibility analysis to constrain initial
!    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      NIMROD diagnostic:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!
!    NOTE: New UM Doc Paper to be produced soon (S.Cusack 5/11/01)
!
!END----------------------------------------------------------------
!
!  Arguments:-------------------------------------------------------
      SUBROUTINE VIS_PRECIP                                             &
     &           (Vis_No_Precip                                         &
                                                      !INPUT
     &           ,LCA,CCA,PCT                                           &
                                                      !INPUT
     &           ,Beta_LS_Rain, Beta_LS_Snow                            &
                                                      !INPUT
     &           ,Beta_C_Rain, Beta_C_Snow                              &
                                                      !INPUT
     &           ,P_FIELD,POINTS,K1STPT                                 &
                                                      !INPUT
     &           ,Vis_overall,Vis_LSP,Vis_CP                            &
                                                      !OUTPUT
     &           ,ERROR)                              !OUTPUT
      IMPLICIT NONE
!---------------------------------------------------------------------
! Workspace usage:----------------------------------------------------
! 3 real arrays of size P_FIELD
!---------------------------------------------------------------------
!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
      INTEGER                                                           &
     & P_FIELD                                                          &
                                        ! IN NO. points in field.
     &,POINTS                                                           &
                ! IN Number of gridpoints being processed.
     &,K1STPT                                                           &
                ! IN First gridpoint processed within complete field.
     &,ERROR    ! OUT Error code
      REAL                                                              &
     & Vis_No_Precip(P_FIELD)                                           &
                                        ! IN Vis outside precip.
     &,LCA(P_FIELD)                                                     &
                                        ! IN Total Layer Cloud.
     &,CCA(P_FIELD)                                                     &
                                        ! IN Convective Cloud.
     &,Beta_LS_Rain(P_FIELD)                                            &
                                        ! IN Scattering in LS Rain.
     &,Beta_LS_Snow(P_FIELD)                                            &
                                        ! IN Scattering in LS Snow.
     &,Beta_C_Rain(P_FIELD)                                             &
                                        ! IN Scattering in Conv Rain
     &,Beta_C_Snow(P_FIELD)             ! IN Scattering in Conv Snow

      LOGICAL                                                           &
     & PCT                              ! IN T:Cloud amounts are in %
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     & Vis_overall(P_FIELD)                                             &
                                       ! OUT Visibility overall
     &,Vis_LSP(P_FIELD)                                                 &
                                       ! OUT Visibility in LS Precip.
     &,Vis_CP(P_FIELD)                 ! OUT Visibility in Conv Precip.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Local varables:------------------------------------------------------
!  Define local variables ---------------------------------------------
      INTEGER I       ! Loop counters: I - horizontal field index;

      REAL                                                              &
     & Beta_No_Precip                                                   &
     &,P_LSP(P_FIELD)                                                   &
     &,P_CP(P_FIELD)
!---------------------------------------------------------------------
!  External subroutine called ----------------------------------------
!---------------------------------------------------------------------
! Local and other physical constants----------------------------------
#include "c_pi.h"
#include "c_visbty.h"
!
      ERROR=0
      IF((K1STPT+POINTS-1) >  P_FIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF
      IF(PCT) THEN
        DO I=K1STPT,K1STPT+POINTS-1
          P_CP(I)=CCA(I)/100.0
          P_LSP(I)=(1.0-P_CP(I))*LCA(I)/100.0
        ENDDO
      ELSE
        DO I=K1STPT,K1STPT+POINTS-1
          P_CP(I)=CCA(I)
          P_LSP(I)=(1.0-P_CP(I))*LCA(I)
        ENDDO
      ENDIF


      DO I=K1STPT,K1STPT+POINTS-1

        Beta_No_Precip=-LnLiminalContrast/Vis_No_Precip(I)

        IF(P_LSP(I)  >   0.0) THEN
          Vis_LSP(I) = -LnLiminalContrast /                             &
     &      (Beta_No_Precip +                                           &
     &       Beta_LS_Rain(I) + Beta_LS_Snow(I))
        ELSE
          Vis_LSP(I)=Vis_No_Precip(I)
        ENDIF

        IF(P_CP(I)  >   0.0) THEN
          Vis_CP(I) = -LnLiminalContrast /                              &
     &      (Beta_No_Precip +                                           &
     &       Beta_C_Rain(I) + Beta_C_Snow(I))
        ELSE
          Vis_CP(I)=Vis_No_Precip(I)
        ENDIF

! Ensure no rounding problems lead to vis > vis in clear air
        Vis_overall(I) = MIN((1.0-P_CP(I)-P_LSP(I))*Vis_No_Precip(I) +  &
     &                   P_LSP(I)*Vis_LSP(I) +P_CP(I)*Vis_CP(I),        &
     &                   Vis_No_Precip(I))

      ENDDO

 9999 Continue

      RETURN
      END SUBROUTINE VIS_PRECIP
#endif
