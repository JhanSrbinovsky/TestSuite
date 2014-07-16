#if defined(A04_3B)
!  SUBROUTINE BETA_PRECIP --------------------------------------------
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
!   6.2   22/08/05   Add spaces to GOTOs. P.Selwood
!   6.2   17/11/05   Add in microphysics variables. Damian Wilson
!
!  6.2   28/02/05  Modifications for STPH_RP
!                  A. Arribas
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
      SUBROUTINE BETA_PRECIP                                            &
     &           (LS_Rain, LS_Snow, C_Rain, C_Snow, qcf                 &
                                                      !INPUT
     &           ,rho, T, Pressure                                      &
                                                      !INPUT
     &           ,LCA,CCA,PCT,AVG                                       &
                                                      !INPUT
     &           ,P_FIELD,POINTS,K1STPT                                 &
                                                      !INPUT
     &           ,x1i,x1ic,x1r,x2r,x4r,l_psd,ai,bi,aic,bic              &
                                                      !INPUT
     &           ,lsp_ei,lsp_fi,lsp_eic,lsp_fic                         &
                                                      !INPUT
     &           ,Beta_LS_Rain, Beta_LS_Snow                            &
                                                      !OUTPUT
     &           ,Beta_C_Rain, Beta_C_Snow,ERROR)     !OUTPUT
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
     & LS_Rain(P_FIELD)                                                 &
                                        ! IN Large scale Rain
     &,LS_Snow(P_FIELD)                                                 &
                                        ! IN Large scale Snow
     &,C_Rain(P_FIELD)                                                  &
                                        ! IN Convective Rain
     &,C_Snow(P_FIELD)                                                  &
                                        ! IN Convective Snow
     &,qcf(p_field)                                                     &
                                        ! IN large-scale ice / kg kg-1
     &,rho(p_field)                                                     &
                                        ! IN Air density / kg m-3
     &,T(p_field)                                                       &
                                        ! IN Temperature / K
     &,Pressure(P_FIELD)                                                &
                                        ! IN Pressure
     &,LCA(P_FIELD)                                                     &
                                        ! IN Total Layer Cloud.
     &,CCA(P_FIELD)                     ! IN Convective Cloud.


        Real                                                            &
                        !, INTENT(IN)
     &  x1i                                                             &
                        ! Intercept in aggregate size distribution
     &, x1ic                                                            &
                        ! Intercept in crystal size distribution
     &, x1r                                                             &
                        ! Intercept in raindrop size distribution
     &, x2r                                                             &
                        ! Scaling parameter in raindrop size distribn.
     &, x4r                                                             &
                        ! Shape parameter in raindrop size distribution
     &, ai, bi, aic, bic                                                &
                        ! Ice mass diameter relationship m(D)=ai D^bi
     &, lsp_ei, lsp_fi, lsp_eic, lsp_fic
                        ! Ice particle Best-Reynolds number 
                        ! relationships: Re(D) =EI Be^FI

      LOGICAL                                                           &
     & PCT                                                              &
                                        ! IN T:Cloud amounts are in %
     &,AVG                                                              &
                                        ! IN T:Precip =local*prob
     &,l_psd
                                        ! Dummy for 3B scheme

!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     & Beta_LS_Rain(P_FIELD)                                            &
                                        ! OUT Scattering in LS Rain.
     &,Beta_LS_Snow(P_FIELD)                                            &
                                        ! OUT Scattering in LS Snow.
     &,Beta_C_Rain(P_FIELD)                                             &
                                        ! OUT Scattering in Conv Rain
     &,Beta_C_Snow(P_FIELD)             ! OUT Scattering in Conv Snow
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!  External subroutine called ----------------------------------------
      EXTERNAL GAMMAF
!---------------------------------------------------------------------
! PI needed to set new constants.
!---------------------------------------------------------------------
#include "c_pi.h"
#include "c_r_cp.h"
#include "c_densty.h"
#include "c_lspdrp.h"

      REAL                                                              &
     & PowerR                                                           &
     &,PowerI                                                           &
     &,FactorR1                                                         &
     &,FactorR2                                                         &
     &,FactorI1                                                         &
     &,FactorI2
      REAL                                                              &
     & Inst_LS_Rain(P_FIELD)                                            &
                                          ! Local Large scale Rain
     &,Inst_LS_Snow(P_FIELD)                                            &
                                          ! Local Large scale Snow
     &,Inst_C_Rain(P_FIELD)                                             &
                                          ! Local Convective Rain
     &,Inst_C_Snow(P_FIELD)               ! Local Convective Snow


! Local varables:-----------------------------------------------------
      REAL                                                              &
     & PFactor(P_FIELD)                                                 &
     &,GammaR                                                           &
     &,GammaI                                                           &
     &,GammaI1                                                          &
     &,SmallValue
      PARAMETER(SmallValue=1.0E-7)

!-----------------------------------------------------------------------
!  Define local variables ----------------------------------------------
      INTEGER I       ! Loop counters: I - horizontal field index;
!
      ERROR=0
      IF((K1STPT+POINTS-1) >  P_FIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF

      PowerR=(3.0-X2R)/(DR+4.0-X2R)
      PowerI=(3.0-X2I+X4I)/(BI+DI+1.0-X2I+X4I)
      FactorR1=PI*X1R
      FactorR2=PI/6.0*RHO_WATER*CR*X1R
      FactorI1=0.5*PI*X1I
      FactorI2=AI*CI*X1I
      DO I=K1STPT,K1STPT+POINTS-1
        PFactor(I)=(PREF/Pressure(I))**0.4
      ENDDO

! DEPENDS ON: gammaf
      Call GAMMAF(DR+4.0,GammaR)
! DEPENDS ON: gammaf
      Call GAMMAF(X4I+BI+DI+1.0,GammaI)
! DEPENDS ON: gammaf
      Call GAMMAF(X4I+3.0,GammaI1)

      IF (AVG) THEN

        IF (PCT) THEN

          DO I=K1STPT,K1STPT+POINTS-1

            IF(LCA(I)  >   SmallValue .and. CCA(I)  <   100.0) THEN
              Inst_LS_Rain(I) = 10000.0*LS_Rain(I)/(100.0-CCA(I))/LCA(I)
              Inst_LS_Snow(I) = 10000.0*LS_Snow(I)/(100.0-CCA(I))/LCA(I)
            ELSE
              Inst_LS_Rain(I) = 0.0
              Inst_LS_Snow(I) = 0.0
            ENDIF
            IF(CCA(I)  >   SmallValue) THEN
              Inst_C_Rain(I) = 100.0*C_Rain(I)/CCA(I)
              Inst_C_Snow(I) = 100.0*C_Snow(I)/CCA(I)
            ELSE
              Inst_C_Rain(I) = 0.0
              Inst_C_Snow(I) = 0.0
            ENDIF

          ENDDO

        ELSE

          DO I=K1STPT,K1STPT+POINTS-1

            IF(LCA(I)  >   SmallValue .and. CCA(I)  <   1.0) THEN
              Inst_LS_Rain(I) = LS_Rain(I)/(1.0-CCA(I))/LCA(I)
              Inst_LS_Snow(I) = LS_Snow(I)/(1.0-CCA(I))/LCA(I)
            ELSE
              Inst_LS_Rain(I) = 0.0
              Inst_LS_Snow(I) = 0.0
            ENDIF
            IF(CCA(I)  >   SmallValue) THEN
              Inst_C_Rain(I) = C_Rain(I)/CCA(I)
              Inst_C_Snow(I) = C_Snow(I)/CCA(I)
            ELSE
              Inst_C_Rain(I) = 0.0
              Inst_C_Snow(I) = 0.0
            ENDIF

          ENDDO

        ENDIF

      ELSE

        DO I=K1STPT,K1STPT+POINTS-1
          Inst_LS_Rain(I) = LS_Rain(I)
          Inst_LS_Snow(I) = LS_Snow(I)
          Inst_C_Rain(I)  = C_Rain(I)
          Inst_C_Snow(I)  = C_Snow(I)
        ENDDO

      ENDIF

      DO I=K1STPT,K1STPT+POINTS-1

        IF(Inst_LS_Rain(I)  >   SmallValue) THEN
          Beta_LS_Rain(I)=FactorR1*(Inst_LS_Rain(I)/                    &
     &                    (PFactor(I)*FactorR2*GammaR))**PowerR
        ELSE
          Beta_LS_Rain(I)=0.0
        ENDIF

        IF(Inst_LS_Snow(I)  >   SmallValue) THEN
          Beta_LS_Snow(I)=FactorI1*GammaI1*(Inst_LS_Snow(I)/            &
     &                    (PFactor(I)*FactorI2*GammaI))**PowerI
        ELSE
          Beta_LS_Snow(I)=0.0
        ENDIF

        IF(Inst_C_Rain(I)  >   SmallValue) THEN
          Beta_C_Rain(I) =FactorR1*(Inst_C_Rain(I) /                    &
     &                    (PFactor(I)*FactorR2*GammaR))**PowerR
        ELSE
          Beta_C_Rain(I) =0.0
        ENDIF

        IF(Inst_C_Snow(I)  >   SmallValue) THEN
          Beta_C_Snow(I) =FactorI1*GammaI1*(Inst_C_Snow(I) /            &
     &                    (PFactor(I)*FactorI2*GammaI))**PowerI
        ELSE
          Beta_C_Snow(I) =0.0
        ENDIF

      ENDDO

 9999 Continue

      RETURN
      END SUBROUTINE BETA_PRECIP
#endif
