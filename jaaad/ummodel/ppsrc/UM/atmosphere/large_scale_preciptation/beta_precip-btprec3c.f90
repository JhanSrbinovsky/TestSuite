
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

       Real                                                             &
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
                        ! Ice mass-diameter relationships m(D) = ai D^bi
     &, lsp_ei, lsp_fi, lsp_eic, lsp_fic     
                  ! Ice particle Best-Reynolds number relationships
                  ! Re(D) =EI Be^FI

      LOGICAL                                                           &
     & PCT                                                              &
                                        ! IN T:Cloud amounts are in %
     &,AVG                                                              &
                                        ! IN T:Precip =local*prob
     &,l_psd
                                        ! IN Use generic ice size distbn

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
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

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
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
! C_LSPDRP start

      ! Microphysics parameters

      ! Drop size distribution for rain: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1R lambda^X2R  and m=X4R
!     REAL, PARAMETER :: X1R is set in the UMUI
!     REAL, PARAMETER :: X2R is set in the UMUI
!     REAL, PARAMETER :: X4R is set in the UMUI

      ! Drop size distribution for graupel: N(D) =  N0 D^m exp(-lambda D)
      ! where N0 = X1G lambda^X2G  and m=X4G
      REAL, PARAMETER :: X1G=5.E25
      REAL, PARAMETER :: X2G=-4.0
      REAL, PARAMETER :: X4G=2.5

      ! Particle size distribution for ice: N(D) = N0 D^m exp(-lambda D)
      ! where N0 = X1I TCG lambda^X2I, m=X4I and TCG=exp(- X3I T[deg C])
!     REAL, PARAMETER :: X1I is set in the UMUI
      REAL, PARAMETER :: X2I=0.0
      REAL, PARAMETER :: X3I=0.1222
      REAL, PARAMETER :: X4I=0.0
!     REAL, PARAMETER :: X1IC is set in the UMUI
      REAL, PARAMETER :: X2IC=0.0
      REAL, PARAMETER :: X3IC=0.1222
      REAL, PARAMETER :: X4IC=0.0

      ! Mass diameter relationship for graupel:  m(D) = AG D^BG
      REAL, PARAMETER :: AG=261.8
      REAL, PARAMETER :: BG=3.0

      ! Mass diameter relationship for ice:  m(D) = AI D^BI
      ! These are set in the UMUI (at 6.6). Values at 6.5 are listed below.
      ! Recommended values for the generic particle size distribution are
      ! from Brown and Francis and are
      ! AI=AIC=1.85E-2, BI=BIC=1.9. If l_calcfall is changed from .true.
      ! then the generic psd values should be set below to ci=cic=8.203
      ! and di=dic=0.2888
      ! REAL, PARAMETER :: AI=0.0444
      ! REAL, PARAMETER :: BI=2.1
      ! REAL, PARAMETER :: AIC=0.587
      ! REAL, PARAMETER :: BIC=2.45

      ! The area diameter relationships are only used if
      ! L_CALCFALL=.TRUE.
      ! Area diameter relationship for ice:  Area(D) = RI D^SI
      REAL, PARAMETER :: RI=0.131
      REAL, PARAMETER :: SI=1.88
      REAL, PARAMETER :: RIC=0.131
      REAL, PARAMETER :: SIC=1.88

      ! The Best/Reynolds relationships are only used if
      ! L_CALCFALL=.TRUE.
      ! Relationship between Best number and Reynolds number:
! Re(D) =LSP_EI(C) Be^LSP_FI(C)
      ! These values are set in the UMUI, but the default values are
      ! listed below. N.B. these have been renamed for VN7.3 with 
      ! 'LSP_' added from the previous versions to avoid conflicts 
      ! later in the code. 
      ! REAL, PARAMETER :: LSP_EI=0.2072  Set in the UMUI
      ! REAL, PARAMETER :: LSP_FI=0.638   Set in the UMUI
      ! REAL, PARAMETER :: LSP_EIC=0.2072 Set in the UMUI
      ! REAL, PARAMETER :: LSP_FIC=0.638  Set in the UMUI


      ! The fall speeds of ice particles are only used if
      ! L_CALCFALL=.FALSE.
      ! Fall speed diameter relationships for ice:
      ! vt(D) = CI D^DI
      REAL, PARAMETER :: CI0=14.3
      REAL, PARAMETER :: DI0=0.416
      REAL, PARAMETER :: CIC0=74.5
      REAL, PARAMETER :: DIC0=0.640

      ! Axial ratio (c-axis divided by a-axis) ESTIMATES. These are not
      ! consistent with those from the area diameter relationships.
      REAL, PARAMETER :: AR=1.0
      REAL, PARAMETER :: ARC=1.0

      ! Fall speed diameter relationship for rain: vt(D) = CR D^DR
      REAL, PARAMETER :: CR=386.8
      REAL, PARAMETER :: DR=0.67

      ! Fall speed diameter relationship for graupel: vt(D) = CG D^DG
      REAL, PARAMETER :: CG=253.0
      REAL, PARAMETER :: DG=0.734

      ! Do we wish to calculate the ice fall velocities?
      ! TRUE if calculate speeds, FALSE if specify speeds
      LOGICAL, PARAMETER :: L_CALCFALL=.TRUE.

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
     &,Inst_C_Snow(P_FIELD)                                             &
                                          ! Local Convective Snow
     &,lcai(p_field)                                                    &
                                          ! 1 / large-scale cloud (lca)
     &,m_si(p_field)
                          ! si moment of ice particle size distribution

! Local varables:-----------------------------------------------------
      REAL                                                              &
     & PFactor(P_FIELD)                                                 &
     &,GammaR                                                           &
     &,GammaR1                                                          &
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

      DO I=K1STPT,K1STPT+POINTS-1
        PFactor(I)=(PREF/Pressure(I))**0.4
      ENDDO

      PowerR=(3.0-X2R+X4R)/(DR+4.0-X2R+X4R)
      PowerI=(1.0+SI-X2I+X4I)/(BI+DI0+1.0-X2I+X4I)
      FactorR1=0.5*PI*X1R
      FactorR2=PI/6.0*RHO_WATER*CR*X1R
      FactorI1=2.0*RI*X1I
      FactorI2=AI*CI0*X1I

! DEPENDS ON: gammaf
      Call GAMMAF(X4R+DR+4.0,GammaR)
! DEPENDS ON: gammaf
      Call GAMMAF(X4R+3.0,GammaR1)
! DEPENDS ON: gammaf
      Call GAMMAF(X4I+BI+DI0+1.0,GammaI)
! DEPENDS ON: gammaf
      Call GAMMAF(X4I+SI+1.0,GammaI1)

      If (l_psd) then
        ! Find inverse of layer cloud amount
        Do i=1, p_field
          lcai(i)=1.0/max(lca(i),0.01)
        End do

        ! Use the generic ice particle size distribution to
        ! calculate the si moment of the ice particle size distribution
! DEPENDS ON: lsp_moments
        Call lsp_moments(p_field,rho,T,qcf,lcai,ai,bi,si,m_si)

      End if  ! l_psd

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
          Beta_LS_Rain(I)=FactorR1*GammaR1*(Inst_LS_Rain(I)/            &
     &                    (PFactor(I)*FactorR2*GammaR))**PowerR
        ELSE
          Beta_LS_Rain(I)=0.0
        ENDIF

        IF(Inst_LS_Snow(I)  >   SmallValue) THEN
          If (l_psd) then
            ! Use the generic ice particle size distribution
            Beta_LS_Snow(i)=2.0*ri*m_si(i)
          Else
            Beta_LS_Snow(I)=FactorI1*GammaI1*(Inst_LS_Snow(I)/          &
     &                      (PFactor(I)*FactorI2*GammaI))**PowerI
          End if  ! l_psd
        ELSE
          Beta_LS_Snow(I)=0.0
        ENDIF

        IF(Inst_C_Rain(I)  >   SmallValue) THEN
          Beta_C_Rain(I) =FactorR1*GammaR1*(Inst_C_Rain(I) /            &
     &                    (PFactor(I)*FactorR2*GammaR))**PowerR
        ELSE
          Beta_C_Rain(I) =0.0
        ENDIF

        IF(Inst_C_Snow(I)  >   SmallValue) THEN
          ! Use the 3C size distribution of snow for the convective
          ! contribution since we do not have the equivalent of qcf
          ! easily available.
          Beta_C_Snow(I) =FactorI1*GammaI1*(Inst_C_Snow(I) /            &
     &                    (PFactor(I)*FactorI2*GammaI))**PowerI
        ELSE
          Beta_C_Snow(I) =0.0
        ENDIF

      ENDDO

 9999 Continue

      RETURN
      END SUBROUTINE BETA_PRECIP
