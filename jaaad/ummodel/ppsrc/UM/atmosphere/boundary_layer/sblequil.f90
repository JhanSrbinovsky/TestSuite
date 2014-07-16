
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE SBLequil-----------------------------------------------
!!!
!!!  Purpose: Estimate full and partial diffusivities
!!!           in an equilibrium Stable Boundary Layer
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  5.3     3/10/01  New deck  A.G. Williams
!!!  5.4    30/08/02  Revised Equilibrium SBL scheme (EqSBL v.2).
!!!                   A.G.Williams
!!!  5.5    26/02/03  Revised Equilibrium SBL scheme (EqSBL v.2a)
!!!                   A.G.Williams/A.Lock
!!!  6.2    22/08/05  Add spaces to DO WHILE. P.Selwood
!  6.2   30/01/06  Revised constants giving more realistic mixing
!                  length.       Bob Beare
!!!
!!!  Programming standard: Unified Model Documentation Paper No 3
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------
!! Abbreviations used in variable descriptions:
!! --------------------------------------------
!! z:      height
!! MO_loc: local Monin-Obukhov length
!! ll:     turbulence mixing length
!! zet:    ll/MO_loc (greek symbol zeta)
!! kappa:  von Karman constant
!! us_loc: local u_star  estimate (|uw|^0.5)
!! WT_loc: local <w*THv> estimate
!! KM/KH:  diffusivities for momentum and heat
!!       [For X=U,THv: <w*X>=-KX*(dX/dz)]
!! KUU/KTT/KUT/KTU: partial diffusivities.
!!       [For X=U,THv: <w*X>=-KXU*(dU/dz)-KXT*(dTHv/dz)]
!! PHI_KX: Normalized diffusivities
!!----------------------------------------------------------------------
      SUBROUTINE SBLequil (                                             &
     &  Zhat                                                            &
     & ,zetkz,PKMzet,PKHzet,PKUU,PKTT,PKUT,PKTU                         &
     & ,PHIe,PHIww,Ri,CE,rpow,rCB,CN,iERR,LTIMER                        &
     &  )

      IMPLICIT NONE
      LOGICAL LTIMER

!  Arguments
      REAL                                                              &
     & Zhat                                                             &
                 ! IN  z/MO_loc (>0.0)
     &,zetkz                                                            &
                 ! OUT l_hat=zeta/(kappa*Zhat)=ll/(kappa*z)
     &,PKMzet                                                           &
                 ! OUT PHI_KM/zet=KM/(ll*us_loc)
     &,PKHzet                                                           &
                 ! OUT PHI_KH/zet=KH/(ll*us_loc)
     &,PKUU                                                             &
                 ! OUT KUU/(ll*us_loc)
     &,PKTT                                                             &
                 ! OUT KTT/(ll*us_loc)
     &,PKUT                                                             &
                 ! OUT KUT/(ll^2*g/Thv)
     &,PKTU                                                             &
                 ! OUT KTU/(ll^2*(g/Thv)*-1.0*(WT_loc/us_loc^2)^2)
     &,PHIe                                                             &
                 ! OUT TKE/stress ratio
     &,PHIww                                                            &
                 ! OUT w-variance/stress ratio
     &,Ri                                                               &
                 ! OUT Equilib. SBL Richardson no.
     &,CE                                                               &
                 ! OUT Turbulent dissipation timescl con
     &,rpow                                                             &
                 ! OUT Copy of parameter npow (see below)
     &,rCB                                                              &
                 ! OUT Copy of parameter CB (see below)
     &,CN        ! OUT Constant in length scale formula

      INTEGER                                                           &
     & iERR      ! OUT Error status. Five-digit integer "MLKJI",
!                !     where each digit has the value 0 (no-error)
!                !     or 1 (error). Errors detected are as follows:
!                !       I: Input Zhat < 0 (set to absolute value)
!                !       J: Outer (zetkz,A1) loop exceeded max iteration
!                !       K: Inner (PHIww) loop exceeded max iterations
!                !          at least once
!                !       L: Error in Newton's method for calculating
!                !          PHIww (under/over-flow or wrong root found).
!                !          Set PHIww to neutral value and continue.
!                !       M: Either or both of PKHzet and PKMzet are too
!                !          big/small, or negative. Return neutral soln.

!  Externals
      EXTERNAL TIMER

!  Symbolic constants
! C_VKMAN start
      REAL,PARAMETER:: VKMAN=0.4 ! Von Karman's constant
! C_VKMAN end

!  Tunable parameters
      REAL npow,gamma,C1,C3,CB
      PARAMETER (                                                       &
     & npow=1.0                                                         &
                         !1/power in mixing-length formula
     &,gamma=4.0                                                        &
                         !Constant in 1st guess mix-len form
     &,C1=-0.002                                                        &
                         !Constant in param. of press-stress cov
     &,C3= -0.5                                                         &
                        !Constant in param. of press-stress cov
     &,CB= 1.05                                                         &
                         !Constant in mixing length formula
     &)

!  Matching conditions for homogeneous shear (HS) layers
      REAL PeHSn,PeHSw,PuHSn,PuHSw,PvHSn,PvHSw,PwHSn,PwHSw
      REAL PrHSn,PrHSw,PtuHSw,PttHSw
      PARAMETER (                                                       &
     & PeHSn= 1.0/0.170                                                 &
                          !PHI_e (no wall)
     &,PeHSw=  1.0/0.0862                                               &
                           !PHI_e (wall)
     &,PuHSn= 0.480*PeHSn                                               &
                          !PHI_u (no wall)
     &,PuHSw=  0.460*PeHSw                                              &
                           !PHI_u (wall)
     &,PvHSn= 0.260*PeHSn                                               &
                          !PHI_v (no wall)
     &,PvHSw=  0.340*PeHSw                                              &
                           !PHI_v (wall)
     &,PwHSn= 0.260*PeHSn                                               &
                          !PHI_w (no wall)
     &,PwHSw=  0.200*PeHSw                                              &
                           !PHI_w (wall)
     &,PrHSn= 0.67                                                      &
                          !Prandtl no. (no wall)
     &,PrHSw=  0.85                                                     &
                           !Prandtl no. (wall)
     &,PtuHSw= -2.1                                                     &
                           !PHI_ut (wall)
     &,PttHSw=  7.7                                                     &
                           !PHI_tt (wall)

     &)

!  Parameters for iteration loops
      REAL epsZet,epsA1,epsPw,vsml,vbig,PwLim
      INTEGER maxit
      PARAMETER (                                                       &
     & epsZet =1.0e-3                                                   &
                          !Tolerance for convergence of zeta
     &,epsA1  =1.0e-4                                                   &
                          !Tolerance for convergence of A1
     &,epsPw  =1.0e-2                                                   &
                          !Tolerance for convergence of PHIww
     &,maxit  =100                                                      &
                          !Maximum permissible number of iterations
     &,vsml   =1.0e-20                                                  &
                          !Limit for small numbers
     &,vbig   =1.0e+20                                                  &
                          !Limit for big numbers
     &,PwLim  =1.0                                                      &
                          !Lower limit in PHIw iteration
                          !(check for convergence towards wrong root)
     &)

!  Local storage

!  Dependent constants
      REAL                                                              &
     & C2                                                               &
                   !Constant in pressure-stress cov
     &,CC                                                               &
                   !Constant in pressure-stress cov
     &,CCn                                                              &
                   !no-wall value for CC
     &,D1                                                               &
                   !Constant in wall-effect terms
     &,D2                                                               &
                   !Constant in wall-effect terms
     &,D1T                                                              &
                   !Constant in wall-effect terms
     &,DD                                                               &
                   !Constant in param. of press-temp cov
     &,A1                                                               &
                   !Constant in param. of press-temp cov
     &,A2                                                               &
                   !Constant in param. of press-temp cov
     &,CT                                                               &
                   !Temperature variance  timescl con
     &,H2                                                               &
                   !Constant in PHI_KM expression
     &,EWS,EUS,EVS                                                      &
                   !Intermediate constants
     &,EWB,EUB,EVB                                                      &
                   !Intermediate constants
     &,EWSn                                                             &
                   !no-wall value for EWS
     &,F1,F2,F3,F4                                                      &
                   !Intermediate constants
     &,F1n,F3n,F4n                                                      &
                   !no-wall values for F1, F3, F4
     &,H1                                                               &
                   !Constant in PHI_KH expression
     &,M1,M2                                                            &
                         !Constants in PHIww equation
     &,N0,N2,N3,N5,N6,N9 !Constants in PHIww equation

! Other local quantities
      REAL                                                              &
     & PHItt,PHIut                                                      &
                     !Norm T-var & horiz heat flx
     &,kc                                                               &
                     !kappa*CE
     &,kz,zeta                                                          &
                     !kappa*Zhat, zeta
     &,fwt                                                              &
                     !Wall-effect weighting fn, f=ll/(kappa*z)
     &,N_zetkz,N_PKHzet,N_PKMzet,N_PHIe,N_PHIww !Neutral quantities

! Temporary variables
      REAL                                                              &
     & rCEsPw                                                           &
                         !1.0/(CE*SQRT(PHIww))
     &,tmp,tmp2                                                         &
                         !temporary computational variables
     &,zetprev,zetchange                                                &
                         !temp values for zetkz
     &,A1prev,A1change                                                  &
                         !temp values for A1
     &,PWprev,PWchange                                                  &
                         !temp values for PHIww
     &,yy,py,dpy         !temp variables in PHIww-iteration

      INTEGER                                                           &
     & iloop,jloop                                                      &
                         !loop counters
     &,NewtErr,KxErr                                                    &
                         !error flags
     &,iE1,iE2,iE3       !temporary error indicators

!-----------------------------------------------------------------------
!! Start of code
!-----------------------------------------------------------------------
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SBLEQUIL ',103)
      ENDIF
      rpow=npow
      rCB =CB
      iERR=0
      iE1 =0
      iE2 =0
      iE3 =0
      IF (Zhat <  0) iERR=1

!-----------------------------------------------------------------------
!! 0. Calculate dependent constants (C2,CC,D1,D2,DD,D1T,A2,CE,CT,H2)
!!    and other useful local quantities
!-----------------------------------------------------------------------
      tmp=2.0*(PwHSn-2.0/(PuHSn-PvHSn))/PeHSn
      C2=1.0-C1/tmp
      CCn=2.0*PwHSn*(1.0-C2)/(PuHSn-PvHSn)
      CC =2.0*PwHSw*(1.0-C2)/(PuHSw-PvHSw)
      tmp=PwHSw/(3.0*PwHSw*PwHSw+2.0)
      D1=tmp*(2.0*(1.0-C2)*(2.0*PwHSw*PuHSw+PwHSw*PvHSw                 &
     &   -3.0*PwHSw*PwHSw-4.0)/(3.0*(PuHSw-PvHSw))-2.0*C1*PeHSw/3.0)
      tmp=C1-(2.0*(1.0-C2)*(PwHSw-2.0/(PuHSw-PvHSw))/PeHSw)             &
     &   +3.0*D1/(PeHSw*PwHSw)
      D2=tmp*PeHSw/(3.0*C2*PwHSw)
      EWSn=CCn/(CCn+2.0*(1.0-C2)/3.0)
      F1n=CCn
      F3n=1.0-C2
      F4n=SQRT((F3n-(3.0*C1/(2.0*EWSn)))/F1n)
      EWS=CC/(CC+2.0*(1.0-C2)/3.0+2.0*(2.0*C2*D2/3.0+D1))
      F1=CC+1.5*D1
      F3=1.0-C2+1.5*D2*C2
      F4=SQRT((F3-(3.0*C1/(2.0*EWS)))/F1)
      DD=PrHSn/(F4n*F4n)
      D1T=PrHSw/(F4*F4)-DD
      F2=DD+D1T
      A2=-1.0-((F2+PtuHSw*PeHSw*DD*EWS/3.0)*(F3-1.5*C1/EWS)/F1)
      CE=(PwHSw)**(-1.5)
      CT=3.0*F2/(PttHSw*PeHSw*EWS)
      H2=(1.0+C3)/DD
      CN=(1.0-((CB*F4/CE)**(-1.0/npow)))**(-npow)

      kc=VKMAN*CE
      kz=VKMAN*Zhat

!-----------------------------------------------------------------------
!! Neutral solution
!-----------------------------------------------------------------------
      N_zetkz = 1.0
      N_PKHzet= 1.0/(F2*CE*SQRT(F4))
      N_PKMzet= 1.0
      N_PHIww = 1.0/F4
      N_PHIe  = 3.0*N_PHIww/EWS

!-----------------------------------------------------------------------
!! 1. Loop to get zetkz (and A1)
!-----------------------------------------------------------------------
      zetprev=-vbig
      A1prev =-vbig
      !first guesses for zetkz and A1
      zetkz=1.0/((1.0+((kz*gamma)**(1.0/npow)))**(npow))
      A1=0.5
      zetchange=ABS(zetprev-zetkz)
      A1change =ABS(A1prev -A1)
      iloop=0
      jloop=0

      DO WHILE (                                                        &
     &          ((zetchange >  epszet).OR.(A1change >  epsA1))          &
     &          .AND.(iloop <  maxit)                                   &
     &        )
        iloop=iloop+1
  !--------------------------------------------------
  !! 1.1 Parametric functions requiring zetkz (and A1)
  !--------------------------------------------------
        H1=(1.0-A1)/CT
        fwt=zetkz
        tmp=CC+2.0*(1.0-C2)/3.0+fwt*2.0*(2.0*C2*D2/3.0+D1)
        EWS=CC/tmp
        EWB=(C2*(1.0-fwt*2.0*D2)-3.0-2.0*C3)/tmp
        EUS=1.0+EWS*2.0*(2.0*(1.0-C2)+fwt*(C2*D2+1.5*D1))/(3.0*CC)
        EUB=(3.0-C2*(2.0-fwt*D2)+C3)/CC
        EVS=1.0+EWS*2.0*(C2-1.0+fwt*(C2*D2+1.5*D1))/(3.0*CC)
        EVB=(C2*(1.0+fwt*D2)+C3)/CC
        F1=CC+fwt*1.5*D1
        F2=DD+fwt*D1T
        F3=1.0-C2+fwt*1.5*D2*C2
        F4=SQRT((F3-(3.0*C1/(2.0*EWS)))/F1)
        zeta=zetkz*kz

  !--------------------------------------------------
  !! 1.2 Calculate PHIw, given zetkz
  !--------------------------------------------------
        !Set temporary constants
        M1=F1*F4*F4
        M2=C1*EWB/EWS-H2*(1.0+a2)
        N0=zeta*zeta*zeta*(-H1*M2/(VKMAN*kc))
        N2=zeta*(CE*(H1*F1-H2*F2))
        N3=zeta*zeta*(((1.0-H1)*M2-H1*M1)/VKMAN)
        N5=-kc*CE*F1
        N6=zeta*(CE*(M2+(1.0-H1)*M1))
        N9=kc*CE*M1

        !Iteration starting with neutral PHIw (Newton's method)
        PHIww=N_PHIww+0.5
        PWprev=-vbig
        PWchange=abs(PWprev-PHIww)
        jloop=0
        NewtErr=0
        DO WHILE ((PWchange >  epsPw).AND.(jloop <  maxit))
          jloop=jloop+1
          PWprev=PHIww
          yy=SQRT(PHIww)
          py= N9*(yy**(9.0))+N6*(yy**(6.0))+N5*(yy**(5.0))              &
     &       +N3*(yy**(3.0))+N2*(yy**(2.0))+N0
          dpy= 9.0*N9*(yy**(8.0))+6.0*N6*(yy**(5.0))                    &
     &        +5.0*N5*(yy**(4.0))+3.0*N3*(yy**(2.0))+2.0*N2*yy
          IF (ABS(dpy) <  vsml) THEN
            PHIww=PWprev
            NewtErr=1
          ELSE
            yy=yy-py/dpy
            PHIww=yy*yy
            IF ((PHIww <  PwLim).AND.(iloop >  1)) THEN
              PHIww=PWprev
              NewtErr=1
            ENDIF
          ENDIF
          PWchange=ABS(PWprev-PHIww)
        ENDDO !WHILE loop to get PHIe
        IF (jloop >= maxit) iE1=1
        IF (NewtErr == 1) THEN
          PHIww=N_PHIww
          iE2=1
        ENDIF

  !--------------------------------------------------
  !! 1.3 Compute other dimensionless functions
  !--------------------------------------------------
        KxErr=0
        PHIe=3.0*(PHIww-zeta*2.0*EWB/(3.0*kc*SQRT(PHIww)))/EWS
        PKHzet=(SQRT(PHIww)-zeta*H1/(kc*PHIww))/(F2*CE)
        IF (ABS(PKHzet) <  vsml) PKHzet=vsml
        IF (ABS(PKHZet) >  vbig) PKHzet=vbig
        IF ((PKHzet <  vsml).AND.(iloop == 1)) PKHzet=1.0

        IF (PKHzet <  vsml) THEN
          KxErr=1
        ELSE
          tmp=kc*F1*F4*F4*(PHIww**(1.5))                                &
     &       +zeta*(C1*EWB/EWS-H2*(1.0+a2))
          PKMzet=tmp/(kc*CE*F1*PHIww+zeta*H2/PKHzet)
          IF (ABS(PKMzet) <  vsml) PKMzet=vsml
          IF (ABS(PKMZet) >  vbig) PKMzet=vbig
          IF ((PKMzet <  vsml).AND.(iloop == 1)) PKMzet=1.0
          IF (PKMzet <  vsml) KxErr=1
        ENDIF

  !--------------------------------------------------
  !! 1.4 Recalculate A1 & zetkz; return to top of loop
  !--------------------------------------------------
        zetprev=zetkz
        A1prev =A1
        IF (KxErr == 0) THEN
          Ri=0.0
          Ri=zeta*PKMzet*PKMzet/(VKMAN*PKHzet)
          A1=0.5+1.5*Ri*Ri-Ri*Ri*Ri
          IF (Ri <= 0.0)  A1=0.5
          IF (Ri >  1.0) THEN
            A1=1.0
            zetkz=1.0/((                                                &
     &            ((CN                             )**(-1.0/npow))      &
     &           +((Zhat/(zetkz*PKHzet*CB*CB*PHIww))**( 0.5/npow))      &
     &                  )**(npow))
          ELSE
            zetkz=1.0/((                                                &
     &            ((CN                         )**(-1.0/npow))          &
     &           +((PKMzet*CB*SQRT(PHIww)*zetkz)**(-1.0/npow))          &
     &                  )**(npow))
          ENDIF
        ELSE
          zetkz=zetprev
          A1   =A1prev
        ENDIF
        zetchange=ABS(zetprev-zetkz)
        A1change =ABS(A1prev -A1)
      ENDDO !WHILE main (zetkz,A1) loop

      IF (KxErr == 1) THEN !return neutral soln
        zetkz =N_zetkz
        PKHzet=N_PKHzet
        PKMzet=N_PKMzet
        PHIe  =N_PHIe
        PHIww =N_PHIww
        iE3=1
      ENDIF

      IF (iloop >= maxit) iERR=iERR+10
      IF (iE1 /= 0)       iERR=iERR+100
      IF (iE2 /= 0)       iERR=iERR+1000
      IF (iE3 /= 0)       iERR=iERR+10000

!-----------------------------------------------------------------------
!! 2. Calculate residual output quantities
!-----------------------------------------------------------------------
      rCEsPw=1.0/(CE*SQRT(PHIww))
      PHItt =rCEsPw/(CT*PKHzet)
      tmp   =1.0/PKHzet+(1.0+A2)/PKMzet
      PHIut =-rCEsPw*tmp/DD

      PKUU  = (F3*PHIww-C1*PHIe/2.0)*rCEsPw/F1
      PKUT  = (1.0+C3)*PHIut*PKHzet *rCEsPw/F1
      PKTT  = PHIww                 *rCEsPw/F2
      PKTU  = (1.0-A1)*PHItt*PKMzet *rCEsPw/F2

!-----------------------------------------------------------------------
!! Finish up
!-----------------------------------------------------------------------
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SBLEQUIL ',104)
      ENDIF

      RETURN
      END SUBROUTINE SBLequil
