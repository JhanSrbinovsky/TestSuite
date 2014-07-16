#if defined(A04_3B)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE LSPCON-----------------------------------------
!   PURPOSE: CALCULATES CONSTANTS USED IN PRECIP
!
!  Current Code Owner:  Jonathan Wilkinson
! ----------------------------------------------------------------
           SUBROUTINE LSPCON(CX,CONSTP,x1i,x1ic,x1r,x2r,x4r,m_ci        &
     &,                         l_psd,ai,bi,aic,bic                     &
     &,                         lsp_ei,lsp_fi,lsp_eic,lsp_fic)
             IMPLICIT NONE
!
      Real                                                              &
                  !, Intent(IN)
     &  x1i                                                             &
                  ! Intercept of aggregate size distribution
     &, x1ic                                                            &
                  ! Intercept of crystal size distribution
     &, x1r                                                             &
                  ! Intercept of raindrop size distribution
     &, x2r                                                             &
                  ! Scaling parameter of raindrop size distribution
     &, x4r                                                             &
                  ! Shape parameter of raindrop size distribution
     &, m_ci                                                            &
                  ! Dummy stochastic physics variable    
     &, ai, bi                                                          &
                  ! Aggregate mass-size information m(D) = ai D^bi
     &, aic, bic                                                        &
                  ! Crystal mass-size information m(D) = aic D^bic
     &, lsp_ei, lsp_fi                                                  &
                  ! Aggregate Best-Reynolds Number relation 
                  ! Re(D) =LSP_EI Be^LSP_FI
     &, lsp_eic, lsp_fic  
                  ! Crystal Best-Reynolds Number relation 
                  ! Re(D) =LSP_EIC Be^LSP_FIC

      Logical                                                           &
                  !, Intent(IN)
     &  l_psd
                  ! Use generic particle size distribution

! LOCAL VARIABLES
             INTEGER I
! Counter to print out the contents of CX and CONSTP
             REAL TEMP                                                  &
! Forms input to the GAMMAF routine which calculates gamma functions.
     &,           GBD1,GB1,GD3,GD52,GDR3,GDR4,GDR52,G1,G2,G3
! Represents the gamma function of BI+DI+1 etc.
!
! PROCEDURE CALL
             EXTERNAL GAMMAF
!
#include "c_lspdrp.h"
! Obtain the size of CONSTP and CX
#include "c_lspsiz.h"
! OUTPUT VARIABLES
!
             CX(1)=NINT(1000.*(2.+X4I-X2I)/(BI+1-X2I+X4I))/1000.
! Used in deposition, evaporation of snow and melting calculations.
!
             CX(2)=NINT(1000.*(5.+DI-2.*X2I+2.*X4I)/                    &
     &                  (2.*(BI+1.-X2I+X4I)))/1000.
! Used in deposition, evaporation of snow and melting calculations.
!
             CX(3)=NINT(1000.*DI/(BI+1-X2I+X4I))/1000.
! Used in fall speed and capture calculations.
!
             CX(4)=NINT(1000.*(3.+DI-X2I+X4I)/(BI+1.-X2I+X4I))/1000.
! Used in riming calculations.
!
             CX(5)=NINT(1000.*DR/(4.+DR-X2R))/1000.
! Used in capture calculations.
!
             CX(6)=NINT(1000./(X2I-X4I-1.-BI))/1000.
! Used in capture calculations.
!
             CX(7)=NINT(1000.*(3.0+DR-X2R)/(4.0+DR-X2R))/1000.
! Used in accretion calculations.
!
             CX(8)=NINT(1000.*X2I)/1000.
! Used in capture calculations.
!
             CX(9)=NINT(1000.*X2R)/1000.
! Used in capture calculations.
!
             CX(10)=NINT(1000./(4.0+DR-X2R))/1000.
! Used in capture and evaporation of rain calculations.
!
             CX(11)=NINT(1000.*((DR+5.0)/2.0-X2R))/1000.
! Used in evaporation of rain calculations.
!
             CX(12)=NINT(1000.*(2.0-X2R))/1000.
! Used in evaporation of rain calculations.
!
             CX(13)=X3I
! Used to define temperature dependence of ice particle distribution.
!
             CX(14)=3.+X4I
! Used in capture calculations.
!
             CX(15)=2.+X4I
! Used in capture calculations.
!
             CX(16)=1.+X4I
! Used in capture calculations.
!
! Define values of GBD1 etc.
             TEMP=BI+DI+1.+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GBD1)
             TEMP=BI+1.+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GB1)
             TEMP=3.+DI+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GD3)
             TEMP=2.5+DI/2.+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GD52)
             TEMP=4.+DR
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDR4)
             TEMP=3.+DR
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDR3)
             TEMP=2.5+DR/2.
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,GDR52)
             TEMP=1.+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,G1)
             TEMP=2.+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,G2)
             TEMP=3.+X4I
! DEPENDS ON: gammaf
             CALL GAMMAF(TEMP,G3)
!
! Define values of CONSTP
!
             CONSTP(1)=1.0/(AI*X1I*GB1)
! Used in fallspeed, deposition, riming, capture, evap of snow and
! melting of snow calculations.
!
             CONSTP(2)=6.2832*X1R
!              6.2832 = 2 * pi
! Used in evaporation of rain calculations.
!
             CONSTP(3)=CI*GBD1/GB1
! Used in fallspeed and capture calculations.
!
             CONSTP(4)=0.7854*X1I*CI*GD3
!              0.7854 = pi / 4
! Used in riming calculations.
!
!
             CONSTP(5)=1.0*6.2832*X1I
! 6.2832 = 2 pi   1.0 represents the relative capacitance of spheres.
! Used in deposition and evaporation of snow calculations.
!
             CONSTP(6)=89.48*SQRT(CI)*GD52
!              89.48 = 0.44Sc**0.333 /dyn viscosity**0.5
! Used in deposition, evaporation of snow and melting calculations.
!
             CONSTP(7)=4.57E-7*X1I
!              4.57E-7 = 2 pi ka / Lf
! Used in melting of snow calculations.
!
             CONSTP(8)=523.6*X1R*GDR4*CR
!              523.6 = pi rho(water)/ 6
! Used in capture, evap of rain and melting of snow calculations.
!
             CONSTP(9)=9869.6*X1I*X1R
!              9869.6 = pi**2 rho(water)
! Used in capture calculations.
!
             CONSTP(10)=0.7854*X1R*CR*GDR3
!              0.7854 = pi / 4
! Used in accretion calculations.
!
             CONSTP(11)=CR*GDR4
! Used in capture calculations.
!
             CONSTP(12)=63.100*GDR52*SQRT(CR)
!              63.1 = 0.31Sc**0.333/dyn viscosity**0.5
! Used in evaporation of rain calculations.
!
             CONSTP(13)=G2
! Used in deposition, evap of snow and melting of snow calculations.
!
             CONSTP(14)=G3*0.25
! Used in capture calulations.
!
             CONSTP(15)=G2*2.
! Used in capture calulations.
!
             CONSTP(16)=G1*5.
! Used in capture calulations.
!
! End the subroutine
             RETURN
           END SUBROUTINE LSPCON
!
!  SUBROUTINE GAMMAF-------------------------------------------
!   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
!   A POLYNOMIAL APPROXIMATION
! ----------------------------------------------------------------
!
#endif
