#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine VGRAV -------------------------------------------------
!
! Purpose: To calculate the gravitational sedimentation velocity of
!          tracer particles according to Stoke's law, including the
!          Cunningham correction factor.
!
! Current owners of code:                 S Woodward, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   4.4    03/10/97   Original code        S Woodward, M Woodage
!   5.5    12/02/03   Updated for vn 5.5   S Woodward
!
!
! Code description:
!  Language: FORTRAN77 + extensions
!  Programming standard: UMDP 3 Vn 6
!
! System components covered:
!
! System task:
!
!Documentation: Ref. Pruppacher & Klett
!                    Microphysics of clouds & ppn    1978,1980 edns.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE VGRAV(                                                 &
     & ROW_LENGTH,ROWS,NLEVS,DIAM,RHOP,P,T,                             &
     & VSTOKES,CCF,ETAA)




!
      implicit none
!
!
      INTEGER ROW_LENGTH         !IN row length
      INTEGER ROWS               !IN number of rows
      INTEGER NLEVS              !IN number of levels
!
      REAL DIAM                  !IN particle diameter
      REAL RHOP                  !IN particles density
      REAL P(ROW_LENGTH,ROWS,NLEVS)!IN pressure
      REAL T(ROW_LENGTH,ROWS,NLEVS)!IN temperature
!
      REAL VSTOKES(ROW_LENGTH,ROWS,NLEVS) !OUT sedimentation velocity
      REAL ETAA(ROW_LENGTH,ROWS,NLEVS)!OUT viscosity of air
      REAL CCF(ROW_LENGTH,ROWS,NLEVS) !OUT cunningham correction factor
!
!
#include "c_g.h"
#include "c_0_dg_c.h"
#include "c_sulchm.h"
#include "c_dustgrav.h"
!
! local variables
!
      INTEGER ILEV               !LOC loop counter for levels
      INTEGER I                  !LOC loop counter
      INTEGER J                  !LOC loop counter
      INTEGER K                  !LOC loop counter
!
      REAL TC(ROW_LENGTH,ROWS)   !LOC temperature in deg C
      REAL LAMDAA(ROW_LENGTH,ROWS)!LOC mean free path of particle
      REAL ALPHACCF(ROW_LENGTH,ROWS)!LOC
!
! Calculate viscosity of air (Pruppacher & Klett p.323)
      DO ILEV=1,NLEVS
        DO J=1,ROWS
          DO I = 1,ROW_LENGTH
           TC(I,J)=T(I,J,ILEV)-ZERODEGC
           IF (TC(I,J)  >=  0.) THEN
            ETAA(I,J,ILEV)=(1.718+0.0049*TC(I,J))*1.E-5
           ELSE
            ETAA(I,J,ILEV)=                                             &
     &            (1.718+0.0049*TC(I,J)-1.2E-5*TC(I,J)*TC(I,J))*1.E-5
           ENDIF
!
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !NLEVS
!
      DO ILEV=1,NLEVS
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
!
! Calculate mean free path of particle (Pruppacher & Klett p.323)
           LAMDAA(I,J)=MFP_REF*PREF_MFP*T(I,J,ILEV)/                    &
     &      (P(I,J,ILEV)*TREF_MFP)
! Calculate Cunningham correction factor(Pruppacher & Klett p.361)
           ALPHACCF(I,J)=ACCF+BCCF*EXP(CCCF*DIAM*.5/LAMDAA(I,J))
           CCF(I,J,ILEV)=(1.+ALPHACCF(I,J)*LAMDAA(I,J)/(.5*DIAM))
! Calculate sedimentation velocity (Pruppacher & Klett p.362)
           VSTOKES(I,J,ILEV)=CCF(I,J,ILEV)*(DIAM*DIAM*G*RHOP)/          &
     &             (18.*ETAA(I,J,ILEV))
!
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !NLEV
!
      RETURN
      END SUBROUTINE VGRAV
#endif
