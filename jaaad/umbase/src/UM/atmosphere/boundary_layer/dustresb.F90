#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine DUSTRESB
!
!
! Purpose:
!   To calculate the surface layer resistance for mineral dust
!
! Called by sfexch
!
! Current owners of code: S.Woodward
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   5.5      12/02/03  Original code   S Woodward
!   6.2      11/05/06  Fix for bit-non-reproducibility   A. Malcolm
!   6.4      15/01/07  Malcolm McVean's portability fix   S.Woodward
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
! Documentation: "Modelling the atmospheric lifecycle..."
!                 Woodward, JGR106, D16, pp18155-18166
!---------------------------------------------------------------------
!
       SUBROUTINE DUSTRESB(                                             &
     &  ROW_LENGTH,ROWS,                                                &
     &  PSTAR,TSTAR,RHOSTAR,ARESIST,VSHR,CD_STD_DUST,                   &
     &  R_B_DUST                                                        &
     &  )
      IMPLICIT NONE

#include "c_dust_ndiv.h"

      INTEGER                                                           &
              !IN
     & ROW_LENGTH                                                       &
                  !IN
     &,ROWS       !IN

      REAL                                                              &
           !IN
     & PSTAR(ROW_LENGTH,ROWS)                                           &
                                    !IN surface pressure
     &,TSTAR(ROW_LENGTH,ROWS)                                           &
                                    !IN surface temperature
     &,RHOSTAR(ROW_LENGTH,ROWS)                                         &
                                    !IN surface air density
     &,ARESIST(ROW_LENGTH,ROWS)                                         &
                                    !IN aerodynamic resistance
     &,VSHR(ROW_LENGTH,ROWS)                                            &
                                    !IN surface to lowest lev windspeed
!                                   !   difference
     &,CD_STD_DUST(ROW_LENGTH,ROWS) !IN surface transfer coeffient for
!                                   !   momentum, excluding orographic
!                                   !   form drag

      REAL                                                              &
           !OUT
     & R_B_DUST(ROW_LENGTH,ROWS,NDIV) !OUT surface layer resistance for
!                                     !    mineral dust

!     local variables

      INTEGER                                                           &
     & IDIV                                                             &
            !loop counter, dust divisions
     &,I                                                                &
            !loop counter
     &,J                                                                &
            !loop counter
     &,LEV1 !number of levels for vstokes calculation

      REAL                                                              &
     & NU(ROW_LENGTH,ROWS)                                              &
                                 !kinematic viscosity
     &,ETAA(ROW_LENGTH,ROWS)                                            &
                                 !dynamic viscosity of air
     &,LAMDAA(ROW_LENGTH,ROWS)                                          &
                                 !mean free path of air molecules
     &,VSTOKES1(ROW_LENGTH,ROWS)                                        &
                                 !gravitational settling velocity, lev1
     &,NSTOKES(ROW_LENGTH,ROWS)                                         &
                                 !stokes number = VstokesVshrVshr/nu g
     &,NSCHMIDT(ROW_LENGTH,ROWS)                                        &
                                 !schmidt number = nu/diffusivit
     &,TC(ROW_LENGTH,ROWS)                                              &
                                 !temperature in deg C
     &,ALPHACCF(ROW_LENGTH,ROWS)                                        &
                                 !alpha in cunningham correction factor
     &,CCF(ROW_LENGTH,ROWS)                                             &
                                 !Cunningham correction factor
     &,FVSQ(ROW_LENGTH,ROWS)                                            &
                                 !friction velocity squared
     &,WORK(ROW_LENGTH,ROWS)                                            &
                                 !workspace
     &,STOKES_EXP                                                       &
                                 !stokes term in R_B_DUST equation
     &,SMALLP                    !small +ve number, negligible compared to 1

#include "c_g.h"
#include "c_pi.h"
#include "c_0_dg_c.h"
#include "c_r_cp.h"
#include "c_dustgen.h"
#include "c_dustgrav.h"
#include "c_sulchm.h"
!
       EXTERNAL VGRAV
!
!... epsilon() is defined as almost negligible, so eps/100 is negligible
!
      SMALLP = epsilon(1.0) / 100.0
!
!...calc stokes number, schmidt number and finally resistance
!
      LEV1=1

      DO IDIV=1,NDIV
! DEPENDS ON: vgrav
        CALL VGRAV(                                                     &
     &  ROW_LENGTH,ROWS,LEV1,DREP(IDIV),RHOP,PSTAR,TSTAR,               &
     &  VSTOKES1,CCF,ETAA                                               &
     &  )

!CDIR NOVECTOR
        DO J = 1,ROWS
          DO I= 1,ROW_LENGTH
            NSCHMIDT(I,J)=3.*PI*ETAA(I,J)*ETAA(I,J)*DREP(IDIV)/         &
     &       (RHOSTAR(I,J)*BOLTZMANN*TSTAR(I,J)*CCF(I,J))
            NSTOKES(I,J)=VSTOKES1(I,J)*CD_STD_DUST(I,J)*RHOSTAR(I,J)*   &
     &       VSHR(I,J)*VSHR(I,J)/(ETAA(I,J)*G)
            ! Avoid underflow in Stokes term by setting to zero if 
            ! negligible compared to Schmidt term, i.e., if NSTOKES
            ! is too small.
            IF ( 3.0 / NSTOKES(I,J) <                                   &
                 - LOG10( SMALLP *NSCHMIDT(I,J)**(-2./3.) ) ) THEN
               STOKES_EXP = 10.**(-3./NSTOKES(I,J))
            ELSE
               STOKES_EXP = 0.0
            ENDIF
            R_B_DUST(I,J,IDIV)=1./( SQRT(CD_STD_DUST(I,J)) *            &
     &       (NSCHMIDT(I,J)**(-2./3.)+STOKES_EXP) )
          ENDDO !ROW_LENGTH
        ENDDO !ROWS
      ENDDO !NDIV
      
      RETURN
      END SUBROUTINE DUSTRESB
#endif
