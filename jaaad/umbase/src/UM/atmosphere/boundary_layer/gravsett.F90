#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine GRAVSETT ----------------------------------------------
!
! Purpose: To perform gravitational settlement of tracer particles
!          down to the lowest layer of the model.
!          This version allows tracers to fall through 1 or 2 layers.
!
! Current owners of code:                 S Woodward, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   4.4    03/10/97   Original code        S Woodward, M Woodage
!
!   5.5    03/01/03   Modified for New Dynamics
!                     Removed if def 17-A, as required whenever dust
!                     code used             S Woodward
! 6.2      15/12/05  Correct dust diagnostics when substepping phys2.
!                                                        M. Diamantakis
!
! Code description:
!  Language: FORTRAN77 + extensions
!  Programming standard: UMDP 3 Vn 6
!
! System components covered:
!
! System task:
!
!Documentation: Not yet available
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GRAVSETT(                                              &
     & ROW_LENGTH,ROWS,NLEVS,TRACFLD,DIAM,RHOP,P_LAYER_CENTRES,         &
     & P_LAYER_BOUNDARIES,T,TIMESTEP,NUM_SUBSTEPS,SUBSTEP_NUMBER,DRYDEP,&
     & L_CAM_DUST)
!
      IMPLICIT NONE
!
      INTEGER ROW_LENGTH           !IN row length
      INTEGER ROWS                 !IN number of rows
      INTEGER NLEVS                !IN number of model levels
      INTEGER NUM_SUBSTEPS         !IN number of phys2 substeps
      INTEGER SUBSTEP_NUMBER       !IN phys2 substep number
!
      REAL DIAM                    !IN tracer particle diameter
      REAL RHOP                    !IN tracer particle density
      REAL P_LAYER_CENTRES(ROW_LENGTH,ROWS,0:NLEVS)    !IN
      REAL P_LAYER_BOUNDARIES(ROW_LENGTH,ROWS,0:NLEVS) !IN
      REAL T(ROW_LENGTH,ROWS,NLEVS)!IN temperature
      REAL TIMESTEP                !IN timestep s
!
      REAL TRACFLD(ROW_LENGTH,ROWS,NLEVS) !IN/OUT tracer field
!
      REAL DRYDEP(ROW_LENGTH,ROWS) !IN/OUT dep flux from
                                   !       layer2(kg m-2 s-1)
      LOGICAL L_CAM_DUST           !IN Use old version of dust_uplift 
                                   !   scheme for use in CAM NWP models
                                   !   If this is the case, we allow
                                   !   tracers to settle directly from level 1
!

! Include COMDECKS
!
#if defined(MPP)
! Parameters and Common blocks
#include "parvars.h"
#endif
!
#include "c_r_cp.h"
#include "c_g.h"
#include "c_sulchm.h"
!
! External subroutines called
      EXTERNAL VGRAV
!
! Local variables
!
      INTEGER K                  !LOC loop counter for levels
      INTEGER J                  !LOC loop counter for points
      INTEGER I                  !LOC loop counter for points
!
      REAL VRHOCTIMESTEP(ROW_LENGTH,ROWS)  !  v*rho*tracer*deltat @lev
      REAL RHOK2(ROW_LENGTH,ROWS)  !  rho(lev+2)
      REAL RHOK1(ROW_LENGTH,ROWS)  !  rho(lev+1)
      REAL RHOK(ROW_LENGTH,ROWS)   ! rho(lev)
      REAL DZK(ROW_LENGTH,ROWS)    ! thickness of layer lev
      REAL DZK1(ROW_LENGTH,ROWS)   ! thickness of layer lev+1
      REAL DZK2(ROW_LENGTH,ROWS)   ! thickness of layer lev+2
      REAL VSTOKES(ROW_LENGTH,ROWS,NLEVS)!deposition velocity
!                                         (vstokes corrected)
      REAL MASSOUT2K2(ROW_LENGTH,ROWS) !flux falling 2 levs from lev k+2
      REAL MASSOUT1K2(ROW_LENGTH,ROWS) !flux falling 1 levs from lev k+2
      REAL MASSOUT2K1(ROW_LENGTH,ROWS) !flux falling 2 levs from lev k+1
      REAL MASSOUT1K1(ROW_LENGTH,ROWS) !flux falling 1 levs from lev k+1
      REAL MASSOUT2K(ROW_LENGTH,ROWS)  !flux falling 2 levs from lev k
      REAL MASSOUT1K(ROW_LENGTH,ROWS)  !flux falling 1 levs from lev k
      REAL DUMMY1(ROW_LENGTH,ROWS,NLEVS) !
      REAL DUMMY2(ROW_LENGTH,ROWS,NLEVS) !
!
!
! Calculate settlement velocity
!
!      CALL VGRAV(PFIELD,NLEVS,DIAM,RHOP,PSTAR,AK,BK,T,V,DUMMY1,DUMMY2,
!     &           FIRST_POINT,LAST_POINT)
! DEPENDS ON: vgrav
       CALL VGRAV(                                                      &
     &  ROW_LENGTH,ROWS,NLEVS,DIAM,RHOP,                                &
     &  P_LAYER_CENTRES(1:ROW_LENGTH,1:ROWS,1:NLEVS),T,                 &
     &  VSTOKES,DUMMY1,DUMMY2)

!
! Calculate new tracer mixing ratios
!
! Initialise deposition flux to zero
      IF ( SUBSTEP_NUMBER == 1 ) THEN
        DO J=1,ROWS
          DO I=1,ROW_LENGTH
            DRYDEP(I,J)=0.
          ENDDO
        ENDDO
      END IF
!
! Level 1 (K at start of loop)
!
      IF (L_CAM_DUST) THEN
!       If we are running with the CAM dust scheme, we allow
!       tracer to settle directly from level 1
        DO J = 1,ROWS
          DO I=1,ROW_LENGTH
            RHOK(I,J)=P_LAYER_CENTRES(I,J,1)/(R*T(I,J,1))
            DZK(I,J)=(P_LAYER_BOUNDARIES(I,J,0)-P_LAYER_BOUNDARIES(I,J,1))&
     &              /(RHOK(I,J)*G)
            MASSOUT2K(I,J)=0.
            MASSOUT1K(I,J)=RHOK(I,J)*TRACFLD(I,J,1)*                   &
     &         VSTOKES(I,J,1)*TIMESTEP 
          ENDDO               !I
        ENDDO                !J
      ELSE
!     If we are not running the CAM dust scheme, we do not allow
!     tracer to settle directly from level 1
        DO J = 1,ROWS
          DO I=1,ROW_LENGTH
            RHOK(I,J)=P_LAYER_CENTRES(I,J,1)/(R*T(I,J,1))
            DZK(I,J)=(P_LAYER_BOUNDARIES(I,J,0)-P_LAYER_BOUNDARIES(I,J,1))&
     &              /(RHOK(I,J)*G)
            MASSOUT2K(I,J)=0.
            MASSOUT1K(I,J)=0.
          ENDDO               !I
        ENDDO                !J
      ENDIF !L_CAM_DUST
!
! Level 2 (K+1 at start of loop)
!   NB  deposit tracer direct to ground from lev 2 if V high enough
!
      DO J = 1,ROWS
        DO I = 1,ROW_LENGTH
!
          RHOK1(I,J)=P_LAYER_CENTRES(I,J,2)/(R*T(I,J,2))
          DZK1(I,J)=(P_LAYER_BOUNDARIES(I,J,1)-                         &
     &     P_LAYER_BOUNDARIES(I,J,2))/(RHOK1(I,J)*G)
!
!   check for deposition :
          IF (VSTOKES(I,J,2)*TIMESTEP  >   DZK(I,J)) THEN
!       some tracer deposited onto ground
!
            IF (VSTOKES(I,J,2)*TIMESTEP  >   DZK1(I,J)+DZK(I,J)) THEN
!         all deposited to ground
              MASSOUT2K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*DZK1(I,J)
              MASSOUT1K1(I,J)=0.
            ELSE IF ( VSTOKES(I,J,2)*TIMESTEP  >   DZK1(I,J) ) THEN
!         some deposited to ground, some to layer 1

              MASSOUT2K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*                &
     &         (VSTOKES(I,J,2)*TIMESTEP-DZK(I,J))
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*                &
     &         DZK1(I,J)-MASSOUT2K1(I,J)
            ELSE
!         some deposited to ground, some to layer1, some left in layer2
              MASSOUT2K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*                &
     &         (VSTOKES(I,J,2)*TIMESTEP-DZK(I,J))
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*DZK(I,J)
            ENDIF
!
            DRYDEP(I,J)=DRYDEP(I,J)                                     &
     &                 +MASSOUT2K1(I,J)/(NUM_SUBSTEPS*TIMESTEP)
!
          ELSE
!         only falls into layer 1
            MASSOUT2K1(I,J)=0.
            IF ( VSTOKES(I,J,2)*TIMESTEP  >   DZK1(I,J)) THEN
!           all falls into layer 1
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*DZK1(I,J)
            ELSE
!           some to layer 1 , some left in layer2
              MASSOUT1K1(I,J)=RHOK1(I,J)*TRACFLD(I,J,2)*VSTOKES(I,J,2)  &
     &         *TIMESTEP
            ENDIF
!
          ENDIF
!
        ENDDO                !END I LOOP
      ENDDO                 !END J LOOP
!
! Main loop through levels, from bottom up
!
      DO K = 1,NLEVS-2
!
        DO J = 1,ROWS
          DO I = 1,ROW_LENGTH
!
            RHOK2(I,J)=P_LAYER_CENTRES(I,J,K+2)/(R*T(I,J,K+2))
            DZK2(I,J)=(P_LAYER_BOUNDARIES(I,J,K+1)-                     &
     &       P_LAYER_BOUNDARIES(I,J,K+2))/(RHOK2(I,J)*G)
!
!       Calculate mass of tracer falling between levels
!
!        limit fall to 2 levs
           IF (VSTOKES(I,J,K+2)*TIMESTEP >  (DZK1(I,J)+DZK(I,J)))       &
     &      VSTOKES(I,J,K+2)=(DZK1(I,J)+DZK(I,J))/TIMESTEP
!
!          check how far tracer falls:
           IF ( VSTOKES(I,J,K+2)*TIMESTEP  >   DZK1(I,J) ) THEN
!          it falls through more than 1 layer
             IF ( VSTOKES(I,J,K+2)*TIMESTEP  >                          &
     &        (DZK2(I,J)+DZK1(I,J)) ) THEN
!             all into layer k
               MASSOUT2K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK2(I,J)
               MASSOUT1K2(I,J)=0.
             ELSE IF ( VSTOKES(I,J,K+2)*TIMESTEP  >   DZK2(I,J) ) THEN
!            some into k+1, some into k
               MASSOUT2K2(I,J)=                                         &
     &          RHOK2(I,J)*TRACFLD(I,J,K+2)*(VSTOKES(I,J,K+2)*          &
     &          TIMESTEP-DZK1(I,J))
               MASSOUT1K2(I,J)=                                         &
     &          RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK2(I,J)-MASSOUT2K2(I,J)
             ELSE
!            some left in k+2, some into k+1, some into k
               MASSOUT2K2(I,J)=                                         &
     &          RHOK2(I,J)*TRACFLD(I,J,K+2)*(VSTOKES(I,J,K+2)*          &
     &          TIMESTEP-DZK1(I,J))
               MASSOUT1K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK1(I,J)
             ENDIF
!
           ELSE
!          falls no more than 1 layer
             MASSOUT2K2(I,J)=0.
             IF (VSTOKES(I,J,K+2)*TIMESTEP  >   DZK2(I,J)) THEN
!            all falls into layer k+1
               MASSOUT1K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*DZK2(I,J)
             ELSE
!            some falls into k+1, some left in k+2
               MASSOUT1K2(I,J)=RHOK2(I,J)*TRACFLD(I,J,K+2)*             &
     &          VSTOKES(I,J,K+2)*TIMESTEP
             ENDIF
!
          ENDIF
!
! Update tracer field
!
            TRACFLD(I,J,K)=TRACFLD(I,J,K)+(MASSOUT2K2(I,J)+             &
     &       MASSOUT1K1(I,J)-MASSOUT2K(I,J)-MASSOUT1K(I,J))/            &
     &       (RHOK(I,J)*DZK(I,J))
!
! Put k+2 vals in k+1's & k+1's in k's
            MASSOUT1K(I,J)=MASSOUT1K1(I,J)
            MASSOUT1K1(I,J)=MASSOUT1K2(I,J)
            MASSOUT2K(I,J)=MASSOUT2K1(I,J)
            MASSOUT2K1(I,J)=MASSOUT2K2(I,J)
            DZK(I,J)=DZK1(I,J)
            DZK1(I,J)=DZK2(I,J)
            RHOK(I,J)=RHOK1(I,J)
            RHOK1(I,J)=RHOK2(I,J)
!
          ENDDO           !END I LOOP
        ENDDO            !END J LOOP
!
      ENDDO              !END K LOOP
!
! Top 2 levels
!
      DO J=1,ROWS
        DO I=1,ROW_LENGTH
!
          TRACFLD(I,J,NLEVS-1)=TRACFLD(I,J,NLEVS-1)+                    &
     &     (MASSOUT1K1(I,J)-MASSOUT2K(I,J)-MASSOUT1K(I,J))/             &
     &     (RHOK(I,J)*DZK(I,J))
          TRACFLD(I,J,NLEVS)=TRACFLD(I,J,NLEVS)-                        &
     &     (MASSOUT2K1(I,J)+MASSOUT1K1(I,J))/                           &
     &                (RHOK1(I,J)*DZK1(I,J))
!
        ENDDO       !I
      ENDDO        !J
!
      RETURN
      END SUBROUTINE GRAVSETT
#endif
