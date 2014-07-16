#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE VERTANL ------------------------------------------------
!LL
!LL  Purpose : High level vertical analysis routine.
!LL            Calls lower level VAN--- routines which are specific
!LL            to particular AC Observation types.
!LL
!LL   For use on Cray Y-MP
!LL   For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  S.Bell     <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    15/7/93 :  Include VANUARS call for 405's; include 306's;
!LL                  :  add MODEL_OBS_TYPE to VANLASS arguments (S Bell)
!LL          15/7/93 :  and switch to VANUTMP for 211's     (R Swinbank)
!LL   3.3    11/11/93:  Pass LAND_MASK and call VANCLD      (N Richards)
!LL   3.3    2/12/93 :Pass 4 RAIN fields and call VANRAIN (B Macpherson)
!LL   3.4    03/08/94 : Pass tracer field and analyse with VANUARS
!LL                                                         (R Swinbank)
!LL   3.4    09/09/94:  Include 1.5m Visibility (P Clark)
!LL   3.4    26/9/94 :  Pass cloud regime boundaries & conv cloud.
!LL                  :  Remove RHCRIT from VANRH arg list
!LL                  :  Add call to VANMOPS for type 406 (B Macpherson)
!    4.2 25/11/96: T3E mods + correct VANUARS call Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!LL  4.5  Mar  98:  Changes for mixed phase cloud microphysics
!LL                                          Bruce Macpherson
!LL  5.2  30/11/00:  amend EXNER dimension,import cloud & pressure
!LL                  remove calls to redundant VAN### routines
!LL                                                       B Macpherson
!    6.0  19/06/03:  Remove non-MPP parts of code. T. White
!    6.2  08/11/05:  Pass through l_eacf logical. Damian Wilson
!LL
!LL   Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL   Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!LL
!LL   Called by: AC2
!*L
      SUBROUTINE VERTANL (PSTGS,                                        &
     &                    KACT,IPASS,LENOBT,NDV,OBS_NO,                 &
     &                    MODEL_OBS_TYPE,OBS_LAT,OBS_LONG,OBS,INOBS,    &
     &                    EXNER,PSTAR,THETA,RH,QCL,QCF,                 &
     &                    CONV_CLD,LAYER_CLOUD,PRESSURE,                &
     &                    LS_RAIN,LS_SNOW,CONV_RAIN,CONV_SNOW,          &
     &                    RHCRIT, L_eacf,                               &
     &                    P_FIELD,WTSFLD,LENWTS,NLEVWT,                 &
     &                    CF1PT,CF2PT,CF3PT,CF4PT,                      &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    OBS_INCR,NORMF,LMISSD,                        &
     &                    P_LEVELS,Q_LEVELS,BL_LEVELS,                  &
     &                    ROW_LENGTH,P_ROWS,                            &
     &                    LENOB,NPTOBT,NO_ANAL_LEVS,NO_ANAL_VAR,        &
     &                    ICODE,CMESSAGE)
!
      IMPLICIT NONE
!
#include "acparm.h"
#include "comacp.h"
#include "comobs.h"
!
      EXTERNAL GETOB3, TIMER, VANRAIN, VANMOPS_MIXED_PHASE
!
      INTEGER P_LEVELS,Q_LEVELS,BL_LEVELS,ROW_LENGTH,P_ROWS
      INTEGER P_FIELD,LENWTS,NLEVWT
      INTEGER LENOBT,NDV,INOBS
      INTEGER LENOB,NO_ANAL_LEVS,NO_ANAL_VAR

      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT)
      INTEGER NP4PT(LENOBT)
      INTEGER OBS_NO(LENOBT)
      INTEGER MODEL_OBS_TYPE(LENOBT)

      REAL                                                              &
     &               EXNER    (P_FIELD,P_LEVELS),                       &
     &               PRESSURE (P_FIELD,P_LEVELS),                       &
     &               LAYER_CLOUD(P_FIELD,Q_LEVELS),                     &
     &               PSTGS    (P_FIELD),                                &
     &               PSTAR    (P_FIELD),                                &
     &               THETA    (P_FIELD,P_LEVELS),                       &
     &               RH       (P_FIELD,Q_LEVELS),                       &
     &               QCL      (P_FIELD,Q_LEVELS),                       &
     &               QCF      (P_FIELD,Q_LEVELS),                       &
     &               CONV_CLD (P_FIELD,Q_LEVELS),                       &
     &               LS_RAIN(P_FIELD),                                  &
     &               LS_SNOW(P_FIELD),                                  &
     &               CONV_RAIN(P_FIELD),                                &
     &               CONV_SNOW(P_FIELD),                                &
     &               WTSFLD   (LENWTS,NLEVWT),                          &
     &               RHCRIT(Q_LEVELS),                                  &
     &               CF1PT(LENOBT),CF2PT(LENOBT),CF3PT(LENOBT),         &
     &               CF4PT(LENOBT),                                     &
     &               OBS_INCR(LENOB+1,NO_ANAL_LEVS,NO_ANAL_VAR),        &
     &               NORMF(LENOB+1,NO_ANAL_LEVS),                       &
     &               OBS_LAT(LENOBT),OBS_LONG(LENOBT),                  &
     &               OBS(INOBS,*)

      LOGICAL LMISSD (LENOB+1,NO_ANAL_LEVS)
      Logical L_eacf      ! Use empirically adjusted cloud fraction
!
      INTEGER KACT, IPASS, KTYPE, NPTOBT, NLEVOB
      INTEGER ICODE
      CHARACTER*256 CMESSAGE
!
!     DYNAMIC ALLOCATION
!
      REAL           OBDATA(LENOBT,NDV)
!*
!-----------------------------------------------------------------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('VERTANL ',3)

!L*** 1       Get obs and error values using GETOB3
      KTYPE  = LACT(KACT)
      NLEVOB = NOBLEV(KACT)

! DEPENDS ON: getob3
      CALL GETOB3 (KACT,OBS,OBDATA,                                     &
#if !defined(GLOBAL)
     &             OBS_LAT,OBS_LONG,                                    &
#endif
     &             LENOBT,INOBS,NDV,OBS_NO,ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999
!
!L*** 2      Call appropriate vertical analysis routine VAN---
!     NB. ICODE is not checked here, because only one VAN routine is
!     called followed by a jump to the end of the routine.
!
      IF (KTYPE == 406) THEN
! for MOPS cloud data

         IF (IPASS == 1) THEN

! DEPENDS ON: vanmops_mixed_phase
             CALL VANMOPS_MIXED_PHASE (                                 &
     &                    KACT,IPASS,THETA,EXNER,CONV_CLD,              &
     &                    LAYER_CLOUD,PRESSURE,                         &
     &                    RH,QCL,QCF,P_FIELD,NO_ANAL_LEVS,RHCRIT,       &
     &                    L_eacf,                                       &
     &                    OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    OBS_INCR,NORMF,OBS_LAT,OBS_LONG,              &
     &                    OBS_NO,LMISSD,NPTOBT,LENOBT,NDV,LENOB,        &
     &                    NO_ANAL_LEVS,NO_ANAL_VAR,                     &
     &                    P_LEVELS,BL_LEVELS,ICODE,CMESSAGE)
         ELSEIF(IPASS == 2) THEN

! DEPENDS ON: vanmops_mixed_phase
             CALL VANMOPS_MIXED_PHASE (                                 &
     &                    KACT,IPASS,THETA,EXNER,CONV_CLD,              &
     &                    LAYER_CLOUD,PRESSURE,                         &
     &                    WTSFLD,QCL,QCF,LENWTS,NLEVWT,RHCRIT,          &
     &                    L_eacf,                                       &
     &                    OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    OBS_INCR,NORMF,OBS_LAT,OBS_LONG,              &
     &                    OBS_NO,LMISSD,NPTOBT,LENOBT,NDV,LENOB,        &
     &                    NO_ANAL_LEVS,NO_ANAL_VAR,                     &
     &                    P_LEVELS,BL_LEVELS,ICODE,CMESSAGE)
         ENDIF


      ELSEIF (KTYPE == 506) THEN
!
        IF (IPASS == 1) THEN

! DEPENDS ON: vanrain
          CALL VANRAIN (KACT,IPASS,LS_RAIN,                             &
     &      LS_SNOW,CONV_RAIN,CONV_SNOW,P_FIELD,                        &
     &      OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,                             &
     &      NP1PT,NP2PT,NP3PT,NP4PT,                                    &
     &      OBS_INCR,NORMF,OBS_LAT,OBS_LONG,                            &
     &      OBS_NO,LMISSD,NPTOBT,                                       &
     &      P_LEVELS,LENOBT,NDV,LENOB,                                  &
     &      NO_ANAL_LEVS,NO_ANAL_VAR,                                   &
     &      ICODE,CMESSAGE)

        ELSEIF(IPASS == 2) THEN

! DEPENDS ON: vanrain
          CALL VANRAIN (KACT,IPASS,WTSFLD,                              &
     &      WTSFLD,WTSFLD,WTSFLD,LENWTS,                                &
     &      OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,                             &
     &      NP1PT,NP2PT,NP3PT,NP4PT,                                    &
     &      OBS_INCR,NORMF,OBS_LAT,OBS_LONG,                            &
     &      OBS_NO,LMISSD,NPTOBT,                                       &
     &      P_LEVELS,LENOBT,NDV,LENOB,                                  &
     &      NO_ANAL_LEVS,NO_ANAL_VAR,                                   &
     &      ICODE,CMESSAGE)

        ENDIF


      ELSE
        ICODE=1
        CMESSAGE = 'VERTANL : Obs Type not known'
        WRITE(6,*) 'VERTANL KTYPE not processed ',KTYPE
      ENDIF
!
 999  CONTINUE

      IF(ICODE >  0) THEN
        WRITE(6,*) 'VERTANL Error code',ICODE
        WRITE(6,*) CMESSAGE
      ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('VERTANL ',4)
      RETURN
      END SUBROUTINE VERTANL
#endif
