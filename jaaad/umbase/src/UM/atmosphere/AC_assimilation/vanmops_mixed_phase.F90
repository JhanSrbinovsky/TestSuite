#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINES VANMOPS_MIXED_PHASE and CC_TO_RHTOT-------------------
!LL
!LL  Purpose : Performs vertical analysis of MOPS cloud data
!LL            for (2A) mixed phase cloud microphysics scheme
!LL     When IPASS=1,
!LL     model field is interpolated to ob locations and increments
!LL     calculated. Preliminary weight normalisation factors are also
!LL     calculated.
!LL
!LL     When IPASS=2,
!LL     data density is interpolated to ob locations and
!LL     final weight normalisation factors are calculated.
!LL
!LL     Documentation in a Working Paper by B Macpherson & D Wilson is
!LL     on Metweb at URL
!LL     http://fr2010/~frff/papers/MOPS_for_new_micro.ps.gz
!LL
!LL  For use on Cray
!LL
!LL
!LL B Macpherson<- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history :
!LL version  Date
!LL   4.5   18/11/97  New routine developed from VANMOPS B Macpherson
!LL   5.2   30/11/00  amend for new dynamics      B Macpherson
!LL                   NB oper 4.5 mods not yet included
!LL   5.3   08/06/01  Remove duplicate declarations.  A van der Wal
!     5.3   17/09/01  fix incorrect changes to
!                     theta and q.  B Macpherson
!     6.2   07/11/05  Include empirically adj. cloud frac. D. Wilson
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND-------------------------------------------------------------
!
!*L  Arguments:---------------------------------------------------
       SUBROUTINE VANMOPS_MIXED_PHASE(                                  &
     &                    KACT,IPASS,THETA,EXNER,CONV_CLD,              &
     &                    CF,PRESSURE,                                  &
     &                    Q,QCL,QCF,LENFLD,NLVFLD,RHCRIT,L_eacf,        &
     &                    OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    QINC,NORMF,OBS_LAT,OBS_LONG,                  &
     &                    OBS_NO,LMISSD,NPTOBT,LENOBT,NDV,LENOB,        &
     &                    NO_ANAL_LEVS,NO_ANAL_VAR,                     &
     &                    P_LEVELS,BL_LEVELS,ICODE,CMESSAGE)

      IMPLICIT NONE
!-----------------------------------------------------------------
!  UM comdecks and functions
!-----------------------------------------------------------------
#include "c_r_cp.h"
#include "c_lheat.h"
#include "c_0_dg_c.h"
!-----------------------------------------------------------------
!  Analysis Correction comdecks
!-----------------------------------------------------------------
#include "acparm.h"
#if defined(MPP)
#include "parvars.h"
#else
      INTEGER mype
      PARAMETER (mype = 0 ) ! always zero in non-MPP code
#endif
#include "comacp.h"
#include "comobs.h"
#include "comacdg.h"
#include "comag.h"
!-----------------------------------------------------------------
      INTEGER LENFLD,NLVFLD,LENOBT,LENOB,NDV,NO_ANAL_LEVS,              &
     &        P_LEVELS,BL_LEVELS,NO_ANAL_VAR
      REAL                                                              &
     &               THETA(LENFLD,P_LEVELS),                            &
     &               EXNER(LENFLD,P_LEVELS),                            &
     &               PRESSURE (LENFLD,P_LEVELS),                        &
     &               CONV_CLD(LENFLD,NLVFLD),                           &
     &               CF(LENFLD,NLVFLD),                                 &
     &               Q   (LENFLD,NLVFLD),                               &
     &               QCL (LENFLD,NLVFLD),                               &
     &               QCF (LENFLD,NLVFLD),                               &
     &               RHCRIT  (NLVFLD),                                  &
     &               OBDATA  (LENOBT,NDV),                              &
     &               QINC(LENOB+1,NO_ANAL_LEVS),                        &
     &               NORMF(LENOB+1,NO_ANAL_LEVS),                       &
     &               OBS_LAT(LENOBT),OBS_LONG(LENOBT),                  &
     &               CF1PT(LENOBT), CF2PT(LENOBT),                      &
     &               CF3PT(LENOBT), CF4PT(LENOBT)
      LOGICAL LMISSD(LENOB+1,NO_ANAL_LEVS)
      LOGICAL L_eacf      ! Use empirically adjusted cloud fraction
      INTEGER OBS_NO(LENOBT)
      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT)
      INTEGER NP4PT(LENOBT)

      INTEGER KACT,IPASS,NPTOBT
      INTEGER ICODE
      CHARACTER*256 CMESSAGE

!-INTENT=IN--------------------------------------------------------
!     KACT      - Observation type
!     IPASS     - Calculate incrs and preliminary norm factors for
!                 IPASS=1 and final norm factors for IPASS=2.
!     LENOBT    - no of obs in type
!     LENOB     - no of obs in group
!     NPTOBT    - pointer to present ob type within group
!     NDV       - no of data values in observation
!     NO_ANAL_LEVS - no of analysis levels
!     NO_ANAL_VAR  - no of analysis variables
!     P_LEVELS     - no of model levels
!     EXNER     - diagnostic variable exner pressure
!     CONV_CLD  - convective cloud amount on model levels
!     CF         - layer cloud amount (liq+ice) on model levels
!     PRESSURE   - p on theta levels
!     Q         - model specific humidity      (IPASS=1)
!               - data density on model grid   (IPASS=2)
!     QCF       - model cloud ice
!     QCL       - model cld liq water
!     THETA     - model potential temperature
!     LENFLD    - length of model field
!     NLVFLD    - no of levels in Q field
!     RHCRIT    - critical rh values for cloud formation (fraction)
!     OBDATA    - observed values
!     OBS_LAT   - ob co-lats
!     OBS_LONG  - ob longitudes
!     OBS_NO    - observation numbers
!     NP1PT,NP2PT,NP3PT,NP4PT,
!               - pointers to nearest model grid points for each ob
!     CF1PT,CF2PT,CF3PT,CF4PT,
!               - weighting of model points in bilinear interpolation
!
!-INTENT=INOUT-----------------------------------------------------
!     QINC        -  ob-model humidity increments
!     NORMF       -  normalisation factors
!-INTENT=OUT-----------------------------------------------------
!     LMISSD          -  logical to indicate missing data
!     ICODE,CMESSAGE  - error code and message
!*
!----------------------------------------------------------------------
!*L   Workspace usage
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION ON CRAY
      REAL           TL(LENFLD,NLVFLD)

      REAL           CLD_INC (LENOB+1,NO_ANAL_LEVS)
      REAL           CLD_NORM(LENOB+1,NO_ANAL_LEVS)

      REAL           CF_LYROB(LENOBT), CF_CONVB(LENOBT)
      REAL           QB(LENOBT), RHTOT_OB(LENOBT)
      REAL           CFB(LENOBT), TLB(LENOBT), QCLB(LENOBT)
      REAL           QCFB(LENOBT), PB(LENOBT)
      REAL           QSAT_ICEB(LENOBT), QSAT_WATB(LENOBT)

!     TL         - model liquid water temperature
!     CLD_INC    - obs - model cloud fraction (layer+convective)
!     CLD_NORM   - normalisation array needed for input to DIAGO call
!     CF_LYROB   - target layer cloud fraction derived from cloud ob
!     CF_CONVB   - model convective cloud fraction interp'd to obs pts
!     QB         - model humidity (or obs density) at obs point
!     RHTOT_OB   - target rhtot derived from CF_LYROB
!     CFB        - CF at obs points
!     TLB        - TL at obs points
!     QCLB,QCFB  - QCL,QCF at obs points
!     PB         - model level pressure at obs pt
!     QSAT_ICEB  - QSAT wrt ice at obs pt
!     QSAT_WATB  - QSAT wrt water at obs pt
!*
!----------------------------------------------------------------------
!*L   External subroutine calls
!-----------------------------------------------------------------------
      EXTERNAL TIMER,HINTMO
      EXTERNAL DIAGO,CC_TO_RHTOT,QSAT,QSAT_WAT
!*
!----------------------------------------------------------------------
!     Define local variables
!----------------------------------------------------------------------
      REAL    TQSATICE,TQSATWAT,DTQSAT
      INTEGER KTYPE,NOBIPT,NOBLPT,NEROPT,JK,JOB,J,ERROR
!     KTYPE    -   Observation type
!     NOBIPT   -   Pointer to data values within OBDATA array
!     NOBLPT   -   Pointer to obs levels within OBDATA array
!     NEROPT   -   Pointer to obs errors within OBDATA array
!     JK       -   Loop counter in loops over levels
!     JOB      -         "      in loops over obs
!     J        -         "      in loops over points
!     ERROR    -   error code
!     TQSATICE -   new cloud likely to be ice below this temp
!                  (see THOMO in C_LSPMIC)
!     TQSATWAT -   new cloud likely to be liquid above this temp
!                  (see TNUC in C_LSPMIC)
!     DTQSAT   -   difference between two previous temperatures
!-----------------------------------------------------------------------

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANMMP  ',3)

!L
!L*** 1. PRELIMINARIES
!L       -------------

      KTYPE = LACT (KACT)

      IF (KTYPE == 406) THEN
        NOBIPT = 1
      ENDIF
      NEROPT = NERLEV1(KACT)-NDVHDR

!     Define temps used in choosing qsat value when creating cloud
!     (These temps are the same as the homogeneous and heterogeneous
!      nucleation thresholds in the mixed phase precip scheme,
!      but need not be, hence separate definition)
      TQSATICE = ZERODEGC - 40.0
      TQSATWAT = ZERODEGC - 10.0
      DTQSAT   = TQSATWAT - TQSATICE

      IF (KTYPE == 406) THEN
!     =================
      IF(IPASS == 1) THEN
!L*** 2.1 Current model cloud fraction is consistent with
!L        latest humidity, since ls_cld called just before AC.
!L        But need to get current TL for use in qsat calcns.


!     2.1.2 calculate TL
         DO JK=1,NLVFLD
           DO J=1,LENFLD
            TL(J,JK) = THETA(J,JK)*EXNER(J,JK)- LC * QCL(J,JK) /CP
           END DO
         END DO

!L*** 2.2  HORIZONTAL INTERP OF FIELDS ON MODEL GRID


        ENDIF ! IPASS == 1

! begin loop over levels
      DO JK=1,NO_ANAL_LEVS

!       PUT Q OR OBS DENSITY INTERPOLATED TO OB POSITIONS IN QB
! DEPENDS ON: hintmo
        CALL HINTMO( Q(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,                   &
     &               NP1PT,NP2PT,NP3PT,NP4PT,                           &
     &               LENFLD,1,LENOBT,QB,ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 9999

!L*** 2.3  CALCULATION OF INCRS AND WEIGHT FACTORS FOR OBS INC

!     2.3.1  FIRST SET UP BIT ARRAY FOR NO OBSERVATIONAL DATA
      DO JOB=1,LENOBT
        LMISSD(NPTOBT+JOB,JK) = OBDATA(JOB,NOBIPT+JK-1) == MISSD .OR.   &
     &                          OBDATA(JOB,NEROPT+JK-1) == MISSD
      ENDDO   !JOB

      IF (IPASS == 1) THEN

!      2.3.2 put conv_cld interp'd to ob locations in CF_CONVB
! DEPENDS ON: hintmo
       CALL HINTMO( CONV_CLD(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,             &
     &              NP1PT,NP2PT,NP3PT,NP4PT,                            &
     &              LENFLD,1,LENOBT,CF_CONVB,ICODE,CMESSAGE)
       IF(ICODE >  0) GO TO 9999

!      2.3.3 calculate target layer cloud cover in CF_LYROB
!            based on total cloud = conv_cld + ls_cld *(1-conv_cld)
!            as assumed in radiation scheme.
        DO JOB=1,LENOBT
          IF(CF_CONVB(JOB) /=  1.) THEN
            CF_LYROB(JOB)=(OBDATA(JOB,NOBIPT+JK-1)-CF_CONVB(JOB))/      &
     &                                           (1.-CF_CONVB(JOB))
          ELSE
!         set target layer cloud of zero where conv cover = 1
            CF_LYROB(JOB)=0.
          ENDIF
!         set target layer cloud of zero where conv cover > MOPS
!         (will apply also where MOPS data is missing data value<0,
!          though obs cloud not used then)
          CF_LYROB(JOB)=MAX(0.,CF_LYROB(JOB) )
        ENDDO

!      2.3.4 convert target layer cloud to rhtot
! DEPENDS ON: cc_to_rhtot
        CALL CC_TO_RHTOT (CF_LYROB,LENOBT,RHCRIT(JK),RHTOT_OB,JK,       &
     &                                             L_eacf,BL_LEVELS)

!      2.3.5 interp model b'ground fields to ob pts for incr calcns
!            calculate derived quantities

!        2.3.5.1 put TL    interp'd to ob pts in TLB
! DEPENDS ON: hintmo
          CALL HINTMO( TL(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,                &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,TLB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.2 put QCL   interp'd to ob pts in QCLB
! DEPENDS ON: hintmo
          CALL HINTMO( QCL(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,QCLB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.3 put QCF   interp'd to ob pts in QCFB
! DEPENDS ON: hintmo
          CALL HINTMO( QCF(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,QCFB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.4 put TOTAL LAYER CLOUD FRACTION interp'd to
!                                                   ob pts in CFB
! DEPENDS ON: hintmo
          CALL HINTMO( CF(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,                &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,CFB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.5 put MODEL LEVEL PRESSURE interp'd to ob pts in PB
! DEPENDS ON: hintmo
      CALL HINTMO( PRESSURE(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,              &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,PB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.6 get QSAT(TL,P) wrt ice in QSAT_ICEB
! DEPENDS ON: qsat
           CALL QSAT (QSAT_ICEB, TLB, PB, LENOBT)

!        2.3.5.7 get QSAT(TL,P) wrt water in QSAT_WATB
! DEPENDS ON: qsat_wat
           CALL QSAT_WAT (QSAT_WATB, TLB, PB, LENOBT)

!        2.3.5.8 update RHTOT_OB where ice present
!                covers eqns 12 & 13 in working paper
           DO JOB=1,LENOBT
             IF (QCFB(JOB)  >   0.0) THEN
               RHTOT_OB(JOB) = RHTOT_OB(JOB) *                          &
     &          (QCLB(JOB) + QCFB(JOB)*                                 &
     &                          QSAT_ICEB(JOB)/QSAT_WATB(JOB) )         &
     &         /( QCLB(JOB) + QCFB(JOB) )
             ENDIF
           END DO

!      2.3.6 calculate increments to Q and normalisation factors
          DO JOB=1,LENOBT
            IF(.NOT.LMISSD(NPTOBT+JOB,JK)) THEN
!             Data is valid
!             Calculate incrs to Q
              IF( CFB(JOB) == CF_LYROB(JOB) ) THEN
!               model and ob cloud equal, so no humidity incr
                QINC(NPTOBT+JOB,JK) = 0.0
              ELSEIF(CFB(JOB) >  0.0) THEN
!               eqns 7 & 8 in working paper
                QINC(NPTOBT+JOB,JK) = RHTOT_OB(JOB)*QSAT_WATB(JOB) -    &
     &                                ( QB(JOB)+QCLB(JOB)+QCFB(JOB) )
              ELSEIF(CFB(JOB) == 0.0 .AND. CF_LYROB(JOB) >  0.0) THEN
!               need to create cloud in model,
!               refer to section 2.5 of working paper
                IF    (TLB(JOB) >  TQSATWAT) THEN
                  QINC(NPTOBT+JOB,JK) = RHTOT_OB(JOB)                   &
     &                                 *QSAT_WATB(JOB) - QB(JOB)
                ELSEIF(TLB(JOB) <  TQSATICE) THEN
                  QINC(NPTOBT+JOB,JK) = RHTOT_OB(JOB)                   &
     &                                 *QSAT_ICEB(JOB) - QB(JOB)
                ELSE
                  QINC(NPTOBT+JOB,JK)  = RHTOT_OB(JOB) *                &
     &              (QSAT_ICEB(JOB)*(1-(TLB(JOB)-TQSATICE)/DTQSAT) +    &
     &               QSAT_WATB(JOB)*(TLB(JOB)-TQSATICE)/DTQSAT )        &
     &              - QB(JOB)
                ENDIF
              ENDIF
!             check humidity incr reasonable, or reset to zero
!             (humidity and cloud incrs must have same sign)
              IF( QINC(NPTOBT+JOB,JK)*                                  &
     &                      (CF_LYROB(JOB)-CFB(JOB)) <  0.0 ) THEN
                QINC(NPTOBT+JOB,JK) = 0.0
              ENDIF
!             Calculate normalisation factor
              NORMF(NPTOBT+JOB,JK) = 1.0 /                              &
     &            ( OBDATA(JOB,NEROPT+JK-1)*OBDATA(JOB,NEROPT+JK-1) )
            ELSE
!             Missing data
              QINC(NPTOBT+JOB,JK)  = 0.0
              NORMF(NPTOBT+JOB,JK) = 0.0
            ENDIF ! test on missing data
          ENDDO  !  JOB

!     2.3.7  store cloud increments (oktas) for diagnostics
      IF(LDIAGAC) THEN
        DO JOB=1,LENOBT
          IF(.NOT.LMISSD(NPTOBT+JOB,JK)) THEN
            CLD_INC (JOB,JK) = 8*(OBDATA(JOB,NOBIPT+JK-1)-CF_CONVB(JOB) &
     &                           -CFB(JOB)*(1. - CF_CONVB(JOB) )  )
            CLD_NORM(JOB,JK) = 1.0
          ELSE
            CLD_INC (JOB,JK) = 0.0
            CLD_NORM(JOB,JK) = 0.0
          ENDIF
        END DO

      ENDIF


      ELSEIF (IPASS == 2) THEN
          DO JOB=1,LENOBT
            IF (.NOT.LMISSD(NPTOBT+JOB,JK)) THEN
              NORMF(NPTOBT+JOB,JK) = NORMF(NPTOBT+JOB,JK) /             &
     &                              ( QB(JOB) + 1.0 )
            ELSE
              NORMF(NPTOBT+JOB,JK) = 0.0
            ENDIF
          ENDDO !  JOB
      ENDIF ! IPASS


      ENDDO   ! JK

!L*** 3.  Diagnostics
!         ===========
!     measure fit of model to MOPS cloud cover
      IF(LDIAGAC .AND. IPASS == 1 ) THEN
          if(mype == 0)                                                 &
     &    PRINT '(/,'' DIAGO called from VANMOPS_MIXED_PHASE            &
     &                 - Obs Type'',I5,                                 &
     &                T50,''No of obs '',I6)', KTYPE,LENOBT

! DEPENDS ON: diago
          CALL DIAGO ('MULTI-LEVEL',KTYPE,6,                            &
     &                CLD_INC,CLD_NORM,OBS_LAT,OBS_LONG,LMISSD,         &
     &                LENOBT,LENOBT,0,NO_ANAL_LEVS,NO_ANAL_VAR)


      ENDIF

      ENDIF  ! KTYPE
!     =====

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANMMP ',4)

 9999 CONTINUE
      RETURN
      END SUBROUTINE VANMOPS_MIXED_PHASE



#endif
