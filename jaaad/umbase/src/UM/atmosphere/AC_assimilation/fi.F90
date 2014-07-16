#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE FI -------------------------------------------------
!LL
!LL  Purpose :
!LL
!LL     Perform the horizontal spreading of observational increments.
!LL     There are three stages:
!LL     1. calculate the horizontal scale to be used
!LL     2. spread the observation weights to get normalization factors
!LL     3. spread the normalized increments.
!LL     Spreading is performed using the adjoint of the model->ob interp
!LL     followed by a filter on the model grid using RF.
!LL
!LL     This is an alternative to the method using HORINF.
!LL
!LL  For use on Cray Y-MP
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL  For default cray code, not using GLOBAL implies ELF
!LL
!LL  Written 30/4/92  by Andrew Lorenc
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.1   8/1/92 : Replace IANWTS by IPASS   (R S Bell)
!LL      16/2/93 : Add option to use adjoint of horizontal interpolation
!LL              : to spread ob incrs to model grid but without further
!LL              : filtering on model grid - NPASS_RF=0 (for MOPS data).
!LL              : Introduce MRAMPFN option for consistency with HORINF.
!LL              :                               (B Macpherson)
!LL 4.2  29/11/96  : Mods for T3e Stuart Bell
!LL 5.2  12/12/00: change A to Earth_Radius, remove refs to
!LL                   DIVFILT,DACO,MMSPTW     B Macpherson
!LL 6.2  05/01/06: Remove use of uninitialised array CSCFACT_H.
!LL                Camilla Mathison
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  Documentation: UM Doc Paper 30 section 6.4
!LL
!LLEND------------------------------------------------------------------

!*L  Arguments:---------------------------------------------------
      SUBROUTINE FI(IPASS,TIMESTEP_NO,TIMESTEP,NAVT,JGROUP,             &
     &  NO_ANAL_VAR,LENOB,NO_ANAL_LEVS,NO_WT_LEVS,FIRST_LEVEL,          &
     &  LENMG,ROW_LENGTH,P_ROWS,OBS_NO,                                 &
     &  OBS_LAT,OBS_LONG,OBS_TIME,NORMF,OBS_INCR,CFPT,NPPT,             &
     &  IP_SCALE,MODEL_INCR,PSTAR,                                      &
     &  P_FIELD,P_LEVELS,                                               &
     &  ICODE,CMESSAGE)
       IMPLICIT NONE

      INTEGER   IPASS,TIMESTEP_NO
      INTEGER   JGROUP,FIRST_LEVEL
      INTEGER   NAVT,NO_ANAL_VAR,NO_ANAL_LEVS,NO_WT_LEVS,IP_SCALE
      INTEGER   LENOB,LENMG,ROW_LENGTH,P_ROWS
      REAL      TIMESTEP
      INTEGER   ICODE  ! 0 if OK
      CHARACTER*256 CMESSAGE

      INTEGER OBS_NO(LENOB),NPPT(LENOB+1,4)
      REAL                                                              &
     &     OBS_LAT(LENOB),OBS_TIME(LENOB),CFPT(LENOB+1,4),              &
     &     OBS_LONG(LENOB),NORMF(LENOB+1,NO_ANAL_LEVS),                 &
     &     OBS_INCR(LENOB+1,NO_ANAL_LEVS,NO_ANAL_VAR),                  &
     &     MODEL_INCR(LENMG,NO_ANAL_LEVS,IP_SCALE),                     &
     &     PSTAR(ROW_LENGTH*P_ROWS)

      INTEGER                                                           &
     & P_FIELD,P_LEVELS !IN PTS ON U GRID;P GRID ,NO OF LEVELS


!-INTENT=IN--------------------------------------------------------
!     IPASS : mode
!            1- do stages 1. calculate scale
!                       & 2. spread NORMF
!            0- do stage  3. spread INCR
!     TIMESTEP_NO: current model timestep (TIMESTEP_NO=1,.....N)
!     NAVT: observation variable type
!     JGROUP : ob group code
!     NO_ANAL_LEVS: number of levels analysed ( >=  NO_WT_LEVS)
!     FIRST_LEVEL : model level corresponding to analysis level 1
! N.B. current AC code assumes this is 1; this could be relaxed later
!     NO_WT_LEVS: number of levels for which NORMF is needed
!     LENOB: number of observations in group
!     LENMG : NUMBER OF POINTS IN FIELD ON model GRID
!     ROW_LENGTH : length of model grid row
!     P_ROWS : no of row on model p grid (assumed > U_ROWS)
!     TIMESTEP   : is model timestep in secs
!     OBS_NO : list of ob numbers for this group within ACOBS file
!     OBS_LAT  :  co-lat of obs
!     OBS_LONG : longitude of obs
!     OBS_TIME : time of obs rel to assm start
!     NORMF    : normalization factors for obs (from routine VERTANL)
!     OBS_INCR : increment (ob-model) at ob pt
!     CFPT : model->ob interpolation coeffs
!     NPPT : model->ob interpolation pointers
!     PSTAR : surface pressure (used in DIVFILT for mass weighting)
!     IP_SCALE : 3rd index to MODEL_INCR, pointing to space for scale
!     IP_SCALE=2 when doing a scalar, =3 when doing a vector analysis.
!-INTENT=OUT---------------------------------------------------------
!    stage 1:
!     MODEL_INCR(,,IP_SCALE)     : horizontal scale for stages 2 & 3
!     ( MODEL_INCR(,,1:2) )      : is used as work space in stage 1
!    stage 2:
!     MODEL_INCR(,,1)            : data density field D for norm.factor
!                                : (only first NO_WT_LEVS levels done)
!    stage 3: (second call to FI)
!     MODEL_INCR(,,1:NO_ANAL_VAR): increments
!-INTENT=IN   (second call to FI)
!     MODEL_INCR(,,IP_SCALE)     : horizontal scale
!*---------------------------------------------------------------------

!*L  Workspace usage:-------------------------------------------------
      REAL COS_LAT(P_ROWS) ! sin of co-latitude of model grid rows
      REAL S_INCR(LENOB,NO_ANAL_LEVS,2) ! incs at obs to be scattered
      REAL TF(LENOB) ! Time Factor for each ob  as in HORINF
      REAL S(LENOB)  ! Scale for each ob        as in HORINF
      REAL WRK1(LENOB)
      INTEGER IWK1(LENOB)
!
!*---------------------------------------------------------------------

!*L  external subroutine calls:---------------------------------------
#if defined(GLOBAL)
      EXTERNAL RFCSG,RFVSG,RFVVG ! global filter routines
      EXTERNAL MMSPT
#else
      EXTERNAL RFCSL,RFVSL,RFVVL ! ltd.area filter routines
#endif
      EXTERNAL TIMER
!*---------------------------------------------------------------------

!----analysis correction comdecks-------------------------------------
#include "acparm.h"
#include "comacp.h"
#include "comag.h"
#include "comacdg.h"
#include "commg.h"
!*---------------------------------------------------------------------

!-- constants comdecks-----------------------------------------------
#include "c_a.h"
#include "c_pi.h"
! ---------------------------------------------------------------------

!----------------------------------------------------------------------
!    Define local variables
      INTEGER JROW,JOB,JOBT,J,JLEV,JC ! loop indices
      INTEGER I,ITYPE
      INTEGER IROWS ! number of rows on model grid for current variable
      LOGICAL ILUV  ! .true. if winds
      LOGICAL FILT  ! .true. if filtering on model grid is reqd.
      INTEGER M_GRID ! 0=Ltd Area, 1=pts at pole, 2=staggd from pole
      REAL    TIMEB,TIMEA,RADINF
      REAL    CSCALE_START,CSCALE_OBTIME,CSCALE_END
      INTEGER FIRST_OBS,LAST_OBS,FIRST_TYPE,LAST_TYPE
      REAL Z,ZZ,ZINTEGRAL,ZIOA
      REAL FI_VAR_FACTOR  ! to change SCALE to match tuning to HORINF
      CHARACTER*10 LABEL
!----------------------------------------------------------------------

!L***     PRELIMINARIES
!         -------------
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('FI      ',3)
      ICODE=0
      FILT = NPASS_RF >  0
      IF (FILT) THEN
      ZIOA=(25.1327-18.2655/NPASS_RF)/(DLATMG*DLONGMG)
!     ZIOA is a constant in the formula for correlation area/box area
!     ZIOA=v (6.5.15) / (dlat*dlong)
!     For NPASS_RF =2, v=16.; NPASS_RF=infinity, v=8*PI.
!     v for other NPASS_RF is got from an approximate interpolation.
!     the COS_LAT(JROW) term in the grid box area is done in the loops
!
!     HORINF truncated the influence area at RADINF corr scales
!     and uses an EXP(-r/SCALE) correlation function for wind -
!     FI_VAR_FACTOR compensates for this so FI gives similar
!     results to HORINF with the default values for scale etc..
      FI_VAR_FACTOR=DEF_FI_VAR_FACTOR(GROUP_INDEX(JGROUP))
      FI_VAR_FACTOR=SQRT(FI_VAR_FACTOR)
      ELSE
        ZIOA=1.
        FI_VAR_FACTOR=1.
      ENDIF  !  FILT

      ILUV=NAVT == 3
      IF(ILUV)THEN  ! model U grid
        Z=ROW1MGUV
        IROWS=P_ROWS-1
        M_GRID=2
      ELSE          ! model P grid
        Z=ROW1MGTH
        IROWS=P_ROWS
        M_GRID=1
      ENDIF
!L    INITIALIZE COS_LAT
      DO      J=1,IROWS
        COS_LAT(J)=SIN(Z)
        Z=Z+DLATMG
      ENDDO ! J
#if defined(GLOBAL)
      IF(M_GRID == 1)THEN   ! set area of triangular gridbox at pole
        COS_LAT(1)=.125*COS_LAT(2)
        COS_LAT(IROWS)=COS_LAT(1)
      ENDIF
#else
      M_GRID=0              ! no pole in ltd area version
#endif

!**** define time-factors and scales for each observation
!     The following code for section 2.4 was copied from HORINF (2.6)
!     (preliminary version 1/6/92) and updated 16/2/93
!     so that FI behaves as similarly as possible to HORINF for testing
!
!L2.4 TF which is time factor R(dT) in 3.18 and 3.21
!L##  AND S(DT) BEING THE CORRELATION SCALE AS IN EQUATION 3.19

!     for synoptic insertion mode treat all obs times as though at
!     the nearest synoptic hour I.E. T(OB) = T+0
      IF(LSYN)THEN
       DO 2400 JOB=1,LENOB
       OBS_TIME(JOB) = OBTIME_NOM
!      value of OBTIME_NOM is relative to start of assimilation
2400   ENDDO ! JOB
      END IF

      FIRST_OBS = 1

!     Get first and last types for group being processed.
      FIRST_TYPE = GROUP_FIRST(JGROUP)
      LAST_TYPE  = GROUP_LAST (JGROUP)

      DO JOBT = FIRST_TYPE,LAST_TYPE   !  Loop over types in group

        IF (LENACT(JOBT) >  0) THEN   !  Any obs for this type

          LAST_OBS = FIRST_OBS + LENACT(JOBT) - 1

!         Get analysis parameters for this type
          ITYPE = TYPE_INDEX(JOBT)

          TIMEB         = DEF_TIMEB(ITYPE)
          TIMEA         = DEF_TIMEA(ITYPE)
          RADINF        = DEF_RADINF(ITYPE)
          CSCALE_START  = DEF_CSCALE_START(ITYPE)
          CSCALE_OBTIME = DEF_CSCALE_OBTIME(ITYPE)
          CSCALE_END    = DEF_CSCALE_END(ITYPE)

      DO JOB = FIRST_OBS,LAST_OBS
!       WRK1= time relative obs time (+ve before and -ve after)
        WRK1(JOB) = OBS_TIME(JOB) - (TIMESTEP_NO-1)*TIMESTEP/60.0
        TF(JOB)=0.0
        S(JOB)=1.0
!       note that these initialised values will not change for an ob
!       with relative time outside the insertion period (as required)
      ENDDO ! JOB

      DO JOB = FIRST_OBS,LAST_OBS
!     calculate ramp functions between start of insertion & ob time
!     data on edge of insertion period not fetched by GETOBS, so may
!     be discounted here too.
      IF (WRK1(JOB) >  0.0 .AND. WRK1(JOB) <= TIMEB) THEN
       IF(MRAMPFN == 2)THEN
! quadratic function
       TF(JOB) =                                                        &
     & WRK1(JOB)*WRK1(JOB)*(TIMEF_START-TIMEF_OBTIME)/(TIMEB*TIMEB)     &
     &  + TIMEF_OBTIME
       ELSEIF(MRAMPFN == 1)THEN
! linear function
       TF(JOB) =                                                        &
     & WRK1(JOB)*(TIMEF_START-TIMEF_OBTIME)/TIMEB + TIMEF_OBTIME
       ENDIF

       S(JOB) =                                                         &
     & WRK1(JOB) * ((CSCALE_START-CSCALE_OBTIME)/TIMEB) + CSCALE_OBTIME
      ENDIF
      ENDDO ! JOB

      DO JOB = FIRST_OBS,LAST_OBS
!     calculate ramp functions between ob time and end of insertion
      IF (WRK1(JOB) <= 0.0 .AND. WRK1(JOB) >= -TIMEA) THEN
       IF(MRAMPFN == 2)THEN
! quadratic function
       TF(JOB) =                                                        &
     & WRK1(JOB)*WRK1(JOB)*(TIMEF_END-TIMEF_OBTIME)/(TIMEA*TIMEA)       &
     &  + TIMEF_OBTIME
       ELSEIF(MRAMPFN == 1)THEN
! linear function
       TF(JOB) =                                                        &
     & WRK1(JOB) * (TIMEF_OBTIME-TIMEF_END)/TIMEA + TIMEF_OBTIME
       ENDIF

       S(JOB) =                                                         &
     & WRK1(JOB) * ((CSCALE_OBTIME-CSCALE_END)/TIMEA) + CSCALE_OBTIME
      ENDIF
      ENDDO ! JOB

!     WRK1 can be reused.

!     IWK1,WRK1 CAN BE REUSED.

      IF (MOD(MHCORFN/16,2) == 1) THEN
        DO JOB = FIRST_OBS,LAST_OBS
!         reduce time factor where scale is large so that the ramp
!         function specified holds for the area integrated weight
          WRK1(JOB) = CSCALE_OBTIME/S(JOB)
          WRK1(JOB) = WRK1(JOB)*WRK1(JOB)
          TF(JOB)   = TF(JOB)*WRK1(JOB)
        ENDDO ! JOB
      ENDIF

      ENDIF   !  Jump here if no obs for this type.

      FIRST_OBS = FIRST_OBS + LENACT(JOBT)
      ENDDO    !  End of JOBT loop over obs types in group.

      IF (LDIAGAC) THEN
       IF (LLDAC(1) .AND. NDAC >  0 .AND. IPASS == 1)THEN
       ENDIF
      ENDIF


      IF(IPASS == 1)THEN ! do stages 1 & 2
!L***                  <<<<    STAGE 1: SCALE    >>>>
!L*** 1.  calculate scale as smoothed average of that of local obs
!         using UM Doc 30 (6.4.17)
!     1.1 Initialize MODEL_INCR
      DO      JC=1,2
      DO      JLEV=1,NO_ANAL_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,JC)=0.
      ENDDO ! J
      ENDDO ! JLEV
      ENDDO ! JC

!     1.2 put values for each ob,level into S_INCR
      DO      JLEV=1,NO_ANAL_LEVS
      DO      JOB=1,LENOB
!       see remark in UM Doc P30 5.1 about the multiple uses of NORMF
!       NORMF here is Q1m of (5.1.2) i.e. 1/E**2
!       TF is used as in (6.1.4)
        S_INCR(JOB,JLEV,1)=NORMF(JOB,JLEV)*TF(JOB)
        S_INCR(JOB,JLEV,2)=NORMF(JOB,JLEV)*TF(JOB)*                     &
     &    S(JOB)*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV) ! in metres
!       FI allows a different scale for each level
!       To get results comparable with HORINF use FI_SCALE_FACTOR=1.
      ENDDO ! JOB
      ENDDO ! JLEV

      IF (FILT) THEN
!     1.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp
!     Loop numbered 130 to prevent fpp problems with compiler version 4
      DO 130    JOB=1,LENOB
!       the 4 passes through the loop over JC can be concurrent as
!       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
!       The version of HINTCF introduced with FI ensures that they are.
!FPP$ NODEPCHK L
        DO      JC=1,4
!       DO      J=1,2               ! the loops on J & JLEV should
!       DO      JLEV=1,NO_ANAL_LEVS !collapse into a single vector loop
                J=1                     ! collapse loops by hand
!       the JLEV loop is too short to multitask
!FPP$ NOCONCUR L
        DO      JLEV=1,NO_ANAL_LEVS*2   ! collapse loops by hand
         MODEL_INCR(NPPT(JOB,JC),JLEV,J)=                               &
     &   MODEL_INCR(NPPT(JOB,JC),JLEV,J)+S_INCR(JOB,JLEV,J)*CFPT(JOB,JC)
        ENDDO ! JLEV
!       ENDDO ! J                       ! collapse loops by hand
        ENDDO ! JC
 130  CONTINUE

!     1.4 divide by grid-box areas, to convert from transpose to adjoint
!         (The grid box area is Px in (6.4.15) ).
!         multiply by area integral of correlation function
      ZINTEGRAL=ZIOA*(FI_SCALE/Earth_Radius)**2 ! convert to radians
      I=0
      DO JROW=1,IROWS
        ZZ=ZINTEGRAL/COS_LAT(JROW)
        DO      JC=1,2
        DO      JLEV=1,NO_ANAL_LEVS
        Z=ZZ*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV)**2
        DO      J=1,ROW_LENGTH
          MODEL_INCR(I+J,JLEV,JC)=MODEL_INCR(I+J,JLEV,JC)*Z
        ENDDO ! J
        ENDDO ! JLEV
        ENDDO ! JC
        I=I+ROW_LENGTH
      ENDDO ! JROW

!     1.5 filter weighted scales, and weights
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFCSloop',3)
!FPP$ CNCALL
!     collapse JC & JLEV loops so that multitasking can be over both
!     DO      JC=1,2
!     DO      JLEV=1,NO_ANAL_LEVS
      DO      J=   0,NO_ANAL_LEVS*2-1 ! J loop replaces JC & JLEV loops
      JC=MOD(J,2)+1                   ! J loop replaces JC & JLEV loops
      JLEV=J/2+1                      ! J loop replaces JC & JLEV loops
      Z=FI_SCALE*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV)/Earth_Radius
!       These calls can be done concurrently
#if defined(GLOBAL)
! DEPENDS ON: rfcsg
        CALL RFCSG(MODEL_INCR(1,JLEV,JC),IROWS,ROW_LENGTH,M_GRID,       &
     &     COS_LAT,DLATMG,DLONGMG,Z,NPASS_RF)
#else
! DEPENDS ON: rfcsl
        CALL RFCSL(MODEL_INCR(1,JLEV,JC),IROWS,ROW_LENGTH,M_GRID,       &
     &   Z,COS_LAT,DLATMG,DLONGMG,Z,NPASS_RF)
#endif
!     ENDDO ! JLEV
!     ENDDO ! JC
      ENDDO ! J                       ! J loop replaces JC & JLEV loops
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFCSloop',4)

!     1.6 divide to give SCALE, using (6.4.15). metre=>radians using A.
      DO      JLEV=1,NO_ANAL_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,IP_SCALE)= (MODEL_INCR(J,JLEV,2) +            &
     &  0.1*FI_SCALE*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV))/              &
     &     ((MODEL_INCR(J,JLEV,1)+0.1)*Earth_Radius*FI_VAR_FACTOR)
      ENDDO ! J
      ENDDO ! JLEV


      ENDIF  !  FILT

!L***                  <<<<    STAGE 2: NORMF    >>>>
!L*** 2.  calculate data-density and hence Q to be used in VERTANL
!         using UM Doc 30 (6.4.15)
!     2.1 Initialize MODEL_INCR
      DO      JLEV=1,NO_WT_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,1)=0.
      ENDDO ! J
      ENDDO ! JLEV

!     2.2 put values for each ob,level into S_INCR
!     This need not be done; the values are already there from step 1.2
!     DO JLEV=1,NO_WT_LEVS
!     DO JOB=1,LENOB
!       S_INCR(JOB,JLEV,1)=NORMF(JOB,JLEV)*TF(JOB)
!     ENDDO ! JOB
!     ENDDO ! JLEV

!     2.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp
!     Loop numbered 230 to prevent fpp problems with compiler version 4
      DO 230 JOB=1,LENOB
!       the 4 passes through the loop over JC could be concurrent as
!       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
!       The HINTCF introducted with FI makes sure that they differ.
!FPP$ NODEPCHK L
        DO      JC=1,4
!       the JLEV loop is too short to multitask
!FPP$ NOCONCUR L
        DO      JLEV=1,NO_WT_LEVS
         MODEL_INCR(NPPT(JOB,JC),JLEV,1)=                               &
     &   MODEL_INCR(NPPT(JOB,JC),JLEV,1)+S_INCR(JOB,JLEV,1)*CFPT(JOB,JC)
        ENDDO ! JLEV
        ENDDO ! JC
 230  CONTINUE

      IF (FILT) THEN
!     2.4 divide by grid-box areas, to convert from transpose to adjoint
!         (The grid box area is Px in (6.4.15) ).
!         multiply by area integral of correlation function
!         Do not add weight given to background; this is done in VERTANL
      I=0
      DO JROW=1,IROWS
        Z=ZIOA/COS_LAT(JROW)
        DO      JLEV=1,NO_WT_LEVS
        DO      J=1,ROW_LENGTH
          MODEL_INCR(I+J,JLEV,1)=MODEL_INCR(I+J,JLEV,1)*Z*              &
     &                           MODEL_INCR(I+J,JLEV,IP_SCALE)**2
        ENDDO ! J
        ENDDO ! JLEV
        I=I+ROW_LENGTH
      ENDDO ! JROW

!     2.5 filter weights to get effective data density
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',3)
!FPP$ CNCALL
      DO      JLEV=1,NO_WT_LEVS
!       These calls can be done concurrently
#if defined(GLOBAL)
! DEPENDS ON: rfvsg
        CALL RFVSG(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,        &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
#else
! DEPENDS ON: rfvsl
        CALL RFVSL(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,0.,     &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
#endif
      ENDDO ! JLEV
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',4)

      ENDIF  !  FILT

      ELSE ! (IPASS == 2)THEN ! do stage  3
!L***                  <<<<    STAGE 3: INCR     >>>>
!L*** 3.  calculate INCR on model grid
!         using UM Doc 30 (6.4.16)
!     3.1 Initialize MODEL_INCR
      DO      JC=1,NO_ANAL_VAR
      DO      JLEV=1,NO_ANAL_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,JC)=0.
      ENDDO ! J
      ENDDO ! JLEV
      ENDDO ! JC

!     3.2 put values for each ob,level into S_INCR
!     Loop numbered 320 to prevent fpp problems with compiler version 4
      DO 320  JC=1,NO_ANAL_VAR
!FPP$ SELECT (CONCUR)
      DO      JLEV=1,NO_ANAL_LEVS
      DO      JOB=1,LENOB
!       NORMF (from VERTANL) is now Q2m of (6.1.3) i.e. (1/E**2) * Q
!       TF is used as in (6.1.1)
        S_INCR(JOB,JLEV,JC)=OBS_INCR(JOB,JLEV,JC)*TF(JOB)**2*           &
     &                         NORMF(JOB,JLEV)
      ENDDO ! JLEV
      ENDDO ! JOB
 320  CONTINUE

!     3.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp
!     Loop numbered 330 to prevent fpp problems with compiler version 4
      DO 330 JOB=1,LENOB
!       the 4 passes through the loop over JC could be concurrent as
!       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
!FPP$ NODEPCHK L
        DO      JC=1,4
!       DO      J=1,NO_ANAL_VAR     ! the loops on J & JLEV should
!       DO      JLEV=1,NO_ANAL_LEVS !collapse into a single vector loop
!       the JLEV*J loop is too short to multitask
                J=1                              ! collapse explicitly
!FPP$ NOCONCUR L
        DO      JLEV=1,NO_ANAL_LEVS*NO_ANAL_VAR  ! collapse explicitly
         MODEL_INCR(NPPT(JOB,JC),JLEV,J)=                               &
     &   MODEL_INCR(NPPT(JOB,JC),JLEV,J)+S_INCR(JOB,JLEV,J)*CFPT(JOB,JC)
        ENDDO ! JLEV  & collapsed J
!       ENDDO ! J                                ! collapse explicitly
        ENDDO ! JC
 330  CONTINUE

      IF (FILT) THEN
!     3.4 divide by grid-box areas, to convert from transpose to adjoint
!         (The grid box area is Px in (6.4.15) ).
!         multiply by area integral of correlation function
!             which is about 16*SCALE**2 for a 2-pass filter.
      I=0
      DO JROW=1,IROWS
        Z=ZIOA/COS_LAT(JROW)
!FPP$ NODEPCHK L
        DO      JC=1,NO_ANAL_VAR ! JC is never equal to IP_SCALE
        DO      JLEV=1,NO_ANAL_LEVS
        DO      J=1,ROW_LENGTH
          MODEL_INCR(I+J,JLEV,JC)=MODEL_INCR(I+J,JLEV,JC)*Z*            &
     &                            MODEL_INCR(I+J,JLEV,IP_SCALE)**2
        ENDDO ! J
        ENDDO ! JLEV
        ENDDO ! JC
        I=I+ROW_LENGTH
      ENDDO ! JROW
      IF (LDIAGAC) THEN
       IF (LLDAC(1).AND.NDAC >  0)THEN
        DO JLEV=1,NO_ANAL_LEVS
         DO JOB=1,LENOB ! interpolate incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,1)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,1)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,1)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,1)*CFPT(JOB, 4)
         ENDDO ! JOB
         IF(ILUV)THEN
         DO JOB=1,LENOB ! interpolate V incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,2)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,2)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,2)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,2)*CFPT(JOB, 4)
         ENDDO ! JOB
         ENDIF ! ILUV
        ENDDO ! JLEV
       ENDIF
      ENDIF

!     3.5 filter increments
      IF(ILUV)THEN  ! filter winds as a vector field
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVVloop',3)
!FPP$ CNCALL
       DO JLEV=1,NO_ANAL_LEVS
!       These calls can be done concurrently
#if defined(GLOBAL)
! DEPENDS ON: rfvvg
        CALL RFVVG(MODEL_INCR(1,JLEV,1),MODEL_INCR(1,JLEV,2),           &
     &   IROWS,ROW_LENGTH,M_GRID,COS_LAT,ROW1MGUV,                      &
     &   DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
#else
! DEPENDS ON: rfvvl
        CALL RFVVL(MODEL_INCR(1,JLEV,1),MODEL_INCR(1,JLEV,2),           &
     &   IROWS,ROW_LENGTH,M_GRID,COS_LAT,ROW1MGUV,                      &
     &   DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
#endif
       ENDDO ! JLEV
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVVloop',4)
      ELSE          ! filter scalar increment field
!FPP$ CNCALL
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',3)
       DO JLEV=1,NO_ANAL_LEVS
!       These calls can be done concurrently
#if defined(GLOBAL)
! DEPENDS ON: rfvsg
        CALL RFVSG(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,        &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
#else
! DEPENDS ON: rfvsl
        CALL RFVSL(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,0.,     &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
#endif
       ENDDO ! JLEV
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',4)
      ENDIF ! ILUV

      ENDIF ! FILT

      IF (LDIAGAC) THEN
       IF (LLDAC(1).AND.NDAC >  0)THEN
        DO JLEV=1,NO_ANAL_LEVS
         DO JOB=1,LENOB ! interpolate incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,1)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,1)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,1)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,1)*CFPT(JOB, 4)
         ENDDO ! JOB
         IF(ILUV)THEN
         DO JOB=1,LENOB ! interpolate V incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,2)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,2)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,2)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,2)*CFPT(JOB, 4)
         ENDDO ! JOB
         ENDIF ! ILUV
        ENDDO ! JLEV
       ENDIF
      ENDIF


      ENDIF ! IPASS

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('FI      ',4)
      RETURN
      END SUBROUTINE FI
#endif
