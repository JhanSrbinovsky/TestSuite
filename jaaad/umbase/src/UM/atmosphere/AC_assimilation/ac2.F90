#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: AC
!
!    Purpose : Main control routine for AC Scheme
!     This makes the analysis correction to model fields at the start
!   of each timestep of model integration.  The corrections are
!   calculated at analysis grid points,then interpolated to the model
!   grid before being added to the model fields. The corrections are
!   weighted sums of increments (ob-model) obtained at analysis grid
!   points. The weights depend on the distance between observation and
!   grid point, Observational error and the timeliness of the
!   observation.  Observations are analysed sequentially by observation
!   type in the order given by the array lact in comdeck comacp. A
!   vertical profile of observational increments and errors at
!   observation points are obtained before the weighted sums are
!   calculated. Horizontal filtering of increments on the model grid may
!   be done,and hydrostatic and geostrophic increments are calculated on
!   the model grid if appropriate.mean and rms diagnostics may be
!   calculated.
!
!    For Cray - Global and Limited area
!             - For Global ; Enable defs GLOBAL
!
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.1  22/12/92: F3 now always passed in AC & AC2 argument lists.
!LL              : WINDBAL arguments revised. New variable added,
!LL              : WB_THETA_INC, to control calculation & use of
!LL              : WINDBAL's potential temperature increments.
!LL              : (Phil Andrews)
!LL              : IANWTS and IANINC replaced by IPASS
!LL              : RELAXC called before ACSTASH
!LL              : ACSTASH args corrected      (R S Bell)
!LL 3.2  23/6/93 : Redimension ZSECG2 for portability
!LL              : and revise WINDBAL theta incrs (using HYDRST
!LL              : for WB Pincs from surface winds)
!LL              : and eliminate QA FORTRAN complaints    S Bell
!LL 3.3  13/9/93 : Correct mistake in setting WKLEN       S Bell
!LL 3.3  11/11/93: Pass OROG, LAND_MASK to VERTANL          N Richards
!LL 3.3  1/12/93 : New RAIN field args passed thro to vertanl   Bruce M
!LL 3.3  14/12/93: Revised WINDBAL control (for LAM option)  P Andrews
!LL 3.4  03/08/94: Tracer assimilation changes              Richard S
!LL 3.4   2/2/94 : Preserve RH on temperature assm           S Bell
!LL 3.4  07/09/94: Remove cloud water/ice from arg list and calls to
!LL              : HMRTORH. Remove RHCRIT from HMRTORH arg list.
!LL              : Pass convective cloud as arg.
!LL              : Include new args and call to HYDROL_INC
!LL              : New cloud regime boundaries in arg list, for
!LL              : passing to VERTANL
!LL              :                                 B Macpherson
!LL  3.4  9/9/94 : Include AEROSOL treatment of visibility.(P Clark)
!LL  4.0  3/5/95  : New call to LHN_INC, and new STASH diagnostics for
!LL               :     LS and LHN theta increments (Chris Jones)
!LL  4.0 10/7/95 : Add variable L_OBS_CHECK to allow bypass of check
!LL              : for no obs to assimilate in non-oper run. (G Bason)
!LL  4.1  6/6/96 : Loop over levels for call to MMSPT with
!LL              :      LHN theta increments (Chris Jones)
!LL  4.1 23/5/96 : Code to cope with single or multi-level hydrology,
!LL                or MOSES/Penman Monteith        (Bruce Macpherson)
!    4.1 18/06/96: Changes to cope with changes in STASH addressing
!                  Author D.M. Goddard.
!    4.2 25/11/96: T3E mods. Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!    4.3 21/1/97: Fix DIAGO if tests and MMSPT prints
!                 change ro use SWAPBOUNDS_shmem. Stuart Bell
!    4.3 17/5/97: Optimisation of ANAL_SUM and preparation
!                 for high resolution.  Deborah Salmond
!    4.4 17/10/97: Change to enable running on odd numbers of PE's
!                                       Deborah Salmond
!    4.5 08/04/98: Bug fix DIAGOPR barriers & DIAGO calls  Stuart Bell
!    4.5     Mar 98: Changes for cloud assimilation with mixed
!                    phase microphysics.      Bruce Macpherson
!    4.5 01/05/98  Restrict murk aerosol calculations to aerosol
!                  levels=boundary levels. P.Clark
!    5.2 30/11/00  changes to prepare for LHN/MOPS functionality
!                  as only valid one.           Bruce Macpherson
!    5.3 16/08/01  Remove CALLs to SWAPBOUNDS_sum.
!                  Remove model_status and zero ob check.
!                                       Adam Maycock
!    5.3 05/12/01  Remove reference to the shmcomm and iovarsac include
!                  files, and to all related variables.       S. Cusack
!    6.0 19/06/03  Remove non-MPP parts of code. T. White
!    6.0 11/09/03  Removed double ? for IBM cpp.              P.Dando
!    6.0 10/10/03  Replace SHMEM with GCOM and remove
!                  timings for SX6. Clive Jones
!    6.2 21/10/05  Replace GSYNCs with SSYNCs. P.Selwood
!    6.2 11/01/06  Replace large hard-wired obs(num)dim
!                  with dynamic allocation. R Barnes
!    6.2 08/11/05  Pass through l_eacf logical. Damian Wilson
!
!LL
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!*
!*L  Arguments and declarations:

!*L  Arguments & declarations:
      SUBROUTINE AC2(P_LEVELS, Q_LEVELS, BL_LEVELS,                     &
     &               ROW_LENGTH, P_ROWS,                                &
     &               P_FIELD,                                           &
     &               TIMESTEP_NO, ITER_NO, TIMESTEP,                    &
     &               OBS, TNDV,                                         &
     &               EXNER, PSTAR, THETA, RH, QCL, QCF,                 &
     &               CONV_CLD, LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,  &
     &               D_THETA_DT_CONV,D_THETA_DT_LS,                     &
     &               LAYER_CLOUD,PRESSURE,                              &
     &               RHCRIT, L_eacf,                                    &
     &               OBS_NO, LENOB, NO_ANAL_LEVS, NO_WT_LEVS,           &
     &               NO_ANAL_VAR,                                       &
     &               LENAG, LENMG, WKLEN, INC_TYPE,                     &
     &               NAVT, JGROUP, LWIND, IACTF, IACTL, DG_ONLY,        &
     &               STINDEX, STLIST, LEN_STLIST, SI, SF,               &
     &               STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,         &
     &               STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,             &
     &               ICODE, CMESSAGE)

      IMPLICIT NONE

! Global parameters:
#include "acparm.h"
#include "cntlatm.h"
#include "parvars.h"
#include "mppac.h"

!  For load balancing
      integer lenob_com,lenob_pos,lenob_average,lenob_len,lenob_total   &
     &        ,lenob_newpos,lenob_rem
      common/load_bal_com/lenob_com,lenob_pos,lenob_average,lenob_total &
     &        ,lenob_newpos


! Import arguments (INTENT=IN):
      INTEGER        P_LEVELS                ! Total number of levels
      INTEGER        Q_LEVELS                ! Number of wet levels
      INTEGER        BL_LEVELS               ! Bdy layer levels
      INTEGER        ROW_LENGTH              ! Number of points on row
      INTEGER        P_ROWS                  ! Number of rows (pstar)
      INTEGER        P_FIELD                 ! Number of points in
                                             ! mass field
      INTEGER        TNDV
      INTEGER        TIMESTEP_NO             ! Current model timestep
      INTEGER        ITER_NO                 ! Iteration number
                                             ! excluding diagnostic
      INTEGER        LENMG                   ! Length of model grid
      INTEGER        LENAG                   ! Length of analysis grid
      INTEGER        WKLEN                   ! Length of 2nd dimension
                                             ! of array for derived incs
                                             ! produced by HTDRST,
                                             ! GEOSTR & WINDBAL.
      INTEGER        INC_TYPE                ! Index for data types in
                                             ! MODEL_INCR
      INTEGER        LENOB                   ! No of obs in group
      INTEGER I,IPROC
      INTEGER        NO_ANAL_LEVS            ! No of analysis levels
                                             ! for group of obs types
      INTEGER        NO_ANAL_VAR             ! No of variables being
                                             ! analysed (2 for winds)
      INTEGER        NO_WT_LEVS              ! No of weight levels
                                             ! for group of obs types
      INTEGER        NAVT                    ! Analysis variable type
                                             ! 1/2/3/4/5 =
                                             ! p*/theta/winds/rh/precip
      INTEGER        JGROUP                  ! Loop counter over groups
      INTEGER        IACTF                   ! First obs type in group
      INTEGER        IACTL                   ! Last  obs type in group
      INTEGER        OBS_NO(LENOB+1)

      LOGICAL        LWIND                   ! Indicator if group
                                             ! has wind obs types
      LOGICAL        DG_ONLY                 ! Indicator if diagnostic
                                             ! only iteration
      LOGICAL        L_eacf                  ! Use emp. adjusted
                                             ! cloud fraction.

      REAL           TIMESTEP                ! Timestep in seconds
      REAL           OBS(TNDV+1)             ! Observation data (set up
                                             ! in RDOBS)
      REAL           CONV_CLD(P_FIELD, Q_LEVELS)  ! conv cld amount
      REAL           LS_RAIN(P_FIELD)        ! large scale rain rate
      REAL           LS_SNOW(P_FIELD)        ! large scale snow rate
      REAL           CONV_RAIN(P_FIELD)      ! convective rain rate
      REAL           CONV_SNOW(P_FIELD)      ! convective snow rate
                                             ! above rates diagnostic
      REAL           D_THETA_DT_CONV(P_FIELD,Q_LEVELS)
                                       ! convective latent heating rate
      REAL           D_THETA_DT_LS(P_FIELD,Q_LEVELS)
                                      ! large scale latent heating rate

      REAL           RHCRIT(Q_LEVELS)        ! Critical RH array
                                             ! for cloud

! Import / export arguments (INTENT=INOUT):
      REAL           EXNER(P_FIELD, P_LEVELS) ! exner on theta levels
      REAL           LAYER_CLOUD(P_FIELD, Q_LEVELS) ! as name says
      REAL           PRESSURE (P_FIELD, P_LEVELS) ! p on theta levels
      REAL           PSTAR(P_FIELD)          ! Prognostic variable P*
      REAL           THETA(P_FIELD, P_LEVELS)! Prognostic variable
                                             ! theta
      REAL           RH(P_FIELD, Q_LEVELS)   ! Prognostic variable HMR
! for 2A cloud microphysics, but otherwise contains rh at this stage
      REAL           QCL(P_FIELD, Q_LEVELS)   ! Prognostic variable QCL
      REAL           QCF(P_FIELD, Q_LEVELS)   ! Prognostic variable QCF

! Export arguments (INTENT=OUT):
      INTEGER        ICODE                   ! Error code (0 = OK)

      CHARACTER*256  CMESSAGE                ! Error message

! Stash variables
      REAL           STASHWORK(*)

      INTEGER        LEN_STLIST
      INTEGER        NUM_STASH_LEVELS
      INTEGER        NUM_STASH_PSEUDO
      INTEGER        STINDEX(2, *)
      INTEGER        STLIST(LEN_STLIST, *)
      INTEGER        SI(*)
      INTEGER        STASH_LEVELS(NUM_STASH_LEVELS +1, *)
      INTEGER        STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO +1, *)

      LOGICAL        SF(*)

!*
!L Single iteration of analysis correction data assimilation
!L for observations or group of observations of variable type
!L defined by NAVT

! More global variables:
#include "comobs.h"
#include "comacp.h"
#include "comag.h"
#include "comacdg.h"

! Local (dynamic) arrays:
      INTEGER        MODEL_OBS_TYPE(LENOB+1)
                                         ! Model observation type No
      INTEGER        NE_AG_PT(LENOB+1)     ! Nearest analysis grid
                                            ! point to obs.
      INTEGER        NPPT(LENOB+1, 4)       ! horizontal model ->
                                            ! obs interp pointers

      REAL           OBS_LAT (LENOB+1)      ! Observation latitude
      REAL           OBS_LONG(LENOB+1)      ! "longitude
      REAL           OBS_TIME(LENOB+1)      ! "time
      REAL           OBS_INCR(LENOB+1, NO_ANAL_LEVS, NO_ANAL_VAR)!
                                            ! Obs - Model increments
                                            ! at obs points.
      REAL           NORMF(LENOB+1, NO_ANAL_LEVS)
                                               ! Normalization factor
      REAL           ANAL_INCR_LOCAL(LENMG*NO_ANAL_LEVS*NO_ANAL_VAR) !
                                             ! Accumulated increments
                                             ! on analysis grid
      REAL           MODEL_INCR(LENMG, NO_ANAL_LEVS, INC_TYPE) !
                                             ! Accumulated increments
                                             ! on model grid
      REAL           DRINCS(P_FIELD, WKLEN)  ! Array to hold derived
                                             ! increment fields
      REAL           CFPT(LENOB+1, 4)
                                             ! obs interp coeffs
      REAL           RELAX_CF(P_ROWS)        ! Relaxation coeffs for
                                             ! current group

      LOGICAL        LMISSD(LENOB+1, NO_ANAL_LEVS)
                                             ! missing data in obs.
      REAL           pstgs(P_FIELD)          ! PSTAR at U points
      integer dum1,dum3
      real dum2
      integer bc,bcount
      integer lev_ave,lev_start,lev_end,num_levs,ntimes,lev_times
      integer lev_rem
      integer num_incr_levs
      integer ic
      integer pe_group,iii,jjj
      integer num_groups,                                               &
     &    g_pe_start(0:3),g_pe_end(0:3),                                &
     &    g_lev_start(0:3),g_lev_end(0:3)
      common/split_pe/num_groups,                                       &
     &     g_pe_start,g_pe_end,g_lev_start,g_lev_end

! Local scalars:
      INTEGER        J                       ! Miscelaneous loop cntr
      INTEGER        JACT,JJACT              ! Loop counter over obs
                                             ! types in group
      INTEGER        JANL                    ! Loop counter over
                                             ! analysis levels
      INTEGER        JLEV,JJLEV              ! Loop counter over model
                                             ! levels
      INTEGER        JL,JJL
      INTEGER        JK,JJK
      INTEGER        JPOINT                  ! Loop counter over points

      INTEGER        IPASS                   ! 1 if weights pass
                                             ! 2 if increment pass
      INTEGER        IAG                     ! Offset for level in
                                             ! ANAL_INCR array.
      INTEGER        IJK                     ! Level number in HYDRST
                                             ! & GEOSTR loops
      INTEGER        INOBS                   ! No of obs for this type
      INTEGER        NPTOBT                  ! Offset to first obs of
                                             ! each type in group
      INTEGER        LENOBT                  ! No of obs for obs type
      INTEGER        NDV                     ! No of data values
      INTEGER        N_ROWS                  ! No of rows
      INTEGER        NO_ANAL_INCR_LEVS       ! No of analysis increment
                                             ! levels.
      INTEGER        MODE_HANAL              ! FI or HORINF for this
                                             ! group
      INTEGER        STASH_ITEM_NO           ! Stash Item Number
      LOGICAL        SURFACE_WIND            ! True if surface wind data
      LOGICAL ILUV


! Subroutine calls:
      EXTERNAL ADDINC, FI, GETOB2
      EXTERNAL RELAXC, SET_RELAX_CF, VERTANL
      EXTERNAL DIAGO, MMSPT
      EXTERNAL HINTCF, TIMER
#if !defined(GLOBAL)

#endif

      EXTERNAL AC_STASH
!- End of header section.

!L Subroutine structure:
!L 1.0 Initialize:
! Mode of Horizontal Analysis for this group
      MODE_HANAL = DEF_MODE_HANAL(GROUP_INDEX(JGROUP))



!L 2.0 Weights calculation starts here

      IF (LDIAGAC.AND.mype == 0) THEN
        PRINT '(/,A,(10I4))',                                           &
     &                       ' START OF WTS ANALYSIS PHASE FOR TYPES ', &
     &                                      (LACT(J), J = IACTF, IACTL)

      END IF

!L 2.1 Setup for interpolation model -> obs & vertical analysis
!L     Get common format data for all obs in group       (getob2)
! Vectorisation over all observation types in group
! Loop over obs types in group
      if(lenob /= 0)then
      NPTOBT = 0

      DO JACT = IACTF, IACTL
      JJACT=JACT
        INOBS  = NOBS(JACT)
        LENOBT = LENACT(JACT)

        IF (INOBS >  0 .AND. LENOBT >  0) THEN
! DEPENDS ON: getob2
          CALL GETOB2 (JJACT, OBS(MDISPOBT(JACT) +1), INOBS,            &
     &                OBS_LAT(NPTOBT +1), OBS_LONG(NPTOBT +1),          &
     &                OBS_TIME(NPTOBT +1), MODEL_OBS_TYPE(NPTOBT +1),   &
     &                OBS_NO(NPTOBT +1), LENOBT, ICODE, CMESSAGE)

         ! Check for error - exit routine if occured
          IF (ICODE  >   0) GOTO 999

        END IF

       ! Move pointer to next obs type in group
        NPTOBT = NPTOBT + LENOBT

      END DO

!L 2.2  Horizontal interpolation coefficients model -> obs  (hintcf)
! DEPENDS ON: hintcf
      CALL HINTCF(LWIND, LENOB, OBS_LAT, OBS_LONG, OBS_NO, ROW_LENGTH,  &
     &           P_ROWS, CFPT(1,1), CFPT(1,2), CFPT(1,3), CFPT(1,4),    &
     &           NPPT(1,1), NPPT(1,2), NPPT(1,3), NPPT(1,4), ICODE,     &
     &           CMESSAGE)

! Check for error - exit routine if occured
      IF (ICODE  >   0) GOTO 999

!L 2.3 Do vertical analysis
!L     Make increment & error factor vectors for each type in this
!L     group. The details of processing method depend on ob type
!L     so vectorization is over LENOBT obs in one type
!L     rather than LENOB obs in the group. NPTOBT gives the
!L     increment to point to each section in turn of those
!L     vectors which go over all obs in group.

! Loop over obs types in group
      endif ! lenob /= 0
      NPTOBT = 0

      DO JACT = IACTF, IACTL
      JJACT=JACT
        INOBS  = NOBS(JACT)
        LENOBT = LENACT(JACT)
!
!
!
      if(inobs == 0.or.lenobt == 0)then
      if(LDIAGAC)then
      if(LACT(JACT) == 506)then
! DEPENDS ON: diagopr
        CALL DIAGOPR (1,dum1,dum2,LENOBT,LENOB,NDV,                     &

     &                LMISSD,dum3,NPTOBT,NO_ANAL_LEVS)
      else
      BCOUNT=2
      if(LACT(JACT) == 101.or.                                          &
     &   (LACT(JACT) >= 202.and.LACT(JACT) <= 204).or.                  &
     &   (LACT(JACT) >= 302.and.LACT(JACT) <= 306).or.                  &
     &   (LACT(JACT) >= 402.and.LACT(JACT) <= 404).or.                  &
     &   LACT(JACT) == 406.or.                                          &
     &   LACT(JACT) == 901)BCOUNT=1
!
! BARRIERS FOR DIAGO
!
      do bc=1,bcount
      do jlev=1,glsize(3,fld_type_p)
      do i=0,8
      s_stat(jlev,i)=0.0
      enddo
      enddo
      if(LACT(JACT) == 406)then
! dummy diago call
! DEPENDS ON: diago
       CALL DIAGO ('MULTI-LEVEL', LACT(JACT),6,OBS_INCR,NORMF,          &
     &              OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT,   &
     &              NO_ANAL_LEVS, NO_ANAL_VAR)
      else
! dummy diago call
! DEPENDS ON: diago
       CALL DIAGO ('VAN?', LACT(JACT), 3+bc-bcount,                     &
     &              OBS_INCR, NORMF,                                    &
     &              OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT,   &
     &              NO_ANAL_LEVS, NO_ANAL_VAR)
      endif ! LACT(JACT) == 406
      enddo
      endif ! LACT(JACT) == 506
      endif ! LDIAGAC
      else !inobs == 0.or.lenobt == 0

        IF (INOBS  >   0 .AND. LENOBT  >   0) THEN
          NDV   = NDATAV(JACT) - NDVHDR
          IPASS = 1

! DEPENDS ON: vertanl
          CALL VERTANL (PSTGS,                                          &
     &          JJACT, IPASS, LENOBT, NDV, OBS_NO(NPTOBT +1),           &
     &          MODEL_OBS_TYPE(NPTOBT +1), OBS_LAT(NPTOBT +1),          &
     &          OBS_LONG(NPTOBT +1), OBS(MDISPOBT(JACT) +1), INOBS,     &
     &          EXNER, PSTAR, THETA, RH, QCL, QCF, CONV_CLD,            &
     &          LAYER_CLOUD,PRESSURE,                                   &
     &          LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,                 &
     &          RHCRIT, L_eacf,                                         &
     &          P_FIELD, MODEL_INCR, LENMG, NO_WT_LEVS,                 &
     &          CFPT(NPTOBT+1, 1), CFPT(NPTOBT+1, 2),                   &
     &          CFPT(NPTOBT+1, 3), CFPT(NPTOBT+1, 4),                   &
     &          NPPT(NPTOBT+1, 1), NPPT(NPTOBT+1, 2),                   &
     &          NPPT(NPTOBT+1, 3), NPPT(NPTOBT+1, 4),                   &
     &          OBS_INCR, NORMF, LMISSD,                                &
     &          P_LEVELS, Q_LEVELS, BL_LEVELS, ROW_LENGTH, P_ROWS,      &
     &          LENOB, NPTOBT, NO_ANAL_LEVS, NO_ANAL_VAR,               &
     &          ICODE, CMESSAGE)

         ! Check for error - exit routine if occured
          IF (ICODE  >   0) GOTO 999

        END IF

       ! Move pointer to next obs type in group
        NPTOBT = NPTOBT + LENOBT

      endif !inobs == 0.or.lenobt == 0
      END DO

! Vectorisation elsewhere over all obs in group
      NPTOBT = 0
      LENOBT = LENOB

      IF (LDIAGAC) THEN
        IF (IACTF  /=  IACTL) THEN
         ! The group has >1 type in it, so group stats are worthwhile.
         ! Print mean & rms stats for all obs in group

         IF (LLDAC(2)) THEN
! DEPENDS ON: diago
           CALL DIAGO ('AC', LACT(IACTF), 1, OBS_INCR, NORMF,           &
     &                OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT, &
     &                NO_ANAL_LEVS, NO_ANAL_VAR)

          END IF
        END IF
      END IF


!L 2.4 Do horizontal analysis (for weights)
          IPASS = 1

      IF (MODE_HANAL  ==  2) THEN   !  Use FI method
      if(lenob /= 0)then

! DEPENDS ON: fi
        CALL FI (IPASS, TIMESTEP_NO, TIMESTEP, NAVT, JGROUP,            &
     &          NO_ANAL_VAR, LENOB, NO_ANAL_LEVS, NO_WT_LEVS, 1, LENMG, &
     &          ROW_LENGTH, P_ROWS, OBS_NO, OBS_LAT, OBS_LONG, OBS_TIME,&
     &          NORMF, OBS_INCR, CFPT, NPPT, INC_TYPE, MODEL_INCR,      &
     &          PSTAR,                                                  &
     &          P_FIELD, P_LEVELS, ICODE,                               &
     &          CMESSAGE)
      else
      do jlev=1,NO_ANAL_LEVS
      do j=1,lenmg
      MODEL_INCR(j,jlev,1)=0.0
      enddo
      enddo
      endif

       ! Check for error - exit from routine if occured
        IF (ICODE  >   0) GOTO 999

      END IF ! End of FI horizontal analysis step

      IF (LDIAGAC) THEN
        IF(mype == 0)PRINT *, 'END OF WTS ANALYSIS PHASE'

      END IF

! Save all levels of wts on model grid

      DO JLEV = 1, NO_WT_LEVS   !  Loop over no of weight levels.
      JJLEV=JLEV


       ! Save weights field

        STASH_ITEM_NO = 200+NAVT
        IF (SF(STASH_ITEM_NO)) THEN
! DEPENDS ON: ac_stash
          CALL AC_STASH (STASH_ITEM_NO,JJLEV, P_LEVELS, JGROUP,         &
     &      N_GROUPS, TIMESTEP_NO, STINDEX, STLIST, LEN_STLIST, SI, SF, &
     &      STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,                  &
     &      STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO, MODEL_INCR(1,JLEV,1),&
     &      LENMG, 'Weights Field', ICODE, CMESSAGE)

        END IF

       ! Get statistics on sum of weights
        IF (LDIAGAC .AND. LLDAC(4)) THEN
         IF (NAVT  ==  4) THEN
! DEPENDS ON: mmspt
            CALL MMSPT (MODEL_INCR(1,JLEV,1), JJLEV, 0,                 &
     &                          'SUM OF WTS - RH ', ROW_LENGTH, P_ROWS)

          ELSE IF (NAVT  ==  5) THEN
! DEPENDS ON: mmspt
            CALL MMSPT (MODEL_INCR(1,JLEV,1), JJLEV, 0,                 &
     &                          'SUM OF WTS - PR ', ROW_LENGTH, P_ROWS)
          ELSE      ! NAVT values other than 4 or 5 excluded
             ICODE = 1
             CMESSAGE = 'AC : Should not reach here (3)'
             GOTO 999   ! RETURN

          END IF
        END IF
      END DO   ! End of loop over no of weight levels.

!L 2.6 Calculate ob normalisation factors
! Vectorization here is over types within groups
      NPTOBT = 0

      DO JACT = IACTF, IACTL
      JJACT=JACT
        INOBS  = NOBS(JACT)
        LENOBT = LENACT(JACT)
!
!
!
      if(inobs == 0.or.lenobt == 0)then
      if(LDIAGAC)then
      if(LACT(JACT) /= 406.and.LACT(JACT) /= 506)then
      BCOUNT=2
      if(LACT(JACT) == 101.or.                                          &
     &   (LACT(JACT) >= 202.and.LACT(JACT) <= 204).or.                  &
     &   (LACT(JACT) >= 302.and.LACT(JACT) <= 306).or.                  &
     &   (LACT(JACT) >= 402.and.LACT(JACT) <= 404).or.                  &
     &   LACT(JACT) == 901)BCOUNT=1
!
! BARRIERS FOR DIAGO
!
      do bc=1,bcount
      do jlev=1,glsize(3,fld_type_p)
      do i=0,8
      s_stat(jlev,i)=0.0
      enddo
      enddo
! dummy diago call
! DEPENDS ON: diago
       CALL DIAGO ('VAN?', LACT(JACT), 5+bc-bcount,                     &
     &              OBS_INCR, NORMF,                                    &
     &              OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT,   &
     &              NO_ANAL_LEVS, NO_ANAL_VAR)
      enddo
      endif ! LACT(JACT) /= 406 (and 506)
      endif ! LDIAGAC
      else ! lenob == 0

        IF (INOBS >  0 .AND. LENOBT >  0) THEN
          NDV   = NDATAV(JACT) - NDVHDR
          IPASS = 2

! DEPENDS ON: vertanl
          CALL VERTANL (PSTGS,                                          &
     &           JJACT, IPASS, LENOBT, NDV, OBS_NO(NPTOBT +1),          &
     &           MODEL_OBS_TYPE(NPTOBT +1), OBS_LAT(NPTOBT +1),         &
     &           OBS_LONG(NPTOBT+1), OBS(MDISPOBT(JACT) +1), INOBS,     &
     &          EXNER, PSTAR, THETA, RH, QCL, QCF, CONV_CLD,            &
     &          LAYER_CLOUD,PRESSURE,                                   &
     &           LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,                &
     &           RHCRIT, L_eacf,                                        &
     &           P_FIELD, MODEL_INCR, LENMG, NO_ANAL_LEVS,              &
     &           CFPT(NPTOBT+1, 1), CFPT(NPTOBT+1, 2),                  &
     &           CFPT(NPTOBT+1, 3), CFPT(NPTOBT+1, 4),                  &
     &           NPPT(NPTOBT+1, 1), NPPT(NPTOBT+1, 2),                  &
     &           NPPT(NPTOBT+1, 3), NPPT(NPTOBT+1, 4),                  &
     &           OBS_INCR, NORMF, LMISSD,                               &
     &          P_LEVELS, Q_LEVELS, BL_LEVELS, ROW_LENGTH, P_ROWS,      &

     &           LENOB, NPTOBT, NO_ANAL_LEVS, NO_ANAL_VAR,              &
     &           ICODE, CMESSAGE)

         ! Check for errors:
          IF (ICODE  >   0) GOTO 999

        END IF

        NPTOBT = NPTOBT + LENOBT
      endif ! lenob == 0

      END DO

! Vectorisation elsewhere over all obs in group
      NPTOBT = 0
      LENOBT = LENOB

!L 3.0 Increments calculation starts here

      IF (LDIAGAC.AND.mype == 0) THEN
        PRINT '(/,A,(10I4))',                                           &
     &           ' START OF INCS ANALYSIS PHASE FOR TYPES  ',           &
     &           (LACT(J),J=IACTF,IACTL)

      END IF

!L 3.1 Horizontal analysis step (increments)
          IPASS = 2

      IF (MODE_HANAL  ==  2) THEN   ! Use FI method

      if(lenob /= 0)then
! DEPENDS ON: fi
        CALL FI (IPASS, TIMESTEP_NO, TIMESTEP, NAVT, JGROUP,            &
     &         NO_ANAL_VAR, LENOB, NO_ANAL_LEVS, NO_WT_LEVS, 1, LENMG,  &
     &         ROW_LENGTH, P_ROWS, OBS_NO, OBS_LAT, OBS_LONG, OBS_TIME, &
     &         NORMF, OBS_INCR, CFPT, NPPT, INC_TYPE, MODEL_INCR,       &
     &         PSTAR,                                                   &
     &         P_FIELD, P_LEVELS, ICODE,                                &
     &         CMESSAGE)
      else
      do jlev=1,NO_ANAL_LEVS
      do j=1,lenmg
      MODEL_INCR(j,jlev,1)=0.0
      enddo
      enddo
      endif

       ! Check for errors:
        IF (ICODE  >   0) GOTO 999

      END IF ! End of FI horizontal analysis step

      IF (LDIAGAC) THEN
        IF(mype == 0)PRINT *, 'END OF INCS ANALYSIS PHASE'

      END IF

!L 3.3 Set up Relaxation Coefficients for this group
        N_ROWS = P_ROWS

! DEPENDS ON: set_relax_cf
      CALL SET_RELAX_CF (JGROUP, N_ROWS, RELAX_CF, LWIND, TIMESTEP,     &
     &     TIMESTEP_NO, ITER_NO, ICODE, CMESSAGE)

! Check for errors:
      IF (ICODE  >   0) GOTO 999

!L 3.4 Interpolation from analysis grid to model grid
!L     Loop over levels to interpolate incrs from analysis grid to
!L     model grid. Add incrs to model.
      DO JLEV = 1, NO_ANAL_LEVS
      JJLEV=JLEV

!L     3.4.1 Scale increment by relaxation coefficent

! DEPENDS ON: relaxc
          CALL RELAXC (MODEL_INCR(1,JLEV,1), LENMG, ROW_LENGTH,         &
     &      P_ROWS, RELAX_CF, ICODE, CMESSAGE)

         ! Check for errors:
          IF (ICODE  >   0) GOTO 999


!L     3.4.2 Save increment fields into output file

          STASH_ITEM_NO = 210+NAVT
          IF (SF(STASH_ITEM_NO)) THEN
! DEPENDS ON: ac_stash
            CALL AC_STASH (STASH_ITEM_NO, JJLEV, P_LEVELS, JGROUP,      &
     &        N_GROUPS, TIMESTEP_NO, STINDEX, STLIST, LEN_STLIST, SI,   &
     &        SF, STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,            &
     &        STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                    &
     &        MODEL_INCR(1,JLEV,1), LENMG, 'Increment Field', ICODE,    &
     &        CMESSAGE)

          END IF


!L     3.4.10 Humidity increments
        IF (NAVT  ==  4) THEN

         ! Add RH increments to RH fields
! DEPENDS ON: addinc
          CALL ADDINC (RH(1,JLEV), MODEL_INCR(1,JLEV,1), LENMG,         &
     &      ROW_LENGTH, P_ROWS, NAVT, ICODE, CMESSAGE)

         ! Check for errors:
          IF (ICODE  >   0) GOTO 999

         ! Statistics of increments on model grid
          IF (LDIAGAC .AND. LLDAC(4)) THEN
! DEPENDS ON: mmspt
            CALL MMSPT (MODEL_INCR(1,JLEV,1), JJLEV, 0,                 &
     &        'RH INCREMENTS   ', ROW_LENGTH, P_ROWS)

          END IF
        END IF   ! NAVT = 4

!L     3.4.11 latent heat nudging
        IF (NAVT  ==  5) THEN
!  Stash LS latent heating theta increments if required (rates in K/s)
            STASH_ITEM_NO = 271

            IF (SF(STASH_ITEM_NO)) THEN
              DO JL = 1 , Q_LEVELS
                JJL = JL
! DEPENDS ON: ac_stash
                CALL AC_STASH (STASH_ITEM_NO, JJL, Q_LEVELS,            &
     &            JGROUP, N_GROUPS, TIMESTEP_NO,                        &
     &            STINDEX, STLIST, LEN_STLIST, SI, SF, STASHWORK,       &
     &            STASH_LEVELS, NUM_STASH_LEVELS,                       &
     &            STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                &
     &            D_THETA_DT_LS(1,JL),P_FIELD,                          &
     &            'LS latent heating rates',                            &
     &            ICODE, CMESSAGE)
              ENDDO    ! JL
            ENDIF

          IF (L_LHN) THEN
!C Call to new subroutine LHN_INC
!
! DEPENDS ON: lhn_inc
      CALL LHN_INC(D_THETA_DT_CONV,D_THETA_DT_LS,                       &
     &             LS_RAIN,LS_SNOW,CONV_RAIN,                           &
     &             CONV_SNOW,P_FIELD,Q_LEVELS,MODEL_INCR,TIMESTEP,      &
     &             P_ROWS,ROW_LENGTH,                                   &
     &             LHN_RANGE,                                           &
     &             DRINCS,LENMG,ICODE,CMESSAGE)
!C Add on increments
            DO J=1,Q_LEVELS
! DEPENDS ON: addinc
              CALL ADDINC(THETA(1,J), DRINCS(1,J), P_FIELD, ROW_LENGTH, &
     &                    P_ROWS, NAVT,ICODE,CMESSAGE)
            ENDDO     !J

!  Stash LHN Theta Increments if required (amounts in K)
            STASH_ITEM_NO = 272

            IF (SF(STASH_ITEM_NO)) THEN
              DO JL = 1 , Q_LEVELS
                JJL = JL
! DEPENDS ON: ac_stash
                CALL AC_STASH (STASH_ITEM_NO, JJL, Q_LEVELS,            &
     &            JGROUP, N_GROUPS, TIMESTEP_NO,                        &
     &            STINDEX, STLIST, LEN_STLIST, SI, SF, STASHWORK,       &
     &            STASH_LEVELS, NUM_STASH_LEVELS,                       &
     &            STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                &
     &            DRINCS(1,JL),P_FIELD,'LHN theta increments',          &
     &            ICODE, CMESSAGE)
              ENDDO    ! JL
            ENDIF

            IF (LDIAGAC .AND. LLDAC(4) ) THEN
            DO JL = 1 , Q_LEVELS
            JJL = JL
! DEPENDS ON: mmspt
              CALL MMSPT(DRINCS(1,JL), JJL, 0,                          &
     &                    'LHN Theta incrs ', ROW_LENGTH, P_ROWS)
            ENDDO     ! JL
            ENDIF
          ENDIF       ! L_LHN

        ENDIF    ! NAVT = 5

      END DO   ! JLEV


!L 5.0 Exit from AC2
 999  CONTINUE
      RETURN
      END SUBROUTINE AC2
#endif
