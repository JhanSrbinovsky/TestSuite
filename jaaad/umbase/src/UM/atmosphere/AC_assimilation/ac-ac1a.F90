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
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!*
!*L  Arguments and declarations:
      SUBROUTINE AC (                                                   &
     &  P_LEVELS, Q_LEVELS, ROW_LENGTH, P_ROWS, BL_LEVELS,              &
     &  TNDVMAX, NOBSMAX, P_FIELD,                                      &
     &  TIMESTEP_NO, TIMESTEP,                                          &
     &  EXNER, PSTAR, PRESSURE,                                         &
     &  THETA, RH, QCL, QCF,                                            &
     &  CONV_CLD, LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,               &
     &  LAYER_CLOUD, D_THETA_DT_CONV,D_THETA_DT_LS, RHCRIT, L_eacf,     &
     & OBS_FLAG,OBS,                                                    &
     &  STINDEX,                                                        &
     &  STLIST, LEN_STLIST, SI, SF,                                     &
     &  STASHWORK, STASH_LEVELS,                                        &
     &  NUM_STASH_LEVELS, STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,        &
#include "argppx.h"
     &  ICODE, CMESSAGE)

      IMPLICIT NONE
#include "acparm.h"
#include "cntlatm.h"

! Imported arguments (INTENT=IN):
      INTEGER        P_LEVELS                ! Total number of levels
      INTEGER        Q_LEVELS                ! Number of wet levels
      INTEGER        BL_LEVELS               ! Bdy layer levels
      INTEGER        ROW_LENGTH              ! Number of points on row
      INTEGER        P_ROWS                  ! Number of rows (pstar)
      INTEGER        TNDVMAX                 ! Total no of data values
                                             ! in all obs files
      INTEGER        NOBSMAX                 ! Total no of obs
                                             ! in all obs files
      INTEGER        P_FIELD                 ! Number of points in
                                             ! mass field
      INTEGER        TIMESTEP_NO             ! Current model timestep
      REAL           CONV_CLD(P_FIELD, Q_LEVELS)  ! conv cld amount
      REAL           LS_RAIN(P_FIELD)        ! large scale rain rate
      REAL           LS_SNOW(P_FIELD)        ! large scale snow rate
      REAL           CONV_RAIN(P_FIELD)      ! convective rain rate
      REAL           CONV_SNOW(P_FIELD)      ! convective snow rate
                                             ! above rates diagnostic
      REAL           TIMESTEP                ! Timestep in seconds
      REAL           RHCRIT(Q_LEVELS)        ! Critical rh array
                                             ! for cloud
      REAL           D_THETA_DT_CONV(P_FIELD,Q_LEVELS)
                                       ! convective latent heating rate
      REAL           D_THETA_DT_LS(P_FIELD,Q_LEVELS)
                                      ! large scale latent heating rate


! Stash variables:
      REAL           STASHWORK(*)

      INTEGER        LEN_STLIST              ! Length of STLIST
      INTEGER        NUM_STASH_LEVELS        ! No of Stash levels lists
      INTEGER        NUM_STASH_PSEUDO        ! No of Stash pseudo lists
      INTEGER        STINDEX(2, *)
      INTEGER        STLIST(LEN_STLIST, *)
      INTEGER        SI(*)                   ! Stash Index
      INTEGER        STASH_LEVELS(NUM_STASH_LEVELS +1, *)
                                             ! Levels lists
      INTEGER        STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO +1, *)
                                             ! Pseudo lists

      LOGICAL        SF(*)                   ! Stash Flags

! Import/export arguments (INTENT=INOUT):
      REAL           EXNER(P_FIELD, P_LEVELS) ! exner on theta levels
      REAL           PRESSURE(P_FIELD, P_LEVELS) ! p  on theta levels
      REAL           LAYER_CLOUD (P_FIELD, Q_LEVELS) ! as name says
      REAL           PSTAR(P_FIELD)          ! Prognostic variable pstar
      REAL           THETA(P_FIELD, P_LEVELS)! Prognostic variable
                                             ! theta
      REAL           RH(P_FIELD, Q_LEVELS)   ! Prognostic variable hmr
      REAL           QCL(P_FIELD, Q_LEVELS)   ! Prognostic variable qcl
      REAL           QCF(P_FIELD, Q_LEVELS)   ! Prognostic variable qcf

! Exported arguments (INTENT=OUT):
      INTEGER        ICODE                   ! Non zero for failure

      CHARACTER*256  CMESSAGE                ! Reason for failure



! Global variables:
#include "comobs.h"
#include "docobs.h"
#include "comacp.h"
#include "docacp.h"
#include "comag.h"
#include "docag.h"
#include "comacdg.h"
#include "docacdg.h"

! Local (dynamic) arrays:
#include "parvars.h"
#include "mppac.h"
!  For load balancing
      integer lenob_com,lenob_pos,lenob_average,lenob_len,lenob_total   &
     &        ,lenob_newpos,lenob_rem
      common/load_bal_com/lenob_com,lenob_pos,lenob_average,lenob_total &
     &        ,lenob_newpos
      integer iproc
      INTEGER :: OBS_FLAG(NOBSMAX)       ! Observation flags
      INTEGER        OBS_NO(NOBSMAX)       ! Observation numbers
      REAL  ::   OBS(TNDVMAX)            ! Observation data
! Local variables:
      INTEGER        TNDV                    ! Total no of data values
                                             ! for obs in assimilation
                                             ! time window
      INTEGER        TNOBS                   ! Total no of obs in
                                             ! assimilation time window
      INTEGER        NAVT                    ! Analysis variable type
                                             ! 1/2/3/4/5/6
                                             ! p*/theta/winds/rh/precip
                                             ! /tracer
      INTEGER        WKLEN                   ! Length of vertical
                                             ! dimension of array for
                                             ! derived increments made
                                             ! by HYDRST, GEOSTR &
                                             ! WINDBAL.
      INTEGER        NPTOBT                  ! Offset to first obs of
                                             ! each type in group
      INTEGER        NO_WT_LEVS              ! No of weight levels
                                             ! for group of obs types
      INTEGER        NO_ANAL_LEVS            ! No of analysis levels
                                             ! for group of obs types
      INTEGER        NO_ANAL_VAR             ! No of variables being
                                             ! analysed (2 for winds)
      INTEGER        LENMG                   ! Length of model grid
      INTEGER        LENAG                   ! Length of analysis grid
      INTEGER        LENOB                   ! No of obs in group
      INTEGER        LENOBT                  ! No of obs for obs type
      INTEGER        ITOTAL                  ! Total no of iterations so
                                             ! far in whole assimilation
      INTEGER        IDIAG                   ! Diagnostic control for
                                             ! this timestep
      INTEGER        TOTAL_NO_ITERS          ! Total no of iterations
                                             ! to be done this timestep
                                             ! (including diagnostic
                                             ! only iterations)
      INTEGER        ITER_NO                 ! Iteration number
                                             ! excluding diagnostic
                                             ! only iterations
      INTEGER        KITER                   ! Iteration no used in
                                             ! diagnostic output.
      INTEGER        IACTF                   ! First obs type in group
      INTEGER        IACTL                   ! Last  obs type in group
      INTEGER        INOBS                   ! No of obs for this type
      INTEGER        INOBSDIM                !  "  " (for dimensioning)
      INTEGER        INC_TYPE                ! Pointer into MODEL_INCR
                                             ! if FI used
      INTEGER        ITNOBS(NOBTYPMX)        ! No of obs in each group
                                             ! assimilated on timestep
      INTEGER        I             !? tracer
      INTEGER        IPT_TRACER    !? tracer ! pointer to one tracer
                                             ! within full array TRACERS
      INTEGER        MODE_HANAL              ! FI or HORINF for this
                                             ! group
      INTEGER        JITER                   ! Loop conter for iteration
      INTEGER        J
      INTEGER        JGROUP,JJGROUP          ! Loop counter over groups
      INTEGER        JACT,JJACT              ! Loop counter over obs
                                             ! types in group
      INTEGER        IGROUP                  ! Group index in group
                                             ! dependent arrays
      INTEGER        HMRMODE                 ! Mode for action in
                                             ! HMRTORH
                                             ! 1 = Convert RH to HMR
                                             ! 2 = Convert HMR to RH
      INTEGER        VMODE                  ! Mode for action in
                                            ! VISTOLV and AEROTOLA
                                            ! 1 = Convert VIS to LOGVIS
                                            ! 2 = Convert LOGVIS to VIS
      INTEGER        ISTAT                  ! status for GCOM
      INTEGER        IMSG                   ! message tag

      REAL           ASSM_TIME               ! Assimilation time
                                             ! Relative to start

      LOGICAL        LWIND                   ! Indicator if group
                                             ! has wind obs types
      LOGICAL        DG_THIS                 ! ) Switches to control
      LOGICAL        DG_BETWEEN              ! ) whether diagnostics
      LOGICAL        DG_END                  ! ) required this timestep
                                             ! ) , between iterations
                                             ! ) and at end of timestep
      LOGICAL        DG_ONLY                 ! Indicator if diagnostic
                                             ! only iteration
      LOGICAL        OLD_LDIAGAC             ! Original value of LDIAGAC
                                             ! stored during timestep
      LOGICAL        L_AC_TIMESTEP           ! Switch to control if
                                             ! iteration done this
                                             ! timestep
      LOGICAL        L_eacf                  ! Use emp. adjusted
                                             ! cloud fraction.

#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cprintst.h"
! Subroutine calls:
      EXTERNAL AC2, GETOBS, HMRTORH, RDOBS

! Extra bits:
      DATA ITOTAL/0/
      SAVE ITOTAL

!- End of header
!*

!L 1.0 Call RDOBS to read in observations either from cache file or
!L     from observation files depending on timestep number.
! DEPENDS ON: rdobs
      CALL RDOBS (NO_OBS_FILES, TIMESTEP_NO, TIMESTEP, OBS, OBS_FLAG,   &
     &               TNDV, TNOBS, P_LEVELS, Q_LEVELS, TNDVMAX, NOBSMAX, &
#include "argppx.h"
     &                                                 ICODE, CMESSAGE)

      IF (MYPE == 0 .AND. PRINTSTATUS >= PRSTATUS_NORMAL) THEN
        WRITE(6,*)'RDOBS done, TIMESTEP_NO = ',TIMESTEP_NO
      END IF

! Check for errors (drop out if one has occured):
      IF (ICODE  >   0) GOTO 999

!L 2.0 Set up Diagnostic Iteration Control from MACDIAG(TIMESTEP_NO)
!     MACDIAG = 0   no diagnostics this timestep
!     MACDIAG = 32  to get diagnostics on each iteration this timestep
!     MACDIAG = +64 for extra diagnostic only iteration between
!                   iterations
!     MACDIAG = + 8 for extra diagnostic only iteration at end
!                   of timestep

      IDIAG = MACDIAG(1 + MOD(TIMESTEP_NO -1, MODEACP))

! Remember LDIAGAC for re-use at end of timestep
      OLD_LDIAGAC = LDIAGAC

! Diagnostics required on each iteration this timestep ?
      DG_THIS    = MOD(IDIAG/32, 2)  ==  1

! Extra diagnostic only iteration required between iterations ?
      DG_BETWEEN = DG_THIS .AND. MOD(IDIAG/64, 2)  ==  1

! Extra diagnostic only iteration required at end of timestep ?
      DG_END     = DG_THIS .AND. MOD(IDIAG/8, 2)  ==  1

! Any diagnostics required this timestep ?
      LDIAGAC = LDIAGAC .AND. (DG_THIS .OR. DG_BETWEEN .OR. DG_END)

!L 3.0 Set up iteration control for ac
!     NO_ITERATIONS : no of iterations for each group.
!     INTERVAL_ITER : interval (in timesteps) between iterations
!                   : for each group.

! Determine total no of iterations.
      TOTAL_NO_ITERS = 0

      DO JGROUP = 1, N_GROUPS
        IGROUP = GROUP_INDEX(JGROUP)

        IF (MOD(TIMESTEP_NO -1, DEF_INTERVAL_ITER(IGROUP))  ==  0) THEN
         ! Do iteration(s) from group JGROUP this timestep.
          TOTAL_NO_ITERS = MAX(TOTAL_NO_ITERS,DEF_NO_ITERATIONS(IGROUP))

        END IF
      END DO

      IF (DG_BETWEEN) TOTAL_NO_ITERS = TOTAL_NO_ITERS*2
      IF (DG_END)     TOTAL_NO_ITERS = TOTAL_NO_ITERS+1

!L 4.0 Start loop iterating ac.
      ITER_NO = 0   !  Iteration No excluding Diagnostic only iters.

! Loop over total no of iterations.
      DO JITER =1, TOTAL_NO_ITERS
       ! Diagnostics only this iteration
        DG_ONLY = (DG_BETWEEN .AND. MOD(JITER,2)  ==  1) .OR.           &
     &            (DG_END     .AND. JITER  ==  TOTAL_NO_ITERS)

        IF (.NOT. DG_ONLY) THEN
          ITOTAL = ITOTAL +1
          ITER_NO = ITER_NO +1

        END IF

       ! KITER is used for printed output.
        KITER = ITER_NO

        IF (DG_ONLY) KITER = 0

        DO J = 1, NOBTYPMX
         ! ITNOBS holds a count of obs used this step of each type
          ITNOBS(J) = 0

        END DO
        NAVT =0

!L      4.1 Loop over groups of ob types
        DO JGROUP = 1, N_GROUPS
        JJGROUP=JGROUP
          IGROUP     = GROUP_INDEX(JGROUP)
          L_AC_TIMESTEP = ITER_NO  <=  DEF_NO_ITERATIONS(IGROUP) .AND.  &
     &                MOD(TIMESTEP_NO-1,DEF_INTERVAL_ITER(IGROUP)) == 0

          IF (DG_ONLY .OR. (.NOT.DG_ONLY .AND. L_AC_TIMESTEP) ) THEN

           ! First obs type in group
            IACTF = GROUP_FIRST(JGROUP)

           ! Last obs type in group
            IACTL = GROUP_LAST (JGROUP)

!L          4.1.1 Define analysis variable indiciator for this
!L                  group (NAVT)
           ! NAVT=4 for humidity (cloud), 5 for precip rate
           ! defined by (1st OB type in group)/100
            NAVT  = LACT(IACTF) / 100
            LWIND = NAVT  ==  3

            IF (LDIAGAC.AND.mype == 0) THEN
              PRINT '(/,A,I3,A,I3,A,I2,A,(20I4),/)',                    &
     &              ' AC STEP =', TIMESTEP_NO, ',', KITER,              &
     &              ' STARTING GROUP NO ', JGROUP,                      &
     &              ' , OBS TYPES ', (LACT(J), J = IACTF, IACTL)

            END IF

           ! Convert Model MIXING RATIO to RELATIVE HUMIDITY
  ! DO for all variables (ensures RH preserved when theta changed)
! (except for humidity with 2A cloud microphysics in use and
!  doing a group of MOPS cloud data
!  -assumes obs type 406 will always be in group on its own)


            IF(NAVT /= 4 ) THEN
              HMRMODE = 1
! DEPENDS ON: hmrtorh
              CALL HMRTORH (HMRMODE, EXNER, PRESSURE, THETA, RH,        &
     &               P_FIELD, P_LEVELS, Q_LEVELS, ICODE, CMESSAGE)
            ENDIF

             ! Check for error - drop out if bad.
              IF (ICODE >  0) GO TO 999

           ! Get Horizontal analysis mode : FI or HORINF?
            MODE_HANAL = DEF_MODE_HANAL(IGROUP)

           ! Get No of Weight/Analysis Levels for this group.
            NO_ANAL_LEVS = DEF_NO_ANAL_LEVS(IGROUP)
            NO_WT_LEVS   = DEF_NO_WT_LEVS  (IGROUP)

           ! Set up work area for derived theta incrs from lhn

            WKLEN = 1

            IF(NAVT == 5 .AND. L_LHN) WKLEN = Q_LEVELS

              NO_ANAL_VAR = 1

           ! Set size of model grid to that used in lower
           ! level routines
              LENMG = P_FIELD

!L          4.1.2  Decide on analysis grid for this group of types.
! parameters for FI method
              LENAG       = 1   ! ANAL_INCR not used in FI method
              INC_TYPE = NO_ANAL_VAR +1


!L          4.1.3 Get a list of obs relevant to this time (getobs)
!L                for each type in the current group of types
           ! ASSM_TIME is time since start in minutes (for GETOBS)
            ASSM_TIME = REAL(TIMESTEP_NO -1) * TIMESTEP / 60.0

            NPTOBT = 0

           ! Loop over observation types in group
            DO JACT = IACTF, IACTL
            JJACT=JACT
             ! LENOBT is no of obs. for type to be used on
             ! this timestep and is determined in GETOBS
              LENOBT = 0

             ! INOBS is total no of observations for type
              INOBS = NOBS(JACT)
              INOBSDIM=MAX(INOBS,1)

               ! Check on available workspace in OBS_NO
                IF (NPTOBT+INOBS  <=  NOBSMAX) THEN
! DEPENDS ON: getobs
                 CALL GETOBS (JJACT, ASSM_TIME, INOBS, OBS_NO(NPTOBT+1),&
     &                   LENOBT, OBS(MDISPOBT(JACT)+1),                 &
     &                   OBS_FLAG(OBS_NO_ST(JACT)+1), TIMESTEP_NO,      &
     &                   DG_ONLY,INOBSDIM,ICODE, CMESSAGE)

      IF (LDIAGAC.AND.mype == 0) THEN
      WRITE(6,*)'AC af GETOBS - LENOBT ',LENOBT
      END IF
                 ! Check for error - drop out if bad
                  IF (ICODE >  0) GO TO 999

                ELSE
                  ICODE = 2
      WRITE(6,*)'TEST FAILED - NPTOBT,INOBS,NOBSMAX ',                  &
     & NPTOBT,INOBS,NOBSMAX

                 ! Drop out of assimilation
                  GOTO 990

                END IF

              LENACT(JACT) = LENOBT
              NPTOBT       = NPTOBT + LENOBT

            END DO   ! End of loop over obs types in group (JACT)

           ! LENOB is now the no of observations in this group
           ! to be used on this timestep
            LENOB = NPTOBT

!
!   HOW MANY OBS (lenob_total) for all pe's?
!
      lenob_com=lenob
      CALL GC_SSYNC(NPROC,ISTAT)
      if(mype == 0)then
        lenob_total=lenob_com
      endif

      Do iproc=1,nproc-1
        IMSG=IPROC
! PE0 sends lenob_pos to each PE in turn
        if(mype == 0) then
          lenob_pos=lenob_total+1
          CALL GC_ISEND(IMSG,1,IPROC,ISTAT,lenob_pos,lenob_pos)
        endif
! PE receives lenob_pos from PE0
        if(mype == IPROC) THEN
          CALL GC_IRECV(IMSG,1,0,ISTAT,lenob_pos,lenob_pos)
        endif
        CALL GC_SSYNC(NPROC,ISTAT)
! Retrieve lenob_com from each PE in turn and adjust lenob_total
        IMSG=IPROC+1000
        IF(mype == 0) THEN
! PE0 receives
          CALL GC_IRECV(IMSG,1,IPROC,ISTAT,lenob_com,lenob_com)
          lenob_total=lenob_total+lenob_com
         ELSEIF(mype == IPROC) THEN
! other PEs send
           CALL GC_ISEND(IMSG,1,0,ISTAT,lenob_com,lenob_com)
        ENDIF
      EndDo
      CALL GC_SSYNC(NPROC,ISTAT)

      if(mype == 0) then
      lenob_pos=1
      lenob_com=lenob
      endif
      CALL GC_SSYNC(NPROC,ISTAT)

! Send lenob_total to all PEs
      IMSG=2000
      CALL GC_IBCAST(IMSG,1,0,NPROC,ISTAT,lenob_total)
      CALL GC_SSYNC(NPROC,ISTAT)

      IF (LDIAGAC.AND.mype == 0) THEN
      WRITE(6,*)'b4 AC2 - lenob_total,TNDV,LENOB ',                     &
     & lenob_total,TNDV,LENOB
      END IF
            if(lenob_total >  0)then !  Skip if no obs for this group
             ! increment ob counter array
              ITNOBS(JGROUP) = ITNOBS(JGROUP) + lenob_total


             ! Lower level routine of ac to begin here now that
             ! all the array dimensions are known - lenob, lenag,
             ! no_anal_levs
! DEPENDS ON: ac2
              CALL AC2(P_LEVELS, Q_LEVELS, BL_LEVELS,                   &
     &            ROW_LENGTH, P_ROWS,                                   &
     &            P_FIELD, TIMESTEP_NO, ITER_NO,                        &
     &            TIMESTEP, OBS, TNDV, EXNER, PSTAR,                    &
     &            THETA, RH, QCL, QCF,                                  &
     &            CONV_CLD, LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,     &
     &               D_THETA_DT_CONV,D_THETA_DT_LS,                     &
     &               LAYER_CLOUD,PRESSURE,                              &
     &            RHCRIT, L_eacf,                                       &
     &            OBS_NO, LENOB, NO_ANAL_LEVS, NO_WT_LEVS, NO_ANAL_VAR, &
     &            LENAG, LENMG, WKLEN, INC_TYPE, NAVT, JJGROUP, LWIND,  &
     &            IACTF, IACTL, DG_ONLY, STINDEX, STLIST, LEN_STLIST,   &
     &            SI, SF, STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,    &
     &            STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                &
     &            ICODE, CMESSAGE)

             ! Check for error - drop out if bad.
              IF (ICODE >  0) GO TO 999

            END IF

           ! Convert RELATIVE HUMIDITY back to MIXING RATIO
          ! (except for MOPS with 2A cld microphys)

            IF(NAVT /= 4 ) THEN
              HMRMODE = 2
! DEPENDS ON: hmrtorh
              CALL HMRTORH (HMRMODE, EXNER, PRESSURE, THETA, RH,        &
     &               P_FIELD, P_LEVELS, Q_LEVELS, ICODE, CMESSAGE)
            ENDIF

             ! Check for error - drop out if bad
              IF (ICODE  >   0) GO TO 999


            IF (LDIAGAC.AND.mype == 0) THEN
              PRINT '(A,I3,A,I3,A,I4)', ' AC STEP',                     &
     &             TIMESTEP_NO, ',', KITER, ' END OF GROUP NO', JGROUP

            END IF
          END IF
        END DO   ! End of loop over groups (JGROUP)
      END DO   ! End of loop over total no of iterations (JITER)

       IF (mype == 0) THEN
        PRINT *, ' '
        PRINT '(A,I3)',  ' End of AC for time step:',TIMESTEP_NO
        PRINT '(A,(10I6))', ' Group No   ',(J,J=1,N_GROUPS)
        PRINT '(A,(10I6))', ' No of obs  ',(ITNOBS(J),J=1,N_GROUPS)
        PRINT *, ' '

      END IF

      LDIAGAC = OLD_LDIAGAC

! Fix to allow vectorization of JACT loop by removing character
! operations.
 990  CONTINUE
      IF (ICODE  ==  2) THEN
        ICODE    = 1
        CMESSAGE = 'AC : Insufficient space in array OBS_NO.'//         &
     &                                               'Increase NOBSMAX'

      END IF

 999  CONTINUE
      RETURN
      END SUBROUTINE AC

!*L  Arguments & declarations:
#endif
