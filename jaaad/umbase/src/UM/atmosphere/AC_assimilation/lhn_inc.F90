#if defined(A18_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LHN_INC ------------------------------------------------
!LL
!LL  Purpose : Latent Heat nudging of model heating profiles
!LL
!LL  For use on Cray C90
!LL  For Cray - Global  ; Enable defs GLOBAL
!LL
!LL             <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL  4.1  6/6/96 : Changes to allow "downward" limiting by a factor
!LL              :    ALPHA_LHN, absolute size limiting, and filtering
!LL              :    of increments. Add "num of pts searched"
!LL              :    diagnostics (Chris Jones)
!    4.2 25/11/96: T3E mods Stuart Bell
!LL
!LL  4.1  4/09/96:  Port to CRAY T3E  Deborah Salmond
!    4.3 21/1/97 : Correct External statements. Stuart Bell
!    4.4 4/7/97 : Correct MPP version on LHN_SEARCH.  Deborah Salmond
!    4.5 9/7/98 : Fix call to timer.  Stuart Bell
!LL   5.2   12/12/00  change A to Earth_Radius
!LL                   -amend glsize dimensions   B Macpherson
!     5.3   24/07/01  remove calls to swapbounds_shmem
!                     and set_sides, revise calls to
!                     shmem_get/put        A. Maycock
!    5.3 05/12/01 Remove reference to the shmcomm and iovarsac include
!                 files, use local dynamic arrays in place of those
!                 arrays previously decared in iovarsac.h.   S. Cusack
!    6.0  10/10/03     Replace SHMEM with GCOM for SX6. Clive JOnes
!    6.2  16/12/05     Replace slow communication pattern with calls to
!                      GATHER_FIELD[_ML]/SCATTER_FIELD[_ML]
!                      Remove MPP def. Replace PRINTs with WRITEs.
!                      R. Johanni (NEC) / P.Selwood
!    6.2  23/02/05     Add option to set negative LH to zero.
!                      Mark Dixon.
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL  Documentation :  FR  WP 171
!LL
!LLEND------------------------------------------------------------------
!L*  ARGUMENTS
      SUBROUTINE LHN_INC(CONV_HEAT,LS_HEAT,LS_RAIN,LS_SNOW,CONV_RAIN,   &
     &                   CONV_SNOW,P_FIELD,Q_LEVELS,PR_INC,TIMESTEP,    &
     &                   IROWS,ICOLS,                                   &
     &                   RANGE,                                         &
     &                   THINCS,LENMG,ICODE,CMESSAGE)

      IMPLICIT NONE

!-----AC common blocks
#include "acparm.h"
#include "parvars.h"
#include "mppac.h"
#include "comacp.h"
#include "comacdg.h"
#include "comag.h"
#include "commg.h"
#include "c_a.h"

!-----DECLARE VARIABLES
      INTEGER       P_FIELD,Q_LEVELS,LENMG
      INTEGER       IROWS,ICOLS,RANGE
      REAL          PR_INC(LENMG),TIMESTEP,                             &
     &              LS_SNOW(P_FIELD),LS_RAIN(P_FIELD),                  &
     &              CONV_SNOW(P_FIELD),CONV_RAIN(P_FIELD),              &
     &              CONV_HEAT(P_FIELD,Q_LEVELS),                        &
     &              LS_HEAT(P_FIELD,Q_LEVELS),                          &
     &              THINCS(P_FIELD,Q_LEVELS)
!*
      INTEGER       ICODE
      CHARACTER*256 CMESSAGE

!-INTENT=IN---------------------------------------------------
!     P_FIELD                    - No. of points
!     Q_LEVELS                   - No. of wet levels
!     IROWS                      - No. of rows in model grid
!     ICOLS                      - No. of pts in a row
!     RANGE                      - Search range in grid points
!     LENMG                      - Length of model grid
!     PR_INC(LENMG)              - Rainfall incrs on model grid
!     TIMESTEP                   - Timestep in seconds
!     LS_SNOW(P_FIELD)          \                                      .
!     LS_RAIN(P_FIELD)           \  Large scale and convective
!     CONV_SNOW(P_FIELD)         /    rain and snow rates
!     CONV_RAIN(P_FIELD)        /           (diagnostic)
!     CONV_HEAT(P_FIELD,Q_LEVELS)- L.H. incrs to theta due to conv'n
!     LS_HEAT(P_FIELD,Q_LEVELS)  - L.H. incrs to theta due to dynamics
!-INTENT=INOUT-----------------------------------------------
!
!-INTENT=OUT-------------------------------------------------
!     THINCS                     - Calculated increments to theta
!     ICODE                      - Non-zero for failure
!     CMESSAGE                   - Reason for failure
!*-----------------------------------------------------------
!*L External subroutine calls
      EXTERNAL TIMER , LHN_SEARCH
      EXTERNAL RFCSL

! Local arrays and variables
      INTEGER       SEARCH(4*RANGE*(RANGE+1),2)
                                             ! Search Template
      INTEGER       JPTS , JLEVS             ! Loop counters over model
                                             !       points and levels
      INTEGER       NEAR(P_FIELD)            ! Nearest point found
                                             !    by LHN_SEARCH
      INTEGER       NO_INCRS                 ! Diagnostic - no. of incrs
      INTEGER       NO_SEARCH                !   ""    - no. of searches
      INTEGER       RADIUS(5)                !   "" - breakdown of
                                             !        search results
      INTEGER       M_GRID                   ! Grid type for filter

      REAL          TOT_PR(P_FIELD)          ! Total model precip rate
      REAL          ANAL_PR(P_FIELD)         ! Obs ppn on model grid
      REAL          TOT_LH(P_FIELD,Q_LEVELS) ! Total LH profile
      REAL          LIMIT                    ! Timestep limit of incrs
      REAL          Z                        ! used to set COS_LAT
      REAL          COS_LAT(P_ROWS_MAX)     ! used in filter routine
      REAL          FI_LHN                   ! filter scale in radians
      LOGICAL       L_FIRST_SEARCH           ! True if first search of
                                             !       the timestep

      REAL    pr_inc_global(Max2DFieldSize)
      REAL    tot_pr_global(Max2DFieldSize)
      REAL    anal_pr_global(Max2DFieldSize)
      INTEGER near_global(Max2DFieldSize)
      REAL    tot_lh_global(glsize(1,fld_type_p)*glsize(2,fld_type_p),  &
     &                      Q_LEVELS/nproc+1)
      REAL    thincs_global(glsize(1,fld_type_p)*glsize(2,fld_type_p),  &
     &                      Q_LEVELS/nproc+1)

      integer row_length_g, p_rows_g, p_field_g, k, istat
      INTEGER                                                           &
     &  MAP(Q_LEVELS)                                                   &
                            ! processor number for level
     &, N_LEVS_ON_PROC(0:nproc-1)                                       &
                                  ! number of levels on each processor
     &, max_levs_per_cpu                                                &
     &, lev_on_gatpe(Q_LEVELS)

!
!  set up global dimensions
!
      row_length_g=glsize(1,fld_type_p)
      p_rows_g=glsize(2,fld_type_p)
      p_field_g= row_length_g*p_rows_g

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_INC ',3)
      IF(LDIAGAC.AND.mype == 0) THEN
        WRITE (6,*) "***** Starting LHN_INC"
      END IF
!
!C 1.0   set parameters and variables
!
      IF (L_LHN_SEARCH .AND. RANGE == 0) THEN
        if(mype == 0)then
          WRITE (6,*) "WARNING : LHN_RANGE=0,"
          WRITE (6,*) "   therefore, setting L_LHN_SEARCH to FALSE"
        endif
        L_LHN_SEARCH = .FALSE.
      ENDIF

      IF (L_LHN_LIMIT) LIMIT = LHN_LIMIT * TIMESTEP / 86400.0

      L_FIRST_SEARCH = .TRUE.
      FI_LHN    = FI_SCALE_LHN / Earth_Radius
      NO_INCRS  = 0
      NO_SEARCH = 0
      RADIUS(1) = 0
      RADIUS(2) = 0
      RADIUS(3) = 0
      RADIUS(4) = 0
      RADIUS(5) = 0

!
!C 2     initial calculations
!
!C 2.1   precipitation variables
!
      DO JPTS = 1 , P_FIELD
        TOT_PR(JPTS) = LS_RAIN(JPTS) + LS_SNOW(JPTS)+                   &
     &                          CONV_RAIN(JPTS) + CONV_SNOW(JPTS)
        IF (TOT_PR(JPTS)  <   0.0) TOT_PR(JPTS)=0.0
        ANAL_PR(JPTS) = TOT_PR(JPTS) + PR_INC(JPTS)
        IF (ANAL_PR(JPTS)  <   0.0) THEN
          ANAL_PR(JPTS) = 0.0
          PR_INC(JPTS)   = ANAL_PR(JPTS) - TOT_PR(JPTS)
        ENDIF
      ENDDO     ! JPTS
!
!C 2.2   heating profiles, (FR Working Paper 171, eq2)
!
      DO JLEVS = 1 , Q_LEVELS
        DO JPTS = 1 , P_FIELD
          TOT_LH(JPTS,JLEVS) = ( CONV_HEAT(JPTS,JLEVS) +                &
     &                     LS_HEAT(JPTS,JLEVS) ) * TIMESTEP
      IF ( REMOVE_NEG_LH .AND. (TOT_LH(JPTS,JLEVS)                      &
     &       <   0.0) ) THEN
          TOT_LH(JPTS,JLEVS) = 0.0
      ENDIF
        ENDDO   ! JPTS
      ENDDO     ! JLEVS
!
!CCCC
!
!C 3      Calculate theta increments
!
!C 3.1    Set up array NEAR to decide where to get profile for scaling
!         NEAR=-1 implies no scaling necessary, 0 implies scale here
!         greater than 0 means scale a nearby profile, or not at all
!
!  PR_INC, TOT_PR and ANAL_PR are needed on all PEs,
!  so gather and broadcast them

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(PR_INC,       PR_INC_global,                    &
     &                  ICOLS,        IROWS,                            &
     &                  row_length_g, p_rows_g,                         &
     &                  fld_type_p,   halo_type_no_halo,                &
     &                  0,            gc_all_proc_group,                &
     &                  ICODE,        CMESSAGE)

      CALL GC_RBCAST(1111, p_field_g, 0, NPROC, ISTAT, PR_INC_global)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(TOT_PR,        TOT_PR_global,                   &
     &                  ICOLS,         IROWS,                           &
     &                  row_length_g,  p_rows_g,                        &
     &                  fld_type_p,    halo_type_no_halo,               &
     &                  0,             gc_all_proc_group,               &
     &                  ICODE,         CMESSAGE)

      CALL GC_RBCAST(2222, p_field_g, 0, NPROC, ISTAT, TOT_PR_global)

! DEPENDS ON: gather_field
      CALL GATHER_FIELD(ANAL_PR,       ANAL_PR_global,                  &
     &                  ICOLS,         IROWS,                           &
     &                  row_length_g,  p_rows_g,                        &
     &                  fld_type_p,    halo_type_no_halo,               &
     &                  0,             gc_all_proc_group,               &
     &                  ICODE,         CMESSAGE)

      CALL GC_RBCAST(3333, p_field_g, 0, NPROC, ISTAT, ANAL_PR_global)

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_SRCH',3)

!  Calculate NEAR_global on all PEs

      DO JPTS = 1 , p_field_g
        NEAR_global(JPTS) = -1
        IF ( PR_INC_global(JPTS)  /=  0.0 ) THEN
          NEAR_global(JPTS) = 0
          NO_INCRS   = NO_INCRS + 1
          IF ( TOT_PR_global(JPTS)  <                                   &
     &          (EPSILON_LHN * ANAL_PR_global(JPTS)) ) THEN
            NEAR_global(JPTS) = JPTS
            NO_SEARCH = NO_SEARCH + 1
            IF (L_LHN_SEARCH) THEN
!  search for suitable profile
! DEPENDS ON: lhn_search
              CALL LHN_SEARCH(JPTS, NEAR_global(JPTS),                  &
     &           RANGE, SEARCH, p_rows_g, row_length_g,                 &
     &           L_FIRST_SEARCH, TOT_PR_global,                         &
     &           ANAL_PR_global(JPTS), p_field_g, RADIUS)
            ENDIF
          ENDIF
        ENDIF
      ENDDO    ! JPTS

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_SRCH',4)
!
!C 3.2   calculate the increments, and scale by the relaxation coeff
!

! Calculate the mapping of which processor each level will go to
      DO K=0,NPROC-1
        N_LEVS_ON_PROC(K)=0
      ENDDO

      DO K=1, Q_LEVELS
       ! assumes first PE is PE 0
        MAP(K)                 = MOD((K-1),NPROC)
        N_LEVS_ON_PROC(MAP(K)) = N_LEVS_ON_PROC(MAP(K))+1
        LEV_ON_GATPE(K)        = N_LEVS_ON_PROC(MAP(K))
      ENDDO

! Distribute TOT_LH_GLOBAL over the processors
      max_levs_per_cpu = Q_LEVELS / nproc + 1

! DEPENDS ON: gather_field_ml
      CALL GATHER_FIELD_ML(                                             &
     &                      TOT_LH, TOT_LH_global,                      &
     &                      ICOLS, IROWS, Q_LEVELS,                     &
     &                      row_length_g, p_rows_g, max_levs_per_cpu,   &
     &                      map, lev_on_gatpe,                          &
     &                      fld_type_p, halo_type_no_halo )

      DO JLEVS = 1 , N_LEVS_ON_PROC(mype)

        DO JPTS = 1 , p_field_g
          IF (NEAR_global(JPTS)  <   0) THEN
!  No scaling necessary
            THINCS_global(JPTS,JLEVS) = 0.0
          ELSEIF (NEAR_global(JPTS)  ==  0) THEN
!  Scale the profile at the point itself
!  (FR WP 171, eq 5)
            THINCS_global(JPTS,JLEVS) =                                 &
     &               RELAX_CF_LHN * TOT_LH_global(JPTS,JLEVS) *         &
     &               PR_INC_global(JPTS) / TOT_PR_global(JPTS)
          ELSEIF (NEAR_global(JPTS)  >   0 .AND.                        &
     &                               NEAR_global(JPTS)  /=  JPTS) THEN
!  Scale a nearby profile
!  (FR WP 171, eq 7)
            THINCS_global(JPTS,JLEVS) = RELAX_CF_LHN *                  &
     &         (ANAL_PR_global(JPTS)/TOT_PR_global(NEAR_global(JPTS)) * &
     &         TOT_LH_global(NEAR_global(JPTS),JLEVS) -                 &
     &         TOT_LH_global(JPTS,JLEVS))
          ELSEIF (NEAR_global(JPTS)  ==  JPTS .AND. L_LHN_SCALE) THEN
!  Scale by 1/EPSILON_LHN if no suitable profile available
!  (FR WP 171, eq 8)
            THINCS_global(JPTS,JLEVS) = RELAX_CF_LHN *                  &
     &         (( 1.0/EPSILON_LHN) - 1.0) * TOT_LH_global(JPTS,JLEVS)
          ELSE
!  If none of the above apply, ignore the ob
            THINCS_global(JPTS,JLEVS) = 0.0
          ENDIF
        ENDDO ! JPTS

      ENDDO ! JLEVS

! Copy calculated THINCS values back to individual PEs

! DEPENDS ON: scatter_field_ml
      CALL SCATTER_FIELD_ML(                                            &
     &                      THINCS, THINCS_global,                      &
     &                      ICOLS, IROWS, Q_LEVELS,                     &
     &                      row_length_g, p_rows_g, max_levs_per_cpu,   &
     &                      map, fld_type_p, halo_type_no_halo )

!
!  Impose limit by factor of alpha

!
      IF (L_LHN_FACT) THEN
        DO JLEVS = 1 , Q_LEVELS
          DO JPTS = 1 , P_FIELD
            IF (PR_INC(JPTS)  <   0.0) THEN
              IF ((PR_INC(JPTS)/TOT_PR(JPTS)) <  (ALPHA_LHN-1.0)) THEN
                THINCS(JPTS,JLEVS) = (ALPHA_LHN - 1.0) *                &
     &                               TOT_LH(JPTS,JLEVS) * RELAX_CF_LHN
              ENDIF
            ENDIF
          ENDDO   ! JPTS
        ENDDO     ! JLEVS
      ENDIF
!
!C 3.3   Recursive filtering of increments
!
      IF(L_LHN_FILT) THEN
!       Initialise COS_LAT
        Z = ROW1MGTH
        DO JPTS = 1, p_rows_g
          COS_LAT(JPTS) = SIN(Z)
          Z = Z + DLATMG
        ENDDO  ! JPTS
#if defined(GLOBAL)
        M_GRID = 1
        COS_LAT(1) = 0.125 * COS_LAT(2)
        COS_LAT(IROWS) = COS_LAT(1)
#else
        M_GRID = 0
#endif

! Calculate the mapping of which processor each level will go to
      DO K=0,NPROC-1
        N_LEVS_ON_PROC(K)=0
      ENDDO

      DO K=1, Q_LEVELS
        ! assumes first PE is PE 0
        MAP(K)                 = MOD((K-1),NPROC)
        N_LEVS_ON_PROC(MAP(K)) = N_LEVS_ON_PROC(MAP(K))+1
        LEV_ON_GATPE(K)        = N_LEVS_ON_PROC(MAP(K))
      ENDDO

! Distribute THINCS_GLOBAL over the processors
      max_levs_per_cpu = Q_LEVELS / nproc + 1

! DEPENDS ON: gather_field_ml
      CALL GATHER_FIELD_ML(                                             &
     &                      THINCS, THINCS_global,                      &
     &                      ICOLS, IROWS, Q_LEVELS,                     &
     &                      row_length_g, p_rows_g, max_levs_per_cpu,   &
     &                      map, lev_on_gatpe,                          &
     &                      fld_type_p, halo_type_no_halo )


      DO JLEVS = 1 , N_LEVS_ON_PROC(mype)
! DEPENDS ON: rfcsl
        CALL RFCSL(THINCS_global(1,JLEVS),p_rows_g,row_length_g,M_GRID, &
     &             0.0,COS_LAT,DLATMG,DLONGMG,FI_LHN,NPASS_RF_LHN)
      ENDDO ! JLEVS

! Copy calculated THINCS values back to individual PEs

! DEPENDS ON: scatter_field_ml
      CALL SCATTER_FIELD_ML(                                            &
     &                      THINCS, THINCS_global,                      &
     &                      ICOLS, IROWS, Q_LEVELS,                     &
     &                      row_length_g, p_rows_g, max_levs_per_cpu,   &
     &                      map, fld_type_p, halo_type_no_halo )

      ENDIF   ! L_LHN_FILT

!
!C Impose absolute limit
!
          IF (L_LHN_LIMIT) THEN
            DO JLEVS = 1 , Q_LEVELS
              DO JPTS = 1 , P_FIELD
                IF (THINCS(JPTS,JLEVS)  <   -1.0 * LIMIT)               &
     &                 THINCS(JPTS,JLEVS) = -1.0 * LIMIT
                IF (THINCS(JPTS,JLEVS)  >   LIMIT)                      &
     &                 THINCS(JPTS,JLEVS) = LIMIT
              ENDDO  ! JPTS
            ENDDO    ! JLEVS
          ENDIF

        CALL GC_SSYNC(NPROC,ISTAT)
!
!C  4  Diagnostics
!
      IF (LDIAGAC .AND. LHN_DIAG) THEN
       if(mype == 0)then
        WRITE (6,*) ' '
        WRITE (6,*) "Latent Heat Nudging Scheme, LHN_INC, diagnostics"
        WRITE (6,*) "    Parameters set   : EPSILON_LHN = ",EPSILON_LHN
        WRITE (6,*) "                     : LHN_RANGE   = ",RANGE
        IF (L_LHN_FACT) THEN
          WRITE (6,*) "    Limit increments to scale down rain rate",   &
     &                " by at most factor of 1/ ALPHA"
          WRITE (6,*) "       ALPHA=",ALPHA_LHN
        ENDIF
        IF (L_LHN_FILT) THEN
          WRITE (6,*) "    Filtering of increments performed"
          WRITE (6,*) "       filter scale = ",FI_SCALE_LHN             &
     &                                                  / 1000.0," Km"
        ENDIF
        IF (L_LHN_LIMIT) THEN
          WRITE (6,*) "    Limiting of increments set to ",LIMIT,       &
     &                " degrees per timestep = ",LHN_LIMIT," degrees",  &
     &                                                     " per day"
        ENDIF
        WRITE (6,*) ' '
        WRITE (6,*) "Number of increments required was ",NO_INCRS
        IF (L_LHN_SEARCH) THEN
          RADIUS(1) = RADIUS(2) + RADIUS(3) + RADIUS(4)
          WRITE (6,*) "LHN_SEARCH was called for ",NO_SEARCH," of them"
          WRITE (6,*) "The search was successful ",RADIUS(1)," times"
          WRITE (6,*) "It found :"
          WRITE (6,*) "     ",RADIUS(2)," pts at 1 point away"
          WRITE (6,*) "     ",RADIUS(3)," pts at 2 points away"
          WRITE (6,*) "     ",RADIUS(4)," pts at 3 or more points away"
          WRITE (6,*) "     ",RADIUS(5)," total pts searched"
        ELSE
          WRITE (6,*) NO_SEARCH," points required the search, but"
          WRITE (6,*) " the search algorithm was disabled"
        ENDIF
        WRITE (6,*) "Scaling at points failing the search was set to: ",&
     &                        L_LHN_SCALE
        WRITE (6,*) ' '
       endif
       CALL GC_SSYNC(NPROC,ISTAT)
      ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('LHN_INC ',4)
      RETURN
      END SUBROUTINE LHN_INC
#endif
