#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routine: ST_MEAN --------------------------------------------------
!LL
!LL  Purpose: Extracts derived diagnostics from climate mean data in the
!LL           D1 array, using the special mean diagnostic sections 21-24
!LL           to define the required diagnostics.
!LL           Designed to be called from within the means subroutine.
!LL           This routine is closely modelled on routines ST_DIAG1 and
!LL           ST_DIAG2, but here the functionality is merged into a
!LL           single routine, although using the same underlying PCRs
!LL           DYN_DIAG and PHY_DIAG.
!LL
!LL  Author:   T C Johns
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 5.1
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: D55
!LL
!LL  Project task: D4
!LL
!LL  External documentation:
!LL    UM Doc Paper C0 - The top-level control system
!LL
!LLEND------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE ST_MEAN (                                              &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "arglndm.h"
#include "argppx.h"
     &             SN,INTMEAN,                                          &
     &             ICODE,CMESSAGE)

      Use level_heights_mod
      Use trignometric_mod, Only : sec_v_latitude, tan_v_latitude,      &
     &                             sec_theta_latitude
      Use dyn_coriolis_mod, Only : f3_at_v
      Use rot_coeff_mod
!
      IMPLICIT NONE

#include "parvars.h"
#include "cmaxsize.h"
#include "csubmodl.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typsts.h"
#include "typptra.h"
#include "typcona.h"
!!! typcona still needed for delta_lambda and delta_phi
#include "typlndm.h"
#include "ppxlook.h"


      INTEGER                                                           &
     &   SN                                                             &
                   ! IN  - Section number for mean diagnostics
     &  ,INTMEAN                                                        &
                   ! IN  - Size of STASHWORK array

     &  ,ICODE     ! OUT - Return code from routine

      CHARACTER*(80) CMESSAGE ! OUT - Return message if failure occurred
!
!*----------------------------------------------------------------------
!  Common blocks
!
#include "chsunits.h"
#include "ccontrol.h"
#include "cphyscon.h"
#include "ctime.h"
#include "c_eta_pmsl.h"
!
!  Dynamically allocated workspace for STASH processing
!
      REAL    STASHWORK(INTMEAN)

! Dummy length for _diag routines
      INTEGER Dummy_len
      PARAMETER (Dummy_len=1)

!  Pressures for DYN_DIAG
      REAL                                                              &
     &        UCOMP_PRESS(NUM_STASH_LEVELS)                             &
     &       ,VCOMP_PRESS(NUM_STASH_LEVELS)                             &
     &       ,CAT_PROB_PRESS(NUM_STASH_LEVELS)                          &
     &       ,PV_THETA(NUM_STASH_LEVELS)                                &
                                          ! requested theta levels
                                          ! for pv.
     &       ,PV_PRESS(NUM_STASH_LEVELS)                                &
                                          ! requested p levels
     &       ,THETA_ON_PV(NUM_STASH_LEVELS)                             &
                                            ! requested pv levels
     &       ,T_PRESS(NUM_STASH_LEVELS)                                 &
     &       ,W_PRESS(NUM_STASH_LEVELS)                                 &
     &       ,Q_PRESS(NUM_STASH_LEVELS)                                 &
     &       ,HEAVY_PRESS(NUM_STASH_LEVELS)                             &
     &       ,Z_PRESS(NUM_STASH_LEVELS)                                 &
     &       ,Dummy_levels(Dummy_len)

!  Dummy array for DYN_DIAG
      REAL                                                              &
     &       PSTAR_OLD(THETA_FIELD_SIZE)

!  Pressures for PHY_DIAG
      REAL                                                              &
     &        T_P_PRESS(NUM_STASH_LEVELS)                               &
     &       ,HTS_PRESS(NUM_STASH_LEVELS)                               &
     &       ,Q_P_ICE_PRESS(NUM_STASH_LEVELS)                           &
     &       ,Q_P_W_PRESS(NUM_STASH_LEVELS)                             &
     &       ,WBPT_PRESS(NUM_STASH_LEVELS)                              &
     &       ,TH_ADV_PRESS(NUM_STASH_LEVELS)                            &
     &       ,DUMMY_TR_PRESS(TR_VARS+1,NUM_STASH_LEVELS)
!  Trig functions
      REAL                                                              &
     &        NMOST_LAT,                                                &
     &        WMOST_LONG,                                               &
     &        EW_SPACE,                                                 &
     &        NS_SPACE,                                                 &
     &        PHI_POLE,                                                 &
     &        LAMBDA_POLE,                                              &
     &        LAT_STEP_INVERSE,                                         &
     &        LONG_STEP_INVERSE
!
!  Local variables
!
      INTEGER                                                           &
     &        I,K,                                                      &
     &        ISL,                                                      &
     &        NI,                                                       &
     &        LEVEL,                                                    &
     &        UCOMP_P_LEVS,                                             &
     &        VCOMP_P_LEVS,                                             &
     &        WCOMP_P_LEVS,                                             &
     &        CAT_PROB_LEVS,                                            &
     &        PV_THETA_LEVS,                                            &
     &        PV_PRESS_LEVS,                                            &
     &        THETA_ON_PV_LEVS,                                         &
     &        DUMMY_TR_LEVS(TR_VARS+1),                                 &
     &        TR_THETA_FIELD_SIZE,                                      &
                                   ! Dummy size for tracer fields
     &        HTS_LEVS,                                                 &
                            !  no. pressure levels for heights
     &        T_P_LEVS,                                                 &
                            !  no. pressure levels for temperature
     &        Q_P_ICE_LEVS,                                             &
                            !  no. pressure levels for humdity wrt ice
     &        Q_P_W_LEVS,                                               &
                            !  no. pressure levels for humdity wrt water
     &        im_ident                                                  &
                            !  Internal Model Identifier
     &       ,im_index      !  Internal Model Index for stash arrays
      INTEGER                                                           &
     &        PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,PT209,    &
     &  PT210,PT211,PT212,PT213,PT214,PT215,PT216,PT217,PT218,PT219,    &
     &  PT220,PT221,PT222,PT223,PT224,PT225,PT226,PT227,PT228,PT229,    &
     &  PT260,PT261,PT262,PT263,PT264,PT265,PT266,                      &
     &  PT256                                                           &
     &  ,PT270,PT271

      LOGICAL ROTATE_UV,ROTATE_MAX_UV                                   &
     &        ,SF_TRACER(TR_VARS+1)
!
! NB: mean P_EXNER has been calculated in MEANCTL1
!
!     Set to atmosphere internal model
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)

! DYN_DIAG and PHY_DIAG do not initialise MPP halos
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INTMEAN
        STASHWORK(I)=0.
      ENDDO

!L----------------------------------------------------------------------
!L 2. Set up levels and pointer information for call to DYN_DIAG
!L    NOTE: Item nos differ from section 15.
!L
      ISL=STINDEX(1,201,SN,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            UCOMP_P_LEVS=STASH_LEVELS(1,NI)
            DO K=1,UCOMP_P_LEVS
              UCOMP_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            ICODE=1
            CMESSAGE='ST_MEAN : STASH_LEVELS not pressure for U_COMP'
            GOTO 999
          ENDIF
        ELSE
          CMESSAGE='ST_MEAN : STASH_LEVELS not a LEVEL list for U_COMP'
          ICODE=1
          GOTO 999
        END IF
      ELSE
        UCOMP_P_LEVS=1
      END IF

      ISL=STINDEX(1,202,SN,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            VCOMP_P_LEVS=STASH_LEVELS(1,NI)
            DO K=1,VCOMP_P_LEVS
              VCOMP_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            ICODE=1
            CMESSAGE='ST_MEAN : STASH_LEVELS not pressure for V_COMP'
            GOTO 999
          ENDIF
        ELSE
          CMESSAGE='ST_MEAN : STASH_LEVELS not a LEVEL list for V_COMP'
          ICODE=1
          GOTO 999
        END IF
      ELSE
        VCOMP_P_LEVS=1
      END IF

      PT201=SI(201,SN,im_index)
      PT202=SI(202,SN,im_index)
      PT205=SI(205,SN,im_index)
      PT206=SI(206,SN,im_index)
      PT207=SI(207,SN,im_index)
      PT208=SI(208,SN,im_index)
      PT209=SI(209,SN,im_index)
      PT260=SI(260,SN,im_index)
      PT261=SI(261,SN,im_index)
      PT262=SI(262,SN,im_index)
      PT263=SI(263,SN,im_index)
      PT264=SI(264,SN,im_index)
      PT265=SI(265,SN,im_index)
      PT266=SI(266,SN,im_index)
      PT270=SI(270,SN,im_index)
      PT271=SI(271,SN,im_index)

!L----------------------------------------------------------------------
!L 3. Call DYN_DIAG for section SN diagnostics
!L
! Set flags to control wind rotation according to whether ELF grid
      IF(A_FIXHD(4) == 3.OR.A_FIXHD(4) == 103) THEN  ! ELF Grid
        ROTATE_UV=.TRUE.
        ROTATE_MAX_UV=.TRUE.
      ELSE
        ROTATE_UV=.FALSE.
        ROTATE_MAX_UV=.FALSE.
      ENDIF
      NMOST_LAT=A_REALHD(3)
      WMOST_LONG=A_REALHD(4)
      NS_SPACE=A_REALHD(2)
      EW_SPACE=A_REALHD(1)
      PHI_POLE=A_REALHD(5)
      LAMBDA_POLE=A_REALHD(6)

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER('DYN_DIAG',3)

!---------------------------------------------------------------------
! Beware! STASH item no.s for meaned Sections are not all identical
!         with item no.s for individual diagnostic sections (15,16) !!
!         section 15  here
!         201         201
!         202         202
!---------------------------------------------------------------------


! DEPENDS ON: dyn_diag
      CALL Dyn_diag(                                                    &
! Primary data: in
     &  D1(jexner_rho_levels(1))                                        &
     & ,D1(JRHO(1)),D1(JU(1)),D1(JV(1)),D1(JW(0))                       &
     & ,D1(jexner_theta_levels(1))                                      &
     &, D1(jtheta(1))                                                   &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,global_rows,global_row_length                                    &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,Model_domain                                                     &
! Grid coordinates: in
     &,delta_lambda,delta_phi                                           &
     &,A_REALHD(rh_deltaEW),A_REALHD(rh_deltaNS)                        &
     &,A_REALHD(rh_baselat),A_REALHD(rh_baselong)                       &
     &,A_REALHD(rh_rotlat),A_REALHD(rh_rotlong)                         &
! Pre-calculated grid associated arrays: in
     &,r_at_u,r_at_v                                                    &
     &, r_theta_levels, r_rho_levels, sec_v_latitude                    &
     &, tan_v_latitude, sec_theta_latitude, f3_at_v                     &
     &,rot_coeff1,rot_coeff2                                            &
! Time information: in
     &,forecast_hrs                                                     &
! Theta levels for output arrays
     &, dummy_levels                                                    &
! Pressure levels for output arrays: in
     &,Dummy_levels,Dummy_levels                                        &
     &,UCOMP_PRESS,VCOMP_PRESS,Dummy_levels                             &
     &,Dummy_levels,Dummy_levels                                        &
     &,Dummy_levels,Dummy_levels,Dummy_levels                           &
     &,Dummy_levels,Dummy_levels,Dummy_levels                           &
! Model levels    for output arrays: in
     &,Dummy_levels,Dummy_levels                                        &
     &,Dummy_levels                                                     &
     &,Dummy_levels,Dummy_levels                                        &
! Flags to request each diagnostic output field: in
! wind related diagnostics
     &,.FALSE.,.FALSE.                                                  &
     &,SF(201,SN),SF(202,SN)                                            &
     &,.FALSE.,.FALSE.,.FALSE.                                          &
     &,.FALSE.,.FALSE.                                                  &
     &,.FALSE.,.FALSE.                                                  &
! PV related diagnostics
     &,.FALSE.,.FALSE.,.FALSE.                                          &
     &,.FALSE.,.FALSE.,.FALSE.                                          &
! test fields
     &,.FALSE.,.FALSE.,.FALSE.,.FALSE.                                  &
! flux diagnostics
     &,SF(260,SN),SF(261,SN),SF(262,SN),SF(263,SN)                      &
     &,SF(264,SN),SF(265,SN),SF(266,SN)                                 &
! height and height level diagnostics
     &,.FALSE.,.FALSE.                                                  &
     &,.FALSE.,.FALSE.,.FALSE.                                          &
     &,.FALSE.,.FALSE.,.FALSE.                                          &
! other diagnostics
     &,SF(270,SN),SF(271,SN)                                            &
! Flags for wind rotation (lam grid): in
     &,.FALSE.                                                          &
! Diagnostics lengths: in
     &,Dummy_len,Dummy_len                                              &
     &,Dummy_len,Dummy_len                                              &
     &,ucomp_p_levs,vcomp_p_levs,Dummy_len                              &
     &,Dummy_len,Dummy_len                                              &
     &,Dummy_len,Dummy_len                                              &
     &,Dummy_len,Dummy_len,Dummy_len,Dummy_len                          &
     &,Dummy_len,Dummy_len,Dummy_len,Dummy_len                          &
! Diagnostic arrays: out
! wind related diagnostics
     &,Stashwork,Stashwork                                              &
     &,Stashwork(pt201),Stashwork(pt202)                                &
     &,Stashwork,Stashwork,Stashwork                                    &
     &,Stashwork,Stashwork                                              &
     &,Stashwork,Stashwork                                              &
! PV related diagnostics
     &,Stashwork,Stashwork,Stashwork                                    &
     &,Stashwork,Stashwork,Stashwork                                    &
! test fields
     &,Stashwork,Stashwork,Stashwork,Stashwork                          &
! flux diagnostics
     &,Stashwork(pt260),Stashwork(pt261),Stashwork(pt262)               &
     &,Stashwork(pt263),Stashwork(pt264),Stashwork(pt265)               &
     &,Stashwork(pt266)                                                 &
! height and height level diagnostics
     &,Stashwork,Stashwork                                              &
     &,Stashwork,Stashwork,Stashwork                                    &
     &,Stashwork,Stashwork,Stashwork                                    &
! other diagnostics
     &,Stashwork(pt270),Stashwork(pt271)                                &
     &)
! DEPENDS ON: timer
                          IF (LTIMER) CALL TIMER('DYN_DIAG',4)
!L----------------------------------------------------------------------
!L 5. Set up levels and pointer information for call to PHY_DIAG
!L    NOTE: Item nos differ from section 16.
!L
      ISL=STINDEX(1,212,SN,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            T_P_LEVS=STASH_LEVELS(1,NI)
            DO K=1,T_P_LEVS
              T_P_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            ICODE=1
            CMESSAGE='ST_MEAN : STASH_LEVELS not pressure for T_P  '
            GOTO 999
          ENDIF
        ELSE
          CMESSAGE='ST_MEAN : STASH_LEVELS not a LEVEL list for T_P  '
          ICODE=1
          GOTO 999
        END IF
      ELSE
        T_P_LEVS=1
      END IF

      ISL=STINDEX(1,211,SN,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            HTS_LEVS=STASH_LEVELS(1,NI)
            DO K=1,HTS_LEVS
              HTS_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            ICODE=1
            CMESSAGE='ST_MEAN : STASH_LEVELS not pressure for HTS  '
            GOTO 999
          ENDIF
        ELSE
          CMESSAGE='ST_MEAN : STASH_LEVELS not a LEVEL list for HTS  '
          ICODE=1
          GOTO 999
        END IF
      ELSE
        HTS_LEVS=1
      END IF

      ISL=STINDEX(1,213,SN,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            Q_P_ICE_LEVS=STASH_LEVELS(1,NI)
            DO K=1,Q_P_ICE_LEVS
              Q_P_ICE_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            ICODE=1
            CMESSAGE='ST_MEAN : STASH_LEVELS not pressure for QPIce'
            GOTO 999
          ENDIF
        ELSE
          CMESSAGE='ST_MEAN : STASH_LEVELS not a LEVEL list for QPIce'
          ICODE=1
          GOTO 999
        END IF
      ELSE
        Q_P_ICE_LEVS=1
      END IF
!L------------ Extract required pressures for R H wrt water ------------
      ISL=STINDEX(1,256,SN,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            Q_P_W_LEVS=STASH_LEVELS(1,NI)
          DO K=1,Q_P_W_LEVS
              Q_P_W_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            ICODE=1
            CMESSAGE='ST_MEAN : STASH_LEVELS not pressure for Q_PWater '
            GOTO 999
          ENDIF
        ELSE
          CMESSAGE='ST_MEAN : STASH_LEVELS not a LEVEL list for Q_PWat'
          ICODE=1
          GOTO 999
        END IF
      ELSE
        Q_P_W_LEVS=1
      END IF

      DO K=1,TR_VARS+1
        SF_TRACER(K)=.FALSE.   ! Set tracer logical indicators to false
      ENDDO

      TR_THETA_FIELD_SIZE=1  ! Set size for tracer arrays to 1
!
      IF(SF(211,SN)) THEN
        SF(210,SN)=.TRUE. !making sure model half heights switched
      ENDIF               !on if heights on pressure surface is reqd

      PT211=SI(211,SN,im_index)
      PT212=SI(212,SN,im_index)
      PT213=SI(213,SN,im_index)
      PT214=SI(214,SN,im_index)
      PT215=SI(215,SN,im_index)
      PT216=SI(216,SN,im_index)
      PT217=SI(217,SN,im_index)
      PT218=SI(218,SN,im_index)
      PT219=SI(219,SN,im_index)
      PT220=SI(220,SN,im_index)
      PT224=SI(224,SN,im_index)
      PT256=SI(256,SN,im_index)

!L----------------------------------------------------------------------
!L 6. Call PHY_DIAG for section SN diagnostics
!L
! DEPENDS ON: timer
                          IF (LTIMER) CALL TIMER('PHY_DIAG',3)

!---------------------------------------------------------------------
! Beware! STASH item no.s for meaned Sections are not identical with
!         item no.s for individual diagnostic sections (15,16) !!
!         section 16  here
!         202         211
!         203         212
!         204         213
!         222         224
!---------------------------------------------------------------------

! DEPENDS ON: phy_diag
      CALL Phy_diag(                                                    &
! Primary data: in
     & D1(JPSTAR),D1(JP(1)),D1(JRHO(1)),D1(JU(1)),D1(JV(1)),D1(JW(0))   &
     &,D1(JTHETA(1)),D1(JQ(1))                                          &
     &,D1(jp_theta_levels(1))                                           &
     &,d1(jexner_rho_levels(1)),d1(jexner_theta_levels(1))              &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,global_row_length, global_rows                                   &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,r_theta_levels,r_rho_levels                                      &
     &,delta_lambda,delta_phi,sec_theta_latitude                        &
! Control information: in
     &,Model_domain,npmsl_height                                        &
! Pressure levels for output arrays: in
     &,HTS_PRESS,T_P_PRESS                                              &
     &,Q_P_ICE_PRESS,Q_P_W_PRESS,Dummy_levels                           &
! Flags to request each diagnostic output field: in
     &,.FALSE.                                                          &
     &,.FALSE.                                                          &
     &,SF(211,SN),SF(212,SN),SF(213,SN),.FALSE.                         &
     &,SF(224,SN)                                                       &
     &,.FALSE.,SF(256,SN)                                               &
! Diagnostics lengths: in
     &,HTS_LEVS,T_P_LEVS,Q_P_ICE_LEVS,Q_P_W_LEVS                        &
     &,Dummy_len                                                        &
! Diagnostic arrays: out
     &,Stashwork                                                        &
     &,Stashwork                                                        &
     &,STASHWORK(PT211),STASHWORK(PT212),STASHWORK(PT213)               &
     &,Stashwork                                                        &
     &,STASHWORK(PT224)                                                 &
     &,Stashwork,STASHWORK(PT256)                                       &
     &)

! DEPENDS ON: timer
                          IF (LTIMER) CALL TIMER('PHY_DIAG',4)
!L----------------------------------------------------------------------
!L 7. Call STASH to perform processing for merged DYN_DIAG/PHY_DIAG
!L    diagnostics for mean section SN
!L
! DEPENDS ON: timer
                          IF (LTIMER) CALL TIMER('STASH   ',3)

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,SN,STASHWORK,                                &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &                                 ICODE,CMESSAGE)

! DEPENDS ON: timer
                          IF (LTIMER) CALL TIMER('STASH   ',4)
 999  CONTINUE
      RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE ST_MEAN
#endif
