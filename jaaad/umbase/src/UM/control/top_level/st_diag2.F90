#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine ST_DIAG2---------------------------------------------
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1   9/02/93  : added comdeck CHSUNITS to define NUNIST for
!LL                   comdeck CCONTROL
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.2  13/04/93  Dynamic allocation of main arrays. PFLD removed
!LL                 - not used. R T H Barnes.
!LL   3.3  15/12/93  Remove hardwired LEVNO_ABOVE_BOUNDARY=5 and
!LL        change name to LEVNO_PMSL_CALC, determined as first model
!LL        level above eta=0.795       C.Wilson
!LL   3.4  21/9/94  Calculates model level geopotential heights in
!LL                 metres             S.A.Woltering
!LL  3.4  29/11/94 Add P_FIELDDA,P_LEVELSDA for portable dyn.allocn.
!LL  4.1  15/05/96 Add code to process tracer data.
!LL                Including TR_VARSDA         D.Podd
!LL   4.3  10/02/97  Added PPX arguments to COPY_DIAG   P.Burton
!LL   4.4  25/09/97  Fix problems when tracers included and no tracer
!LL                  diagnostics required.  D. Podd.
!LL   4.4  09/09/97  Remove calculation of LEVNO_PMSL_CALC. D. Robinson
!LL   4.5 20/04/98   Initialise STASHWORK so that PHY_DIAG does not
!LL                  need to initialise halos S.D.Mullerworth
!LL   4.5  05/06/98  New arguments L_VINT_TP, L_LSPICE for PHY_DIAG.
!ll                  D. Robinson
!LL   5.0  11/05/99  Change STASH argument list to generic form.
!LL                  R. Rawlins
!LL   5.0  24/05/99  Changed variable names; remove *DA arguments
!LL                  P.Burton.
!LL                  Change interface to Phy_diag. R Rawlins.
!LL   5.1  25/02/00  Remove redundant initialisation of trig variables.
!LL                  R Rawlins
!LL   5.1  17/02/00  Add diagnostic for T on model levels. [And replace
!LL                  non-activated, untested height diagnostic
!LL                  (16,225).] R Rawlins
!LL   5.2  19/01/01  Add diagnostic for wet bulb potential temperature
!LL                  on pressure levels. R Rawlins
!LL   5.3  05/12/01  Add diagnostics for geopotential height on
!                    model theta and rho levels. D.M. Goddard
!     5.4  12/04/02  Add diagnostics for relative humidity wrt
!                    water 16,256. D.M. Goddard
!     5.4  21/03/02  Remove comment on same line as #include
!                                                 S. Carroll
!     6.0  26/09/03  Prevents overwriting of other diagnostics
!                    when geopotential height requested not
!                    for all levels. D.M. Goddard, P. Selwood
!LL
!LL Purpose : To provide the interface for PHY_DIAG
!LL
!LL Control routine for CRAY YMP
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL                       version no. 1, dated 15/01/90
!LL
!LL Logical components covered : D4
!LL
!LL System task : P0
!LL
!LL Documentation : Unified Model Documentation Paper No P0
!LL
!LLEND---------------------------------------------------------------
!*L Arguments

      SUBROUTINE ST_DIAG2( INT16,                                       &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &                    ICODE,CMESSAGE)

Use level_heights_mod, Only : eta_theta_levels, eta_rho_levels,         &
                               r_theta_levels, r_rho_levels
Use trignometric_mod,  Only : sec_theta_latitude

      IMPLICIT NONE
!*L
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typsts.h"
#include "ctracera.h"
#include "typ_atm_fields.h"
#include "typcona.h"
#include "typlndm.h"
#include "ppxlook.h"

      INTEGER                                                           &
     &        INT16                                                     &
                                ! STASHWORK size = STASH_MAXLEN(16)
     &       ,ICODE             ! Out return code : 0 Normal exit
!                               !                 :>0 Error exit

      CHARACTER*(80)                                                    &
     &        CMESSAGE          ! Out error message if ICODE > 0

#include "chsunits.h"
#include "ccontrol.h"
#include "c_r_cp.h"
! for radians
#include "c_pi.h"
#include "c_eta_pmsl.h"

!L External subroutines called

      EXTERNAL                                                          &
     &        STASH,                                                    &
     &        PHY_DIAG,                                                 &
     &        TIMER                                                     &
     &       ,copydiag_3D

!L Locally dynamically allocated work area

      REAL                                                              &
     &        STASHWORK(INT16),                                         &
     &        T_P_PRESS(NUM_STASH_LEVELS),                              &
     &        HTS_PRESS(NUM_STASH_LEVELS),                              &
     &        RHice_press(NUM_STASH_LEVELS),                            &
     &        RHwat_press(NUM_STASH_LEVELS),                            &
     &        WBPT_PRESS(NUM_STASH_LEVELS),                             &
     &        TH_ADV_PRESS(NUM_STASH_LEVELS),                           &
     &        PRESS_LEVS(NUM_STASH_LEVELS),                             &
     &        TR_PRESS(TR_VARS+1,NUM_STASH_LEVELS),                     &
     &        T_m(THETA_FIELD_SIZE,MODEL_LEVELS)                        &
                                                 ! T on model levels
     &       ,hts_theta(theta_field_size,model_levels)                  &
                                                 ! hts on theta levels
     &       ,hts_rho(theta_field_size,model_levels)
                                                 ! hts on rho levels

! Local variables

      INTEGER                                                           &
     &        I,J,                                                      &
     &        NI,                                                       &
     &        K,                                                        &
     &        ISL,                                                      &
     &        BL,                                                       &
     &        TL,                                                       &
     &        LEVEL,                                                    &
     &        LAST_POINT,                                               &
     &        FIRST_POINT                                               &
     &       ,ITR                                                       &
     &       ,STASH_TR_FIRST                                            &
     &       ,STASH_TR_LAST                                             &
     &       ,im_ident                                                  &
                            !  Internal Model Identifier
     &       ,im_index                                                  &
                            !  Internal Model Index for Stash arrays
     &      ,item                                                       &
                            ! STASH item
     &      ,sect           ! STASH section
      PARAMETER( sect = 16 ) ! for 'physics' end of timestep diagnostics

      INTEGER                                                           &
     &        T_P_LEVS,                                                 &
     &        HTS_LEVS, H2_P_LEVS,                                      &
     &        RHice_levs,RHwat_levs,                                    &
     &        WBPT_LEVS,                                                &
     &        TH_ADV_P_LEVS                                             &
     &       ,TR_PRESS_LEVS(TR_VARS+1)                                  &
     &       ,TR_THETA_FIELD_SIZE   ! size for DA of tracer fields

      INTEGER                                                           &
     &        H2_IND(NUM_STASH_LEVELS)

      INTEGER                                                           &
     &        PT_TRACER(TR_VARS+1)
      INTEGER                                                           &
     &        PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,PT209,    &
     &  PT210,PT211,PT212,PT213,PT214,PT215,PT216,PT217,PT218,PT219,    &
     &  PT220,PT221,PT222,PT223,PT224,PT225,PT255,PT256

      REAL                                                              &
     &        EW_SPACE,                                                 &
     &        NS_SPACE

      LOGICAL                                                           &
     &        SF_TRACER(TR_VARS+1)

!L Initialisation
        FIRST_POINT = 1
        LAST_POINT  = THETA_FIELD_SIZE

!L Internal Structure:

!     Set to atmosphere internal model
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)


!L----- Calculate additional diagnostic quantities------------------

!L------------ Extract required pressures for T_P ----------------------
      NS_SPACE=A_REALHD(2)
      EW_SPACE=A_REALHD(1)

      ISL=STINDEX(1,203,16,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            T_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,T_P_LEVS
              T_P_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
        T_P_LEVS=1
      END IF

!L------------ Extract required pressures for Heights ------------------

      ISL=STINDEX(1,202,16,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            HTS_LEVS=STASH_LEVELS(1,NI)
            DO K =1,HTS_LEVS
              HTS_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
        HTS_LEVS=1
      END IF


!L------------ Extract required pressures for R H wrt ice --------------

      ISL=STINDEX(1,204,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            RHice_levs=STASH_LEVELS(1,NI)
            DO K =1,RHIce_levs
              RHice_press(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
             CMESSAGE='ST_DIAG2 : Level not pressure for RHice'
            ICODE=1
            RETURN
          END IF
        ELSE
           CMESSAGE='ST_DIAG2 : Level not a levels list for RHice'
          ICODE=1
          RETURN
        END IF
      ELSE
        RHice_Levs=1
      END IF

!L------------ Extract required pressures for R H wrt water ------------
      ISL=STINDEX(1,256,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            RHwat_levs=STASH_LEVELS(1,NI)
            DO K =1,RHwat_levs
              RHwat_press(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
            CMESSAGE='ST_DIAG2 : Level not pressure for RHwat'
            ICODE=1
            RETURN
          END IF
        ELSE
          CMESSAGE='ST_DIAG2 : Level not a levels list for RHwat'
          ICODE=1
          RETURN
        END IF
      ELSE
        RHwat_levs=1
      END IF
!L------------ Extract required pressures for Wet bulb pot temp---------

      ISL=STINDEX(1,205,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            WBPT_LEVS=STASH_LEVELS(1,NI)
            DO K =1,WBPT_LEVS
              WBPT_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
             CMESSAGE='ST_DIAG2 : Level not pressure for WBPT'
            ICODE=1
            RETURN
          END IF
        ELSE
           CMESSAGE='ST_DIAG2 : Level not a levels list for WBPT'
          ICODE=1
          RETURN
        END IF
      ELSE
        WBPT_LEVS=1
      END IF
!L------------ Extract required pressures for Thermal advection --------

      ISL=STINDEX(1,219,16,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI=-STLIST(10,ISL)
            TH_ADV_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,TH_ADV_P_LEVS
              TH_ADV_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
          ELSE
             CMESSAGE='ST_DIAG2 : Level not pressure for THADV'
            ICODE=1
            RETURN
          END IF
        ELSE
           CMESSAGE='ST_DIAG2 : Level not a levels list for THADV'
          ICODE=1
          RETURN
        END IF
      ELSE
        TH_ADV_P_LEVS=1
      END IF

!L------------ Extract required pressures for Tracers ------------------

      STASH_TR_FIRST=226
      STASH_TR_LAST=254
!     ITR is a count of tracers found to be using this diagnostic
      ITR=0
      IF (TR_VARS >  0) THEN

!     Initialize PT_TRACER and SF_TRACER arrays
        DO I=1,TR_VARS
          PT_TRACER(I)=1
          SF_TRACER(I)=.FALSE.
        END DO

      DO J=STASH_TR_FIRST,STASH_TR_LAST
        ISL=STINDEX(1,J,16,im_index)
        IF(ISL >  0) THEN
          IF(STLIST(10,ISL) <  0) THEN
            IF(STLIST(11,ISL) == 2) THEN
                ITR=ITR+1
              NI=-STLIST(10,ISL)
                TR_PRESS_LEVS(ITR)=STASH_LEVELS(1,NI)
                DO K=1,TR_PRESS_LEVS(ITR)
                  TR_PRESS(ITR,K)=STASH_LEVELS(K+1,NI)/1000.0
              ENDDO
                PT_TRACER(ITR)=SI(J,16,im_index)
                SF_TRACER(ITR)=SF(J,16)
            ELSE
               CMESSAGE='ST_DIAG2 : Level not pressure for Tracers'
              ICODE=1
              RETURN
            END IF
          ELSE
             CMESSAGE='ST_DIAG2 : Level not a levels list for Tracers'
            ICODE=1
            RETURN
          END IF
        END IF
      END DO
      END IF

!     Set last (or only) values in tracer pointer arrays
      PT_TRACER(TR_VARS+1)=1
      SF_TRACER(TR_VARS+1)=.FALSE.

!     Set size of TR_P_FIELD_DA depending on whether any tracers
!     are using this diagnostic. (Used for dynamic allocation of
!     tracer arrays in phy_diag).
      IF (ITR >  0) THEN
        TR_THETA_FIELD_SIZE=THETA_FIELD_SIZE
      ELSE
        TR_THETA_FIELD_SIZE=1

      END IF

!L----- check height available for calculation of height**2 -----------


      IF (SF(224,16)) THEN
        IF (.NOT.SF(202,16)) THEN
          CMESSAGE='ST_DIAG2 : ERROR h**2 requires H at same timestep'
          ICODE=1
          GOTO 999
         ELSE
      ISL=STINDEX(1,224,16,im_index)
          IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            H2_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,H2_P_LEVS
              PRESS_LEVS(K)=STASH_LEVELS(K+1,NI)/1000.0
              DO I=1,HTS_LEVS
                IF(PRESS_LEVS(K) == HTS_PRESS(I)) THEN
                  H2_IND(k)=I
                ENDIF
              ENDDO
            ENDDO
          ELSE
            H2_P_LEVS=1
          END IF
         ENDIF
      ELSE
        H2_P_LEVS=1
      END IF


      PT201=SI(201,16,im_index)
      PT202=SI(202,16,im_index)
      PT203=SI(203,16,im_index)
      PT204=SI(204,16,im_index)
      PT205=SI(205,16,im_index)
      PT206=SI(206,16,im_index)
      PT207=SI(207,16,im_index)
      PT208=SI(208,16,im_index)
      PT209=SI(209,16,im_index)
      PT210=SI(210,16,im_index)
      PT211=SI(211,16,im_index)
      PT212=SI(212,16,im_index)
      PT213=SI(213,16,im_index)
      PT214=SI(214,16,im_index)
      PT215=SI(215,16,im_index)
      PT216=SI(216,16,im_index)
      PT217=SI(217,16,im_index)
      PT218=SI(218,16,im_index)
      PT219=SI(219,16,im_index)
      PT220=SI(220,16,im_index)
      PT221=SI(221,16,im_index)
      PT222=SI(222,16,im_index)
      PT223=SI(223,16,im_index)
      PT224=SI(224,16,im_index)
      PT225=SI(225,16,im_index)
      PT255=SI(255,16,im_index)
      PT256=SI(256,16,im_index)

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('PHY_DIAG',3)
      END IF

! Initialise STASHWORK as PHY_DIAG avoids MPP halo calculations
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INT16
        STASHWORK(I)=0.
      ENDDO

! DEPENDS ON: phy_diag
      CALL Phy_diag(                                                    &
! Primary data: in
     & PSTAR,P,RHO,U,V,W                                                &
     &,THETA,Q                                                          &
     &,p_theta_levels                                                   &
     &,exner_rho_levels,exner_theta_levels                              &
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
     &,RHice_press,RHwat_press,WBPT_PRESS                               &
! Flags to request each diagnostic output field: in
     &,SF(  4,16),SF(201,16),SF(202,16),SF(203,16),SF(204,16)           &
     &,SF(205,16),SF(222,16),SF(255,16),SF(256,16)                      &
! Diagnostics lengths: in
     &,HTS_LEVS,T_P_LEVS,RHice_levs,RHwat_levs,WBPT_LEVS                &
! Diagnostic arrays: out
     &,T_m                                                              &
     &,hts_theta,STASHWORK(PT202),STASHWORK(PT203)                      &
     &,STASHWORK(PT204),STASHWORK(PT205),STASHWORK(PT222)               &
     &,hts_rho,STASHWORK(PT256)                                         &
     &)



      item = 4      !  (4,16) is T on model levels
      IF( SF(item,sect) ) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),T_m,         &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

      item = 201    !  (201,16) is geopotential height on theta levels

      IF( SF(item,sect) ) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),hts_theta,   &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)

      item = 255    !  (255,16) is geopotential height on rho levels

      IF( SF(item,sect) ) THEN
! DEPENDS ON: copydiag_3d
        CALL copydiag_3D(STASHWORK(SI(item,sect,im_index)),hts_rho,     &
     &        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        im_ident,sect,item,                                       &
     &        icode,cmessage)

      ENDIF   ! SF(  item,sect)


      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('PHY_DIAG',4)
      END IF

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',3)
      END IF

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,16,STASHWORK,                                &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &           ICODE,CMESSAGE)

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',4)
      END IF

  999 CONTINUE
      RETURN
      END SUBROUTINE ST_DIAG2
#endif
