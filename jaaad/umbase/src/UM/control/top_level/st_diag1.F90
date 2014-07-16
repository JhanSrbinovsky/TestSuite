#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine ST_DIAG1 ----------------------------------------------
!LL
!LL Purpose: Calculates STASH output diagnostics from 'dynamical'
!LL    fields, ie the wind fields.
!LL    Called at timestep 0 and at the end of timesteps.
!LL
!LL Control routine for CRAY YMP
!LL
!LL TJ, RS      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.1    9/02/93 : added comdeck CHSUNITS to define NUNITS for
!LL                    comdeck CCONTROL.
!LL 3.1   25/01/93  Change arguments to DYN_DIAG to include extra test
!LL diagnostics, items 231,232,233,234. R. Rawlins
!LL 3.1      14/01/93 Include code to output potential vorticity on
!LL                   pressure surfaces and theta on a PV surface.
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL  3.2  13/04/93  Dynamic allocation of main arrays. R T H Barnes.
!LL 3.3       01/11/93 Correct calculations of LAT_STEP_INVERSE
!LL                    and LONG_STEP_INVERSE.   A.S.Lawless
!LL   3.4   26/05/94  LOGICAL LLINTS passed to DYN_DIAG
!LL                                               S.J.Swarbrick
!LL   3.5   10/04/95  Sub-model changes : Timestep length removed
!LL                   from Atmos dump header. D Robinson.
!LL  4.4  10/04/97 : Add new diagnostics wq, Heavyside function and
!LL                  total column KE. Nos 235, 236, 237  R.A.Stratton
!LL       30/07/97 : Also Z, uZ, vZ, nos 238, 239, 240 where Z is
!LL                  geopotential height on u grid. R A Stratton.
!LL       19/08/97 : 241 mountain torque added. R A Stratton.
!LL   4.4   03/10/97  Pass LEVNO_PMSL_CALC to DYN_DIAG. D. Robinson
!LL   4.5 20/04/98   Initialise STASHWORK so that DYN_DIAG does not
!LL                  need to initialise halos S.D.Mullerworth
!LL   5.0 11/05/99   Change STASH argument list to generic form.
!LL                  R. Rawlins
!LL   5.0 24/05/99   Changed variable names and interface to dyn_diag,
!LL                  remove *DA references, for
!LL                  'C-P C dynamics upgrade' project. Rick Rawlins.
!LL   5.1 25/01/00   (1) Change u,v diagnostics from 'C' to 'B' grid.
!LL                  Retain old diagnostics but with new STASHcodes.
!LL                  (2) Remove 'product' diagnostic references (now
!LL                  calculated in section 30).
!LL                  (3) Remove references to diagnostics not yet
!LL                  re-enabled from vn4.5.
!LL                  (4) Introduce u,v on model levels interpolated
!LL                  from 'C' to 'B' grid.
!LL                  R Rawlins
!LL   5.2   19/03/01 Changed argument list for Dyn_diag
!LL   5.2 19/02/01   Introduce rotation of winds for a subset of
!LL                  diagnostics in lam. R Rawlins
!LL   5.2 15/11/00   Add in stash variables for pv on theta levels
!LL                  diagnostic.   Zoe Gardner
!     5.3 06/09/01   Set up stash variables levels required for pv on
!                    pressure levels diagnostic.  D.M. Goddard
!  5.4  16/07/02   Introduce new diagnostics: u,v,w,theta,rho,p on
!                  geometric height levels and height (in meters) of
!                  theta, rho model levels from sea level.
!                           M. Diamantakis
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!     5.4 14/05/02   Add a theta on pv=2 surface diagnostic and theta
!                    at pv points and pv on model levels T.J. Hinton
!  5.5  28/02/03   Add new mass flux diagnostics       Carol Roadnight
!  5.5  26/02/03   Get averaged spectra of vertical velocity field
!                                                       Carol Roadnight
!  6.0  18/07/03   Merging 15215 and 15218 T.J.Hinton
!  6.1  17/05/04   delta_lamda and delta_phi added to Dyn_diag
!  6.1             arguments. Adam Clayton
!  6.4  02/08/06   Get out true unscaled density   Andy Malcolm
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL                       version no. 1, dated 15/01/90
!LL
!LL Logical components covered:
!LL
!LL System task : P0
!LL
!LL Documentation : Unified Model Documentation Paper No P0
!LL                version number 11 dated 26/11/90
!LL             and Unified Model documentation paper No C4
!LL                version number 11 dated 23/11/90
!LL
!LLEND---------------------------------------------------------------
!*L Arguments

      SUBROUTINE ST_DIAG1( INT15,                                       &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &                    ICODE,CMESSAGE)

      Use level_heights_mod
      Use trignometric_mod, Only : sec_v_latitude, tan_v_latitude,      &
     &                             sec_theta_latitude
      Use dyn_coriolis_mod, Only : f3_at_v
      Use rot_coeff_mod
!*

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
     &        INT15,                                                    &
                                ! Dummy for STASH_MAXLEN(15)
     &        ICODE              ! Out return code : 0 Normal exit
!                               !                 : >0 Error exit

      CHARACTER*(80)                                                    &
     &        CMESSAGE          ! Out error message if ICODE > 0

      CHARACTER*(*) RoutineName
      PARAMETER ( RoutineName = 'ST_DIAG1')

#include "chsunits.h"
#include "ccontrol.h"
#include "cphyscon.h"
#include "ctime.h"

!L External subroutines called

      EXTERNAL                                                          &
     &        STASH,                                                    &
     &        TIMER,                                                    &
     &        DYN_DIAG                                                  &
     &       ,Ereport

! Integer function
      INTEGER exppxi

!L Dynamically allocated workspace for stash processing
      INTEGER                                                           &
     &        ucomB_model(NUM_STASH_LEVELS)                             & 
                                                 ! 'B' grid model levs
     &       ,vcomB_model(NUM_STASH_LEVELS)                             & 
                                                 ! 'B' grid model levs
     &       ,Htheta_model(num_stash_levels)                            &
     &       ,Hrho_model(num_stash_levels)


      REAL                                                              &
     &        STASHWORK(INT15)                                          &
     &       ,ucomB_press(NUM_STASH_LEVELS)                             & 
                                                 ! 'B' grid press levs
     &       ,vcomB_press(NUM_STASH_LEVELS)                             & 
                                                 ! 'B' grid press levs
     &       ,UCOMP_PRESS(NUM_STASH_LEVELS)                             &
     &       ,VCOMP_PRESS(NUM_STASH_LEVELS)                             &
     &       ,WCOMP_PRESS(NUM_STASH_LEVELS)                             &
     &       ,TESTD_PRESS(NUM_STASH_LEVELS)                             &
     &       ,TESTD_MODEL(NUM_STASH_LEVELS)                             &
     &       ,PV_THETA(NUM_STASH_LEVELS)                                &
     &       ,p_height(num_stash_levels)                                &
     &       ,theta_height(num_stash_levels)                            &
     &       ,rho_height(num_stash_levels)                              &
     &       ,w_height(num_stash_levels)                                &
     &       ,u_height(num_stash_levels)                                &
     &       ,v_height(num_stash_levels)                                &
     &       ,PV_PRESS(NUM_STASH_LEVELS)

! Local variables

      INTEGER                                                           &
     &        I,                                                        &
     &        NI,                                                       &
     &        K,                                                        &
     &        ISL,                                                      &
     &        BL,                                                       &
     &        TL,                                                       &
     &        LEVEL,                                                    &
     &        irotu,irotv,                                              &
     &        ucomB_m_levs,                                             &
     &        vcomB_m_levs,                                             &
     &        ucomB_p_levs,                                             &
     &        vcomB_p_levs,                                             &
     &        UCOMP_P_LEVS,                                             &
     &        Htheta_m_levs,Hrho_m_levs,                                &
     &        p_h_levs,theta_h_levs,rho_h_levs,                         &
     &        w_h_levs,u_h_levs,v_h_levs,                               &
     &        VCOMP_P_LEVS,                                             &
     &        WCOMP_P_LEVS,                                             &
     &        pv_theta_levs, pv_press_levs,                             &
     &        TESTD_P_LEVS,TESTD_M_LEVS                                 &
     &       ,im_ident                                                  &
                            !  Internal Model Identifier
     &       ,im_index                                                  &
                            !  Internal Model Index for Stash Arrays
     &       ,item                                                      &
                            !  STASH item no.
     &       ,sect          !  STASH section

      Parameter( sect = 15 ) ! for 'dynamical' related diagnostics

      INTEGER                                                           &
     &        PT002,PT003,                                              &
     &        PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,PT209,    &
     &        PT210,PT211,PT212,PT213,PT214,PT215,PT216,PT217,PT218,    &
     &        PT219,PT220,PT221,PT222,PT223,PT224,PT225,PT226,PT227,    &
     &        PT228,PT229,PT230,                                        &
     &        PT231,PT232,PT233,PT234,PT235,PT236,PT237,                &
     &        PT238,PT239,PT240,PT241                                   &
     &       ,PT242,PT243,PT244,PT245,PT246                             &
     &       ,PT260,PT261,PT262,PT263,PT264,PT265,PT266                 &
     &       ,pt101,pt102,pt108,pt119                                   &
     &       ,pt127,pt142,pt143,pt144                                   &
     &       ,PT270,PT271

      LOGICAL                                                           &
! Flags for wind rotation (lam grid)
     & rot_uvcomB_press


!L Internal Structure:

!     Set to atmosphere internal model
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)
      ICODE = 0
      CMESSAGE=''


!L  Section 15  Dynamics diagnostics
!L
!L  Local workspace definitions
!L ---------------------------------------------------------------------
!L      call DYN_DIAG to calculate dynamical diagnostics and
!L      call STASH to process output
!L ---------------------------------------------------------------------
!L
!L This section of code contains numbers used to
!L check on the type of level which determines
!L the interpolation,
!L the codes are as follows (all for stashlist entry 11)
!L 1 -- model levels
!L 2 -- Pressure Levels
!L 3 -- Height Levels
!L 4 -- Theta Levels
!L 5 -- Potential Vorticity Levels

!L-------------------Extract Reqd Model levs for u on B grid (15,002)--

      item = 2    ! u component on B grid, model levels
      ISL=STINDEX(1,item,sect,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            ucomB_m_levs=STASH_LEVELS(1,NI)
            DO K =1,ucomB_m_levs
              ucomB_model(K)=STASH_LEVELS(K+1,NI)  ! integer levels
            ENDDO
      ELSE
            ucomB_m_levs=0
      END IF

!L-------------------Extract Reqd Model levs for v on B grid (15,003)--

      item = 3    ! v component on B grid, model levels
      ISL=STINDEX(1,item,sect,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            vcomB_m_levs=STASH_LEVELS(1,NI)
            DO K =1,ucomB_m_levs
              vcomB_model(K)=STASH_LEVELS(K+1,NI)  ! integer levels
            ENDDO
      ELSE
            vcomB_m_levs=0
      END IF

!L-------------------Extract Reqd Pressures for u_comB_p-------------

      ISL=STINDEX(1,201,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            ucomB_p_levs=STASH_LEVELS(1,NI)
            DO K =1,ucomB_p_levs
              ucomB_press(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
         ucomB_p_levs=1
      END IF

!L-------------------Extract Reqd Pressures for v_comB_p-------------

      ISL=STINDEX(1,202,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            vcomB_p_levs=STASH_LEVELS(1,NI)
            DO K =1,vcomB_p_levs
              vcomB_press(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
         vcomB_p_levs=1
      END IF
!L-------------------Extract Reqd Pressures for U_COMP_P-------------

      ISL=STINDEX(1,243,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            UCOMP_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,UCOMP_P_LEVS
              UCOMP_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
         UCOMP_P_LEVS=1
      END IF

!L-------------------Extract Reqd Pressures for V_COMP_P-------------

      ISL=STINDEX(1,244,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            VCOMP_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,VCOMP_P_LEVS
              VCOMP_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
         VCOMP_P_LEVS=1
      END IF

!L-------------------Extract Reqd Pressures for W_COMP_P-------------

      ISL=STINDEX(1,242,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            WCOMP_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,WCOMP_P_LEVS
              WCOMP_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
         WCOMP_P_LEVS=1
      END IF
!L----------Extract required thetas for Potn_vort on theta ----

      ISL=STINDEX(1,214,15,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 4) THEN
            NI = -STLIST(10,ISL)
            PV_THETA_LEVS = STASH_LEVELS(1,NI)
            DO K = 1,PV_THETA_LEVS
              PV_THETA(K) = STASH_LEVELS(K+1,NI)/1000.
! ***** levels are stored as integers so divide by a thousand **
            ENDDO
          ELSE
            CMESSAGE =                                                  &
     &       ' ST_DIAG1 level not theta for pv_theta: exit routine'
            WRITE(6,*) RoutineName,CMESSAGE
            ICODE = 1
            RETURN
          END IF
        ELSE
          CMESSAGE =                                                    &
     &   ' ST_DIAG1 level not a LEVELS list for PV_Theta: exit routine'
            WRITE(6,*) RoutineName,CMESSAGE
          ICODE = 1
          RETURN
        END IF
      ELSE
        PV_THETA_LEVS = 1
      END IF

!L----------Extract required pressures for Potn_vort on press ----

      ISL=STINDEX(1,229,15,im_index)
      IF(ISL >  0) THEN
        IF(STLIST(10,ISL) <  0) THEN
          IF(STLIST(11,ISL) == 2) THEN
            NI = -STLIST(10,ISL)
            PV_PRESS_LEVS = STASH_LEVELS(1,NI)
            DO K = 1,PV_PRESS_LEVS
              PV_PRESS(K) = STASH_LEVELS(K+1,NI)/1000.0
! ***** levels are stored as integers so divide by a thousand **
            ENDDO
          ELSE
            CMESSAGE =                                                  &
     &       ' ST_DIAG1 level not pressure for pv_press: exit routine'
            WRITE(6,*) RoutineName,CMESSAGE
            ICODE = 1
            RETURN
          END IF
        ELSE
          CMESSAGE =                                                    &
     &   ' ST_DIAG1 level not a LEVELS list for PV_press exit routine'
            WRITE(6,*) RoutineName,CMESSAGE
          ICODE = 1
          RETURN
        END IF
      ELSE
        PV_PRESS_LEVS = 1
      END IF

!L----------Extract required PVs for Theta on pv -----------------

      ISL=STINDEX(1,230,15,im_index)
      IF(ISL >  0) THEN
         CMESSAGE =                                                     &
     & ' PV theta_on_pv (230,15) not enabled: exit diagnostic routine '
         WRITE(6,*) RoutineName,CMESSAGE
         ICODE = -1
      END IF

!L----------Extract required pressures for CAT_PROB_SINGLE--------------

      ISL=STINDEX(1,205,15,im_index)
      IF(ISL >  0) THEN
         CMESSAGE =                                                     &
     & 'CAT_PROB_SINGLE (205,15) not enabled: exit diagnostic routine '
         WRITE(6,*) RoutineName,CMESSAGE
         ICODE = -1
      END IF


!L-------------------Extract Reqd Pressures for Test Diagnostic 233--

      ISL=STINDEX(1,233,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            TESTD_P_LEVS=STASH_LEVELS(1,NI)
            DO K =1,TESTD_P_LEVS
              TESTD_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
            ENDDO
      ELSE
            TESTD_P_LEVS=1
      END IF

!L-------------------Extract Reqd Model levs for Test Diagnostic 234--

      ISL=STINDEX(1,234,15,im_index)
      IF(ISL >  0) THEN
            NI=-STLIST(10,ISL)
            TESTD_M_LEVS=STASH_LEVELS(1,NI)
            DO K =1,TESTD_M_LEVS
              TESTD_MODEL(K)=STASH_LEVELS(K+1,NI)  ! Converts to real
            ENDDO
      ELSE
            TESTD_M_LEVS=1
      END IF

!------ Extract Reqd model levels for Height on theta levels --------

      isl = stindex(1,101,15,im_index)
      If ( isl  >   0 ) Then
            ni = -stlist(10,isl)
            Htheta_m_levs = stash_levels(1,ni)
            Do k = 1, Htheta_m_levs
              Htheta_model(k) = stash_levels(k+1,ni)  ! integer levels
            End do
      Else
            Htheta_m_levs = 0
      End if

!------ Extract Reqd model levels for Height on rho levels --------

      isl = stindex(1,102,15,im_index)
      If ( isl  >   0 ) Then
            ni = -stlist(10,isl)
            Hrho_m_levs = stash_levels(1,ni)
            Do k = 1, Hrho_m_levs
              Hrho_model(k) = stash_levels(k+1,ni)  ! integer levels
            End do
      Else
            Hrho_m_levs = 0
      End if

!-------------------Extract Reqd Heights for p -------------

      If ( sf(108,15) ) Then
         isl=stindex(1,108,15,im_index)
         If ( isl  >   0 ) Then
            ni       = -stlist(10,isl)
            p_h_levs = stash_levels(1,ni)
            Do k=1,p_h_levs
               p_height(k) = stash_levels(k+1,ni)/1000.
            Enddo
         Else
            p_h_levs = 1
         End if
      End if

!-------------------Extract Reqd Heights for theta -------------

      isl=stindex(1,119,15,im_index)
      If ( isl  >   0 ) Then
         ni           = -stlist(10,isl)
         theta_h_levs = stash_levels(1,ni)
         Do k=1,theta_h_levs
            theta_height(k) = stash_levels(k+1,ni)/1000.
         End do
      Else
         theta_h_levs = 1
      End if


!-------------------Extract Reqd Heights for rho -------------

      isl=stindex(1,127,15,im_index)
      If ( isl  >   0 ) Then
         ni         = -stlist(10,isl)
         rho_h_levs = stash_levels(1,ni)
         Do k=1,rho_h_levs
            rho_height(k) = stash_levels(k+1,ni)/1000.
         End do
      Else
         rho_h_levs = 1
      End if

!-------------------Extract Reqd Heights for w -------------

      isl=stindex(1,142,15,im_index)
      If ( isl  >   0 ) Then
         ni       = -stlist(10,isl)
         w_h_levs = stash_levels(1,ni)
         Do k=1,w_h_levs
            w_height(k) = stash_levels(k+1,ni)/1000.
         End do
      Else
         w_h_levs = 1
      End if

!-------------------Extract Reqd Heights for u -------------

      isl=stindex(1,143,15,im_index)
      If ( isl  >   0 ) Then
         ni         = -stlist(10,isl)
         u_h_levs = stash_levels(1,ni)
         Do k=1,u_h_levs
            u_height(k) = stash_levels(k+1,ni)/1000.
         End do
      Else
         u_h_levs = 1
      End if

!-------------------Extract Reqd Heights for v -------------

      isl=stindex(1,144,15,im_index)
      If ( isl  >   0 ) Then
         ni         = -stlist(10,isl)
         v_h_levs = stash_levels(1,ni)
         Do k=1,v_h_levs
            v_height(k) = stash_levels(k+1,ni)/1000.
         End do
      Else
         v_h_levs = 1
      End if


      IF(ICODE == 0) THEN        ! Check error code

!----------------------------------------------------------------------
! Set up flags for rotating winds for a subset of lam diagnostics. This
! converts u,v of native equatorial lat-long lam grid back to u,v with
! respect to standard lat-long grid. Note that rotation of winds
! is not performed generically by STASH, but code needs to be
! introduced specifically for each (u,v) pair. The rotation flag in
! STASHmaster needs to be consistent and hence the flag is checked
! before performing rotation in dyn_diag. [Hence also, rotation can be
! suppressed by a STASHmaster change.]
! Note no explicit error trapping here.

      rot_uvcomB_press = .false.               ! default

      IF (model_domain == mt_lam .OR.                                   &
     &    model_domain == mt_cyclic_lam .or.                            &
     &    model_domain == mt_bi_cyclic_lam) THEN

! winds (B grid) on pressure levels =201/202
        IF(SF(201,sect).AND.SF(202,sect)) THEN

! DEPENDS ON: exppxi
           irotu=exppxi(im_ident,sect,201,ppx_rotate_code,              &
#include "argppx.h"
     &                  icode,cmessage)
! DEPENDS ON: exppxi
           irotv=exppxi(im_ident,sect,202,ppx_rotate_code,              &
#include "argppx.h"
     &                  icode,cmessage)

           IF( irotu == ppx_elf_rotated .AND.                           &
     &         irotv == ppx_elf_rotated )     THEN
              rot_uvcomB_press = .true.
           ENDIF

        ENDIF     ! 201/202

      ENDIF ! lam only
!L------------------Set up Pointers for STASHWORK -------------------

      PT002=SI(  2,sect,im_index)
      PT003=SI(  3,sect,im_index)
      PT201=SI(201,15,im_index)
      PT202=SI(202,15,im_index)
      PT203=SI(203,15,im_index)
      PT204=SI(204,15,im_index)
      PT205=SI(205,15,im_index)
      PT206=SI(206,15,im_index)
      PT207=SI(207,15,im_index)
      PT208=SI(208,15,im_index)
      PT209=SI(209,15,im_index)
      PT210=SI(210,15,im_index)
      PT211=SI(211,15,im_index)
      PT212=SI(212,15,im_index)
      PT213=SI(213,15,im_index)
      PT214=SI(214,15,im_index)
      PT215=SI(215,15,im_index)
      PT216=SI(216,15,im_index)
      PT217=SI(217,15,im_index)
      PT218=SI(218,15,im_index)
      PT219=SI(219,15,im_index)
      PT220=SI(220,15,im_index)
      PT221=SI(221,15,im_index)
      PT222=SI(222,15,im_index)
      PT223=SI(223,15,im_index)
      PT224=SI(224,15,im_index)
      PT225=SI(225,15,im_index)
      PT226=SI(226,15,im_index)
      PT227=SI(227,15,im_index)
      PT228=SI(228,15,im_index)
      PT229=SI(229,15,im_index)
      PT230=SI(230,15,im_index)
      PT231=SI(231,15,im_index)
      PT232=SI(232,15,im_index)
      PT233=SI(233,15,im_index)
      PT234=SI(234,15,im_index)
      PT235=SI(235,15,im_index)
      PT236=SI(236,15,im_index)
      PT237=SI(237,15,im_index)
      PT238=SI(238,15,im_index)
      PT239=SI(239,15,im_index)
      PT240=SI(240,15,im_index)
      PT241=SI(241,15,im_index)
      PT242=SI(242,15,im_index)
      PT243=SI(243,sect,im_index)
      PT244=SI(244,sect,im_index)
      PT245=SI(245,sect,im_index)
      PT246=SI(246,sect,im_index)
      PT260=SI(260,sect,im_index)
      PT261=SI(261,sect,im_index)
      PT262=SI(262,sect,im_index)
      PT263=SI(263,sect,im_index)
      PT264=SI(264,sect,im_index)
      PT265=SI(265,sect,im_index)
      PT266=SI(266,sect,im_index)
      pt101 = si(101,sect,im_index)
      pt102 = si(102,sect,im_index)
      pt108 = si(108,sect,im_index)
      pt119 = si(119,sect,im_index)
      pt127 = si(127,sect,im_index)
      pt142 = si(142,sect,im_index)
      pt143 = si(143,sect,im_index)
      pt144 = si(144,sect,im_index)
      PT270=SI(270,sect,im_index)
      PT271=SI(271,sect,im_index)

! Initialise STASHWORK array because DYN_DIAG does not initialise halos
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INT15
        STASHWORK(I)=0.
      ENDDO

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('DYN_DIAG',3)
      END IF

! DEPENDS ON: dyn_diag
      CALL Dyn_diag(                                                    &
! Primary data: in
     &  exner_rho_levels                                                &
     & ,RHO,U,V,W                                                       &
     & ,exner_theta_levels                                              &
     &,theta                                                            &
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
     &,pv_theta                                                         &
! Pressure levels for output arrays: in
     &,ucomB_press,vcomB_press                                          &
     &,UCOMP_PRESS,VCOMP_PRESS,WCOMP_PRESS                              &
     &,testd_press                                                      &
     &,p_height,theta_height,rho_height                                 &
     &,w_height,u_height,v_height                                       &
     &,pv_press                                                         &
! Model levels    for output arrays: in
     &,ucomB_model,vcomB_model                                          &
     &,testd_model                                                      &
     &,Htheta_model,Hrho_model                                          &
! Flags to request each diagnostic output field: in
! wind related diagnostics
     &,SF(  2,sect),SF(  3,sect)                                        &
     &,SF(201,sect),SF(202,sect)                                        &
     &,SF(243,sect),SF(244,sect),SF(242,sect)                           &
     &,SF(212,sect),SF(213,sect)                                        &
     &,SF(245,sect),SF(246,sect)                                        &
! PV related diagnostics
     &,SF(214,sect),SF(215,sect),SF(216,sect)                           &
     &,SF(217,sect),SF(218,sect),SF(229,sect)                           &
! test fields
     &,SF(231,sect),SF(232,sect),SF(233,sect),SF(234,sect)              &
! flux diagnostics
     &,SF(260,sect),SF(261,sect),SF(262,sect),SF(263,sect)              &
     &,SF(264,sect),SF(265,sect),SF(266,sect)                           &
! height and height level diagnostics
     &,SF(101,sect),SF(102,sect)                                        &
     &,SF(108,sect),SF(119,sect),SF(127,sect)                           &
     &,SF(143,sect),SF(144,sect),SF(142,sect)                           &
! other diagnostics
     &,SF(270,sect),SF(271,sect)                                        &
! Flags for wind rotation (lam grid): in
     &,rot_uvcomB_press                                                 &
! Diagnostics lengths: in
     &,ucomB_m_levs,vcomB_m_levs                                        &
     &,ucomB_p_levs,vcomB_p_levs                                        &
     &,ucomp_p_levs,vcomp_p_levs,wcomp_p_levs                           &
     &,pv_theta_levs, pv_press_levs                                     &
     &,testd_p_levs,testd_m_levs                                        &
     &,Htheta_m_levs,Hrho_m_levs                                        &
     &,p_h_levs,theta_h_levs,rho_h_levs,w_h_levs,u_h_levs,v_h_levs      &
! Diagnostic arrays: out
! wind related diagnostics
     &,Stashwork(pt002),Stashwork(pt003)                                &
     &,Stashwork(pt201),Stashwork(pt202)                                &
     &,Stashwork(pt243),Stashwork(pt244),Stashwork(pt242)               &
     &,Stashwork(pt212),Stashwork(pt213)                                &
     &,Stashwork(pt245),Stashwork(pt246)                                &
! PV related diagnostics
     &,Stashwork(pt214),Stashwork(pt215),Stashwork(pt216)               &
     &,Stashwork(pt217),Stashwork(pt218),Stashwork(pt229)               &
! test fields
     &,Stashwork(pt231),Stashwork(pt232)                                &
     &,Stashwork(pt233),Stashwork(pt234)                                &
! flux diagnostics
     &,Stashwork(pt260),Stashwork(pt261),Stashwork(pt262)               &
     &,Stashwork(pt263),Stashwork(pt264),Stashwork(pt265)               &
     &,Stashwork(pt266)                                                 &
! height and height level diagnostics
     &,Stashwork(pt101),Stashwork(pt102)                                &
     &,Stashwork(pt108),Stashwork(pt119),Stashwork(pt127)               &
     &,Stashwork(pt143),Stashwork(pt144),Stashwork(pt142)               &
! other diagnostics
     &,Stashwork(pt270),Stashwork(pt271)                                &
     &)

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('DYN_DIAG',4)
      END IF


      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',3)
      END IF

      IF(ICODE /= 0) THEN
        RETURN
      ENDIF

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,15,STASHWORK,                                &
#include "argd1.h"
#include "argdumga.h"
#include "argsts.h"
#include "argppx.h"
     &           ICODE,CMESSAGE)

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',4)
      END IF

      ENDIF   ! Check on error code ICODE

! Error or warning exit
      IF(ICODE /= 0) THEN
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ICODE,CMESSAGE)
      ENDIF

      RETURN
      END SUBROUTINE ST_DIAG1
#endif
