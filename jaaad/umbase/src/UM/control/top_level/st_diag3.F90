#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ End of Timestep climate diagnostics control routine
!
! Subroutine Interface:
      SUBROUTINE ST_DIAG3( STASHWORK,INT30,                             &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arg_atm_fields.h"
#include "argcona.h"
#include "argppx.h"
     &    energy_corr_now,                                              &
     &    inc_u, inc_v, inc_w, inc_t,                                   &
     &    inc_q, inc_qcl, inc_qcf, inc_rho,                             &
     &                    ICODE,CMESSAGE)

      Use level_heights_mod
      Use trignometric_mod,  Only : sin_theta_latitude,                 &
     &     cos_theta_latitude, sin_theta_longitude, cos_theta_longitude
      Use rad_mask_trop_mod, Only : max_trop_level,min_trop_level

      IMPLICIT NONE
! Description:
! Calculates end of timestep diagnostics (held in STASH section 30)
! of which there are three types:
! 1) Model level fields and products
! 2) Pressure level fields and products ie U V T UU VV TT UT VT
! 3) Single level fields eg mountain torque and column integrations
!    eq atmospheric mass
!
! Method:
!   This routine sets controls the EoT diagnostics, and sets up
!   array sizes etc for a call to EOT_DIAG which calculates the
!   fields. STASHWORK is filled with the required data, and
!   passed up to atm_step2, ready for a call to stash.
!   The stash items are split into centuries, as follows:
!
!   000s   Standard Model level fields and products
!   100s   Other Model level fields
!   200s   Standard Pressure level fields and products
!   300s   Other Pressure level fields
!   400s   Surface fields and Column integrals
!
!  Standard diagnostics are assigned a number, which are:
!  1 U wind, 2 V wind, 3 W wind, 4 Temperature, 5 Q, 6 RH, 7 Z, 8 Omega
!
!
!  Further user and technical documentation can be found under:
!  http://www-hc/~hadsu/
!
! Current Code Owner: Simon Wilson
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.1     03/03/2000 Original code. SS Wilson
!
! 5.2     15/08/2000 Pstar added. SS Wilson
!
! 5.2     22/09/2000 Extra diagnostics and bug fixes. SS Wilson
!
! 5.3     08/08/2001 Tropopause diagnostics added. DGH Tan
! 5.4     24/04/2002 Add RH on model levels, wbig and the square of
!                    increments. R A Stratton
!         04/07/2002 Add column integral of cvT and gr. R A Stratton.
! 5.5     13/12/2002 Add wbig=0.1, Correct an extend model level
!                    products. R A Stratton
!         18/12/2002 Add more column integrals. R A Stratton.
!         31/12/2002 Dry mass weighting on model levels. R A Stratton
! 5.5     23/01/2003 Correct mountain torque diagnostic and add
!                    surface pressure drag diagnostics. S. Webster
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
! Global variables (*CALLed COMDECKs etc...):

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
#include "ppxlook.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "cphyscon.h"
#include "ctime.h"
#include "c_eta_pmsl.h"

      Integer, parameter ::                                             &
     &  npress_diags=8                                                  &
                                 ! No. of diags on pressure levels
     &, nmodel_diags=7           ! No. of diags on modellevels

!   Scalar  arguments with intent(in):

!   Array  arguments with intent(in):
      real energy_corr_now                                              &
     &,  inc_u(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &        model_levels)                                             &
     &, inc_v(1-offx:row_length+offx, 1-offy:n_rows+offy,               &
     &        model_levels)                                             &
     &, inc_w(row_length, rows, model_levels)                           &
     &, inc_t(1-offx:row_length+offx,                                   &
     &               1-offy:rows+offy, model_levels)                    &
     &, inc_q(1-halo_i:row_length+halo_i,                               &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_qcl(1-halo_i:row_length+halo_i,                             &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_qcf(1-halo_i:row_length+halo_i,                             &
     &               1-halo_j:rows+halo_j,wet_levels)                   &
     &, inc_rho(1-offx:row_length+offx,                                 &
     &             1-offy:rows+offy, model_levels)

!   Scalar arguments with intent(InOut):
      INTEGER                                                           &
     &        INT30,                                                    &
                                ! Dummy for STASH_MAXLEN(30)
     &        ICODE              ! Out return code : 0 Normal exit
                                !                 : >0 Error exit

      CHARACTER*256                                                     &
     &        CMESSAGE          ! Out error message if ICODE > 0

!   Array  arguments with intent(InOut):
      REAL                                                              &
     &  STASHWORK(*)

!   Scalar arguments with intent(out):
      INTEGER                                                           &
     &  U_M_LEVS,                                                       &
     &  V_M_LEVS,                                                       &
     &  W_M_LEVS,                                                       &
     &  T_M_LEVS,                                                       &
     &  Q_M_LEVS,                                                       &
     &  Z_M_LEVS,                                                       &
     &  KE_M_LEVS,                                                      &
     &  wbig_m_levs,                                                    &
     &  wbig2_m_levs,                                                   &
     &  dry_mass_m_levs,                                                &
     &  Rh_m_levs,                                                      &
     &  U_MM_LEVS,                                                      &
     &  V_MM_LEVS,                                                      &
     &  W_MM_LEVS,                                                      &
     &  T_MM_LEVS,                                                      &
     &  Q_MM_LEVS,                                                      &
     &  Z_MM_LEVS,                                                      &
     &  KE_MM_LEVS,                                                     &
     &  TEOT_M_LEVS,                                                    &
     &  U_INC_LEVS,                                                     &
     &  V_INC_LEVS,                                                     &
     &  W_INC_LEVS,                                                     &
     &  T_INC_LEVS,                                                     &
     &  Q_INC_LEVS,                                                     &
     &  QCL_INC_LEVS,                                                   &
     &  QCF_INC_LEVS,                                                   &
     &  RHO_INC_LEVS,                                                   &
     &  U2_INC_LEVS,                                                    &
     &  V2_INC_LEVS,                                                    &
     &  W2_INC_LEVS,                                                    &
     &  T2_INC_LEVS,                                                    &
     &  Q2_INC_LEVS,                                                    &
     &  QCL2_INC_LEVS,                                                  &
     &  QCF2_INC_LEVS,                                                  &
     &  RHO2_INC_LEVS,                                                  &
     &  U_P_LEVS,                                                       &
     &  V_P_LEVS,                                                       &
     &  W_P_LEVS,                                                       &
     &  T_P_LEVS,                                                       &
     &  Q_P_LEVS,                                                       &
     &  RH_P_LEVS,                                                      &
     &  Z_P_LEVS,                                                       &
     &  OM_P_LEVS,                                                      &
     &  HEAVY_P_LEVS,                                                   &
     &  TV_P_LEVS,                                                      &
     &  TVOM_P_LEVS,                                                    &
     &  FIELD1_P_LEVS,                                                  &
     &  FIELD2_P_LEVS,                                                  &
     &  n_levels                                                        &
     &  ,im_ident                                                       &
                                !  Internal Model Identifier
     &  ,im_index               !  Internal Model Index for Stash Arrays

      INTEGER                                                           &
     &  PT001,PT002,PT003,PT004,PT005,PT006,PT007,                      &
     &  PT101,PT102,PT103,PT104,PT105,PT106,PT107,                      &
     &  PT111,PT112,PT113,PT114,PT115,                                  &
     &  PT181,PT182,PT183,PT184,PT185,PT186,PT187,PT188,                &
     &  PT171,PT172,PT173,PT174,PT175,PT176,PT177,PT178,                &
     &  PT201,PT202,PT203,PT204,PT205,PT206,PT207,PT208,                &
     &  PT301,PT302,PT303,                                              &
     &  PT401,PT402,PT403,PT404,PT405,PT406,PT407,PT408,PT409,PT410,    &
     &  PT411,PT412,PT413,PT414,PT415,PT416,PT417,PT418,PT419,PT420,    &
     &  PT421,PT422,PT423,PT424,PT425,PT426,PT427,PT428,PT429,PT430,    &
     &  PT431,PT432,PT433,PT434,PT435,PT436,PT437,PT438,PT439           &
     & ,PT440,PT441,PT451,PT452,PT453,PT454

!   Array  arguments with intent(out):
      REAL                                                              &
     &  U_PRESS(NUM_STASH_LEVELS)                                       &
     &  ,V_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,W_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,T_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,Q_PRESS(NUM_STASH_LEVELS)                                      &
     &  ,RH_PRESS(NUM_STASH_LEVELS)                                     &
     &  ,Z_PRESS(NUM_STASH_LEVELS)                                      &
                                   ! Z pressure levels
     &  ,OM_PRESS(NUM_STASH_LEVELS)                                     &
     &  ,HEAVY_PRESS(NUM_STASH_LEVELS)                                  &
     &  ,PRESS_LEVS(NUM_STASH_LEVELS)

      LOGICAL U_M_LIST(MODEL_LEVELS)                                    &
     &  ,V_M_LIST(MODEL_LEVELS)                                         &
     &  ,W_M_LIST(MODEL_LEVELS)                                         &
     &  ,T_M_LIST(MODEL_LEVELS)                                         &
     &  ,Q_M_LIST(MODEL_LEVELS)                                         &
     &  ,Z_M_LIST(MODEL_LEVELS)                                         &
     &  ,KE_M_LIST(MODEL_LEVELS)                                        &
     &  ,wbig_m_list(MODEL_LEVELS)                                      &
     &  ,wbig2_m_list(MODEL_LEVELS)                                     &
     &  ,RH_m_list(MODEL_LEVELS)                                        &
     &  ,DRY_MASS_M_LIST(MODEL_LEVELS)                                  &
     &  ,U_MM_LIST(MODEL_LEVELS)                                        &
     &  ,V_MM_LIST(MODEL_LEVELS)                                        &
     &  ,W_MM_LIST(MODEL_LEVELS)                                        &
     &  ,T_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Q_MM_LIST(MODEL_LEVELS)                                        &
     &  ,Z_MM_LIST(MODEL_LEVELS)                                        &
     &  ,KE_MM_LIST(MODEL_LEVELS)                                       &
     &  ,TEOT_M_LIST(MODEL_LEVELS)                                      &
     &  ,U_INC_LIST(MODEL_LEVELS)                                       &
     &  ,V_INC_LIST(MODEL_LEVELS)                                       &
     &  ,W_INC_LIST(MODEL_LEVELS)                                       &
     &  ,T_INC_LIST(MODEL_LEVELS)                                       &
     &  ,Q_INC_LIST(MODEL_LEVELS)                                       &
     &  ,QCL_INC_LIST(MODEL_LEVELS)                                     &
     &  ,QCF_INC_LIST(MODEL_LEVELS)                                     &
     &  ,RHO_INC_LIST(MODEL_LEVELS)                                     &
     &  ,U2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,V2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,W2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,T2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,Q2_INC_LIST(MODEL_LEVELS)                                      &
     &  ,QCL2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,QCF2_INC_LIST(MODEL_LEVELS)                                    &
     &  ,RHO2_INC_LIST(MODEL_LEVELS)

      logical prod_m_list(model_levels,nmodel_diags,nmodel_diags)

      INTEGER                                                           &
     &  prod_IND(NUM_STASH_LEVELS*2,npress_diags,npress_diags),         &
     &  prod_p_levs(npress_diags,npress_diags)

! Local parameters:

! Local scalars:
      INTEGER                                                           &
     &  I,                                                              &
     &  NI,                                                             &
     &  K,                                                              &
     &  ISL,                                                            &
     &  BL,                                                             &
     &  TL,                                                             &
     &  LEVEL,                                                          &
     &  code30,nii,nij,j

      character*2 cfield1,cfield2

! Local dynamic arrays:
      REAL                                                              &
     &  FIELD1_PRESS(NUM_STASH_LEVELS)                                  &
     &  ,FIELD2_PRESS(NUM_STASH_LEVELS)

! Function & Subroutine calls:

      EXTERNAL                                                          &
     &  STASH,                                                          &
     &  TIMER,                                                          &
     &  EOT_DIAG,                                                       &
     &  SET_LEVELS_LIST,                                                &
     &  CHECK_PROD_LEVS

!- End of header

!     Set to atmosphere internal model
      im_ident = atmos_im
      im_index = internal_model_index(im_ident)



!--------------------Extract Reqd Model Levels for U-------------
      ISL=STINDEX(1,001,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_M_LIST(I)) U_M_LEVS=U_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_M_LIST(I)=.FALSE.
        ENDDO
        U_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for V-------------
      ISL=STINDEX(1,002,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_M_LIST(I)) V_M_LEVS=V_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_M_LIST(I)=.FALSE.
        ENDDO
        V_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W-------------
      ISL=STINDEX(1,003,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_M_LIST(I)) W_M_LEVS=W_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_M_LIST(I)=.FALSE.
        ENDDO
        W_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T-------------
      ISL=STINDEX(1,004,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_M_LIST(I)) T_M_LEVS=T_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_M_LIST(I)=.FALSE.
        ENDDO
        T_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q-------------
      ISL=STINDEX(1,005,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_M_LIST(I)) Q_M_LEVS=Q_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_M_LIST(I)=.FALSE.
        ENDDO
        Q_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Z-------------
      ISL=STINDEX(1,006,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Z_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Z_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Z_M_LIST(I)) Z_M_LEVS=Z_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Z_M_LIST(I)=.FALSE.
        ENDDO
        Z_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for KE-------------
      ISL=STINDEX(1,007,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    KE_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        KE_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (KE_M_LIST(I)) KE_M_LEVS=KE_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          KE_M_LIST(I)=.FALSE.
        ENDDO
        KE_M_LEVS=1
      ENDIF
!--------------Extract Reqd Model Levels for Dry mass weighting----
      ISL=STINDEX(1,115,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &  dry_mass_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        DRY_MASS_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (DRY_MASS_M_LIST(I)) DRY_MASS_M_LEVS=DRY_MASS_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          DRY_MASS_M_LIST(I)=.FALSE.
        ENDDO
        DRY_MASS_M_LEVS=1
      ENDIF

!-Extract Reqd model levels for products and store in prod_m_list

      do i=1,nmodel_diags
        do j=1,nmodel_diags
          code30=i*10+j
          IF (SF(code30,30)) THEN
            ISL=STINDEX(1,code30,30,im_index)
! DEPENDS ON: set_levels_list
            CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),     &
     &        PROD_M_LIST(1,i,j),                                       &
     &        STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
          ELSE
            DO k=1,MODEL_LEVELS
              prod_M_LIST(k,i,j)=.FALSE.
            ENDDO
          ENDIF
        enddo
      enddo


!--------------------Extract Reqd Model Levels for wbig----------
      ISL=STINDEX(1,112,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Wbig_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Wbig_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Wbig_M_LIST(I)) Wbig_M_LEVS=Wbig_M_LEVS+1
        ENDDO
      ELSE
        Wbig_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for wbig2----------
      ISL=STINDEX(1,114,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Wbig2_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Wbig2_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Wbig2_M_LIST(I)) Wbig2_M_LEVS=Wbig2_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Wbig2_M_LIST(I)=.FALSE.
        ENDDO
         Wbig2_M_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for RH ----------
      ISL=STINDEX(1,113,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RH_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RH_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RH_M_LIST(I)) RH_M_LEVS=RH_M_LEVS+1
        ENDDO
      ELSE
        RH_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U mass weighted
      ISL=STINDEX(1,101,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_MM_LIST(I)) U_MM_LEVS=U_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_MM_LIST(I)=.FALSE.
        ENDDO
        U_MM_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for V mass weighted
      ISL=STINDEX(1,102,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_MM_LIST(I)) V_MM_LEVS=V_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_MM_LIST(I)=.FALSE.
        ENDDO
        V_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W mass weighted
      ISL=STINDEX(1,103,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_MM_LIST(I)) W_MM_LEVS=W_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_MM_LIST(I)=.FALSE.
        ENDDO
        W_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T mass weighted
      ISL=STINDEX(1,104,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_MM_LIST(I)) T_MM_LEVS=T_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_MM_LIST(I)=.FALSE.
        ENDDO
        T_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q mass weighted
      ISL=STINDEX(1,105,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_MM_LIST(I)) Q_MM_LEVS=Q_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_MM_LIST(I)=.FALSE.
        ENDDO
        Q_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Z mass weighted
      ISL=STINDEX(1,106,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Z_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Z_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Z_MM_LIST(I)) Z_MM_LEVS=Z_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Z_MM_LIST(I)=.FALSE.
        ENDDO
        Z_MM_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for KE mass weighted
      ISL=STINDEX(1,107,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    KE_MM_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        KE_MM_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (KE_MM_LIST(I)) KE_MM_LEVS=KE_MM_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          KE_MM_LIST(I)=.FALSE.
        ENDDO
        KE_MM_LEVS=1
      ENDIF


!--------------------Extract Reqd Model Levels for T at EOT
      ISL=STINDEX(1,111,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    TEOT_M_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        TEOT_M_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (TEOT_M_LIST(I)) TEOT_M_LEVS=TEOT_M_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          TEOT_M_LIST(I)=.FALSE.
        ENDDO
        TEOT_M_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T increment
      ISL=STINDEX(1,181,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T_INC_LIST(I)) T_INC_LEVS=T_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          T_INC_LIST(I)=.FALSE.
        ENDDO
        T_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q increment
      ISL=STINDEX(1,182,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q_INC_LIST(I)) Q_INC_LEVS=Q_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          Q_INC_LIST(I)=.FALSE.
        ENDDO
        Q_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCL increment
      ISL=STINDEX(1,183,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCL_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCL_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCL_INC_LIST(I)) QCL_INC_LEVS=QCL_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          QCL_INC_LIST(I)=.FALSE.
        ENDDO
        QCL_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCF increment
      ISL=STINDEX(1,184,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCF_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCF_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCF_INC_LIST(I)) QCF_INC_LEVS=QCF_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          QCF_INC_LIST(I)=.FALSE.
        ENDDO
        QCF_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U increment
      ISL=STINDEX(1,185,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U_INC_LIST(I)) U_INC_LEVS=U_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          U_INC_LIST(I)=.FALSE.
        ENDDO
        U_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for V increment
      ISL=STINDEX(1,186,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V_INC_LIST(I)) V_INC_LEVS=V_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          V_INC_LIST(I)=.FALSE.
        ENDDO
        V_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for W increment
      ISL=STINDEX(1,187,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W_INC_LIST(I)) W_INC_LEVS=W_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          W_INC_LIST(I)=.FALSE.
        ENDDO
        W_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for RHO increment
      ISL=STINDEX(1,188,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RHO_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RHO_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RHO_INC_LIST(I)) RHO_INC_LEVS=RHO_INC_LEVS+1
        ENDDO
      ELSE
        DO I=1,MODEL_LEVELS
          RHO_INC_LIST(I)=.FALSE.
        ENDDO
        RHO_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for T increment **2
      ISL=STINDEX(1,171,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    T2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        T2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (T2_INC_LIST(I)) T2_INC_LEVS=T2_INC_LEVS+1
        ENDDO
      ELSE
        T2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for Q increment **2
      ISL=STINDEX(1,172,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    Q2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        Q2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (Q2_INC_LIST(I)) Q2_INC_LEVS=Q2_INC_LEVS+1
        ENDDO
      ELSE
        Q2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCL increment **2
      ISL=STINDEX(1,173,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCL2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCL2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCL2_INC_LIST(I)) QCL2_INC_LEVS=QCL2_INC_LEVS+1
        ENDDO
      ELSE
        QCL2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for QCF increment **2
      ISL=STINDEX(1,174,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    QCF2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        QCF2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (QCF2_INC_LIST(I)) QCF2_INC_LEVS=QCF2_INC_LEVS+1
        ENDDO
      ELSE
        QCF2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for U increment **2
      ISL=STINDEX(1,175,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    U2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        U2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (U2_INC_LIST(I)) U2_INC_LEVS=U2_INC_LEVS+1
        ENDDO
      ELSE
        U2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for V increment **2
      ISL=STINDEX(1,176,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    V2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        V2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (V2_INC_LIST(I)) V2_INC_LEVS=V2_INC_LEVS+1
        ENDDO
      ELSE
        V2_INC_LEVS=1
      ENDIF
!--------------------Extract Reqd Model Levels for W increment **2
      ISL=STINDEX(1,177,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    W2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        W2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (W2_INC_LIST(I)) W2_INC_LEVS=W2_INC_LEVS+1
        ENDDO
      ELSE
        W2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Model Levels for RHO increment **2
      ISL=STINDEX(1,178,30,im_index)
      IF(ISL >  0) THEN
! DEPENDS ON: set_levels_list
        CALL SET_LEVELS_LIST(MODEL_LEVELS,NITEMS,STLIST(1,ISL),         &
     &    RHO2_INC_LIST,STASH_LEVELS,NUM_STASH_LEVELS+1,ICODE,CMESSAGE)
        RHO2_INC_LEVS=0
        DO I=1,MODEL_LEVELS
          IF (RHO2_INC_LIST(I)) RHO2_INC_LEVS=RHO2_INC_LEVS+1
        ENDDO
      ELSE
        RHO2_INC_LEVS=1
      ENDIF

!--------------------Extract Reqd Pressures for U_COMP_P-------------

      ISL=STINDEX(1,201,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        U_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,U_P_LEVS
          U_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        U_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for V_COMP_P-------------

      ISL=STINDEX(1,202,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        V_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,V_P_LEVS
          V_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        V_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for W_COMP_P-------------

      ISL=STINDEX(1,203,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        W_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,W_P_LEVS
          W_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        W_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for T on wind grid-------

      ISL=STINDEX(1,204,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        T_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,T_P_LEVS
          T_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        T_P_LEVS=1
      END IF


!--------------------Extract Reqd Pressures for q on wind grid-------

      ISL=STINDEX(1,205,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        Q_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,Q_P_LEVS
          Q_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        Q_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for RH-----

      ISL=STINDEX(1,206,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        RH_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,RH_P_LEVS
          RH_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        RH_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for geopotential height-----

      ISL=STINDEX(1,207,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        Z_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,Z_P_LEVS
          Z_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        Z_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for Omega---

      ISL=STINDEX(1,208,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        OM_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,OM_P_LEVS
          OM_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        OM_P_LEVS=1
      END IF

!--------------------Extract Reqd Pressures for products---


      do i=1,npress_diags
        do j=1,npress_diags
          if (i  ==  1) then
            cfield1='u'
            field1_press=u_press
            field1_p_levs=u_p_levs
          else if (i  ==  2) then
            cfield1='v'
            field1_press=v_press
            field1_p_levs=v_p_levs
          else if (i  ==  3) then
            cfield1='w'
            field1_press=w_press
            field1_p_levs=w_p_levs
          else if (i  ==  4) then
            cfield1='t'
            field1_press=t_press
            field1_p_levs=t_p_levs
          else if (i  ==  5) then
            cfield1='q'
            field1_press=q_press
            field1_p_levs=q_p_levs
          else if (i  ==  6) then
            cfield1='rh'
            field1_press=rh_press
            field1_p_levs=rh_p_levs
          else if (i  ==  7) then
            cfield1='z'
            field1_press=z_press
            field1_p_levs=z_p_levs
          else if (i  ==  8) then
            cfield1='om'
            field1_press=om_press
            field1_p_levs=om_p_levs
          endif

          if (j  ==  1) then
            cfield2='u'
            field2_press=u_press
            field2_p_levs=u_p_levs
          else if (j  ==  2) then
            cfield2='v'
            field2_press=v_press
            field2_p_levs=v_p_levs
          else if (j  ==  3) then
            cfield2='w'
            field2_press=w_press
            field2_p_levs=w_p_levs
          else if (j  ==  4) then
            cfield2='t'
            field2_press=t_press
            field2_p_levs=t_p_levs
          else if (j  ==  5) then
            cfield2='q'
            field2_press=q_press
            field2_p_levs=q_p_levs
          else if (j  ==  6) then
            cfield2='rh'
            field2_press=rh_press
            field2_p_levs=rh_p_levs
          else if (j  ==  7) then
            cfield2='z'
            field2_press=z_press
            field2_p_levs=z_p_levs
          else if (j  ==  8) then
            cfield2='om'
            field2_press=om_press
            field2_p_levs=om_p_levs
          endif
          code30=200+i*10+j
          IF (SF(code30,30)) THEN
            IF ((.NOT.SF(200+i,30)).OR.(.NOT.SF(200+j,30))) THEN
              CMESSAGE=' '
              CMESSAGE='ST_DIAG3 : '//cfield1//'*'//cfield2//           &
     &          ' error '//cfield1//' and '//cfield2//                  &
     &          ' must be requested'
              ICODE=1
              GOTO 999
            ELSE
              NI=-STLIST(10,STINDEX(1,code30,30,im_index))
              NII=-STLIST(10,STINDEX(1,200+i,30,im_index))
              NIJ=-STLIST(10,STINDEX(1,200+j,30,im_index))

! DEPENDS ON: check_prod_levs
              CALL CHECK_PROD_LEVS(                                     &
     &          NUM_STASH_LEVELS,STASH_LEVELS,NI,                       &
     &          FIELD1_P_LEVS,FIELD1_PRESS,                             &
     &          FIELD2_P_LEVS,FIELD2_PRESS,                             &
     &          prod_P_LEVS(i,j),prod_IND(1,i,j))
            END IF
          ELSE
            prod_P_LEVS(i,j)=1
          ENDIF
        enddo
      enddo

!     Do checks for other pressure fields

      if (sf(302,30).and..not.                                          &
     &  (sf(245,30).and.sf(204,30).and.sf(205,30))) then
        CMESSAGE='ST_DIAG3 : Virtual temperature'//                     &
     &          ' error T*Q must also be requested'
        ICODE=1
        GOTO 999
      endif

      if (sf(303,30).and..not.                                          &
     &  (sf(245,30).and.sf(204,30).and.                                 &
     &  sf(205,30).and.sf(208,30))) then
        CMESSAGE='ST_DIAG3 : Virtual temperature*omega'//               &
     &          ' error T*Q and omega must also be requested'
        ICODE=1
        GOTO 999
      endif



!--------------------Extract Reqd Pressures for Heavyside function---

      ISL=STINDEX(1,301,30,im_index)
      IF(ISL >  0) THEN
        NI=-STLIST(10,ISL)
        HEAVY_P_LEVS=STASH_LEVELS(1,NI)
        DO K =1,HEAVY_P_LEVS
          HEAVY_PRESS(K)=STASH_LEVELS(K+1,NI)/1000.0
        ENDDO
      ELSE
        HEAVY_P_LEVS=1
      END IF

!Lset up level information for tv
      if (sf(302,30)) then
        if (sf(204,30).and.sf(205,30)) then
          tv_p_levs=prod_p_levs(4,5)
        else
          tv_p_levs=1
          CMESSAGE='ST_DIAG3 : T and Q must be selected to calculate '//&
     &      'virtual temperature'
          ICODE=1
          GOTO 999
        endif
      else
        tv_p_levs=1
      endif

!Lset up level information for tv*om
      if (sf(303,30)) then
        if (sf(204,30).and.sf(205,30).and.sf(208,30)) then
          tvom_p_levs=prod_p_levs(4,5)
        else
          tvom_p_levs=1
          CMESSAGE='ST_DIAG3 : T,Q and Omega'//                         &
     &      'must be selected to calculate '//                          &
     &      'virtual temperature'
          ICODE=1
          GOTO 999
        endif
      else
        tvom_p_levs=1
      endif



!-------------------Set up Pointers for STASHWORK -------------------
      PT001=SI(001,30,im_index)
      PT002=SI(002,30,im_index)
      PT003=SI(003,30,im_index)
      PT004=SI(004,30,im_index)
      PT005=SI(005,30,im_index)
      PT006=SI(006,30,im_index)
      PT007=SI(007,30,im_index)
      PT101=SI(101,30,im_index)
      PT102=SI(102,30,im_index)
      PT103=SI(103,30,im_index)
      PT104=SI(104,30,im_index)
      PT105=SI(105,30,im_index)
      PT106=SI(106,30,im_index)
      PT107=SI(107,30,im_index)
      PT111=SI(111,30,im_index)
      PT112=SI(112,30,im_index)
      PT113=SI(113,30,im_index)
      PT114=SI(114,30,im_index)
      PT115=SI(115,30,im_index)
      PT171=SI(171,30,im_index)
      PT172=SI(172,30,im_index)
      PT173=SI(173,30,im_index)
      PT174=SI(174,30,im_index)
      PT175=SI(175,30,im_index)
      PT176=SI(176,30,im_index)
      PT177=SI(177,30,im_index)
      PT178=SI(178,30,im_index)
      PT181=SI(181,30,im_index)
      PT182=SI(182,30,im_index)
      PT183=SI(183,30,im_index)
      PT184=SI(184,30,im_index)
      PT185=SI(185,30,im_index)
      PT186=SI(186,30,im_index)
      PT187=SI(187,30,im_index)
      PT188=SI(188,30,im_index)
      PT201=SI(201,30,im_index)
      PT202=SI(202,30,im_index)
      PT203=SI(203,30,im_index)
      PT204=SI(204,30,im_index)
      PT205=SI(205,30,im_index)
      PT206=SI(206,30,im_index)
      PT207=SI(207,30,im_index)
      PT208=SI(208,30,im_index)
      PT301=SI(301,30,im_index)
      PT302=SI(302,30,im_index)
      PT303=SI(303,30,im_index)
      PT401=SI(401,30,im_index)
      PT402=SI(402,30,im_index)
      PT403=SI(403,30,im_index)
      PT404=SI(404,30,im_index)
      PT405=SI(405,30,im_index)
      PT406=SI(406,30,im_index)
      PT407=SI(407,30,im_index)
      PT408=SI(408,30,im_index)
      PT409=SI(409,30,im_index)
      PT410=SI(410,30,im_index)
      PT411=SI(411,30,im_index)
      PT412=SI(412,30,im_index)
      PT413=SI(413,30,im_index)
      PT414=SI(414,30,im_index)
      PT415=SI(415,30,im_index)
      PT416=SI(416,30,im_index)
      PT417=SI(417,30,im_index)
      PT418=SI(418,30,im_index)
      PT419=SI(419,30,im_index)
      PT420=SI(420,30,im_index)
      PT421=SI(421,30,im_index)
      PT422=SI(422,30,im_index)
      PT423=SI(423,30,im_index)
      PT424=SI(424,30,im_index)
      PT425=SI(425,30,im_index)
      PT426=SI(426,30,im_index)
      PT427=SI(427,30,im_index)
      PT428=SI(428,30,im_index)
      PT429=SI(429,30,im_index)
      PT430=SI(430,30,im_index)
      PT431=SI(431,30,im_index)
      PT432=SI(432,30,im_index)
      PT433=SI(433,30,im_index)
      PT434=SI(434,30,im_index)
      PT435=SI(435,30,im_index)
      PT436=SI(436,30,im_index)
      PT437=SI(437,30,im_index)
      PT438=SI(438,30,im_index)
      PT439=SI(439,30,im_index)
      PT440=SI(440,30,im_index)
      PT441=SI(441,30,im_index)
      PT451=SI(451,30,im_index)
      PT452=SI(452,30,im_index)
      PT453=SI(453,30,im_index)
      PT454=SI(454,30,im_index)

! Initialise STASHWORK array because EOT_DIAG does not initialise halos
!* DIR$ CACHE_BYPASS STASHWORK
      DO I=1,INT30
        STASHWORK(I)=0.
      ENDDO

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EOT_DIAG',3)
      END IF

! DEPENDS ON: eot_diag
      CALL Eot_diag(                                                    &
! Primary data: in
     &  PSTAR,P,RHO,U,V,W                                               &
     &  ,THETA,Q,QCL,QCF                                                &
     &  ,p_theta_levels ,exner_rho_levels                               &
     &  ,exner_theta_levels,energy_corr_now                             &
     &  ,inc_u, inc_v, inc_w, inc_t ,inc_q, inc_qcl, inc_qcf, inc_rho   &
! Grid sizes and definition: in
     &  ,rows,n_rows,row_length,model_levels,wet_levels,bl_levels       &
     &  ,theta_field_size,u_field_size,v_field_size                     &
     &  ,eta_theta_levels,eta_rho_levels,r_theta_levels,r_rho_levels    &
     &  ,sin_theta_latitude,cos_theta_latitude                          &
     &  ,sin_theta_longitude,cos_theta_longitude                        &
     &  ,Model_domain ,delta_lambda ,delta_phi,global_row_length        &
! Grid coordinates: in
     &  ,A_REALHD(rh_deltaEW),A_REALHD(rh_deltaNS)                      &
     &  ,A_REALHD(rh_baselat),A_REALHD(rh_baselong)                     &
     &  ,A_REALHD(rh_rotlat),A_REALHD(rh_rotlong)                       &
! Pre-calculated grid associated arrays: in
     &  ,r_at_u,r_at_v                                                  &
! Time information: in
     &  ,forecast_hrs                                                   &
! Pressure levels for output arrays: in
     &  ,U_PRESS,V_PRESS,W_PRESS,T_PRESS,Q_PRESS,RH_PRESS               &
     &  ,Z_PRESS,OM_PRESS,HEAVY_PRESS                                   &
! Flags to request each diagnostic output field: in
     &  ,SF(001,30),SF(002,30),SF(003,30),SF(004,30),SF(005,30)         &
     &  ,SF(006,30),SF(007,30)                                          &
     &  ,SF(101,30),SF(102,30),SF(103,30),SF(104,30),SF(105,30)         &
     &  ,SF(106,30),SF(107,30)                                          &
     &  ,SF(111,30),SF(112,30),SF(114,30),SF(113,30),SF(115,30)         &
     &  ,SF(181,30),SF(182,30),SF(183,30),SF(184,30),SF(185,30)         &
     &  ,SF(186,30),SF(187,30),SF(188,30)                               &
     &  ,SF(171,30),SF(172,30),SF(173,30),SF(174,30),SF(175,30)         &
     &  ,SF(176,30),SF(177,30),SF(178,30)                               &
     &  ,SF(201,30),SF(202,30),SF(203,30)                               &
     &  ,SF(204,30),SF(205,30),SF(206,30),SF(207,30),SF(208,30)         &
     &  ,SF(301,30),SF(302,30),SF(303,30)                               &
     &  ,SF(401,30),SF(402,30),SF(403,30),SF(404,30),SF(405,30)         &
     &  ,SF(406,30),SF(407,30),SF(408,30),SF(409,30),SF(410,30)         &
     &  ,SF(411,30),SF(412,30),SF(413,30),SF(414,30),SF(415,30)         &
     &  ,SF(416,30),SF(417,30),SF(418,30),SF(419,30),SF(420,30)         &
     &  ,SF(421,30),SF(422,30),SF(423,30),SF(424,30),SF(425,30)         &
     &  ,SF(426,30),SF(427,30),SF(428,30),SF(429,30),SF(430,30)         &
     &  ,SF(431,30),SF(432,30),SF(433,30),SF(434,30),SF(435,30)         &
     &  ,SF(436,30),SF(437,30),SF(438,30),SF(439,30)                    &
     &  ,SF(440,30),SF(441,30)                                          &
     &  ,SF(451,30),SF(452,30),SF(453,30),SF(454,30)                    &
! Diagnostics lengths: in
     &  ,u_m_levs,v_m_levs,w_m_levs,t_m_levs,q_m_levs,z_m_levs          &
     &  ,ke_m_levs,u_mm_levs,v_mm_levs,w_mm_levs,t_mm_levs,q_mm_levs    &
     &  ,z_mm_levs,ke_mm_levs,teot_m_levs,wbig_m_levs,wbig2_m_levs      &
     &  ,RH_m_levs,dry_mass_m_levs                                      &
     &  ,t_inc_levs,q_inc_levs,qcl_inc_levs,qcf_inc_levs                &
     &  ,u_inc_levs,v_inc_levs,w_inc_levs,rho_inc_levs                  &
     &  ,t2_inc_levs,q2_inc_levs,qcl2_inc_levs,qcf2_inc_levs            &
     &  ,u2_inc_levs,v2_inc_levs,w2_inc_levs,rho2_inc_levs              &
     &  ,u_p_levs,v_p_levs,w_p_levs, t_p_levs,q_p_levs,rh_p_levs        &
     &  ,z_p_levs,om_p_levs,heavy_p_levs,tv_p_levs,tvom_p_levs          &
     &  ,prod_p_levs                                                    &
                     ! pressure levels required for each product
     &  ,npress_diags,nmodel_diags                                      &
! Tropopause index bounds
     &  ,min_trop_level,max_trop_level                                  &
! Model levels indicies: in
     &  ,u_m_list,v_m_list,w_m_list,t_m_list,q_m_list,z_m_list,ke_m_list&
     &  ,u_mm_list,v_mm_list,w_mm_list,t_mm_list,q_mm_list,z_mm_list    &
     &  ,ke_mm_list,teot_m_list,wbig_m_list,wbig2_m_list,RH_m_list      &
     &  ,dry_mass_m_list                                                &
     &  ,t_inc_list,q_inc_list,qcl_inc_list,qcf_inc_list                &
     &  ,u_inc_list,v_inc_list,w_inc_list,rho_inc_list                  &
     &  ,t2_inc_list,q2_inc_list,qcl2_inc_list,qcf2_inc_list            &
     &  ,u2_inc_list,v2_inc_list,w2_inc_list,rho2_inc_list              &
     &  ,prod_m_list                                                    &
                     ! index of model levels required for product
! Product levels: in
     &  ,prod_ind                                                       &
                  ! index of pressure levels required for product
! Diagnostic arrays: out
     &  ,SI(1,30,im_index) ,Stashwork                                   &
     &  ,Stashwork(pt001),Stashwork(pt002),Stashwork(pt003)             &
     &  ,Stashwork(pt004),Stashwork(pt005),Stashwork(pt006)             &
     &  ,Stashwork(pt007)                                               &
     &  ,Stashwork(pt101),Stashwork(pt102),Stashwork(pt103)             &
     &  ,Stashwork(pt104),Stashwork(pt105),Stashwork(pt106)             &
     &  ,Stashwork(pt107),Stashwork(pt111),Stashwork(pt112)             &
     &  ,Stashwork(pt114),Stashwork(pt113),Stashwork(pt115)             &
     &  ,Stashwork(pt181),Stashwork(pt182),Stashwork(pt183)             &
     &  ,Stashwork(pt184),Stashwork(pt185),Stashwork(pt186)             &
     &  ,Stashwork(pt187),Stashwork(pt188)                              &
     &  ,Stashwork(pt171),Stashwork(pt172),Stashwork(pt173)             &
     &  ,Stashwork(pt174),Stashwork(pt175),Stashwork(pt176)             &
     &  ,Stashwork(pt177),Stashwork(pt178)                              &
     &  ,Stashwork(pt201),Stashwork(pt202),Stashwork(pt203)             &
     &  ,Stashwork(pt204),Stashwork(pt205),Stashwork(pt206)             &
     &  ,Stashwork(pt207),Stashwork(pt208)                              &
     &  ,Stashwork(pt301),Stashwork(pt302),Stashwork(pt303)             &
     &  ,Stashwork(pt401),Stashwork(pt402),Stashwork(pt403)             &
     &  ,Stashwork(pt404),Stashwork(pt405),Stashwork(pt406)             &
     &  ,Stashwork(pt407),Stashwork(pt408),Stashwork(pt409)             &
     &  ,Stashwork(pt410),Stashwork(pt411),Stashwork(pt412)             &
     &  ,Stashwork(pt413),Stashwork(pt414),Stashwork(pt415)             &
     &  ,Stashwork(pt416),Stashwork(pt417),Stashwork(pt418)             &
     &  ,Stashwork(pt419),stashwork(pt420),Stashwork(pt421)             &
     &  ,Stashwork(pt422),stashwork(pt423),Stashwork(pt424)             &
     &  ,Stashwork(pt425),stashwork(pt426),Stashwork(pt427)             &
     &  ,Stashwork(pt428),stashwork(pt429),Stashwork(pt430)             &
     &  ,Stashwork(pt431),stashwork(pt432),Stashwork(pt433)             &
     &  ,Stashwork(pt434),stashwork(pt435),Stashwork(pt436)             &
     &  ,Stashwork(pt437),stashwork(pt438),Stashwork(pt439)             &
     &  ,Stashwork(pt440),Stashwork(pt441)                              &

     &  ,Stashwork(pt451),Stashwork(pt452),Stashwork(pt453)             &
     &  ,Stashwork(pt454)  )


      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('EOT_DIAG',4)
      END IF


      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('STASH   ',3)
      END IF

      IF(ICODE /= 0) THEN
        RETURN
      ENDIF

! DEPENDS ON: stash
      CALL STASH(a_sm,a_im,30,STASHWORK,                                &
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
      END SUBROUTINE ST_DIAG3
#endif
