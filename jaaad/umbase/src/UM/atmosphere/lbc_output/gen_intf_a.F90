#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Generates Atmosphere Lateral Boundary Conditions
!
! Subroutine Interface:

      SUBROUTINE GEN_INTF_A (                                           &
#include "argduma.h"
#include "arginfa.h"
#include "argptra.h"
#include "argsts.h"
#include "argppx.h"
     &  JINTF,NFTOUT,                                                   &
     &  d1_array,                                                       &
     &  INTERNAL_MODEL                                                  &
     & )

      IMPLICIT NONE

!
! Description:
!   <Say what this routine does>
!
! Method:
!   <Say how it does it: include references to external documentation>
!   <If this routine is very complex, then include a "pseudo code"
!    description of it to make its structure and method clear>
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!   5.3    07/06/01   Remove duplicate declaration.  A van der Wal
!   5.3    22/10/01   Remove Intf_RimW_Orog. Dave Robinson
!   5.3    16/10/01   New logical, Source_Cyclic. D Robinson.
!   5.5    03/02/03   Include code for qcf2,qrain,qgraup.  R.M.Forbes
!   5.5    06/08/02   New logical, Source_Rotated. Copy model fields
!                     into local arrays before calling make_lbcs. If
!                     source grid is rotated, unrotate model winds
!                     before calling make_lbcs. D.Robinson.
!   6.0    29/07/03   Include cloud fraction lbcs. Damian Wilson
!   6.1    01/09/04   Pass ltimer to lbc_writflds.  R.Barnes
!   6.1    18/08/04   Remove repeated question mark.   P.Dando
!   6.2    15/08/05   Free format fixes. P.Selwood.
!   6.2    01/10/04   Include code for murk aerosol lbcs.  R.M.Forbes
!   6.2    06/03/06   Re-calculate weights and indexes for horizontal
!                     interpolation if required. D.Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

#include "cmaxsize.h"
#include "cintfa.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typinfa.h"
#include "typptra.h"
#include "csubmodl.h"
#include "typsts.h"

      Integer  ::  jintf           !  Index to interface area
      Integer  ::  nftout          !  Unit number for interface area
      Integer  ::  internal_model  !  Internal Model Number

      Real     ::  d1_array(*)     !  D1 array with prognostic data

      CHARACTER*(80) CMESSAGE ! Error message if ICODE>0

      Integer :: i,j,var   ! Loop indices
      Integer :: lbc_num
      Integer :: lookup_start
      Integer :: len_io
      Integer :: ntime
      Integer :: len_ppname
      Integer :: im_ident      !   Internal model identifier
      Integer :: len_data
      Integer :: len_data_p
      Integer :: len_data_u
      Integer :: len_data_v
      Integer :: len_data_max

      Real A_IO


#include "chsunits.h"
#include "chistory.h"
#include "ccontrol.h"
#include "clookadd.h"
#include "ctracera.h"
#include "ppxlook.h"

!*L External subroutines called :
      EXTERNAL SETPOS, BUFFOUT

!*---------------------------------------------------------------------
       integer ifld, ihalo, iside
       integer intf_halosize(2,Nhalo_max)
       integer lbc_row_len, lbc_nrows
       integer len_data_sect
       Integer var1
       Integer rimwidth

       Integer Item_Prog(intf_lookupsa)
       Integer IPt_D1(intf_lookupsa)
       Integer lbc_halo_type(intf_lookupsa)
       Integer lbc_fld_type (intf_lookupsa)
       Integer lbc_rim_type (intf_lookupsa)
       Integer lbc_level_type(intf_lookupsa)
       Integer lbc_first_level(intf_lookupsa)
       Integer lbc_last_level(intf_lookupsa)
       Integer lbc_levels   (intf_lookupsa)

       Integer lbc_src_halo_type (intf_lookupsa)
       Integer lbc_src_fld_type  (intf_lookupsa)
       Integer lbc_src_levels    (intf_lookupsa)
       Integer lbc_src_level_type(intf_lookupsa)
       Integer lbc_src_first_level(intf_lookupsa)
       Integer lbc_src_last_level(intf_lookupsa)

       Integer lbc_fld_type_interp
       Integer   ::  ErrorStatus
       Character (Len=*), Parameter :: RoutineName= 'Gen_Intf_A'
       Real, Dimension (:,:),    Allocatable :: LBC_Data
       Real, Dimension (:),      Allocatable :: LBC_Data_u
       Real, Dimension (:),      Allocatable :: LBC_Data_v
       Real, Allocatable :: LBC_Coeff1 (:)
       Real, Allocatable :: LBC_Coeff2 (:)
       Real, Allocatable :: Lambda_p_in(:)
       Real, Allocatable :: Phi_p_in(:)
       Real, Allocatable :: Lambda_u_in(:)
       Real, Allocatable :: Phi_v_in(:)

       Logical u_or_v
!      -------------------------
       Integer local_row_length
       Integer local_rows
       Integer local_row_length_v
       Integer local_rows_v
       Integer source_halo_x
       Integer source_halo_y
       Integer source_fld_type
       Integer source_halo_type
       Integer source_levels

       Real    source_delta_lat
       Real    source_delta_long
       Real    source_first_lat
       Real    source_first_long
       Real    source_pole_lat
       Real    source_pole_long
! -----------------------------------------
       Integer orog_global_row_length
       Integer orog_global_rows
       Integer orog_local_row_length
       Integer orog_local_rows
       Integer orog_halo_x
       Integer orog_halo_y
       Integer orog_fld_type
       Integer orog_halo_type

       Real    orog_first_lat
       Real    orog_first_long
! -----------------------------------------
       Integer lbc_halo_x
       Integer lbc_halo_y

       Real    lbc_delta_lat
       Real    lbc_delta_long
       Real    lbc_first_lat
       Real    lbc_first_long
       Real    lbc_pole_lat
       Real    lbc_pole_long

       Real , allocatable :: source_field   (:,:,:)
       Real , allocatable :: source_field_v (:,:,:)
       Integer :: lev, jrow
       Integer :: field_size, field_size_v

! Indexes for bottom left & right points in bi-linear interpolation

       Integer, dimension (:), allocatable :: lbc_index_bl
       Integer, dimension (:), allocatable :: lbc_index_br

! Weights for 4 points in bi-linear interpolation

       Real,    dimension (:), allocatable :: lbc_weights_tr
       Real,    dimension (:), allocatable :: lbc_weights_br
       Real,    dimension (:), allocatable :: lbc_weights_tl
       Real,    dimension (:), allocatable :: lbc_weights_bl

       Integer lbc_size, lbc_size_u, lbc_size_v, lbc_size_interp
       Integer max_levels_per_pe
       Integer gather_pe
       Integer i_uv
       Integer halo_x,halo_y,lbc_rows
       Integer row,pt,ipt,jvar
       Logical l_lbc_u, l_lbc_v, l_lbc_winds
       Logical :: l_calc_lbc_wts
#include "parlbcs.h"
       Integer :: lbc_rim_size(Nrima_max)
       logical l_vi
       integer n_segs
       integer max_seg_size
       Integer :: prev_src_fld_type
       Integer :: prev_lbc_fld_type
       Integer :: prev_src_halo_type
       Integer :: prev_lbc_halo_type
       integer :: pt_lbc 
       integer :: LBCrow_len 
       integer :: LBCrows      
       Real    :: dlam_wk
       Real    :: dphi_wk 

       Logical :: Source_Cyclic   ! T : Source Grid is cyclic
       Logical :: Source_Rotated  ! T : Source Grid is rotated

!      ------------------------
      CHARACTER*80 STRING         ! work array
      CHARACTER*14 PPNAME         ! boundary output filename

!*---------------------------------------------------------------------
!     Stash item numbers for interface fields
!     Any change to code generating and testing ITEM_INTFA should also
!     consider the corresponding use of ITEM_BOUNDA in INBOUND1/CHKLKBA1
      INTEGER ITEM_INTFA (INTF_LOOKUPSA)
!*---------------------------------------------------------------------
!*
!L Internal structure:

       ErrorStatus = 0
       CMESSAGE=' '

       Source_Cyclic = (A_FIXHD(4) < 3)  !  If Global or NH or SH
       Source_Rotated= (A_REALHD(5)/= 90.0 .or. A_REALHD(6) /= 0.0)

! Set Data Time

       LBC_DT_Year  = A_FIXHD(21)
       LBC_DT_Month = A_FIXHD(22)
       LBC_DT_Day   = A_FIXHD(23)
       LBC_DT_Hour  = A_FIXHD(24)
       LBC_DT_Min   = A_FIXHD(25)
       LBC_DT_DayNo = A_FIXHD(27)

!      =======================================================
       Do Ifld = 1, NFld_Max
         If     (Ifld == fld_type_p) Then   ! p-grid
           Src_Grid (Ifld,1) = a_realhd(rh_baselat)
           Src_Grid (Ifld,2) = a_realhd(rh_baselong)
         Elseif (Ifld == fld_type_u) Then   ! u-grid
           Src_Grid (Ifld,1) = a_realhd(rh_baselat)
           Src_Grid (Ifld,2) = a_realhd(rh_baselong) +                  &
     &                         0.5 * a_realhd(rh_deltaEW)
         Elseif (Ifld == fld_type_v) Then   ! v-grid
           Src_Grid (Ifld,1) = a_realhd(rh_baselat) +                   &
     &                         0.5 * a_realhd(rh_deltaNS)
           Src_Grid (Ifld,2) = a_realhd(rh_baselong)
         Endif
       Enddo

!      =======================================================

       Do Ifld = 1, NFld_Max
         If     (Ifld == fld_type_p) Then   ! p-grid
           LBC_Grid (Ifld,1) = Intf_FirstLat(jintf)
           LBC_Grid (Ifld,2) = Intf_FirstLong(jintf)
           LBC_Grid (Ifld,3) = Intf_Row_Length(jintf)
           LBC_Grid (Ifld,4) = Intf_P_Rows(jintf)
         Elseif (Ifld == fld_type_u) Then   ! u-grid
           LBC_Grid (Ifld,1) = Intf_FirstLat(jintf)
           LBC_Grid (Ifld,2) = Intf_FirstLong(jintf) +                  &
     &                         0.5 * Intf_ewspace(jintf)
           LBC_Grid (Ifld,3) = Intf_Row_Length(jintf) - 1
           LBC_Grid (Ifld,4) = Intf_P_Rows(jintf)
         Elseif (Ifld == fld_type_v) Then   ! v-grid
           LBC_Grid (Ifld,1) = Intf_FirstLat(jintf) +                   &
     &                         0.5 * Intf_nsspace(jintf)
           LBC_Grid (Ifld,2) = Intf_FirstLong(jintf)
           LBC_Grid (Ifld,3) = Intf_Row_Length(jintf)
           LBC_Grid (Ifld,4) = Intf_P_Rows(jintf) - 1
         Endif
       Enddo

!      =======================================================

!      Set up intf_halosize for this area
       intf_halosize(1,1)=1
       intf_halosize(2,1)=1

       intf_halosize(1,2)=intf_exthalo_ew(jintf)
       intf_halosize(2,2)=intf_exthalo_ns(jintf)

       intf_halosize(1,3)=0
       intf_halosize(2,3)=0
!      =======================================================      

       LBCrow_len = intf_row_length(jintf)
       LBCrows    = intf_p_rows(jintf)
       lbc_halo_x = intf_exthalo_ew(jintf)
       lbc_halo_y = intf_exthalo_ns(jintf)
       
       allocate (Lambda_p_in(1-lbc_halo_x: LBCrow_len + lbc_halo_x)) 
       allocate (Phi_p_in(1-lbc_halo_y: LBCrows + lbc_halo_y))
       allocate (Lambda_u_in(1-lbc_halo_x: LBCrow_len + lbc_halo_x)) 
       allocate (Phi_v_in(1-lbc_halo_y: LBCrows + lbc_halo_y))
     
       do pt_lbc = 1, LBCrow_len
         Lambda_u_in(pt_lbc) = Lambda_intf_u(pt_lbc,jintf)
         Lambda_p_in(pt_lbc) = Lambda_intf_p(pt_lbc,jintf)
       end do
       do pt_lbc = 1, LBCrows
         Phi_p_in(pt_lbc) = Phi_intf_p(pt_lbc,jintf)
         Phi_v_in(pt_lbc) = Phi_intf_v(pt_lbc,jintf)
       end do 
             
! compute lambda and phi intervals dlam_wk and dphi_wk in the 
! halo region

       If ( intf_l_var_lbc(jintf) ) Then
         dlam_wk = Lambda_p_in(LBCrow_len)-Lambda_p_in(LBCrow_len-1) 
         dphi_wk = Phi_p_in(LBCrows)-Phi_p_in(LBCrows-1)
       Else       !regular grid
         dlam_wk =  Intf_ewspace(jintf)
         dphi_wk =  Intf_nsspace(jintf) 
       End if

! compute lambda and phi in the halo region
                       
       do pt_lbc = 0, 1 - lbc_halo_x,  - 1
         Lambda_p_in(pt_lbc) = Lambda_p_in(1) - (1 - pt_lbc)*dlam_wk                
         Lambda_u_in(pt_lbc) = Lambda_u_in(1) - (1 - pt_lbc)*dlam_wk
       end do 
         
       do pt_lbc = 0, 1 - lbc_halo_y,  - 1
         Phi_p_in(pt_lbc) = Phi_p_in(1) - (1 - pt_lbc)*dphi_wk
         Phi_v_in(pt_lbc) = Phi_v_in(1) - (1 - pt_lbc)*dphi_wk
       end do 
          
       do pt_lbc = LBCrow_len + 1, LBCrow_len + lbc_halo_x  
         Lambda_p_in(pt_lbc) = Lambda_p_in(pt_lbc-1) + dlam_wk
         Lambda_u_in(pt_lbc) = Lambda_u_in(pt_lbc-1) + dlam_wk
       end do
    
       do pt_lbc =  LBCrows + 1,  LBCrows + lbc_halo_y  
         Phi_p_in(pt_lbc) = Phi_p_in(pt_lbc-1) + dphi_wk
         Phi_v_in(pt_lbc) = Phi_v_in(pt_lbc-1) + dphi_wk 
       end do
! ==============================================================

! DEPENDS ON: lbc_grid_sizes
      Call LBC_Grid_Sizes (jintf)

! ==============================================================

       im_ident = internal_model

! Logical to indicate if model grid is rotated

!      Stash codes of variables to be put in LBC file.
       ITEM_INTFA( 1) = lbc_stashcode_orog
       ITEM_INTFA( 2) = lbc_stashcode_u
       ITEM_INTFA( 3) = lbc_stashcode_v
       ITEM_INTFA( 4) = lbc_stashcode_w
       ITEM_INTFA( 5) = lbc_stashcode_density
       ITEM_INTFA( 6) = lbc_stashcode_theta
       ITEM_INTFA( 7) = lbc_stashcode_q
       ITEM_INTFA( 8) = lbc_stashcode_qcl
       ITEM_INTFA( 9) = lbc_stashcode_qcf
       ITEM_INTFA(10) = lbc_stashcode_exner
       ITEM_INTFA(11) = lbc_stashcode_u_adv
       ITEM_INTFA(12) = lbc_stashcode_v_adv
       ITEM_INTFA(13) = lbc_stashcode_w_adv
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_qcf2
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_qrain
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_qgraup
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       If (L_pc2) Then
         lbc_num = lbc_num+1
         ITEM_INTFA(lbc_num) = lbc_stashcode_cf_bulk
         lbc_num = lbc_num+1
         ITEM_INTFA(lbc_num) = lbc_stashcode_cf_liquid
         lbc_num = lbc_num+1
         ITEM_INTFA(lbc_num) = lbc_stashcode_cf_frozen
       EndIf
       ! Setup for murk aerosol lbcs if active
       If (L_murk) Then
         lbc_num = lbc_num + 1
         ITEM_INTFA(lbc_num) = lbc_stashcode_murk
       EndIf

!      Stash Codes for prognostics corresponding to above variables.
       ITEM_PROG( 1) = 33
       ITEM_PROG( 2) = 2
       ITEM_PROG( 3) = 3
       ITEM_PROG( 4) = 150
       ITEM_PROG( 5) = 253
       ITEM_PROG( 6) = 4
       ITEM_PROG( 7) = 10
       ITEM_PROG( 8) = 254
       ITEM_PROG( 9) = 12
       ITEM_PROG(10) = 255
       ITEM_PROG(11) = 256
       ITEM_PROG(12) = 257
       ITEM_PROG(13) = 258
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         ITEM_PROG(lbc_num) = 271
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         ITEM_PROG(lbc_num) = 272
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         ITEM_PROG(lbc_num) = 273
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       If (L_pc2) Then
         lbc_num = lbc_num+1
         ITEM_PROG(lbc_num) = 266
         lbc_num = lbc_num+1
         ITEM_PROG(lbc_num) = 267
         lbc_num = lbc_num+1
         ITEM_PROG(lbc_num) = 268
       EndIf
       ! Setup for murk aerosol lbcs if active
       If (L_murk) Then
         lbc_num = lbc_num + 1
         ITEM_PROG(lbc_num) = 90
       EndIf

!      Addresses in D1 pointing to above prognostics.
       IPT_D1( 1) = JOROG
       IPT_D1( 2) = JU(1)
       IPT_D1( 3) = JV(1)
       IPT_D1( 4) = JW(0)
       IPT_D1( 5) = JRHO(1)
       IPT_D1( 6) = JTHETA(1)
       IPT_D1( 7) = JQ(1)
       IPT_D1( 8) = JQCL(1)
       IPT_D1( 9) = JQCF(1)
       IPT_D1(10) = JEXNER_RHO_LEVELS(1)
       IPT_D1(11) = JU_ADV(1)
       IPT_D1(12) = JV_ADV(1)
       IPT_D1(13) = JW_ADV(0)
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JQCF2(1)
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JQRAIN(1)
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JQGRAUP(1)
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       If (L_pc2) Then
         lbc_num = lbc_num+1
         IPT_D1(lbc_num) = JCF_BULK(1)
         lbc_num = lbc_num+1
         IPT_D1(lbc_num) = JCF_LIQUID(1)
         lbc_num = lbc_num+1
         IPT_D1(lbc_num) = JCF_FROZEN(1)
       EndIf
       ! Setup for murk aerosol lbcs if active
       If (L_murk) Then
         lbc_num = lbc_num + 1
         IPT_D1(lbc_num) = JMURK(1)
       EndIf

       lbc_rim_size(rima_type_norm) = IntfWidthA(jintf)
       lbc_rim_size(rima_type_orog) = IntfWidthA(jintf)

!L 2.0 Update Information in Headers

!L     Open boundary output file if reinitialised during run

      IF (FT_STEPS(NFTOUT) >  0) THEN
        STRING=MODEL_FT_UNIT(NFTOUT)
        PPNAME=STRING(18:31)
        LEN_PPNAME=LEN(PPNAME)
! DEPENDS ON: file_open
        CALL FILE_OPEN(NFTOUT,PPNAME,LEN_PPNAME,1,1,ErrorStatus)
        IF (ErrorStatus /= 0) THEN
          Write (CMessage,*) 'Error opening preassigned boundary file'
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        ENDIF

!      Determine position where to Buffer out data to

       NTIME=FT_LASTFIELD(NFTOUT)+1
      ELSE
       NTIME=FT_LASTFIELD(NFTOUT)+1

      ENDIF

       if (ntime == 1) Then
         Var1 = 1  ! Include orography
       else
         Var1 = 2  ! Do not include orography
       endif

!L 2.1 Fixed length header
!     FIXHD_INTFA(152,JINTF) =  INTF_LOOKUPSA * NTIME
      FIXHD_INTFA(152,JINTF) =  INTF_LOOKUPSA +                         &
     &                         (INTF_LOOKUPSA-VAR1+1)*(NTIME-1)

!     FIXHD_INTFA(161,JINTF) = LEN_INTF_DATA*NTIME

!L 2.2 Integer Constants
      INTHD_INTFA(3,JINTF) = NTIME

!L 2.3 LOOKUP Table

!     Determine position in LOOKUP table
      LOOKUP_START=FIXHD_INTFA(150,JINTF) +                             &
     &             FIXHD_INTFA(151,JINTF) *                             &
     &            (FIXHD_INTFA(152,JINTF) - (INTF_LOOKUPSA-VAR1+1) ) - 1

!     Check that there is enough space for this entry in LOOKUP table
      IF (FIXHD_INTFA(150,JINTF)+                                       &
     &    FIXHD_INTFA(151,JINTF)*FIXHD_INTFA(152,JINTF) >               &
     &   FIXHD_INTFA(160,JINTF)) THEN
        Write (CMessage,*)                                              &
     &  'Insufficient space for headers in boundary dataset'
        ErrorStatus = 1
! DEPENDS ON: ereport
        Call Ereport ( RoutineName, ErrorStatus, CMessage )
      ENDIF

! =======================================================

! DEPENDS ON: lbc_setup
      Call LBC_SetUp (                                                  &
#include "argsts.h"
#include "argppx.h"
     &     item_intfa,                                                  &
     &     lbc_fld_type,                                                &
     &     lbc_halo_type,                                               &
     &     lbc_rim_type,                                                &
     &     lbc_level_type,                                              &
     &     lbc_first_level,                                             &
     &     lbc_last_level,                                              &
     &     lbc_levels,                                                  &
     &     intf_lookupsa,                                               &
     &     jintf,                                                       &
     &     im_ident                                                     &
     & )

! =============================================================
!      Hardwire level type for now
!      level type for orog lbcs = 5, so ok
!      level type for rest      = 0, needs to be 1 or 2
!      so we know which prognostics are on rho or theta levels
       lbc_level_type (2)  = 1
       lbc_level_type (3)  = 1
       lbc_level_type (4)  = 2
       lbc_level_type (5)  = 1
       lbc_level_type (6)  = 2
       lbc_level_type (7)  = 2
       lbc_level_type (8)  = 2
       lbc_level_type (9)  = 2
       lbc_level_type (10) = 1
       lbc_level_type (11) = 1
       lbc_level_type (12) = 1
       lbc_level_type (13) = 2
       ! Setup for additional microphysics lbcs if active
       lbc_num = 13
       IF (L_mcr_qcf2) THEN  ! qcf2 lbcs active
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       ENDIF
       IF (L_mcr_qrain) THEN  ! qrain lbcs active
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       ENDIF
       IF (L_mcr_qgraup) THEN  ! qgraup lbcs active
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       ENDIF
       ! Setup for additional cloud fraction lbcs if active
       If (L_pc2) Then
         lbc_num = lbc_num+1
         lbc_level_type(lbc_num) = 2
         lbc_num = lbc_num+1
         lbc_level_type(lbc_num) = 2
         lbc_num = lbc_num+1
         lbc_level_type(lbc_num) = 2
       EndIf
       ! Setup for murk aerosol lbcs if active
       If (L_murk) Then
         lbc_num = lbc_num + 1
         lbc_level_type (lbc_num) = 2
       EndIf

! ===============================================================

! DEPENDS ON: lbc_src_setup
      CALL LBC_Src_Setup (                                              &
#include "argsts.h"
#include "argppx.h"
     &     item_prog,                                                   &
     &     lbc_src_fld_type,                                            &
     &     lbc_src_halo_type,                                           &
     &     lbc_src_level_type,                                          &
     &     lbc_src_first_level,                                         &
     &     lbc_src_last_level,                                          &
     &     lbc_src_levels,                                              &
     &     intf_lookupsa,                                               &
     &     im_ident                                                     &
     & )

!  ==================================================

! DEPENDS ON: lbc_validity_time
      Call LBC_Validity_Time (LCAL360)

!  ==================================================

! DEPENDS ON: lbc_setup_lookup
      Call LBC_Setup_LookUp (                                           &
     &     lookup_intfa(1,1,jintf)                                      &
     &,    Len1_lookup                                                  &
     &,    intf_lookupsa                                                &
     &,    Item_intfa                                                   &
     &,    jintf                                                        &
     &,    ntime                                                        &
     &,    fixhd_intfa(1,jintf),                                        &
#include "argppx.h"
     &     lbc_rim_size                                                 &
     &,    lbc_rim_type                                                 &
     &,    lbc_halo_type                                                &
     &,    lbc_fld_type                                                 &
     &,    intf_halosize                                                &
     &,    lbc_levels                                                   &
     & )

!  ==================================================

!L 3.0 Main loop to generate LBCs

       prev_src_fld_type  = -1
       prev_lbc_fld_type  = -1
       prev_src_halo_type = -1
       prev_lbc_halo_type = -1

       lbc_fld_type_interp = 1

       do j=var1,intf_lookupsa

!      Test if weights need to be calculated

         l_calc_lbc_wts =                                               &
     &     ( lbc_src_fld_type(j)  /= prev_src_fld_type ) .or.           &
     &     ( lbc_fld_type(j)      /= prev_lbc_fld_type ) .or.           &
!    &     ( lbc_src_halo_type(j) /= prev_src_halo_type) .or.
     &     ( lbc_halo_type(j)     /= prev_lbc_halo_type)

         l_lbc_u     = lbc_fld_type(j) == fld_type_u
         l_lbc_v     = lbc_fld_type(j) == fld_type_v
         l_lbc_winds = l_lbc_u .or. l_lbc_v

         len_data      = lookup_intfa(lblrec,j,jintf)
         len_data_sect = lookup_intfa(lbnrec,j,jintf)

         lbc_size = len_data / lbc_levels(j)

         lbc_size_interp =                                              &
     &   lbc_interp_lenrima(lbc_fld_type(j),lbc_halo_type(j))

         if (l_lbc_u) then

! allocate actual size for u and v in lbc_data

           allocate ( lbc_data (lbc_size_interp*lbc_levels(j),2) )
           jvar = 1
           lbc_size_u = lbc_size
         elseif (l_lbc_v) then
!          space allocated with u-component
           jvar = 2
           lbc_size_v = lbc_size
         else

! need to cater for rounding to sector boundary

           len_data_max  = max ( len_data, len_data_sect )

!          allocate ( lbc_data  (lbc_size_interp*lbc_levels(j),1) )
           allocate ( lbc_data (len_data_max,1) )
           jvar = 1
         endif

         lbc_data(:,jvar) = 0.0

         source_fld_type = lbc_src_fld_type (j)

         global_row_length = glsize(1,source_fld_type)
         global_rows       = glsize(2,source_fld_type)

         local_row_length  = blsize(1,source_fld_type)
         local_rows        = blsize(2,source_fld_type)

         local_row_length_v = blsize(1,fld_type_v)
         local_rows_v       = blsize(2,fld_type_v)

         source_halo_type = lbc_src_halo_type(j)
         source_halo_x    = halosize(1,source_halo_type)
         source_halo_y    = halosize(2,source_halo_type)
         source_levels    = lbc_src_levels(j)

         allocate ( lbc_coeff1(lbc_size_interp) )
         allocate ( lbc_coeff2(lbc_size_interp) )

         max_levels_per_pe = ((source_levels-1)/nproc)+1

         gather_pe = 0

         l_vi = (intf_vert_interp(jintf)  .and. j /= 1)

         if (l_lbc_winds) then
           i_uv = 1
         else
           i_uv = 0
         endif

         source_delta_lat  = a_realhd(rh_deltaNS)
         source_delta_long = a_realhd(rh_deltaEW)
         source_first_lat  = Src_Grid(Source_fld_type,1)
         source_first_long = Src_Grid(Source_fld_type,2)
         source_pole_lat   = a_realhd(rh_rotlat)
         source_pole_long  = a_realhd(rh_rotlong)

         lbc_row_len    = LBC_Grid(lbc_fld_type_interp,3)
         lbc_rows       = LBC_Grid(lbc_fld_type_interp,4)
         lbc_first_lat  = LBC_Grid(lbc_fld_type_interp,1)
         lbc_first_long = LBC_Grid(lbc_fld_type_interp,2)

         lbc_delta_lat  = Intf_nsspace(jintf)
         lbc_delta_long = Intf_ewspace(jintf)
         lbc_pole_lat   = Intf_polelat(jintf)
         lbc_pole_long  = Intf_polelong(jintf)
         lbc_halo_x     = Intf_halosize(1,lbc_halo_type(j))
         lbc_halo_y     = Intf_halosize(2,lbc_halo_type(j))

         rimwidth = lbc_rim_size(lbc_rim_type(j))

! Orography field
         orog_fld_type          = lbc_src_fld_type(1)
         orog_local_row_length  = blsize(1,orog_fld_type)
         orog_local_rows        = blsize(2,orog_fld_type)
         orog_global_row_length = glsize(1,orog_fld_type)
         orog_global_rows       = glsize(2,orog_fld_type)
         orog_first_lat         = src_grid(orog_fld_type,1)
         orog_first_long        = src_grid(orog_fld_type,2)

         orog_halo_type = lbc_src_halo_type(1)
         orog_halo_x    = halosize(1,orog_halo_type)
         orog_halo_y    = halosize(2,orog_halo_type)
! vi
         If (l_vi) Then
           max_seg_size = ((lbc_size_interp-1)/nproc) + 1
           n_segs       = ((lbc_size_interp-1)/max_seg_size) + 1
         Else
           max_seg_size = 1
           n_segs       = 1
         End If

         If (l_calc_lbc_wts) Then

           If ( allocated(lbc_index_bl)   ) deallocate (lbc_index_bl)
           If ( allocated(lbc_index_br)   ) deallocate (lbc_index_br)
           If ( allocated(lbc_weights_tr) ) deallocate (lbc_weights_tr)
           If ( allocated(lbc_weights_br) ) deallocate (lbc_weights_br)
           If ( allocated(lbc_weights_bl) ) deallocate (lbc_weights_bl)
           If ( allocated(lbc_weights_tl) ) deallocate (lbc_weights_tl)

           allocate ( lbc_index_bl  (lbc_size_interp) )
           allocate ( lbc_index_br  (lbc_size_interp) )
           allocate ( lbc_weights_tr(lbc_size_interp) )
           allocate ( lbc_weights_br(lbc_size_interp) )
           allocate ( lbc_weights_bl(lbc_size_interp) )
           allocate ( lbc_weights_tl(lbc_size_interp) )

          End If

! Copy prognostic into local space.
! Only required if source grid is LAM as model winds need to be
! un-rotated before generating the LBCs.
! Note that copy is done for all variables ; makes it easier to
! generalise code

         Allocate (source_field                                         &
     &            (1-source_halo_x:local_row_length+source_halo_x,      &
     &             1-source_halo_y:local_rows+source_halo_y,            &
     &             source_levels) )

         If (.not. l_lbc_v) Then

!          Copy prognostic from d1_array into local space - for
!          winds, the u-component is copied here.

           ipt = 0
           Do lev = 1, source_levels
             Do jrow = 1-source_halo_y, local_rows+source_halo_y
               Do i = 1-source_halo_x, local_row_length+source_halo_x
                 ipt = ipt+1
                 source_field(i,jrow,lev) = d1_array(ipt_d1(j)+ipt-1)
               End Do
             End Do
           End Do

!          Check length of data copied.

           field_size = lasize(1,source_fld_type,source_halo_type) *    &
     &                  lasize(2,source_fld_type,source_halo_type)

!          Why just u ? Checking out.
           If (l_lbc_u) Then
             If (ipt /= field_size*source_levels) Then
               Write (CMessage,*) ' ipt /= field_size ?'
               ErrorStatus = 1
! DEPENDS ON: ereport
               Call Ereport ( RoutineName, ErrorStatus, CMessage )
             End If
           End If

         End If

         If (l_lbc_u) Then

!          This pass is for a u-component so copy the v-component
!          as well. Store in source_field_v.

           Allocate (source_field_v                                     &
     &              (1-source_halo_x:local_row_length_v+source_halo_x,  &
     &               1-source_halo_y:local_rows_v+source_halo_y,        &
     &               source_levels) )

           ipt = 0
           Do lev = 1, source_levels
             Do jrow = 1-source_halo_y, local_rows_v+source_halo_y
               Do i = 1-source_halo_x, local_row_length_v+source_halo_x
                 ipt = ipt+1
                 source_field_v(i,jrow,lev) =                           &
     &           d1_array(ipt_d1(j+1)+ipt-1)
               End Do
             End Do
           End Do

!          Check length of v field copied.

           field_size_v = lasize(1,fld_type_v,source_halo_type) *       &
     &                    lasize(2,fld_type_v,source_halo_type)

           If (ipt /= field_size_v * source_levels) Then
             write (6,*) ' ipt ',ipt
             write (6,*) ' field_size_v ',field_size_v
             write (6,*) ' field_size_v * src_levs ',                   &
     &                     field_size_v*source_levels
             Write (CMessage,*) ' ipt /= field_size_v ?'
             ErrorStatus = 2
! DEPENDS ON: ereport
             Call Ereport ( RoutineName, ErrorStatus, CMessage )
           End If

         End If

         If (source_rotated .and. l_lbc_u) Then

!          The model grid is a LAM grid. On the pass through for the
!          u-component, unrotate the model winds.

! DEPENDS ON: lbc_unrotate_model_winds
           Call LBC_Unrotate_Model_Winds (                              &
     &          source_field                                            &
     &,         source_field_v                                          &
     &,         field_size                                              &
     &,         field_size_v                                            &
     &,         row_length                                              &
     &,         rows                                                    &
     &,         n_rows                                                  &
     &,         source_levels                                           &
     &,         source_halo_type                                        &
     &,         source_pole_lat                                         &
     &,         source_pole_long                                        &
     &,         src_grid(fld_type_p,1)                                  &
     &,         src_grid(fld_type_p,2)                                  &
     &,         source_delta_lat                                        &
     &,         source_delta_long                                       &
     & )

         End If

         If (l_lbc_v) Then

!          The winds are un-rotated during the pass through for the
!          u-component and the v-component is stored in source_field_v.
!          Before calling make_lbcs for the v-component, copy
!          source_field_v into source_field.

           ipt = 0
           Do lev = 1, source_levels
             Do jrow = 1-source_halo_y, local_rows+source_halo_y
               Do i = 1-source_halo_x, local_row_length+source_halo_x
                 ipt = ipt+1
                 source_field(i,jrow,lev) = source_field_v(i,jrow,lev)
               End Do
             End Do
           End Do

         End If

! DEPENDS ON: make_lbcs
         call make_lbcs (                                               &
! Prognostic
     &        source_field                                              &
     &,       local_row_length                                          &
     &,       local_rows                                                &
     &,       global_row_length                                         &
     &,       global_rows                                               &
     &,       source_halo_x                                             &
     &,       source_halo_y                                             &
     &,       source_levels                                             &
     &,       source_fld_type                                           &
     &,       source_halo_type                                          &
     &,       source_delta_lat                                          &
     &,       source_delta_long                                         &
     &,       source_first_lat                                          &
     &,       source_first_long                                         &
     &,       source_pole_lat                                           &
     &,       source_pole_long                                          &
     &,       source_cyclic                                             &
     &,       source_rotated                                            &
! Orography Field
     &,       d1_array(ipt_d1(1))                                       &
     &,       orog_local_row_length                                     &
     &,       orog_local_rows                                           &
     &,       orog_global_row_length                                    &
     &,       orog_global_rows                                          &
     &,       orog_fld_type                                             &
     &,       orog_halo_type                                            &
     &,       orog_halo_x                                               &
     &,       orog_halo_y                                               &
     &,       orog_first_lat                                            &
     &,       orog_first_long                                           &
! LBC field
     &,       lbc_row_len                                               &
     &,       lbc_rows                                                  &
     &,       lbc_levels(j)                                             &
     &,       lbc_delta_lat                                             &
     &,       lbc_delta_long                                            &
     &,       lbc_first_lat                                             &
     &,       lbc_first_long                                            &
     &,       lbc_pole_lat                                              &
     &,       lbc_pole_long                                             &
     &,       lbc_halo_x                                                &
     &,       lbc_halo_y                                                &
     &,       rimwidth                                                  &
     &,       lbc_data(1,jvar)                                          &
     &,       lbc_coeff1                                                &
     &,       lbc_coeff2                                                &
     &,       lbc_size_interp                                           &
     &,       lbc_src_level_type(j)                                     &
     &,       intf_l_var_lbc(jintf)                                     &
     &,       Lambda_p_in                                               &
     &,       Phi_p_in                                                  &
! hi - indexes and weights
     &,       l_calc_lbc_wts                                            &
     &,       lbc_index_bl  , lbc_index_br                              &
     &,       lbc_weights_tr, lbc_weights_br                            &
     &,       lbc_weights_bl, lbc_weights_tl                            &
! vi
     &,       max_seg_size                                              &
     &,       n_segs                                                    &
     &,       max_levels_per_pe                                         &
     &,       gather_pe                                                 &
     &,       l_vi                                                      &
! vi - src
     &,       model_levels                                              &
     &,       a_inthd(ih_height_gen)                                    &
     &,       a_inthd(ih_1_c_rho_level)                                 &
     &,       a_realhd(rh_z_top_theta)                                  &
     &,       a_levdepc(jetatheta)                                      &
     &,       a_levdepc(jetarho)                                        &
     &,       lbc_src_first_level(j)                                    &
     &,       lbc_src_last_level(j)                                     &
! vi - lbc
     &,       intf_p_levels(jintf)                                      &
     &,       inthd_intfa(ih_height_gen,jintf)                          &
     &,       lbc_first_r_rho(jintf)                                    &
     &,       lbc_z_top_model(jintf)                                    &
     &,       lbc_eta_theta(1,jintf)                                    &
     &,       lbc_eta_rho(1,jintf)                                      &
     &,       lbc_first_level(j)                                        &
     &,       lbc_last_level(j)                                         &
     &,       i_uv                                                      &
     & )

        deallocate (source_field)
        If (l_lbc_v) Then
          deallocate (source_field_v)
        Endif

        If (.not. l_lbc_v) Then
          deallocate (lbc_coeff1)
          deallocate (lbc_coeff2)
        Endif

      if (l_lbc_v) Then

!     LBCs now available for both u and v

! Input to rotate_winds : u and v are on A grid, ie, u and v on P grid

! DEPENDS ON: lbc_rotate_winds
          call lbc_rotate_winds (                                       &
     &         lbc_data(1,1)                                            &
                                !  u lbcs on p grid
     &,        lbc_data(1,2)                                            &
                                !  v lbcs on p grid
     &,        lbc_size_interp                                          &
                                !  size of u & v lbcs on A grid
     &,        lbc_levels(j)                                            &
                                !  no of lbc levels
     &,        lbc_coeff1                                               &
                                !  \ coefficients to
     &,        lbc_coeff2                                               &
                                !  / rotate winds
     & )

! Output from rotate_winds :
! u and v still on A grid. Winds are now w.r.t the rotated grid.

          deallocate (lbc_coeff1)
          deallocate (lbc_coeff2)

! Interpolate u from A grid to C u-grid

          len_data_u    = lbc_size_u*lbc_levels(j)
          len_data_sect = lookup_intfa(lbnrec,j-1,jintf)
          len_data_max  = max ( len_data_u, len_data_sect )

          allocate ( lbc_data_u(len_data_max) )

! DEPENDS ON: lbc_u_a_to_c
          call lbc_u_a_to_c (                                           &
     &         lbc_data(1,1)                                            &
                                !  u lbcs on p grid
     &,        lbc_data_u                                               &
                                !  u lbcs on u grid
     &,        lbc_size_interp                                          &
                                !  field size on p grid
     &,        lbc_size_u                                               &
                                !  field size on u grid
     &,        lbc_levels(j)                                            &
                                !  no of levels
     &,        intf_row_length(jintf)                                   &
     &,        intf_p_rows(jintf)                                       &
     &,        lbc_rim_size(lbc_rim_type(j))                            &
     &,        intf_halosize(1,2)                                       &
     &,        intf_halosize(2,2)                                       &
     &,        intf_l_var_lbc(jintf)                                    &
     &,        Lambda_p_in                                              &
     &,        Lambda_u_in                                              &
     & )

! DEPENDS ON: lbc_writflds
        call lbc_writflds (nftout,lbc_data_u,len_data_max,              &
     &                     lookup_intfa(1,j-1,jintf),                   &
     &                     fixhd_intfa(1,jintf), ltimer )

        deallocate (lbc_data_u)

! Interpolate v from A grid to C v-grid

          len_data_v    = lbc_size_v*lbc_levels(j)
          len_data_sect = lookup_intfa(lbnrec,j,jintf)
          len_data_max  = max ( len_data_v, len_data_sect )

          allocate ( lbc_data_v(len_data_max) )

! DEPENDS ON: lbc_v_a_to_c
          call lbc_v_a_to_c (                                           &
     &         lbc_data(1,2)                                            &
                                !  v lbcs on p grid
     &,        lbc_data_v                                               &
                                !  v lbcs on v grid
     &,        lbc_size_interp                                          &
                                !  field size on p grid
     &,        lbc_size_v                                               &
                                !  field size on v grid
     &,        lbc_levels(j)                                            &
                                !  no of levels
     &,        intf_row_length(jintf)                                   &
     &,        intf_p_rows(jintf)                                       &
     &,        lbc_rim_size(lbc_rim_type(j+1))                          &
     &,        intf_halosize(1,2)                                       &
     &,        intf_halosize(2,2)                                       &
     &,        intf_l_var_lbc(jintf)                                    &
     &,        Phi_p_in                                                 &
     &,        Phi_v_in                                                 &
     & )

! DEPENDS ON: lbc_writflds
          call lbc_writflds (nftout,lbc_data_v,len_data_max,            &
     &                       lookup_intfa(1,j,jintf),                   &
     &                       fixhd_intfa(1,jintf), ltimer )

          deallocate (lbc_data_v)
          deallocate (lbc_data)

        endif   !  If (l_lbc_v)

        If ( .not. l_lbc_winds ) Then

! DEPENDS ON: lbc_post_interp_transf
          call lbc_post_interp_transf (                                 &
     &         lbc_size,                                                &
     &         lbc_first_level(j),                                      &
     &         lbc_last_level(j),                                       &
     &         item_intfa(j),                                           &
     &         lbc_data )

        End If

        if ( .not. l_lbc_winds ) then

! DEPENDS ON: lbc_writflds
          call lbc_writflds (nftout,lbc_data,len_data_max,              &
     &                       lookup_intfa(1,j,jintf),                   &
     &                       fixhd_intfa(1,jintf), ltimer )

          deallocate (lbc_data)

        endif

        prev_src_fld_type  = lbc_src_fld_type(j)
        prev_lbc_fld_type  = lbc_fld_type(j)
        prev_src_halo_type = lbc_src_halo_type(j)
        prev_lbc_halo_type = lbc_halo_type(j)

       enddo

       If ( allocated(lbc_index_bl)   ) deallocate (lbc_index_bl)
       If ( allocated(lbc_index_br)   ) deallocate (lbc_index_br)
       If ( allocated(lbc_weights_tr) ) deallocate (lbc_weights_tr)
       If ( allocated(lbc_weights_br) ) deallocate (lbc_weights_br)
       If ( allocated(lbc_weights_bl) ) deallocate (lbc_weights_bl)
       If ( allocated(lbc_weights_tl) ) deallocate (lbc_weights_tl)
        
        deallocate (Lambda_p_in)
        deallocate (Phi_p_in)
        deallocate (Lambda_u_in)
        deallocate (Phi_v_in)
        
!L 4.0 Write out headers/data

!L 4.1 Fixed length header

! DEPENDS ON: setpos
        CALL SETPOS (NFTOUT,0,ErrorStatus)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,FIXHD_INTFA(1,JINTF),LEN_FIXHD,LEN_IO,A_IO)

        If (A_IO /= -1.0 .or. LEN_IO /= LEN_FIXHD) THEN

          Write (6,*) ' Return Code from BUFFOUT    ',A_IO
          Write (6,*) ' Length of data transferred  ',LEN_IO
          Write (6,*) ' Expected transferred length ',LEN_FIXHD

          Write (CMessage,*) 'Error in BUFFOUT of fixed length header.'
          ErrorStatus = 2
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        END IF

!L 4.2 Integer constants

! DEPENDS ON: buffout
        CALL BUFFOUT (NFTOUT,INTHD_INTFA(1,JINTF),                      &
     &                PP_LEN_INTHD,LEN_IO,A_IO)

        If (A_IO /= -1.0 .or. LEN_IO /= PP_LEN_INTHD) Then

          Write (6,*) ' Return Code from BUFFOUT    ',A_IO
          Write (6,*) ' Length of data transferred  ',LEN_IO
          Write (6,*) ' Expected transferred length ',PP_LEN_INTHD

          Write (CMessage,*) 'Error in BUFFOUT of integer header.'
          ErrorStatus = 3
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        END IF

!L 4.3 PP headers in LOOKUP table
! DEPENDS ON: setpos
        CALL SETPOS(NFTOUT,LOOKUP_START,ErrorStatus)
! DEPENDS ON: buffout
        CALL BUFFOUT(NFTOUT,LOOKUP_INTFA(1,VAR1,JINTF),                 &
     &               LEN1_LOOKUP*(INTF_LOOKUPSA-VAR1+1),LEN_IO,A_IO)

        IF(A_IO /= -1.0.OR.                                             &
     &     LEN_IO /= LEN1_LOOKUP*(INTF_LOOKUPSA-VAR1+1)) THEN

          Write (6,*) ' Return Code from BUFFOUT    ',A_IO
          Write (6,*) ' Length of data transferred  ',LEN_IO
          Write (6,*) ' Expected transferred length ',                  &
     &    LEN1_LOOKUP*(INTF_LOOKUPSA-VAR1+1)

          Write (CMessage,*) 'Error in BUFFOUT of lookup headers.'
          ErrorStatus = 4
! DEPENDS ON: ereport
          Call Ereport ( RoutineName, ErrorStatus, CMessage )
        END IF

!L     Close boundary output file if reinitialised during run
      IF (FT_STEPS(NFTOUT) >  0) THEN
        LEN_PPNAME=LEN(PPNAME)
! DEPENDS ON: file_close
        CALL FILE_CLOSE(NFTOUT,PPNAME,LEN_PPNAME,1,0,ErrorStatus)
      END IF

!L     Update FT_LASTFIELD
      FT_LASTFIELD(NFTOUT) = FT_LASTFIELD(NFTOUT) + 1

      RETURN
      END SUBROUTINE GEN_INTF_A
#endif
