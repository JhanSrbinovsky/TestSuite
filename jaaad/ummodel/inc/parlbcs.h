! Include file : parlbcs.h
!
! Must be called after parvars.h
!
! Description:
!   Contains variables in connection with generating LBCs.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.2   13/11/00  Original code. D.Robinson
!   5.5   03/02/03  Include qcf2,qrain,qgraup lbc stashcodes. R.Forbes
!   6.0   30/07/03  Include pc2 lbc stashcodes. Damian Wilson
!   6.2   01/10/05  Include murk aerosol lbc stashcodes.  R.M.Forbes
!
! -----------------------------------------------------------
! Stash Codes for LBCs
!
      Integer, Parameter :: lbc_stashcode_orog    = 32001
      Integer, Parameter :: lbc_stashcode_u       = 32002
      Integer, Parameter :: lbc_stashcode_v       = 32003
      Integer, Parameter :: lbc_stashcode_w       = 32004
      Integer, Parameter :: lbc_stashcode_density = 32005
      Integer, Parameter :: lbc_stashcode_theta   = 32006
      Integer, Parameter :: lbc_stashcode_q       = 32007
      Integer, Parameter :: lbc_stashcode_qcl     = 32008
      Integer, Parameter :: lbc_stashcode_qcf     = 32009
      Integer, Parameter :: lbc_stashcode_exner   = 32010
      Integer, Parameter :: lbc_stashcode_u_adv   = 32011
      Integer, Parameter :: lbc_stashcode_v_adv   = 32012
      Integer, Parameter :: lbc_stashcode_w_adv   = 32013
      Integer, Parameter :: lbc_stashcode_qcf2    = 32014
      Integer, Parameter :: lbc_stashcode_qrain   = 32015
      Integer, Parameter :: lbc_stashcode_qgraup  = 32016
      Integer, Parameter :: lbc_stashcode_cf_bulk = 32017
      Integer, Parameter :: lbc_stashcode_cf_liquid=32018
      Integer, Parameter :: lbc_stashcode_cf_frozen=32019
      Integer, Parameter :: lbc_stashcode_murk    = 32020

! -----------------------------------------------------------
!     Data Time for LBC data
      Integer :: LBC_DT_Year
      Integer :: LBC_DT_Month
      Integer :: LBC_DT_Day
      Integer :: LBC_DT_Hour
      Integer :: LBC_DT_Min
      Integer :: LBC_DT_DayNo

      COMMON /LBC_DT/ LBC_DT_Year, LBC_DT_Month, LBC_DT_Day,            &
     &                LBC_DT_Hour, LBC_DT_Min,   LBC_DT_DayNo

! -----------------------------------------------------------

!     Validity Time for LBC data
      Integer :: LBC_VT_Year
      Integer :: LBC_VT_Month
      Integer :: LBC_VT_Day
      Integer :: LBC_VT_Hour
      Integer :: LBC_VT_Min
      Integer :: LBC_VT_DayNo

      COMMON /LBC_VT/ LBC_VT_Year, LBC_VT_Month, LBC_VT_Day,            &
     &                LBC_VT_Hour, LBC_VT_Min,   LBC_VT_DayNo

! -----------------------------------------------------------

      Integer, Parameter :: P_Src_Grid = 2
      Integer, Parameter :: P_LBC_Grid = 4

!     1 : Start Latitude
!     2 : Start Longitude
!     3 : Row Length
!     4 : Rows

      Real :: Src_Grid (Nfld_max, P_Src_Grid)
      Real :: LBC_Grid (Nfld_max, P_LBC_Grid)

      COMMON /LBC_Grids/ Src_Grid, LBC_Grid

! -------------------------------------------------------------

      Integer :: LBC_Global_LenRimA (Nfld_max, Nhalo_max)
      Integer :: LBC_Interp_LenRimA (Nfld_max, Nhalo_max)

      COMMON /LBC_Sizes/ LBC_Global_LenRimA, LBC_Interp_LenRimA

! -------------------------------------------------------------
