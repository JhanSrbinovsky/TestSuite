! ----------------------- include file: AMAXSACM -----------------------
! Description: Quick fix replacement of AMAXSIZE reference in ACPARM.
!              Contains sizes superseded at 5.0, but still embedded
!              throughout AC assimilation code. This allows AC
!              routines to compile without having to re-analyse and
!              correct every routine at this stage. This file should
!              only be of transient use as an interim solution.
!
! Current Code Owner: R. Rawlins
!         5.2   30/11/00 remove ROW_LENGTH_MAX and HORIZ_DIM_MAX
!                        already in amaxsize (now called in ac_ctl)
!                                          B Macpherson
!         6.1   31/08/04 Allow up to 100 levels.  R.Barnes
!         6.2   25/11/05 Set p_rows_max for N320.  R.Barnes

      INTEGER,PARAMETER:: P_ROWS_MAX = 481 ! Max number of rows
      INTEGER,PARAMETER:: P_LEVELS_MAX = 100 ! Max no. of total levels
      INTEGER,PARAMETER:: Q_LEVELS_MAX = 100 ! Max no. of wet levels

! AMAXSACP end
