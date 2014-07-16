#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE WRITEHDR(month,year,daym0,daymas,filename,flheader,    &
     &                    fixhd12,umstepno,z_top_of_model,              &
     &                    first_constant_r_rho_level)
!-----------------------------------------------------------------------
!-
!-   Purpose and Methods : Write out Fixed Length Header to file.
!-   See UM Documentation Paper F3 for details of header
!-
!-   Inputs  :
!-             MONTH,YEAR,      - time
!-             OUTDAY,PERIOD,   - first day and days in month
!-
!-   Outputs : NONE
!-   Controls:
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.5    28/07/99  Created.  C.E. Johnson
!  5.1    04/09/00  360/365 day calendar switch. C.E. Johnson
!  5.2    11/12/01  Changed format for vn5.2. C.E. Johnson
!  6.1    21/10/04  No Change.
!
!-
!VVV  V2.6  WRITEHDR  5/IX/00  - Validity now DAYM
!-----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_INTF
      IMPLICIT NONE
!-----------------------------------------------------------------------
#include "typsize.h"
      INTEGER,                 INTENT(IN) :: month
      INTEGER,                 INTENT(IN) :: year
      INTEGER,                 INTENT(IN) :: fixhd12
      INTEGER,                 INTENT(IN) :: umstepno
      INTEGER,                 INTENT(IN) :: first_constant_r_rho_level

      REAL,                    INTENT(IN) :: daym0
      REAL,                    INTENT(IN) :: daymas
      REAL,                    INTENT(IN) :: z_top_of_model
      INTEGER, DIMENSION(Len_Fixhd), INTENT(OUT) :: flheader
      CHARACTER(LEN=14),             INTENT(OUT) :: filename

      INTEGER                         :: icode
      INTEGER, DIMENSION(a_len_inthd) :: i_const
      REAL, DIMENSION(a_len_realhd)   :: r_const
      REAL, DIMENSION(nlev+1,a_len2_levdepc) :: l_d_const
      REAL, DIMENSION(a_len1_rowdepc,a_len2_rowdepc) :: lenrow
      REAL, DIMENSION(a_len1_coldepc,a_len2_coldepc) :: lencol
      
      CHARACTER(LEN=80)               :: cmessage

#include "c_mdi.h"
! for nunits in cntlall
#include "chsunits.h"
! for lcal360
#include "cntlall.h"
#include "iheadapm.h"
#include "rheadapm.h"

! Specify the Fixed Length Header
      flheader=0
      flheader(1)=20               ! for MASS
      flheader(2)=1                ! Sub-model Indicator (1 = Atmosphere
      flheader(3)=1                ! Hybrid grid
      flheader(4)=0                ! Global data
      flheader(5)=3                ! Fieldsfile
      IF (lcal360) THEN
        flheader(8)=2              ! 360 day calendar
      ELSE
        flheader(8)=1              ! 365 day calendar
      ENDIF
      flheader(9)=3                ! Grid stagger 3 = C grid
      flheader(12)=fixhd12         ! UM version
      flheader(21)=year            ! | First
      flheader(22)=month           ! | Validity
      flheader(23)=INT(daym0)      ! | Time
      flheader(24)=INT((daym0-INT(daym0))*24)
      flheader(25)=MOD(INT(daym0*24*60),60)
      flheader(28)=year            ! + Last
      flheader(29)=month           ! + Validity
      flheader(30)=INT(daymas)     ! + Time
      flheader(31)=INT((daymas-INT(daymas))*24)
      flheader(32)=MOD(INT(daymas*24*60),60)
      flheader(100)=Len_Fixhd+1   ! Start Integer Consts.
      flheader(101)=a_len_inthd    ! Dimension of Integer Consts.
      flheader(105)=flheader(100)+flheader(101) ! Start Real Consts.
      flheader(106)=a_len_realhd    ! Dimension of Real Consts.
      flheader(110)=flheader(105)+flheader(106) ! Start Lev Dep Consts.
      flheader(111)=nlev+1         ! 1st Dimension of Lev. Dep. Consts.
      flheader(112)=a_len2_levdepc  ! 2nd Dimension of Lev. Dep. Consts.
      flheader(115)=flheader(110)+flheader(111)*flheader(112) 
!                                   ! Start Row Dep. Consts.
      flheader(116)=a_len1_rowdepc  ! 1st Dimension of Row Dep. Consts.
      flheader(117)=a_len2_rowdepc  ! 2nd Dimension of Row Dep. Consts.
      flheader(120)=flheader(115)+flheader(116)*flheader(117)            
!                                  ! Start Row Dep. Consts.
      flheader(121)=a_len1_coldepc ! 1st Dimension of Col. Dep. Consts.
      flheader(122)=a_len2_coldepc ! 2nd Dimension of Col. Dep. Consts.
      flheader(125)=imdi            ! Start Fields of Consts.
      flheader(126)=1              ! 1st Dimension of Fields of Consts.
      flheader(127)=1              ! 2nd Dimension of Fields of Consts.
      flheader(130)=imdi            ! Start Extra Consts.
      flheader(131)=1              ! 1st Dimension of Extra Consts.
      flheader(135)=imdi            ! Start History File.cd
      flheader(136)=1              ! 1st Dimension of History File.
      flheader(140)=imdi            ! Start C.F.1
      flheader(141)=1              ! Length
      flheader(142)=imdi            ! Start C.F.1
      flheader(143)=1              ! Length
      flheader(144)=imdi            ! Start C.F.1
      flheader(145)=1              ! Length
                                   ! Start Lookup Table
      flheader(150)=flheader(110)+flheader(111)*flheader(112)
      flheader(151)=Len1_lookup   ! 1st Dimension of Lookup Table
! 2nd Dimension of Lookup Table
      flheader(152)=nlev*(num3dflux+1+numchem)
! Start of Data
      flheader(160)=flheader(150)-1+flheader(151)*flheader(152)
      flheader(161)=nlong*mnlat    ! Dimension of Data

! Specify Integer Constants
      i_const=0
      i_const(ih_a_step)=umstepno   ! No of timesteps since start of run
      i_const(2)=nint((daym0-daymas)*24.0)    ! Meaning interval in hour
! Number of dumps used:
      i_const(3)=nint((daym0-daymas)*daysec/stochem_advection_step)
      i_const(ih_rowlength)=nlong   ! No of points in E-W direction
      i_const(ih_rows)=mnlat        ! No of points in N-S direction
      i_const(ih_model_levels)=nlev ! No of points in vertical direction
      i_const(ih_n_types)=num3dflux+1+numchem  ! No of different field t
      i_const(ih_height_gen)=2      ! Smooth transition of height levels
      i_const(ih_1_c_rho_level)=first_constant_r_rho_level

! Specify real constants
      r_const(rh_deltaEW)=360.0/nlong             ! Grid spacing E-W
      r_const(rh_deltaNS)=180.0/mnlat             ! Grid spacing N-S
      r_const(rh_baselat)=r_const(rh_deltaNS)/2.  ! First PTR row latitu
      r_const(rh_baselong)=r_const(rh_deltaEW)/2. ! First PTR row longti
      r_const(rh_rotlat)=90.0                     ! Pseudo N Pole Latitu
      r_const(rh_rotlong)=-90.0                   ! Pseudo S Pole Latitu
      r_const(rh_z_top_theta)=Z_Top_of_Model      ! Height of top level
      r_const(29)=rmdi                            ! Missing data indicat

! Specify Level Dependant Constants
      l_d_const(:,1)=eta_stochem
      l_d_const(:,3)=rmdi               !RHCrit ? use rmdi
      l_d_const(:,4)=rmdi               !Soil thick
      l_d_const(:,5)=Zsea_stochem
      l_d_const(:,6)=Ck_stochem

      l_d_const(1:nlev,2)=(Eta_stochem(1:nlev)+                         &
     &                     Eta_stochem(0:nlev-1))/2.0
      l_d_const(1:nlev,7)=Zsea_stochem_half
      l_d_const(1:nlev,8)=Ck_stochem_half

      l_d_const(nlev+1,2)=rmdi
      l_d_const(nlev+1,7)=rmdi
      l_d_const(nlev+1,8)=rmdi

! DEPENDS ON: findname
      CALL FINDNAME('p',filetype2,'c',0,0,filename,                     &
     &  month,year)

! DEPENDS ON: file_close
      CALL FILE_CLOSE(20,filename,14,1,0,icode)
      WRITE(6,*)'WRITEHDR: Opening new file ',filename,' on unit 20'
! DEPENDS ON: file_open
      CALL FILE_OPEN(20,filename,14,1,1,icode)
      IF (icode /= 0) THEN
        cmessage='Error opening new PPfile'
! DEPENDS ON: ereport
        CALL EREPORT('WRITEHDR',1,cmessage)
      ENDIF

! DEPENDS ON: init_pp
      CALL INIT_PP(20,'p',flheader(151),flheader(152),                  &
     &             flheader,i_const,r_const,l_d_const,lenrow,lencol,    &
     &             Len_Fixhd,a_len_inthd,a_len_realhd,nlev+1,           &
     &             a_len2_levdepc,a_len1_rowdepc,a_len2_rowdepc,        &
     &             a_len1_coldepc,a_len2_coldepc,                       &
     &             a_len_inthd,a_len_realhd,                            &
     &             icode,cmessage)
      IF (Icode /= 0) THEN
        cmessage='Error initialising new PPfile'
        write(6,*) cmessage
! DEPENDS ON: ereport
        CALL EREPORT('WRITEHDR',1,cmessage)
      ENDIF

      END SUBROUTINE WRITEHDR
#endif
