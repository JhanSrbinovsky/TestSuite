
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+Decode the STASH level code

Module Rcf_Level_Code_Mod

!  Subroutine Rcf_Level_Code_Mod - decodes the STASH level code
!
! Description:
!   Sets ILOUT to an appropriate level size according to the value of
!   ILIN.
!
! Method:
!******************************************************************
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
! Any Changes to this routine must be accompanied with equivalent
! changes to the deck LEVCOD1
!******************************************************************
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

SUBROUTINE Rcf_Level_Code( ILIN, ILOUT, Grid )

Use Rcf_Model_Mod, Only :  &
    STLevGWDrag,           &
    BotVDiffLev,           &
    TopVDiffLev,           &
    Oaslev

Use Rcf_CntlAtm_Mod, Only : &
    H_SWBands,              &
    H_LWBands

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Ereport_Mod, Only : &
    Ereport

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)    :: ILIN           ! Model level code
Integer, Intent(Out)   :: ILOUT          ! The actual level

Type (Grid_Type), Intent(In) :: Grid     ! The grid that decoded
                                         ! values should correspond to


! Local variables
Integer                      :: ErrorStatus
Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName = 'Rcf_Level_Code'


Select Case ( ILIN )
  Case ( 1 )                            ! First atmos level
    ILOUT=1

  Case ( 2 )                            ! Top atmos level
    ILOUT= Grid % MODEL_LEVELS

  Case ( 3 )                            ! Top wet level
    ILOUT= Grid % WET_LEVELS

  Case ( 4 )
    ILOUT= Grid % MODEL_LEVELS - 1

  Case ( 5 )                            ! First boundary layer level
    ILOUT=MIN(1,Grid % BL_LEVELS)

  Case ( 6 )                            ! Last boundary layer level
    ILOUT=Grid % BL_LEVELS

  Case ( 7 )
    ILOUT= Grid % BL_LEVELS+1

  Case ( 8 )                            ! First soil level
    ILOUT=MIN(1, Grid % ST_LEVELS)

  Case ( 9 )                            ! Last soil level
    ILOUT= Grid % ST_LEVELS

  Case ( 10 )                           ! First tracer level
    ILOUT= Grid % MODEL_LEVELS - Grid % TR_LEVELS+1

  Case ( 11 )                           ! Last tracer level
    ILOUT= Grid % MODEL_LEVELS

  Case ( 12 )
    ILOUT= Grid % MODEL_LEVELS+1

  Case ( 13 )                           ! First gravity wave drag level
    ILOUT=StLevGWdrag

  Case ( 14 )                           ! Last gravity wave drag level
    ILOUT= Grid % MODEL_LEVELS

  Case ( 15 )
    ILOUT=BotVDiffLev

  Case ( 16 )
    ILOUT=TopVDiffLev-1

  Case ( 17 )
    ILOUT=TopVDiffLev

  Case ( 18 )
    ILOUT= Grid % BL_LEVELS-1

  Case ( 19 )
    ILOUT= Grid % MODEL_LEVELS+1

  Case ( 20 )
    ILOUT=MIN(2, Grid % ST_LEVELS)

  ! Ocean removed at vn7.0 so this is redundant
  Case ( 21 )
    ILOUT=1
  
  ! Ocean removed at vn7.0 so this is redundant
  Case ( 22 )
    ILOUT=0

  Case ( 23 )
    ILOUT= Grid % OZONE_LEVELS

  Case ( 24 )
    ILOUT= Grid % MODEL_LEVELS*H_SWBANDS

  Case ( 25 )
    ILOUT=( Grid % MODEL_LEVELS+1)*H_SWBANDS

  Case ( 26 )
    ILOUT= Grid % WET_LEVELS*H_SWBANDS

  Case ( 27 )
    ILOUT= Grid % MODEL_LEVELS*H_LWBANDS

  Case ( 28 )
    ILOUT=( Grid % MODEL_LEVELS+1)*H_LWBANDS

  Case ( 29 )
    ILOUT= Grid % WET_LEVELS*H_LWBANDS

  Case ( 30 )
    ILOUT=2

  Case ( 32 )
    ILOUT=H_SWBANDS

  Case ( 33 )
    ILOUT=H_LWBANDS

  Case ( 34 )
    ILOUT= Grid % SM_LEVELS

  Case ( 35 )
    ILOUT= Grid % CLOUD_LEVELS

  Case ( 36 )                       ! Wave model first level (direction)
    ILOUT=1

  Case ( 37 )                       ! Wave model last level (direction)
                                    ! No wave model in rcf
!    ILOUT=NANG
     ILOUT=0

  Case ( 38 )                       ! Surace theta level
    ILOUT=0

  Case ( 41 )
! Allow room for expansion of ocean assimilation groups.
    ILOUT=OASLEV(1)

  Case ( 42 )
! Allow room for expansion of ocean assimilation groups.
    ILOUT=OASLEV(2)

  Case ( 43 )
! Allow room for expansion of ocean assimilation groups.
    ILOUT=OASLEV(3)

  Case ( 44 )
! Allow room for expansion of ocean assimilation groups.
    ILOUT=OASLEV(4)

  Case ( 45 )
! Allow room for expansion of ocean assimilation groups.
    ILOUT=OASLEV(5)

  Case ( 46 )
! Allow room for expansion of ocean assimilation groups.
    ILOUT=OASLEV(6)

  Case Default
    WRITE(Cmessage,*) 'LEVCOD: IMPOSSIBLE LEVEL CODE FOUND ',ILIN
    ErrorStatus=1
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

RETURN
END Subroutine Rcf_Level_Code
End Module Rcf_Level_Code_Mod


