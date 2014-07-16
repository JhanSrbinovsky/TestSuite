#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert GRIB LSM to Logical UM one

Module Rcf_Grib_Spcl_LSM_Mod

! SUBROUTINE Rcf_Grib_Spcl_LSM
!
! Description: Converts LSM from 'real' type to 'logical'
!
! Method: flip flop the logical based on data > 0.9999 for .True.
!
! Current Code Owner: Roddy Sharp
!
! History:
! Version   Date      Comment
! -------  --------   -------------------------
!  5.4     15/07/02   Original code. Roddy Sharp (frtz)
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_LSM(Current,FpData,LgData,UM_Hdr,Marker)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     &
    LenArrayMax,     &
    Grb_Data_Real,   &
    Grb_Data_Log

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Diag

Use Rcf_UMhead_Mod, Only :  &
    um_header_type            ! Derived containing UM header info

Use Rcf_HeadAddress_Mod, Only : &
    IC_NumLandPoints

Use EReport_Mod, Only :     &
    EReport

Use Rcf_Parvars_mod, Only : &
    mype

Implicit None

! Declarations:
! Global variables (#include statements etc):
#include "clookadd.h"
! contains LBLREC (amongst others)

! Subroutine arguments

!< Scalar arguments with intent(In):>
Integer                             :: Marker  ! position within lookups

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer          :: Current
Real, Intent(InOut)                 :: FpData(LenArrayMax)
Logical, Intent(InOut)              :: LgData(LenArrayMax)
Type (Um_Header_type)               :: UM_Hdr

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_LSM'

! Local variables

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

Integer                          :: Count,I       ! used to count no.
                                                  ! of land points

!=======================================================================
!  Routine Code Start :
!=======================================================================

Count = 0

! Double Check the datatype is still 'real'
If (Current % Data_Type == Grb_Data_Real ) Then

  If ( PrintStatus >= PrStatus_Diag  ) Then
    If ( mype == 0 ) Then
      Write (6,'(A)') 'Converting Real LSM to Logical LSM'
    End If
  End If

  Where (FpData >0.99999)
    LgData = .True.  ! It's a land point
  ElseWhere
    LgData = .False. ! It's not
  EndWhere

  Current % Data_Type = Grb_Data_Log

Else   ! It wasn't real data I recieved.
  Cmessage(1) = 'Couldnt convert data type recieved to logical'
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )

End If

Do I = 1, UM_Hdr % Lookup (LBLREC,Marker)
  If (LgData(I)) Then
    Count =Count + 1
  End If
End Do

! record the no. of land points in the header
UM_Hdr % IntC( IC_NumLandPoints ) = Count

Return

End Subroutine Rcf_Grib_Spcl_LSM
End Module Rcf_Grib_Spcl_LSM_Mod
#endif
