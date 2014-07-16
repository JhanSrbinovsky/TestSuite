
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reconfiguration Ancillary Processing

Module Rcf_Ancil_Mod

!  Subroutine Rcf_Ancil  - Reconfiguration Ancillary Processing
!
! Description:
!    Controls ancillary processing in the reconfiguration
!
! Method:
!    Reserves an unit number for the ancillary files and
!    calls Rcf_Ancil_Atmos for atmosphere ancillary processing.
!    Ocean ancillary processing to be added later.
!
! Current Code Owner: D. Robinson
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Ancil ( Hdr_In, Hdr_Out,              &
                       Fields_In, Field_Count_In,    &
                       Fields_Out, Field_Count_Out )

Use rcf_ancil_atmos_mod, Only : &
    rcf_ancil_atmos

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    Rcf_Free_Unit

Use Ancil_mod, Only : &
    AncF_UnitNo

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Implicit None

! Arguments
Type (um_header_type), Intent(In)   :: Hdr_In
Type (um_header_type), Intent(In)   :: Hdr_Out
Type (field_type), Pointer          :: Fields_In (:)
Type (field_type), Pointer          :: Fields_Out (:)
Integer, Intent(In)                 :: Field_Count_In
Integer, Intent(Out)                :: Field_Count_Out

External Timer

! DEPENDS ON: timer
If (LTimer) Call Timer( 'Rcf_Ancil', 3)
!-------------------------------------------------------------------
! Get a unit number for Ancillary files
!-------------------------------------------------------------------
Call Rcf_Get_Unit (AncF_UnitNo)

!-------------------------------------------------------------------
! Do processing for Atmosphere Ancillaries
!-------------------------------------------------------------------
Call Rcf_Ancil_Atmos ( Hdr_In, Hdr_Out,                &
                         Fields_In, Field_Count_In,      &
                         Fields_Out, Field_Count_Out )

!-------------------------------------------------------------------
! Free the unit number for Ancillary files
!-------------------------------------------------------------------
Call Rcf_Free_Unit (AncF_UnitNo)

! DEPENDS ON: timer
If (LTimer) Call Timer( 'Rcf_Ancil', 4)

Return
End Subroutine Rcf_Ancil

End Module Rcf_Ancil_Mod
