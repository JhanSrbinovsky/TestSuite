#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Wrapper for the READFLDS routine.

Module Rcf_Read_Field_Mod

!  Subroutine Rcf_Read_Field - read in a field
!
! Description:
!   A wrappr for the READFLDS routine - thus reads a field into the
!   field data-type
!
! Method:
!   Selects real,int or logical data pointer into which to read data.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Read_Field( field, hdr, decomp )

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Ereport_mod, Only : &
    Ereport

use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_PrintStatus_Mod, Only : &
    LTimer

Use Rcf_Parvars_Mod, Only : &
    current_decomp_type


Implicit None

! Arguments
Type( field_type ), Intent( InOut )    :: field
Type( um_header_type), Intent( In )    :: hdr
Integer, Intent( In )                  :: decomp

! Comdecks
#include "cppxref.h"

! Local vars
Character (Len=*), Parameter :: RoutineName='Rcf_Read_Field'
Character (Len=80)           :: Cmessage
Integer                      :: ErrorStatus
Integer                      :: orig_decomp     ! decomposition on
                                                ! entry to routine

External Timer
!----------------------------------------------------------------

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)
!------------------------------------------------------------------
! Change decomposition if required not same as current
!------------------------------------------------------------------
orig_decomp = current_decomp_type
If (orig_decomp /= decomp) Then
  Call Rcf_Change_Decomposition( decomp )
End If

!------------------------------------------------------------------
! Read in the data - need 3 cases real, logical, integer
!------------------------------------------------------------------
Select Case( field % stashmaster % data_type )
  Case (ppx_type_real)
! DEPENDS ON: rcf_readflds
    Call Rcf_ReadFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                       hdr % Lookup, hdr % Len1Lookup, field % Data,   &
                       field % level_size, hdr % FixHd,                &
                       ErrorStatus, Cmessage )

  Case (ppx_type_int)
! DEPENDS ON: rcf_readflds
    Call Rcf_ReadFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                      hdr % Lookup, hdr % Len1Lookup, field % Data_Int,&
                       field % level_size, hdr % FixHd,                &
                       ErrorStatus, Cmessage )

  Case (ppx_type_log)
! DEPENDS ON: rcf_readflds
    Call Rcf_ReadFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                      hdr % Lookup, hdr % Len1Lookup, field % Data_Log,&
                       field % level_size, hdr % FixHd,                &
                       ErrorStatus, Cmessage )

  Case Default
    ErrorStatus = -10
    Cmessage = 'Unable to read in field as datatype cannot be &
               &determined'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

!-------------------------------------------------------------------
! Check the returned error conditions
!-------------------------------------------------------------------
If ( ErrorStatus /= 0 ) Then
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!--------------------------------------------------------------------
! Change the decomposition back if required
!---------------------------------------------------------------------
If (current_decomp_type /= orig_decomp) Then
  Call Rcf_Change_Decomposition( orig_decomp )
End If

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return

End Subroutine Rcf_Read_Field
End Module Rcf_Read_Field_Mod
#endif
