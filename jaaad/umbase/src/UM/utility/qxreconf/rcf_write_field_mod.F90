#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper for WRITFLDS

Module Rcf_Write_Field_Mod

!  Subroutine Rcf_Write_Field
!
! Description:
!    Wrapper for writflds - utilising fields data types
!
! Method:
!    Deals seperately with real, integer and logical data.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.4   12/06/02   GRIB PRep. Added Arg to call to WritFlds. Rod.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Write_Field( field, hdr, decomp, free )

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Dealloc_Field

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Ereport_mod, Only : &
    Ereport

use Rcf_UMhead_Mod, ONly : &
    um_header_type

use Rcf_PrintStatus_Mod, Only : &
    LTimer

Use Rcf_Parvars_Mod, Only : &
    current_decomp_type

Implicit None

! Arguments
Type( field_type ), Intent( InOut )    :: field
Type( um_header_type), Intent( In )    :: hdr
Integer, Intent( In )                  :: decomp
Logical, Optional, Intent( In )        :: free

! Comdecks
#include "cppxref.h"

! Local vars
Character (Len=*), Parameter :: RoutineName='Rcf_Write_Field'
Character (Len=80)           :: Cmessage
Integer                      :: ErrorStatus
Logical                      :: free_data
Integer                      :: orig_decomp ! decomposition when
                                            ! routine is entered

External Rcf_WritFlds, Timer

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
! First set the free_data logical appropriately
!------------------------------------------------------------------
If ( Present( free ) ) Then
  free_data = free
Else
  free_data = .FALSE.
End If

!------------------------------------------------------------------
! Now write out the data - need 3 cases real, logical, integer
!------------------------------------------------------------------
Select Case( field % stashmaster % data_type )
  Case (ppx_type_real)
! DEPENDS ON: rcf_writflds
    Call Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                   hdr % Lookup, hdr % Len1Lookup, field % Data,       &
                   field % level_size, hdr % FixHd,                    &
                   ErrorStatus, Cmessage , .True.)

  Case (ppx_type_int)
! DEPENDS ON: rcf_writflds
    Call Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                   hdr % Lookup, hdr % Len1Lookup, field % Data_Int,   &
                   field % level_size, hdr % FixHd,                    &
                   ErrorStatus, Cmessage , .True.)

  Case (ppx_type_log)
! DEPENDS ON: rcf_writflds
    Call Rcf_WritFlds( hdr % UnitNum, field % levels, field % dump_pos,&
                   hdr % Lookup, hdr % Len1Lookup, field % Data_Log,   &
                   field % level_size, hdr % FixHd,                    &
                   ErrorStatus, Cmessage , .True.)

  Case Default
    ErrorStatus = -10
    Cmessage = 'Unable to write out field as datatype cannot be &
               &determined'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

!-------------------------------------------------------------------
! Check the returned error conditions
!-------------------------------------------------------------------
If ( ErrorStatus /= 0 ) Then
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-------------------------------------------------------------------
! Free the allocated data-space if required to
!-------------------------------------------------------------------
If (free_data) Then
  Call Rcf_Dealloc_Field( field )
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

End Subroutine Rcf_Write_Field
End Module Rcf_Write_Field_Mod
#endif
