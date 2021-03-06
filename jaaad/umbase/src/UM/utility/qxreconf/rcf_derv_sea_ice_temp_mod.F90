#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Derive sea ice temp from surface temp.

Module Rcf_Derv_Sea_Ice_Temp_Mod

!  Subroutine Rcf_Derv_Sea_Ice_Temp
!
! Description:
!   Copies T* to sea ice temp.
!
! Method:
!   A basic copy with no calculation.
!
! Current Code Owner: R.Sharp
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   06/02/03   Original code.  R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Derv_Sea_Ice_Temp( fields_out, field_count_out,        &
                                  sea_ice_temp, hdr_out )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_tstar,           &
    stashcode_prog_sec

Implicit None

! Arguments
Integer,               Intent(In)    :: field_count_out
Type( field_type ),    Target        :: fields_out( field_count_out )
Type( field_type ),    Intent(InOut) :: sea_ice_temp
Type( um_header_type), Intent(In)    :: hdr_out

! Comdecks
#include "c_mdi.h"

! Local variables
Integer                        :: pos
Integer                        :: i
Type( field_type ), Pointer    :: t_star

!-----------------------------------------------------------------
! Write out our action if appropriate
!-----------------------------------------------------------------
If ( mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
  Write (6,*) 'Copying T* to Sea Ice Temp'
End If

!-----------------------------------------------------------------
! Get T* from output dump
!-----------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_tstar,                &
                 fields_out, field_count_out, pos )
t_star => fields_out( pos )
Call Rcf_Alloc_Field( t_star )
Call Rcf_Read_Field( t_star, Hdr_Out, decomp_rcf_output )

!------------------------------------------------------------------
! Do the copy...
!------------------------------------------------------------------
sea_ice_temp % Data( : , 1) = t_star % Data( : , 1)

!-------------------------------------------------------------------
! clean up
!-------------------------------------------------------------------
Call Rcf_DeAlloc_Field( t_star )

Return
End Subroutine Rcf_Derv_Sea_Ice_Temp
End Module Rcf_Derv_Sea_Ice_Temp_Mod
#endif
