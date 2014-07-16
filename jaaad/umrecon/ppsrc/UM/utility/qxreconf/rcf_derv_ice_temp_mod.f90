
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reinitialises sea ice surface temperature

Module Rcf_Derv_Ice_Temp_Mod

! Description:
!     This subroutine is the initialises the sea ice temperature
!     on ice catagories when these are called for in stash
!     It sets all ice category pseudo levels to have temperatures
!     equal to the single level sea ice temperature (stash=49)
!
!
! Current Code Owner: J.Ridley
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.5   12/03/03   Original code.  J.Ridley
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Derv_Ice_Temp( fields_out, field_count_out, hdr_out,  &
                            ice_temp_cat )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_sea_ice_temp,    &
    stashcode_prog_sec

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field


Implicit None

! Arguments
Type( field_type ), Pointer       :: fields_out(:)
Type( um_header_type), Intent(In) :: hdr_out
Type( field_type ), Intent(InOut) :: ice_temp_cat
Integer, Intent(In)               :: field_count_out

! Internal variables
Type( field_type ), Pointer       ::  ice_temp


Integer                           ::  pos   ! position in array
Integer                           ::  i     ! loop index



!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Reinitialising ice surface temperature'
  End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------
! GBM ice temperature; will abort if ice_temp not found

  Call Rcf_Locate( stashcode_prog_sec, stashcode_sea_ice_temp,       &
                   fields_out,field_count_out,pos)

  ice_temp => fields_out(pos)
  Call Rcf_Alloc_Field( ice_temp )
  Call Rcf_Read_Field( ice_temp, hdr_out, decomp_rcf_output )


!----------------------------------------------------------------------
! Loop through ice_temp_cat
!----------------------------------------------------------------------

  Do i = 1,ice_temp_cat % levels
    ice_temp_cat % data(:,i) = ice_temp % data(:,1)
  EndDo

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------
  Call Rcf_Dealloc_Field( ice_temp )

Return
End Subroutine Rcf_Derv_Ice_Temp
End Module Rcf_Derv_Ice_Temp_Mod
