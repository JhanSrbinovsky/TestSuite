#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT 2000, Met Office, All Rights Reserved.
! Please refer to file $UMDIR/vn$VN/copyright.txt for further details
! *****************************COPYRIGHT*******************************

Module Rcf_smc_stress_Mod

!  Subroutine Rcf_smc_stress - calculate the soil moisture stress
!
! Description:
! This module calculate soil moisture stress fields,
! and then interpolate smc_stress into the new domain.
! This smc_stress will be converted into the soil moisture
! in the post process routine.
! 
! Method:
!  Read input of smc, smc at wilting and critical point.
!  Calculate smc stress using these inputs and then interpolate
!  smc_stress to the new domain. The smc in the new domain will be
!  then calculated from smc_stress in the post process routine.
!
!
! Current Code Owner: Changgui Wang
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_smc_stress(smc_field, fields_in, field_count_in, hdr_in)

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Parvars_Mod, Only : &
    mype

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid
    
Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_input

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field 

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_vol_smc_wilt,&
    stashcode_vol_smc_cri, &
    stashcode_prog_sec

Implicit None

! Arguments
Type( field_type), Intent( InOut )   :: smc_field

Type( field_type ), Pointer          :: fields_in(:)

Type( field_type ), Pointer          :: smc_wlt
Type( field_type ), Pointer          :: smc_cri

Type( um_header_type ), Intent(In)   :: hdr_in
Integer, Intent(In)                  :: field_count_in

#include "c_mdi.h" 
#include "c_densty.h"

! Local variables
Integer                              :: i
Integer                              :: j
Integer                              :: k
Integer                              :: pos
Integer                              :: decomp_old   ! old decomposition
Integer                              :: stat         ! gcom status

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Processing soil moisture (stashcode 9) '
End If

! Convert smc to VOL SMC:
  Do j = 1, Input_Grid % sm_levels
     smc_field % data(:,j) = smc_field % data(:,j)/ &
      ( Input_Grid % soil_depths(j)*RHO_WATER )               
  End do

!------------------------------------------------------------------
! find and read input smc_wlt
!------------------------------------------------------------------
Call Rcf_Locate(stashcode_prog_sec, stashcode_vol_smc_wilt, fields_in, field_count_in, pos )
smc_wlt  => fields_in(pos)
Call Rcf_Alloc_Field( smc_wlt )
Call Rcf_Read_Field(  smc_wlt, Hdr_In, decomp_rcf_input )

!------------------------------------------------------------------
! find and read input smc_cri
!------------------------------------------------------------------
Call Rcf_Locate(stashcode_prog_sec, stashcode_vol_smc_cri, fields_in, field_count_in, pos )
smc_cri  => fields_in(pos)
Call Rcf_Alloc_Field( smc_cri )
Call Rcf_Read_Field( smc_cri, Hdr_In, decomp_rcf_input )

! Calculate smc_stress and use the new data for smc:
Do j = 1, smc_field % levels

    Do k=1, smc_field % level_size

      If (smc_cri%data(k,1)-smc_wlt%data(k,1) == 0 .or. &
        smc_cri%data(k,1)== RMDI.or.                  &
        smc_wlt%data(k,1)== RMDI.or.                  &
        smc_field % data(k,j) == RMDI)then
        smc_field % data(k,j) = RMDI
      Else
        smc_field % data(k,j)= &
        (smc_field % data(k,j)-smc_wlt % data(k,1))/ &
        (smc_cri % data(k,1)-smc_wlt % data(k,1))
      End if
    End do            
Enddo
  If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
    write (6,*)'Input Soil moisture has been converted to soil stress.'                 
  End If

  Call Rcf_DeAlloc_Field( smc_cri )
  Call Rcf_DeAlloc_Field( smc_wlt )
Return
End Subroutine Rcf_smc_stress
End Module Rcf_smc_stress_Mod
#endif
