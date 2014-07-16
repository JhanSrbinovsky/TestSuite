
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisation of STASH related variables

Module Rcf_Stash_Init_Mod

!  Subroutine Rcf_Stash_Init - initialisation of STASH variables
!
! Description:
!   Reads headers for stashmasters and calls rcf_stash_proc to
!   read stashmaster and set up output dump addressing.
!
! Method:
!   Calls rcf_hdppxrf and rcf_stash_proc
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

Subroutine Rcf_Stash_Init( )

Use Rcf_Parvars_mod, Only : &
    mype

Use Rcf_hdppxrf_mod, Only : &
    Rcf_hdppxrf

Use Rcf_Ppx_Info_Mod, Only: &
    ppxrecs

Use Rcf_Stash_Proc_mod, Only : &
    Rcf_Stash_Proc

Use Rcf_Submodel_Mod, Only : &
    N_Internal_Model,    &
    Internal_Model_List, &
    Atmos_IM

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Implicit None

! Local variables
Integer                :: i

!-------------------------------------------------------------------
! Read in STASHmaster headers
!-------------------------------------------------------------------
Do i = 1, N_Internal_Model
  If ( Internal_Model_List(i) == Atmos_IM ) Then
    Call Rcf_hdppxrf( 'STASHmaster_A','STASHMSTR' )
  End If
End Do

! User STASHmaster header
Call Rcf_hdppxrf( '             ','         ', .TRUE. )

If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
  Write (6,*) 'Rcf_Stash_Init : Total No of STASHmaster records ', &
               ppxRecs
End If

!---------------------------------------------------------------------
! Read in the STASHmasters, and do the output dump addressing
!---------------------------------------------------------------------
Call Rcf_Stash_Proc( 'STASHMSTR', 'USTSHMSTR' )


End Subroutine Rcf_Stash_Init

End Module Rcf_Stash_Init_Mod
