
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top level program for reconfiguration

Program Reconfigure

! Description:
!  Top level program for reconfiguration
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   25/09/01   Added call to Rcf_Finalise. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Use Rcf_Initialise_Mod, Only : &
    Rcf_Initialise

Use Rcf_Finalise_Mod, Only : &
    Rcf_Finalise

Use Rcf_Read_Namelists_Mod, Only : &
    Rcf_Read_Namelists

Use Rcf_Control_Mod, Only : &
    Rcf_Control

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Parvars_mod, Only : &
    mype,               &
    nproc

Use Ereport_mod, Only :&
    Ereport

Implicit None

Type (um_header_type)        :: hdr_in       ! header from input dump
Type (um_header_type)        :: hdr_out      ! header from output dump
Integer                      :: ErrorStatus
Integer                      :: info         ! GCOM dummy
Character (Len=*), Parameter :: RoutineName = 'Reconfigure'
Character (Len=80)           :: Cmessage

External gc_init, gc_exit

!------------------------------------------------------------------
! Initialise GCOM
!------------------------------------------------------------------
Call Gc_Init( ' ', mype, nproc )
If ( nproc < 0 ) Then
  ErrorStatus = 20
  Cmessage = 'Parallel Initialisation Failed!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-------------------------------------------------------------------
! Perform initialisation
!-------------------------------------------------------------------
Call Rcf_Initialise( hdr_in, hdr_out )

!-------------------------------------------------------------------
! Do the real work
!-------------------------------------------------------------------
Call Rcf_Control( hdr_in, hdr_out )

!------------------------------------------------------------------
! Tidy Up
!------------------------------------------------------------------
Call Rcf_Finalise( hdr_in, hdr_out )

!------------------------------------------------------------------
! Kill off gcom
!------------------------------------------------------------------
Write (6,*) ' End of rcf program reached. PE ',mype
Call gc_exit()

End Program Reconfigure
