#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Interface to STASH processing

Module Rcf_Stash_Proc_Mod

!  Subroutine Rcf_Stash_Proc - STASH processing
!
! Description:
! This routine is an interface to most of the stash processing -
! reading in the  stashmaster files and setting up the addressing
! of the output dump
!
! Method:
!   Does some initialisation, then calls rcf_getppx, submodl and
!   rcf_address.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains


Subroutine Rcf_Stash_Proc( EnvDir, UsmEnvDir )

Use Ereport_Mod, Only:   &
    Ereport

Use Rcf_Getppx_Mod, Only    :     &
    Rcf_Getppx

Use Rcf_Address_Mod, Only :       &
    Rcf_Address

Use Rcf_Submodel_Mod, Only :      &
    a_im,                     &
    n_submodel_partition_max, &
    n_internal_model_max,     &
    Internal_Model_Index

Use Rcf_Model_Mod, Only :         &
    H_Vers

Use Rcf_NRecon_Mod, Only :    &
    DumpProgLevs,             &
    PrimDataLen

Use Rcf_PPX_Info_Mod, Only : &
    STM_Record

Implicit None

! Subroutine arguments
Character (Len=9), Intent(In)   :: EnvDir
                                    ! Env. Variable for STASHmaster
                                    ! directory
Character (Len=9), Intent(In)   :: UsmEnvDir
                                    ! Env. Var for user STASHmaster
                                    ! files

! Local Data
Character (Len=*), Parameter    :: RoutineName = 'Rcf_StashProc'
Character (Len=80)              :: CMESSAGE     ! Error return message
Integer                         :: ErrorStatus  ! +ve = fatal error

!--------------------------------------------------------------------

! Initialisation
ErrorStatus =0

! Initialisation of data length arrays
H_VERS( :, : )  = 0

DumpProgLevs(:) = 0
PrimDataLen (:) = 0

! Nullify STM_Record array
Nullify( STM_Record )

! Read STASHmaster files into look-up arrays PPXI, PPXC
If (INTERNAL_MODEL_INDEX(A_IM) > 0) Then
  Call Rcf_Getppx( 'STASHmaster_A', EnvDir )
End If

!Read user STASHmaster files (which may be empty)
Call Rcf_Getppx('             ', UsmEnvDir, .TRUE. )

! Define submodel and section/version configuration
! DEPENDS ON: setmodl
Call SETMODL(ErrorStatus,CMESSAGE)
If (ErrorStatus  /=  0) Then
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Output dump addressing
Call Rcf_Address( )

Return
End Subroutine Rcf_Stash_Proc
End Module Rcf_Stash_Proc_Mod
#endif
