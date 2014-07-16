#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper for WritHead

Module Rcf_WriteUMhdr_mod

!  Subroutine Rcf_WriteUMhdr - wrapper for UM routine
!
! Description:
!   Simple wrapper routine for WritHead
!
! Method:
!   Wrapper routine for WritHead - setting the addressing correctly
!   and converting from the GRID data-type.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   24/09/01   Correction of Intent for UMhdr. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_WriteUMhdr( UMhdr )

Use Rcf_UMhead_Mod, Only : &
    UM_header_type,    &
    LenFixHd

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Parvars_Mod, Only : &
    mype

Implicit None

! Arguments
Type (UM_header_type), Intent( InOut ) :: UMhdr

! Local constant
Character (Len=*), Parameter :: RoutineName = 'Rcf_WriteUMhdr'

! Local variables
Character (Len=80)           :: Cmessage
Integer                      :: ErrorStatus
Integer   ::  WordAddress        ! Position on file, used in SETPOS
Integer   ::  number_of_data_words_in_memory
Integer   ::  number_of_data_words_on_disk
Integer   ::  disk_address       ! Current rounded disk address
                                 ! and final data length

Integer   :: ppxi_dummy         ! dummy for subroutine call
Integer   :: ppxrecs_dummy      ! dummy for subroutine call
Character :: ppxc_dummy         ! dummy for subroutine call
External Writhead, SetPos, Set_Dumpfile_Address

!----------------------------------------------------------------

  WordAddress = 0
  ErrorStatus = 0
! DEPENDS ON: setpos
  Call SetPos (UMhdr % UnitNum, WordAddress, ErrorStatus)

  If ( ErrorStatus /= 0) Then
    Cmessage = 'SETPOS failure, repositioning to start of dump!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

! DEPENDS ON: set_dumpfile_address
  Call Set_Dumpfile_Address( UMhdr % FixHd, LenFixHd, UMhdr % Lookup, &
                             UMhdr % Len1Lookup,  UMhdr % Len2Lookup, &
                             number_of_data_words_in_memory, &
                             number_of_data_words_on_disk, disk_address)

! DEPENDS ON: writhead
  Call WritHead(  UMhdr % UnitNum,                        &
                  UMhdr % FixHd,     LenFixHd,            &
                  UMhdr % IntC,      UMhdr % LenIntC,     &
                  UMhdr % RealC,     UMhdr % LenRealC,    &
                  UMhdr % LevDepC,   UMhdr % Len1LevDepC, &
                                     UMhdr % Len2LevDepC, &
                  UMhdr % RowDepC,   UMhdr % Len1RowDepC, &
                                     UMhdr % Len2RowDepC, &
                  UMhdr % ColDepC,   UMhdr % Len1ColDepC, &
                                     UMhdr % Len2ColDepC, &
                  UMhdr % FldsOfC,   UMhdr % Len1FldsOfC, &
                                     UMhdr % Len2FldsOfC, &
                  UMhdr % ExtraC,    UMhdr % LenExtraC,   &
                  UMhdr % HistFile,  UMhdr % LenHistFile, &
                  UMhdr % CompFldI1, UMhdr % LenCompFldI1, &
                  UMhdr % CompFldI2, UMhdr % LenCompFldI2, &
                  UMhdr % CompFldI3, UMhdr % LenCompFldI3, &
                  UMhdr % Lookup,    UMhdr % Len1Lookup,   &
                                     UMhdr % Len2Lookup,   &
                  UMhdr % LenData,                         &
                  ppxi_dummy, ppxc_dummy, ppxrecs_dummy,   &
                  UMhdr % StartData,                       &
                  ErrorStatus,        CMessage )

  If (ErrorStatus /= 0) Then
    Cmessage='Error in Rcf_WritHead!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If


Return

End Subroutine Rcf_WriteUMhdr

End Module Rcf_WriteUMhdr_Mod
#endif
