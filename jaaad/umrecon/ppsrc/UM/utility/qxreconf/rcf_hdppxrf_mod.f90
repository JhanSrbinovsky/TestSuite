
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reads sizes of stashmaster for later space allocation.

Module Rcf_Hdppxrf_Mod

!  Subroutine Rcf_hdppxrf  - finds sizes of stashmaster
!
! Description:
!   Counts the number of records in a stashmaster or reads the
!   header from a user stashmaster
!
! Method:
!   Relies on the fixed formatting of a STASHmaster file to count the
!   the number of records.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   11/10/01   Alter to use Max_Filename_Len parameter. R.Sharp
!   5.3   20/01/02   Using function get_umversion to get um version
!                    number from environment variable $VN
!                    D.M. Goddard
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

SUBROUTINE Rcf_Hdppxrf( StmsrNam , EnvDir , User  )

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_FortranIO_Mod, Only   : &
    Rcf_Get_Unit,               &
    Rcf_Free_Unit,              &
    Max_Filename_Len

Use Rcf_Ppx_Info_Mod, Only: &
    ppxrecs,            &
    nrecs_ustash

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Min

Use Rcf_Parvars_Mod, Only : &
    mype,               &
    nproc

IMPLICIT NONE

! Arguments
Character (Len=13), Intent(In) :: StmsrNam   ! Name of STASHmaster
                                             ! file
Character (Len=9), Intent(In)  :: EnvDir     ! Name of Env. Var for
                                             ! directory for STASHmaster
Logical, Optional, Intent(In)  :: User       ! Flag for use of
                                             ! user-stmaster.

! Local parameters:
Integer, PARAMETER :: no_version = -9999 ! Value returned by function
                                         ! get_umversion if environment
                                         ! variable VN not set

! Local variables
Character (Len=*), Parameter :: RoutineName = 'Rcf_Hdppxrf'
Integer            :: icode
Integer            :: Int_Model_No
Integer            :: ErrorStatus   ! Error return code
Integer            :: FirstBlank
Integer            :: IOStatus
Integer            :: old_ppxRecs   ! Stores old value of ppxRecs
Integer            :: nft           ! Fortran unit number
Integer            :: msg           ! GCOM message id
Integer            :: info          ! GCOM info flag
Integer            :: um_version,   & !Version of UM
                      um_revision     !Revision of UM
Integer            :: stm_version,  & !Version of STASHmaster file
                      stm_revision    !Revision of STASHmaster file
Logical            :: UUsrSTM

Logical            :: found_version   !Indicates presence of STM version
Logical            :: l_exist         !T if STASHmaster file exists

Character (Len=80) :: cmessage
Character (Len=1)  :: char1
Character (Len=8)  :: c_um_version   !UM version as string
Character (Len=8)  :: c_stm_version  !STASHmaster version string

Character (Len=Max_Filename_Len) :: STASH_MSTR
                                    ! Do. STASH master files

! Function & Subroutine calls:
Integer            :: get_um_version

External Fort_Get_Env

!-------------------------------------------------------------------
! Most of the processing here only needs to be done by pe0 as it
! involves file manipulations; the results will be broadcast at the end.
!-------------------------------------------------------------------
pe0: If (mype == 0 ) Then

  IOStatus = 0

! Assign optional argument to no user STASHmaster by default
  If ( PRESENT( User) ) Then
    UUsrSTM = User
  Else
    UUsrSTM = .FALSE.
  Endif

! Store current value of ppxRecs
  old_ppxRecs = ppxRecs

  IF ( .NOT. UUsrSTM ) THEN
! Get a unit number for STASHmaster file
    Call Rcf_Get_Unit( nft )

! Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
    Call Fort_Get_Env( EnvDir, 9, STASH_MSTR, Max_Filename_Len,    &
                       icode)

    If ( icode /= 0 ) Then
      cmessage = 'Cannot get STASH directory from Environment Vars'
      Call Ereport ( RoutineName, icode, cmessage)
    Endif

    FirstBlank = Len_Trim( STASH_MSTR ) + 1
    STASH_MSTR = Trim( STASH_MSTR )

    STASH_MSTR( FirstBlank:FirstBlank ) = '/'
    STASH_MSTR( FirstBlank+1:FirstBlank+13 ) = StmsrNam

    !   Check that the STASHmaster file exists
    INQUIRE ( file=STASH_MSTR, exist=l_exist, iostat=IOStatus )

    If ( .not. l_exist ) then
      Write (6,*) 'STASHmaster File does not exist.'
      Write (6,*) 'FILE=',STASH_MSTR
      ErrorStatus=100
      CMESSAGE=' HDPPXRF: STASHmaster file does not exist.'

      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    !    Open the STASHmaster file

    Open( Unit=nft, File=STASH_MSTR, Iostat=IOStatus )
    If ( PrintStatus >= PrStatus_Normal ) Then
      Write(6,*) '!!!! STASH_MSTR ',STASH_MSTR
    End If

    IF (IOStatus /= 0) THEN
      CMESSAGE= &
        'Error opening STASHmaster file, routine HDPPXRF'
      If ( PrintStatus >= PrStatus_Min ) Then
        Write (6,*) &
          'HDPPXRF: Fortran Error Response = ',IOStatus, &
          ' Opening STASHmaster file ',StmsrNam
      End If

      Call Ereport ( RoutineName, icode, cmessage)
    End If

!    Get the UM version from the environment variable $VN.
! DEPENDS ON: get_um_version
      um_version = get_um_version()

      IF ( um_version .ne. no_version ) THEN

!     Now check through the header section of the STASHmaster
!     file looking for H3
      found_version = .false.
      Read( nft, '(A1)') char1

      Do While( char1 .eq. 'H' .or. char1 .eq. '#')

        If (char1 .eq. 'H') Then
          Backspace nft
          Read (nft, '(1X, A1)') char1

          If (char1 .eq. '3') Then
!  This line starts with H3 and should indicate the
!  STASHmaster version. The line should look like
!     H3| UM_VERSION=4.3
            found_version = .true.
            Backspace nft
            Read (nft, '(15x,a8)') c_stm_version
            Read(c_stm_version, '(i1,1x,i1)' ) &
                   stm_version, stm_revision
            stm_version = stm_version*100 + stm_revision
!     Now perform the check against the UM version

            If (stm_version .ne. um_version) then
              Write (cmessage,*) &
                  'HDPPXRF : UM version and STASHmaster version differ'
              icode = -90

              If ( PrintStatus >= PrStatus_Min ) Then
                Write (6,*) 'Version of STASHmaster file (' &
                  ,stm_version, &
                  ') does not match UM version (' &
                  ,um_version,') in file ',StmsrNam
               End If

              Call Ereport ( RoutineName, icode, cmessage)
            End If  ! version check

          End If  ! char1 == '3'

        End If                 ! char1 == 'H'
        Read (nft, '(a1)') char1

      End Do


      If (.not. found_version) Then
        If ( PrintStatus >= PrStatus_Min ) Then
          Write (6,*) &
             'HDPPXRF : No STASHmaster version available; Unable to'
          Write (6,*) &
             'check against UM version for file ',StmsrNam
        End If

        cmessage = 'HDPPXRF : No STASHmaster version available'
        icode = -100

        Call Ereport( RoutineName, icode, cmessage )
      End If

!     For safety, rewind to the start of the STASH file.
      Rewind (nft)

    End If

! Count records - ppxRecs is counter
    Int_Model_No = 0
    Do While (Int_Model_No /= -1 )
      Read( nft, '(A1)' ) char1
      If ( char1 == '1' ) Then
        Backspace nft
        Read( nft, '(2X,I5)' ) Int_Model_No
        If ( Int_Model_No /= -1 ) Then
          ppxRecs = ppxRecs + 1
        End If
      End If
    End Do

    If (PrintStatus >= PrStatus_Normal ) Then
      Write (6,*) 'Hdppxrf : No of STASHmaster records in ', &
                   StmsrNam,' ',ppxRecs - old_ppxRecs
    End If

    Close( Unit = nft )
    Call Rcf_Free_Unit( nft )

  Else   ! User STASHmaster

    ppxRecs = ppxRecs + NRECS_USTASH

    If (PrintStatus >= PrStatus_Normal ) Then
      Write (6,*) 'Hdppxrf : No of User STASHmaster records ', &
                  NRECS_USTASH
    End If

  End If

End If pe0

!-------------------------------------------------------------------
! Now to broadcast the the variables calculated
!-------------------------------------------------------------------
msg = 7001
Call gc_ibcast( msg, 1, 0, nproc, info, ppxRecs )

Return
End Subroutine Rcf_hdppxrf

End Module Rcf_Hdppxrf_Mod
