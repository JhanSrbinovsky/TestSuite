#if defined(RECON) || defined(VAROPSVER)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Error/Warning reporting code

Module Ereport_Mod

!  Subroutine Ereport - handles errors/warnings
!
! Description:
!   The standard reconfiguration error and warning reporting routine.
!   Prints the relevant error and exits if required.
!
! Method:
!   ErrorStatus > 0 is a hard error
!   ErrorStatus < 0 is a warning
!   ErrorStatus = 0 is hunky-dory
!   If on T3E, Cray abort is called if so required - otherwise,
!   GC_Abort is used.
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

  Subroutine Ereport ( RoutineName, ErrorStatus, Cmessage )

  Use Rcf_Parvars_Mod, Only : &
      mype,               &
      nproc

  Implicit None

! Arguments
  Character (Len=*), Intent( In )   :: RoutineName
  Character (Len=*), Intent( In )   :: Cmessage
  Integer, Intent( InOut )          :: ErrorStatus

! Local scalars
  Integer, Parameter            :: unset = -99
  Character (Len=*), Parameter  :: astline = '************************&
  &*****************************************************'
  Character (Len=*), Parameter  :: msg ='Job Aborted from Ereport'

  If ( ErrorStatus > 0 ) Then
    If (mype == 0) Then    ! don't get mangled-up messages in STDERR
    Write (0,*) astline
    Write (0,*) 'ERROR!!! in reconfiguration in routine ', &
                 RoutineName( 1 : Len_Trim(RoutineName) )
    Write (0,*) 'Error Code:- ',ErrorStatus
    Write (0,*) 'Error Message:- ',Cmessage( 1 : Len_Trim(Cmessage) )
    Write (0,*) 'Error generated from processor ', mype
    Write (0,*) astline
    End If

    ! Send errors to STDOUT as well
    Write (6,*) astline
    Write (6,*) 'ERROR!!! in reconfiguration in routine ', &
                 RoutineName( 1 : Len_Trim(RoutineName) )
    Write (6,*) 'Error Code:- ',ErrorStatus
    Write (6,*) 'Error Message:- ',Cmessage( 1 : Len_Trim(Cmessage) )
    Write (6,*) 'Error generated from processor ', mype
    Write (6,*) astline

! DEPENDS ON: um_fort_flush
    Call Um_Fort_Flush(6, unset)
    Call Um_Fort_Flush(0, unset)

    ! On T3E use Cray abort
#if defined (T3E)
    Call Abort( msg )
#else
    Call GC_Abort( mype, nproc, msg )
    STOP
#endif


  Else If ( ErrorStatus < 0) Then
    Write (6,*) astline
    Write (6,*) 'WARNING in reconfiguration in routine ', &
                 RoutineName( 1 : Len_Trim(RoutineName) )
    Write (6,*) 'Warning Code:- ',ErrorStatus
    Write (6,*) 'Warning Message:- ',Cmessage( 1 : Len_Trim(Cmessage) )
    Write (6,*) 'Warning generated from processor ', mype
    Write (6,*) astline
  Endif

! Reset ErrorStatus
  ErrorStatus = 0

  Return
  End Subroutine Ereport

End Module Ereport_Mod
#endif
