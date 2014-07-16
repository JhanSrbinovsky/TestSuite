#if defined(SCMA)
!
! MODULE SCM_UTILS--------------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT***************************************


MODULE scm_utils

!-------------------------------------------------------------------------------
! Description:
!   SCM UTILITY MODULE; For use with SCM in routines with "USE scm_utils"
!   Contains subroutines:
!   SCM_MESSAGE: to send messages to stdout or unit number with timestep
!                information. Calls with no arguments will output the current 
!                model timestep.
!-------------------------------------------------------------------------------
!
! Current code owner: R. Wong
!
! Version    Date       Comment
!-------------------------------------------------------------------------------
!     6.5    26/04/07   New Module Created                               R. Wong
!-------------------------------------------------------------------------------
!

  Implicit None
  
  Integer :: scm_timestep_count ! Number of current timestep
  Real    :: scm_timestep       ! SCM run timestep

  Integer      , Private :: Nan_list_count = 0 ! Updated by scm trap nan
  Character(50), Private :: Nan_list(50)       ! Updated by scm trap nan

  Logical :: cv_run_mess = .false.  ! Display runtime convective messages

Contains  

  Subroutine scm_message(string, unit)

    CHARACTER(*), Optional, Intent(IN) :: string
    INTEGER,      Optional, Intent(IN) :: unit

    CHARACTER(10)                      :: ts_str
    INTEGER                            :: unitno 

    If (present(unit)) Then 
      unitno = unit
    Else
      unitno = 6
    End If
     
    Write(ts_str,'(A4,I4,A2)') '[TS ', scm_timestep_count,'] '

    If (present(string)) Then
      Write(unitno,*) ts_str, TRIM(ADJUSTL(string))  
    Else
      Write(unitno,*) ts_str  
    End If

    Return
  End Subroutine scm_message

!-----------------------------------------------------------------------------

  Subroutine scm_trap_nan(string, routine)

    Character(*), Intent(IN)                 :: string  ! Diagnostic with NaN
    Character(*), Intent(IN)                 :: routine ! Calling routine
    Integer                                  :: i

    If (nan_list_count <= 50) Then 

      Do i=1, nan_list_count
        If (nan_list(i) == string) Then
          Return
        End If
      End Do

      If (nan_list_count == 50) Then
        call scm_message ('50+ Diagnostics contain NaNs')
        nan_list_count = nan_list_count + 1
      Else  
        call scm_message ('NaN detected: '//string//' from '//routine)
        nan_list(nan_list_count+1) = string
        nan_list_count = nan_list_count + 1
      End If
    End If

    Return
  End Subroutine scm_trap_nan
!=============================================================================
END MODULE scm_utils
#endif
