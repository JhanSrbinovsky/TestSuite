#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Write a message to the operator for operational runs
!
! Subroutine Interface:
      SUBROUTINE OperatorMessage(nproc)

      IMPLICIT NONE
!
! Description: OperatorMessage writes a message to the operator for
! operational runs to indicate that jobs have started correctly.
!
!
! Method: Use SHELL routine to fork a process without copying a memory
! image. This allows unix command msgi to write a messsage to the
! operator. Only activated if environment variable UM_OPER_MESS="true".
! Note: fort_get_env fails if contents of environment variable is longer
! than the available storage, so length of UM_OPER_MESS (or RUNID)
! should not exceed env_char_len (=8 here).
!
! Current Code Owner: Rick Rawlins
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.5   16/02/98  Original code. Rick Rawlins
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:

! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER  nproc                  ! Number of processors

! Local parameters:

! Local scalars:
      CHARACTER*80 shell_arg     ! text for argument of shell call
      CHARACTER*60 message
      CHARACTER*8  env_char      ! text of env. variable

      INTEGER shell_arg_len                                             &
                                 ! length of shell argument
     &       ,command_len                                               &
                                 ! length of shell command
     &       ,icode                                                     &
                                 ! error return code
     &       ,RUNID_char_len                                            &
                                 ! length of RUNID contents
     &       ,env_char_len       ! length of env_char array

      Parameter(                                                        &
     &           RUNID_char_len=5)

! Function & Subroutine calls:
      External fort_get_env,shell

!- End of header
      message='"xxxx UM successfully running on xxxx PEs "' ! template
      shell_arg='msgi'            ! shell command
      command_len=4               ! length of shell command
      env_char_len=len(env_char)

! Check that env.var. UM_OPER_MESS="true" in top level unix script
      env_char='false'
      CALL fort_get_env('UM_OPER_MESS',12,env_char,env_char_len,icode)

      IF(icode == 0 .AND. env_char(1:4) == 'true') THEN

! Construct message: place RUNID and no. of PEs in message
        CALL fort_get_env('RUNID',5,env_char,env_char_len,icode)
        IF(icode /= 0) THEN
          write(6,*) 'Routine OperatorMessage: Warning: problem in ',   &
     &               'reading environment variable RUNID, icode=',icode
        ENDIF

        message(2:RUNID_char_len+1)=env_char  ! substitute RUNID
        write(message(34:37),'(I4)') nproc    ! substitute no. of PEs
        shell_arg(command_len+2:)=message

! Send message to operator
        shell_arg_len=len(shell_arg)
        CALL shell(shell_arg,shell_arg_len)
        write(6,*) shell_arg                ! and standard output

      ENDIF                               ! Test on UM_OPER_MESS

      Return
      END SUBROUTINE OperatorMessage
#endif
