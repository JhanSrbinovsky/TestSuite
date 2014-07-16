
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
program flumeMain




  use flumerun  ! Declares Flume_run switch - required by Flume and UM

  implicit none




  integer          :: err                 ! used in FORT_GET_ENV call
  character(len=1) :: c_flume_run         ! used in FORT_GET_ENV call
  character(len=1) :: c_flume_print_level ! used in FORT_GET_ENV call

  call FORT_GET_ENV('FLUME_RUN',9,c_flume_run,1,err)

  IF (err  /=  0) THEN
    Flume_run=.FALSE.
  ELSE
    IF (c_flume_run=='T') THEN
      Flume_run=.TRUE.
    ELSE
      Flume_run=.FALSE.
    END IF
  END IF

! DEPENDS ON: um_shell
    call UM_SHELL(.FALSE.)

end program flumeMain
