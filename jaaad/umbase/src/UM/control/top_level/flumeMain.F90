#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
program flumeMain

#if defined(FLUME)
  use MFlm, only : MFlm_run
#endif
  use flumerun  ! Declares Flume_run switch - required by Flume and UM

  implicit none

#if defined(FLUME)
  integer          :: debugLevel
#endif
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

#if defined(FLUME)
  IF (Flume_run) THEN
    call FORT_GET_ENV('FLUME_print_level',17,c_flume_print_level,1,err)

    debugLevel=0
    IF (err /= 0) THEN
      write(6,*) 'Flume: failed to get print level, setting to 1'
      debugLevel=1
    ELSE
      select case (c_flume_print_level)
      case ('1') 
        debugLevel=1
      case ('2') 
        debugLevel=2
      case ('3') 
        debugLevel=3
      case ('4') 
        debugLevel=4
      case ('5') 
        debugLevel=5
      case ('6') 
        debugLevel=6
      case ('7') 
        debugLevel=7
      case ('8') 
        debugLevel=8
      case default
        debugLevel=1
      end select
    END IF

    write(6,*) 'FLUME UM RUN: print level', debugLevel 
    call MFlm_run(debugLevel)
  ELSE
! DEPENDS ON: um_shell
    call UM_SHELL(.FALSE.)
  END IF
#else
! DEPENDS ON: um_shell
    call UM_SHELL(.FALSE.)
#endif

end program flumeMain
#endif
