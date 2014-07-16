! Description: COMDECK containing uv horizontal types
!  for use in idealised problems
!
! Author : T. Davies
! History:
! Version  Date      Comment.
! 5.3      15/11/01  New code
! 5.5      03/02/03  Add uv horizontal and vertical options. R.Forbes

! u,v wind horizontal variation options
      INTEGER, PARAMETER :: uv_horiz_const   = 0 ! constant in horiz.
      INTEGER, PARAMETER :: uv_horiz_ramp    = 1 ! symmetric ramp fn
      INTEGER, PARAMETER :: uv_horiz_balance = 2 ! balance from pressure
      INTEGER, PARAMETER :: uv_horiz_deform  = 3 ! deformation field

! u,v wind vertical profile options
      INTEGER, PARAMETER :: uv_vert_const    = 0  ! constant with hght
      INTEGER, PARAMETER :: uv_vert_interp   = 1  ! interp using 4 vals
      INTEGER, PARAMETER :: uv_vert_namelist = 9  ! read from namelist
      INTEGER, PARAMETER :: uv_vert_dump     = 10 ! read from dump
