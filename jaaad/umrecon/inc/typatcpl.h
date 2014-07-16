#if defined(ATMOS) || defined (A34_1A)
! TYPATCPL start
! Description: Gridline coordinates for interpolation and area-averaging
! between atmosphere and river-routing grids (Part of TYPAOCPL.hdk)
! Author: C.Bunton 28.02.03
!
! History:
! Version  Date    Comment
!  5.5  28/02/03   Original code. C. Bunton
!
!
      REAL :: XPA(AOCPL_ROW_LENGTH+1)  ! Atmosphere TP longitude coordina
      REAL :: XUA(0:AOCPL_ROW_LENGTH)  ! Atmosphere U longitude coordinat
      REAL :: XVA(AOCPL_ROW_LENGTH+1)  ! Atmosphere V longitude coordinat
      REAL :: YPA(AOCPL_P_ROWS)        ! Atmosphere TP latitude coordinat
      REAL :: YUA(AOCPL_P_ROWS)        ! Atmosphere U latitude coordinate
      REAL :: YVA(0:AOCPL_P_ROWS)      ! Atmosphere V latitude coordinate
! TYPATCPL end
#endif
