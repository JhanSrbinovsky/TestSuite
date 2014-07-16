!----------------------------------------------------------------------
! comdeck: CFLDDIMA
! Purpose: declares dimensions of fields.
!          This deck is linked to AFLDDIMA and CFLDDIMS
! History:
! version  date         change
! 5.3      21/09/00     New comdeck
! Author:  M. J. Bell
!----------------------------------------------------------------------
! declarations:

! atmosphere grid
      integer ncols     !  number of columns
      integer nrowst    !  number of rows (tracer grid)
      integer nrowsu    !  number of rows (u grid)
      integer nrowsv    !  number of rows (v grid)
      integer nrowsuv   !  number of rows on grid with coincident velys
                        !  (at B grid locations both for B & C grids)

!----------------------------------------------------------------------
