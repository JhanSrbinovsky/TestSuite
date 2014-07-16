#if defined(A26_2A)
! C_GRID2GRID start
! Description: Parameters for LAM river routing
! Author: Vicky Bell, CEH Wallingford
!
! History:
! Version  Date    Comment
!  5.5  28/02/03   Original code. Vicky Bell
!  6.0  17/09/03   Change def from A20 to A26.  D. Robinson.
!  6.0  12/09/03   Remove g2g_timestep. V.A.Bell

      REAL, PARAMETER :: cland  = 0.2  ! land wave speed (m/s)
      REAL, PARAMETER :: criver = 0.2  ! subsurf river wave speed (m/s)
      REAL, PARAMETER :: cbland  = 0.18! subsurf land wave speed (m/s)
      REAL, PARAMETER :: cbriver = 0.18! subsurf river wave speed (m/s)
      REAL, PARAMETER :: runoff_factor = 0.7
!                                      ! runoff volume factor
      REAL, PARAMETER :: retl = 0.     ! return flow (land squares) (<1)
      REAL, PARAMETER :: retr = 0.15   ! return flow (river squares) (<1)
      REAL, PARAMETER :: slfac = 0.    ! slope factor (not used yet)
      INTEGER, PARAMETER :: a_thresh = 1 ! threshold area
! END C_GRID2GRID
#endif
!---------------------------------------------------------------------
