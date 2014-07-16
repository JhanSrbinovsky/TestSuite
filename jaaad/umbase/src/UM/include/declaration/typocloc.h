! TYPOCLOC start
!
! Description:
!    Define variables which are local to a particular block of rows
!
! Current Code Owner: S. Foreman
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   3.4   01.09.94  Original code. R. Hill
!   5.3   10.10.01  Include ENGINT for ocean restructuring. R. Hill
!   6.0   22/08/03  Minor update for SX6 porting. R. Hill
! ====================================================================
      REAL :: ENGINT(8) ! OUT  Internal energy components
      REAL :: DTABS(NT) ! OUT  Absolute change in tracers
      REAL :: TVAR(NT)  ! OUT  Variance of tracer change
      REAL :: TTDTOT(6,NT) ! OUT Integral diagnostics of tracers

      INTEGER :: KMT(IMT)  !      No of levels at T points
      INTEGER :: KMTP(IMT) !      and at row to north
      INTEGER :: KMU(IMT)  !      No of levels at U points
      INTEGER :: KMUP(IMT) !      and at row to north
      INTEGER :: KMTPP(IMT)
      INTEGER :: KMUPP(IMT_BIH)
      INTEGER :: KMTJM(IMT)      ! No of levels at T points at row J-1
! TYPOCLOC end
