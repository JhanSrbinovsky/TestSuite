! ----------------------- include file: INTERPOR -----------------------
! Description: Hold magic numbers for order of vertical interpolation,
!              where prime use is for atmosphere physics/dynamics
!              diagnostics.
!              interp_order should be set to a _default value where
!              used generally in the code, but can be set explicitly
!              locally where required.
!
! Current Code Owner: R. Rawlins
! History:
! Version  Date      Comment.
!  5.0 19/06/99  Original version. R. Rawlins.

      INTEGER,PARAMETER:: interp_order_linear  = 1
      INTEGER,PARAMETER:: interp_order_cubic   = 3
      INTEGER,PARAMETER:: interp_order_quintic = 5
      INTEGER,PARAMETER:: interp_order_default = interp_order_cubic

! INTERPOR end
