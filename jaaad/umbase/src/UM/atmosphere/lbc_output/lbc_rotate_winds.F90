#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Rotates LBC winds to rotated grid
!
! Subroutine Interface:

      Subroutine LBC_Rotate_Winds (                                     &
     &           u_p                                                    &
     &,          v_p                                                    &
     &,          lbc_size_p                                             &
     &,          lbc_levels                                             &
     &,          coeff1                                                 &
     &,          coeff2                                                 &
     & )

      IMPLICIT NONE
!
! Description:
!   Rotate LBC wind components for a rotated grid.
!
! Method:
!   Rotation Coefficients calculated in LBC_Interp_Coeffs are used to
!   transform u and v to rotated u and v.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:

! Subroutine arguments

      Integer :: lbc_size_p                      !  size of lbc grid
      Integer :: lbc_levels                      !  no of lbc levels

      Real    :: u_p ( lbc_size_p, lbc_levels )  ! u lbcs on p grid
      Real    :: v_p ( lbc_size_p, lbc_levels )  ! v lbcs on p grid
      Real    :: coeff1 (lbc_size_p)             ! \ rotation
      Real    :: coeff2 (lbc_size_p)             ! / coefficients

! Local scalars:

      Integer :: i      ! loop index for points
      Integer :: level  ! loop index for levels
      Real    :: u, v   ! u and v component

!- End of header

! ---------------------------------
! Perform the rotation of the winds
! ---------------------------------

!     w_lltoeq not used here.
!     use loop from w_lltoeq and modify to reduce arrays required.

      Do level = 1, lbc_levels

        Do i = 1, lbc_size_p

          u = u_p(i,level)
          v = v_p(i,level)

          u_p(i,level) = Coeff1(i) * u - Coeff2(i) * v
          v_p(i,level) = Coeff1(i) * v + Coeff2(i) * u

        Enddo

      Enddo

!     u_p and v_p now contain winds for the rotated grid

      Return
      END SUBROUTINE LBC_Rotate_Winds
#endif
