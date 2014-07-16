#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Friction_Suarez_Held

      Subroutine IDL_Friction_Suarez_Held                               &
     &                        (row_length, rows, n_rows                 &
     &,                        model_levels, timestep                   &
     &,                        model_domain                             &
     &,                        off_x, off_y, at_extremity               &
     &,                        friction_level                           &
     &,                        base_frictional_timescale                &
     &,                        SuHe_sigma_cutoff, SuHe_fric             &
     &,                        p, p_star, u, v, R_u, R_v )

! Purpose:
!          Provides temperature relaxation for Suarez Held test problem.
!
! Original Progammer: Terry Davies
! Current code owner: Andrew J Malcolm
!
! History:
! Version   Date     Comment
! ----     -------   -------
! 5.3      25/11/00  This deck introduced     Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      Implicit None

! Parallel setup variables
      Integer                                                           &
     &  off_x                                                           &
                   ! Size of small halo in i
     &, off_y      ! Size of small halo in j.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a u field
     &, n_rows                                                          &
                           ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, model_domain     ! type of domain

      Real                                                              &
     &  timestep

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &         model_levels)                                            &
     &, p_star(row_length, rows)                                        &
     &, u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &        model_levels)                                             &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &        model_levels)                                             &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &        model_levels)                                             &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)

! Suarez Held variables
      Real                                                              &
     &  friction_level(model_levels)                                    &
     &, base_frictional_timescale                                       &
     &, SuHe_sigma_cutoff

      Integer                                                           &
     &  SuHe_fric       ! Switch to choose form of friction

! local variables
      Integer                                                           &
     &  i, j, k                                                         &
     &, j0, j1                                                          &
     &, levels_h

      Real                                                              &
     &  temp

       Real                                                             &
     &  work_u_to_p(row_length, rows, model_levels)                     &
     &, work_v_to_p(row_length, rows, model_levels)                     &
     &, work_p_to_v(row_length, n_rows, model_levels)                   &
     &, work_p_halo(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &                 model_levels)

#include "fldtype.h"

! No External routines

! ----------------------------------------------------------------------
! Section 0. Initialise number of working levels levels_h
! ----------------------------------------------------------------------

!  Expect boundary -layer levels to be no more than half model_levels
      levels_h = model_levels / 2

! ----------------------------------------------------------------------
! Section 1.  Calculate friction terms
! ----------------------------------------------------------------------
      If(SuHe_fric  ==  1)then
! ----------------------------------------------------------------------
! Section 1.1 Original simple friction applied at u,v points
! ----------------------------------------------------------------------
        j0 = 1
        j1 = rows
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows - 1

! add on to increment field
        Do k = 1, levels_h
          Do j = j0, j1
            Do i = 1, row_length
              R_u(i,j,k) = R_u(i,j,k) - timestep *                      &
     &                          friction_level(k) * u(i,j,k)
            End Do
          End Do
        End Do

! add on to increment field
        Do k = 1, levels_h
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) -  timestep *                     &
     &                           friction_level(k) * v(i,j,k)
            End Do
          End Do
        End Do

      else If(SuHe_fric  ==  2)then
! ----------------------------------------------------------------------
! Section 1.2 Simple friction evaluated at p points then
!             interpolated to u,v points
! ----------------------------------------------------------------------

! DEPENDS ON: u_to_p
        Call u_to_p (u, row_length, rows, levels_h,                     &
     &                 off_x, off_y, model_domain,                      &
     &                 at_extremity, work_u_to_p)

! DEPENDS ON: v_to_p
        Call v_to_p (v, row_length, rows, n_rows, levels_h,             &
     &                 off_x, off_y, model_domain,                      &
     &                 at_extremity, work_v_to_p)

! set polar values to zero.
        If (model_domain  ==  1 .and. at_extremity(PSouth) ) Then
          Do k = 1, levels_h
            Do i = 1, row_length
              work_u_to_p(i,1,k) = 0.
              work_v_to_p(i,1,k) = 0.
            End Do
          End Do
        End If
        If (model_domain  ==  1 .and. at_extremity(PNorth) ) Then
          Do k = 1, levels_h
            Do i = 1, row_length
              work_u_to_p(i,rows,k) = 0.
              work_v_to_p(i,rows,k) = 0.
            End Do
          End Do
        End If

! Calculate k_v*u and k_v*v
        Do k = 1, levels_h
          Do j = 1, rows
            Do i = 1, row_length
              temp = (p(i,j,k)/p_star(i,j) - SuHe_sigma_cutoff) /       &
     &                          (1.0 - SuHe_sigma_cutoff)
              work_u_to_p(i,j,k) = base_frictional_timescale *          &
     &                            max(0.0, temp) * work_u_to_p(i,j,k)
              work_v_to_p(i,j,k) = base_frictional_timescale *          &
     &                            max(0.0, temp) * work_v_to_p(i,j,k)
            End Do
          End Do
        End Do

!
! copy u-tendencies into haloed arrays
        Do k = 1, levels_h
          Do j = 1, rows
            Do i = 1, row_length
              work_p_halo(i,j,k) = work_u_to_p(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(work_p_halo, row_length, rows,                 &
     &           levels_h, off_x, off_y, fld_type_p, .false.  )
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_p_halo, row_length, rows,         &
     &                                      levels_h, off_x, off_y )
! interpolate to u grid
! DEPENDS ON: p_to_u
        Call p_to_u (work_p_halo, row_length, rows, levels_h,           &
     &                 off_x, off_y, work_u_to_p)

! add on to increment field
        Do k = 1, levels_h
          Do j = 1, rows
            Do i = 1, row_length
              R_u(i,j,k) = R_u(i,j,k) - work_u_to_p(i,j,k) * timestep
            End Do
          End Do
        End Do

! copy v-tendencies into haloed arrays
        Do k = 1, levels_h
          Do j = 1, rows
            Do i = 1, row_length
              work_p_halo(i,j,k) = work_v_to_p(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(work_p_halo, row_length, rows,                 &
     &           levels_h, off_x, off_y, fld_type_p, .false.  )
! DEPENDS ON: fill_external_halos
        Call FILL_EXTERNAL_HALOS(work_p_halo, row_length, rows,         &
     &                                      levels_h, off_x, off_y )
! interpolate to v grid
! DEPENDS ON: p_to_v
        Call p_to_v (work_p_halo, row_length, rows, n_rows,             &
     &                      levels_h, off_x, off_y, work_p_to_v)

! add on to increment field
        Do k = 1, levels_h
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) - work_p_to_v(i,j,k) * timestep
            End Do
          End Do
        End Do

      else
        print*,'SuHe_fric = ',SuHe_fric,' not supported'
        print*,'Put correct value in IDEALISED NAMELIST - 1 or 2'
        stop
      end If !  SuHe_fric  ==  1

      Return
      END SUBROUTINE IDL_Friction_Suarez_Held

#endif
