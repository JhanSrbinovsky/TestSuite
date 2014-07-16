
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_Vector_Wind

      Subroutine Polar_Vector_Wind_n(                                   &
     &                             v_comp, sin_longitude,               &
     &                             cos_longitude, row_length,           &
     &                             rows, model_levels,                  &
     &                             mag_vector_np, dir_vector_np,        &
     &                             mag_vector_sp, dir_vector_sp,        &
     &                             halo_i, halo_j, g_row_length,        &
     &                             proc_row_group, at_extremity)



! Purpose:
!          Calculates vector wind at each pole and returns magnitude
!          and direction.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
!     5.3  19/10/01   Use appropriate gcg routines.   S. Cusack
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.


      Integer                                                           &
     &  row_length                                                      &
                       ! Dimension of v component field in i direction.
     &, rows                                                            &
                       ! Dimension of v component field in j direction
                       ! plus one extra row.
     &, model_levels                                                    &
                       ! Dimension of v component field in k direction.
     &, halo_i                                                          &
                       ! Size of halo in i direction.
     &, halo_j                                                          &
                       ! Size of halo in j direction.
     &, g_row_length                                                    &
                       ! global number of points on a row
     &, proc_row_group ! Group id for processors on the same row

      Real                                                              &
     &  v_comp (1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,        &
     &          model_levels)                                           &
                                                 ! v component of wind
                                                 ! on C grid.
     &, cos_longitude(row_length,rows)                                  &
                                       ! cosine of longitude
     &, sin_longitude(row_length,rows) ! sine of longitude

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  mag_vector_np(model_levels)                                     &
                                    ! mag. of the vector wind at N pole
     &, dir_vector_np(model_levels)                                     &
                                    ! dirn. of the vector wind at N pole
     &, mag_vector_sp(model_levels)                                     &
                                    ! mag. of the vector wind at S pole
     &, dir_vector_sp(model_levels) ! dirn. of the vector wind at S pole

! Local Variables.

      Integer                                                           &
     &  i, k                                                            &
                 ! Loop indices
     &, info  ! Status variable

      Real                                                              &
     &  a_np                                                            &
     &, b_np                                                            &
     &, a_sp                                                            &
     &, b_sp

      Real                                                              &
     &  wrk(2*model_levels), rwrk(row_length,2*model_levels)

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
     &   NoDomain      ! Value in neighbor array if the domain has
                       !  no neighbor in this direction. Otherwise
                       !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! External Routines: None.

! ----------------------------------------------------------------------
! Section 1.   Calculate magnitude and direction of polar vector wind
!              from v component of wind on row around pole.
! ----------------------------------------------------------------------

      If (at_extremity(PNorth)) then

         Do k = 1, 2*model_levels
            wrk(k) = 0.
         End Do

!         Do k = 1, model_levels
!            Do i = 1, row_length
!               wrk(2*(k-1)+1) = wrk(2*(k-1)+1)
!     &              + v_comp(i,rows,k) * cos_longitude(i,rows)
!               wrk(2*k) = wrk(2*k)
!     &              + v_comp(i,rows,k) * sin_longitude(i,rows)
!            End Do
!         End Do
!
!         call gcg_rsum(2*model_levels, proc_row_group, info, wrk)


         Do k = 1, model_levels
            Do i = 1, row_length
               rwrk(i,2*(k-1)+1) = v_comp(i,rows,k) *                   &
     &              cos_longitude(i,rows)
               rwrk(i,2*k) = v_comp(i,rows,k) *                         &
     &              sin_longitude(i,rows)
            End Do
         End Do


         call gcg_rvecsumr(row_length, row_length, 1, 2*model_levels,   &
     &        rwrk, proc_row_group, info, wrk)





         Do k = 1, model_levels

            a_np = 2. * wrk(2*(k-1)+1) / g_row_length
            b_np = 2. * wrk(2*k) / g_row_length

            mag_vector_np(k) = sqrt (a_np*a_np + b_np*b_np)
            If (a_np  ==  0. .and. b_np  ==  0.) Then
               dir_vector_np(k) = 0.
            Else
               dir_vector_np(k) = atan2 (b_np, a_np)
            End If

         End Do

!         write(*,100) a_np, b_np, dir_vector_np

      End If

      If (at_extremity(PSouth)) then

         Do k = 1, 2*model_levels
            wrk(k) = 0.
         End Do

!         Do k = 1, model_levels
!            Do i = 1, row_length
!               wrk(2*(k-1)+1) = wrk(2*(k-1)+1)
!     &              + v_comp(i,1,k) * cos_longitude(i,1)
!               wrk(2*k) = wrk(2*k)
!     &              + v_comp(i,1,k) * sin_longitude(i,1)
!            End Do
!         End Do
!
!         call gcg_rsum(2*model_levels, proc_row_group, info, wrk)

         Do k = 1, model_levels
            Do i = 1, row_length
               rwrk(i,2*(k-1)+1) = v_comp(i,1,k) * cos_longitude(i,1)
               rwrk(i,2*k) = v_comp(i,1,k) * sin_longitude(i,1)
            End Do
         End Do


         call gcg_rvecsumr(row_length, row_length, 1, 2*model_levels,   &
     &        rwrk, proc_row_group, info, wrk)





         Do k = 1, model_levels

            a_sp = 2. * wrk(2*(k-1)+1) / g_row_length
            b_sp = 2. * wrk(2*k) / g_row_length

            mag_vector_sp(k) = sqrt (a_sp*a_sp + b_sp*b_sp)
            If (a_sp  ==  0. .and. b_sp  ==  0.) Then
               dir_vector_sp(k) = 0.
            Else
               dir_vector_sp(k) = atan2 (b_sp, a_sp)
            End If

         End Do

!         write(*,100) a_sp, b_sp, dir_vector_sp

      End If

 100  format( 3z20 )

! End of routine.
      return
      END SUBROUTINE Polar_Vector_Wind_n

