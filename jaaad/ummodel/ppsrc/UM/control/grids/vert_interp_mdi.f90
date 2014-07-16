
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine vert_interp_mdi

      Subroutine vert_interp_mdi(                                       &
     &                           data_in, row_length, data_rows,        &
     &                           data_levels, desired_r,                &
     &                           halo_x1, halo_y1,                      &
     &                           halo_x2, halo_y2,                      &
     &                           r_at_data, interp_order,               &
     &                           mdi, data_out )

! Purpose:
!          Performs vertical cubic interpolation of a field to a
!          desired r surface given the r value at each of
!          the data points. Where the desired surface is below/above
!          the bottom/top data point a missing data indicator is
!          returned.
!          r can be any quantity that is monotonic increasing with
!          model level, such as height, or potential temperature
!          (assuming that the model is statically stable).
!
! Method:
!          Cubic Lagrange interpolation.
!
! Original Programmer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, data_rows                                                       &
                        ! number of rows of data
     &, data_levels                                                     &
                        ! number of levels of data
     &, interp_order                                                    &
                        ! 1 = linear, 3 = cubic, 5=quintic
     &, halo_x1                                                         &
     &, halo_y1                                                         &
     &, halo_x2                                                         &
     &, halo_y2

      Real                                                              &
     &  desired_r                                                       &
                        ! desired value to which data should be
                        ! interpolated to.
     &, mdi             ! missing data indicator.

      Real                                                              &
     &  data_in (1-halo_x1:row_length+halo_x1,                          &
     &           1-halo_y1:data_rows+halo_y1, data_levels)              &
     &, r_at_data (1-halo_x2:row_length+halo_x2,                        &
     &           1-halo_y2:data_rows+halo_y2, data_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  data_out (row_length, data_rows)

! Local variables

      Integer                                                           &
     & i,j,k

      Integer                                                           &
     &  level_below(row_length, data_rows)

      Real                                                              &
     &  r_here                                                          &
     &, r_here_plus                                                     &
     &, r_here_plus2                                                    &
     &, r_here_minus                                                    &
     &, r_here_plus3                                                    &
     &, r_here_minus2

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
! ----------------------------------------------------------------------

      Do j = 1, data_rows
        Do i = 1, row_length
          level_below(i,j) = 0
        End Do
      End Do

      Do k = 1, data_levels - 1
        Do j = 1, data_rows
          Do i = 1, row_length
            If ( r_at_data(i,j,k)  <=  desired_r ) Then
              level_below(i,j) = k
            End If
          End Do
        End Do
      End Do

! if requested level is above top of model, set to zero, which will
! be converted to missing data indicator.

      Do j = 1, data_rows
        Do i = 1, row_length
          If ( desired_r  >   r_at_data(i,j,data_levels) ) Then
            level_below(i,j) = 0
          End If
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

      Do j = 1, data_rows
        Do i = 1, row_length

          If (level_below(i,j)  ==  0) Then
            data_out(i,j) = mdi

          Else if (level_below(i,j)  ==  1 .or.                         &
     &             level_below(i,j)  ==  data_levels - 1                &
     &             .or. interp_order  ==  1 ) Then
! linearly interpolate
            data_out (i,j) = ( (desired_r -                             &
     &                          r_at_data(i,j,level_below(i,j)) )       &
     &                          * data_in (i,j,level_below(i,j)+1)      &
     &                         -(desired_r -                            &
     &                           r_at_data(i,j,level_below(i,j)+1)) *   &
     &                           data_in (i,j,level_below(i,j)) ) /     &
     &                        ( r_at_data(i,j,level_below(i,j)+1) -     &
     &                          r_at_data(i,j,level_below(i,j)) )

          Else if (level_below(i,j)  ==  2 .or.                         &
     &             level_below(i,j)  ==  data_levels - 2                &
     &             .or. interp_order  ==  3 ) Then

! cubicly interpolate

            r_here_minus = r_at_data(i,j,level_below(i,j)-1)
            r_here = r_at_data(i,j,level_below(i,j))
            r_here_plus = r_at_data(i,j,level_below(i,j)+1)
            r_here_plus2 = r_at_data(i,j,level_below(i,j)+2)

            data_out (i,j) = ( (desired_r - r_here) *                   &
     &                                (desired_r - r_here_plus )*       &
     &                                (desired_r - r_here_plus2 ) ) /   &
     &                              ( (r_here_minus - r_here) *         &
     &                                (r_here_minus - r_here_plus )*    &
     &                                (r_here_minus - r_here_plus2 ) ) *&
     &                              data_in (i,j,level_below(i,j)-1) +  &
     &                              ( (desired_r - r_here_minus) *      &
     &                                (desired_r - r_here_plus )*       &
     &                                (desired_r - r_here_plus2 ) ) /   &
     &                              ( (r_here - r_here_minus) *         &
     &                                (r_here - r_here_plus )*          &
     &                                (r_here - r_here_plus2 ) ) *      &
     &                              data_in (i,j,level_below(i,j)) +    &
     &                              ( (desired_r - r_here_minus) *      &
     &                                (desired_r - r_here )*            &
     &                                (desired_r - r_here_plus2 ) ) /   &
     &                              ( (r_here_plus - r_here_minus) *    &
     &                                (r_here_plus - r_here )*          &
     &                                (r_here_plus - r_here_plus2 ) ) * &
     &                              data_in (i,j,level_below(i,j)+1) +  &
     &                              ( (desired_r - r_here_minus) *      &
     &                                (desired_r - r_here )*            &
     &                                (desired_r - r_here_plus ) ) /    &
     &                              ( (r_here_plus2 - r_here_minus) *   &
     &                                (r_here_plus2 - r_here )*         &
     &                                (r_here_plus2 - r_here_plus ) ) * &
     &                              data_in (i,j,level_below(i,j)+2)

          Else
! interpolate quinticly

            r_here_minus2 = r_at_data(i,j,level_below(i,j)-2)
            r_here_minus = r_at_data(i,j,level_below(i,j)-1)
            r_here = r_at_data(i,j,level_below(i,j))
            r_here_plus = r_at_data(i,j,level_below(i,j)+1)
            r_here_plus2 = r_at_data(i,j,level_below(i,j)+2)
            r_here_plus3 = r_at_data(i,j,level_below(i,j)+3)

            Data_out (i,j) = ((desired_r - r_here_minus) *              &
     &                              (desired_r - r_here )*              &
     &                              (desired_r - r_here_plus )*         &
     &                              (desired_r - r_here_plus2 )*        &
     &                              (desired_r - r_here_plus3 ))/       &
     &                            ( (r_here_minus2 - r_here_minus) *    &
     &                              (r_here_minus2 - r_here )*          &
     &                              (r_here_minus2 - r_here_plus )*     &
     &                              (r_here_minus2 - r_here_plus2 )*    &
     &                              (r_here_minus2 - r_here_plus3 ) ) * &
     &                              data_in (i,j,level_below(i,j)-2) +  &
     &                            ((desired_r - r_here_minus2) *        &
     &                              (desired_r - r_here )*              &
     &                              (desired_r - r_here_plus )*         &
     &                              (desired_r - r_here_plus2 )*        &
     &                              (desired_r - r_here_plus3 ))/       &
     &                            ( (r_here_minus - r_here_minus2) *    &
     &                              (r_here_minus - r_here )*           &
     &                              (r_here_minus - r_here_plus )*      &
     &                              (r_here_minus - r_here_plus2 )*     &
     &                              (r_here_minus - r_here_plus3 ) ) *  &
     &                              data_in (i,j,level_below(i,j)-1) +  &
     &                            ((desired_r - r_here_minus2) *        &
     &                              (desired_r - r_here_minus )*        &
     &                              (desired_r - r_here_plus )*         &
     &                              (desired_r - r_here_plus2 )*        &
     &                              (desired_r - r_here_plus3 ))/       &
     &                            ( (r_here - r_here_minus2) *          &
     &                              (r_here - r_here_minus )*           &
     &                              (r_here - r_here_plus )*            &
     &                              (r_here - r_here_plus2 )*           &
     &                              (r_here - r_here_plus3 ) ) *        &
     &                              data_in (i,j,level_below(i,j)) +    &
     &                            ((desired_r - r_here_minus2) *        &
     &                              (desired_r - r_here_minus )*        &
     &                              (desired_r - r_here )*              &
     &                              (desired_r - r_here_plus2 )*        &
     &                              (desired_r - r_here_plus3 ))/       &
     &                            ( (r_here_plus - r_here_minus2) *     &
     &                              (r_here_plus - r_here_minus )*      &
     &                              (r_here_plus - r_here )*            &
     &                              (r_here_plus - r_here_plus2 )*      &
     &                              (r_here_plus - r_here_plus3 ) ) *   &
     &                              data_in (i,j,level_below(i,j)+1) +  &
     &                            ((desired_r - r_here_minus2) *        &
     &                              (desired_r - r_here_minus )*        &
     &                              (desired_r - r_here )*              &
     &                              (desired_r - r_here_plus )*         &
     &                              (desired_r - r_here_plus3 ))/       &
     &                            ( (r_here_plus2 - r_here_minus2) *    &
     &                              (r_here_plus2 - r_here_minus )*     &
     &                              (r_here_plus2 - r_here )*           &
     &                              (r_here_plus2 - r_here_plus )*      &
     &                              (r_here_plus2 - r_here_plus3 ) ) *  &
     &                              data_in (i,j,level_below(i,j)+2) +  &
     &                            ((desired_r - r_here_minus2) *        &
     &                              (desired_r - r_here_minus )*        &
     &                              (desired_r - r_here )*              &
     &                              (desired_r - r_here_plus )*         &
     &                              (desired_r - r_here_plus2 ))/       &
     &                            ( (r_here_plus3 - r_here_minus2) *    &
     &                              (r_here_plus3 - r_here_minus )*     &
     &                              (r_here_plus3 - r_here )*           &
     &                              (r_here_plus3 - r_here_plus )*      &
     &                              (r_here_plus3 - r_here_plus2 ) ) *  &
     &                              data_in (i,j,level_below(i,j)+3)

          End If

        End Do
      End Do

! end of routine

      Return
      END SUBROUTINE vert_interp_mdi

