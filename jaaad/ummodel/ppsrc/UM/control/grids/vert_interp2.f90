
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine vert_interp2

        Subroutine vert_interp2(                                        &
     &                           data_in, row_length, data_rows,        &
     &                           data_levels, desired_p,                &
     &                           halo_x1, halo_y1,                      &
     &                           halo_x2, halo_y2,                      &
     &                           p_at_data, interp_order,               &
     &                           data_out )

! Purpose:
!          Performs vertical cubic interpolation of a field to a
!          desired p surface given the value of p at each of
!          the data points. Where the desired surface is below/above
!          the bottom/top data point the lowest/highest level data
!          value is returned.
!          p can be any quantity that is monotonic decreasing with
!          model level, such as pressure, or Exner pressure.
!          It is usual, and recommended, that the p surface be exner.
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
! 19/03/01 5.2      Removed mdi for points above top level. C. Wilson
! 05/12/05 6.2      NEC optimisation. JC Rioual/P Selwood
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
     &, halo_x2                                                         &
     &, halo_y1                                                         &
     &, halo_y2

      Real                                                              &
     &  desired_p       ! desired value to which data should be
                        ! interpolated to.

      Real                                                              &
     &  data_in (1-halo_x1:row_length+halo_x1,                          &
     &           1-halo_y1:data_rows+halo_y1, data_levels)              &
     &, p_at_data (1-halo_x2:row_length+halo_x2,                        &
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
     &  p_here                                                          &
     &, p_here_plus                                                     &
     &, p_here_plus2                                                    &
     &, p_here_minus                                                    &
     &, p_here_plus3                                                    &
     &, p_here_minus2

! ----------------------------------------------------------------------
! Section 1. Find level below which desired surface is.
!   (if requested level is above top of model, set to -1, which will
!    use value from top level)
! ----------------------------------------------------------------------

      Do j = 1, data_rows
        Do i = 1, row_length
          level_below(i,j) = -1
        End Do
      End Do


      Do k = 1, data_levels
        Do j = 1, data_rows
          Do i = 1, row_length
            If ((p_at_data(i,j,k) < desired_p) .AND.                    &
     &          (level_below(i,j) == -1)) Then
              level_below(i,j) = k-1
            End If
          End Do
        End Do
      End Do


! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

      Do j = 1, data_rows
        Do i = 1, row_length

          If (level_below(i,j)  ==  -1) Then
            data_out(i,j) = data_in(i,j,data_levels)

          Else If (level_below(i,j)  ==  0) Then
            data_out(i,j) = data_in(i,j,1)

          Else if (level_below(i,j)  ==  1 .or.                         &
     &             level_below(i,j)  ==  data_levels - 1                &
     &             .or. interp_order  ==  1 ) Then
! linearly interpolate
            data_out (i,j) = ( (desired_p -                             &
     &                          p_at_data(i,j,level_below(i,j)) )       &
     &                          * data_in (i,j,level_below(i,j)+1)      &
     &                         -(desired_p -                            &
     &                           p_at_data(i,j,level_below(i,j)+1)) *   &
     &                           data_in (i,j,level_below(i,j)) ) /     &
     &                        ( p_at_data(i,j,level_below(i,j)+1) -     &
     &                          p_at_data(i,j,level_below(i,j)) )

          Else if (level_below(i,j)  ==  2 .or.                         &
     &             level_below(i,j)  ==  data_levels - 2                &
     &             .or. interp_order  ==  3 ) Then
! cubicly interpolate

            p_here_minus = p_at_data(i,j,level_below(i,j)-1)
            p_here = p_at_data(i,j,level_below(i,j))
            p_here_plus = p_at_data(i,j,level_below(i,j)+1)
            p_here_plus2 = p_at_data(i,j,level_below(i,j)+2)

            data_out (i,j) = ( (desired_p - p_here) *                   &
     &                                (desired_p - p_here_plus )*       &
     &                                (desired_p - p_here_plus2 ) ) /   &
     &                              ( (p_here_minus - p_here) *         &
     &                                (p_here_minus - p_here_plus )*    &
     &                                (p_here_minus - p_here_plus2 ) ) *&
     &                              data_in (i,j,level_below(i,j)-1) +  &
     &                              ( (desired_p - p_here_minus) *      &
     &                                (desired_p - p_here_plus )*       &
     &                                (desired_p - p_here_plus2 ) ) /   &
     &                              ( (p_here - p_here_minus) *         &
     &                                (p_here - p_here_plus )*          &
     &                                (p_here - p_here_plus2 ) ) *      &
     &                              data_in (i,j,level_below(i,j)) +    &
     &                              ( (desired_p - p_here_minus) *      &
     &                                (desired_p - p_here )*            &
     &                                (desired_p - p_here_plus2 ) ) /   &
     &                              ( (p_here_plus - p_here_minus) *    &
     &                                (p_here_plus - p_here )*          &
     &                                (p_here_plus - p_here_plus2 ) ) * &
     &                              data_in (i,j,level_below(i,j)+1) +  &
     &                              ( (desired_p - p_here_minus) *      &
     &                                (desired_p - p_here )*            &
     &                                (desired_p - p_here_plus ) ) /    &
     &                              ( (p_here_plus2 - p_here_minus) *   &
     &                                (p_here_plus2 - p_here )*         &
     &                                (p_here_plus2 - p_here_plus ) ) * &
     &                              data_in (i,j,level_below(i,j)+2)


          Else
! interpolate quinticly

            p_here_minus2 = p_at_data(i,j,level_below(i,j)-2)
            p_here_minus = p_at_data(i,j,level_below(i,j)-1)
            p_here = p_at_data(i,j,level_below(i,j))
            p_here_plus = p_at_data(i,j,level_below(i,j)+1)
            p_here_plus2 = p_at_data(i,j,level_below(i,j)+2)
            p_here_plus3 = p_at_data(i,j,level_below(i,j)+3)

            Data_out (i,j) = ((desired_p - p_here_minus) *              &
     &                              (desired_p - p_here )*              &
     &                              (desired_p - p_here_plus )*         &
     &                              (desired_p - p_here_plus2 )*        &
     &                              (desired_p - p_here_plus3 ))/       &
     &                            ( (p_here_minus2 - p_here_minus) *    &
     &                              (p_here_minus2 - p_here )*          &
     &                              (p_here_minus2 - p_here_plus )*     &
     &                              (p_here_minus2 - p_here_plus2 )*    &
     &                              (p_here_minus2 - p_here_plus3 ) ) * &
     &                              data_in (i,j,level_below(i,j)-2) +  &
     &                            ((desired_p - p_here_minus2) *        &
     &                              (desired_p - p_here )*              &
     &                              (desired_p - p_here_plus )*         &
     &                              (desired_p - p_here_plus2 )*        &
     &                              (desired_p - p_here_plus3 ))/       &
     &                            ( (p_here_minus - p_here_minus2) *    &
     &                              (p_here_minus - p_here )*           &
     &                              (p_here_minus - p_here_plus )*      &
     &                              (p_here_minus - p_here_plus2 )*     &
     &                              (p_here_minus - p_here_plus3 ) ) *  &
     &                              data_in (i,j,level_below(i,j)-1) +  &
     &                            ((desired_p - p_here_minus2) *        &
     &                              (desired_p - p_here_minus )*        &
     &                              (desired_p - p_here_plus )*         &
     &                              (desired_p - p_here_plus2 )*        &
     &                              (desired_p - p_here_plus3 ))/       &
     &                            ( (p_here - p_here_minus2) *          &
     &                              (p_here - p_here_minus )*           &
     &                              (p_here - p_here_plus )*            &
     &                              (p_here - p_here_plus2 )*           &
     &                              (p_here - p_here_plus3 ) ) *        &
     &                              data_in (i,j,level_below(i,j)) +    &
     &                            ((desired_p - p_here_minus2) *        &
     &                              (desired_p - p_here_minus )*        &
     &                              (desired_p - p_here )*              &
     &                              (desired_p - p_here_plus2 )*        &
     &                              (desired_p - p_here_plus3 ))/       &
     &                            ( (p_here_plus - p_here_minus2) *     &
     &                              (p_here_plus - p_here_minus )*      &
     &                              (p_here_plus - p_here )*            &
     &                              (p_here_plus - p_here_plus2 )*      &
     &                              (p_here_plus - p_here_plus3 ) ) *   &
     &                              data_in (i,j,level_below(i,j)+1) +  &
     &                            ((desired_p - p_here_minus2) *        &
     &                              (desired_p - p_here_minus )*        &
     &                              (desired_p - p_here )*              &
     &                              (desired_p - p_here_plus )*         &
     &                              (desired_p - p_here_plus3 ))/       &
     &                            ( (p_here_plus2 - p_here_minus2) *    &
     &                              (p_here_plus2 - p_here_minus )*     &
     &                              (p_here_plus2 - p_here )*           &
     &                              (p_here_plus2 - p_here_plus )*      &
     &                              (p_here_plus2 - p_here_plus3 ) ) *  &
     &                              data_in (i,j,level_below(i,j)+2) +  &
     &                            ((desired_p - p_here_minus2) *        &
     &                              (desired_p - p_here_minus )*        &
     &                              (desired_p - p_here )*              &
     &                              (desired_p - p_here_plus )*         &
     &                              (desired_p - p_here_plus2 ))/       &
     &                            ( (p_here_plus3 - p_here_minus2) *    &
     &                              (p_here_plus3 - p_here_minus )*     &
     &                              (p_here_plus3 - p_here )*           &
     &                              (p_here_plus3 - p_here_plus )*      &
     &                              (p_here_plus3 - p_here_plus2 ) ) *  &
     &                              data_in (i,j,level_below(i,j)+3)

          End If

        End Do
      End Do

! end of routine

      Return
      END SUBROUTINE vert_interp2

