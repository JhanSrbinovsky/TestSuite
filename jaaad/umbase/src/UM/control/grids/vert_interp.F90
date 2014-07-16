#if defined(RECON) || defined(C92_2A) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:
      Subroutine vert_interp                                            &
     &                      (data_in, data_points,                      &
     &                       data_levels, desired_r,                    &
     &                       r_at_data, interp_order,                   &
     &                       data_out )

! Purpose:
!          Performs vertical interpolation of a field to a
!          desired r surface given the r value at each of
!          the data points. Where the desired surface is below/above
!          the top/bottom point, then a value is assigned from
!          the top/bottom data point, except when
!          interp_order=2 is chosen, when linear extrapolation is
!          performed.
!          r can be any quantity that is monotonic increasing with
!          model level, such as height, or potential temperature
!          (assuming that the model is statically stable).
!
! Method:
!          Cubic Lagrange interpolation.
!
! Original Progammer: Mark H. Mawson
! Current code owner: I Edmond
!
! History:
! Date     Version     Comment
! ----     -------     -------
!
! 22/5/96     4.1     Linear extrapolation added to routine to
!                     obtain data at points above/below
!                     top/bottom levels.                Ian Edmond
! 29/07/98    4.5     Optimisation changes for T3E
!                     Author D.M. Goddard
! 04/05/00    5.1     Addition of interp_order == 2 option and bug-fix
!                     for equal top levels. P.Selwood.
! 10/08/00    5.2     Addition of data copy for equal input and
!                     output heights.  P. Selwood.
! 08/09/03    6.0     Added new def for use with makebc. R.Sempers
! 06/12/05    6.2     Simplify loop structure and remove goto to
!                     allow for vectorisation
!                     Jean-Christophe Rioual/R Sempers
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
! System component covered: ??
! System Task:              ??
!
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  data_points                                                     &
                        ! number of rows of data
     &, data_levels                                                     &
                        ! number of levels of data
     &, interp_order    ! 1 = linear, 3 = cubic, 5=quintic
                        ! 2 = linear with no extrapolation at top
                        !     or bottom of model.

      Real                                                              &
     &  data_in (data_points, data_levels)                              &
     &, r_at_data (data_points, data_levels)                            &
     &, desired_r(data_points) ! desired value to which
                                         ! data should be interpolated t

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  data_out (data_points)

! Local variables

      Integer                                                           &
     & j,k

      Integer last
      Integer                                                           &
     &  level_below(data_points)

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

      Do j=1, data_points
         level_below(j)=data_levels
      End Do

      Do k=1, data_levels
         Do j = 1, data_points
            If((r_at_data(j,k)  >   desired_r(j))                       &
     &               .and. (level_below(j)  ==  data_levels)) then
               level_below(j)=k-1
            End If
         End Do
      End Do

      if (interp_order  /=  2) then

      Do j = 1, data_points

      ! If requested level is above top of model, do linear
      ! extrapolation using data on top and second top levels.
          If ( desired_r(j)  >   r_at_data(j,data_levels) ) Then
            data_out(j) = data_in(j,data_levels) + (desired_r(j)        &
     &      - r_at_data(j,data_levels)) * (data_in(j,data_levels)       &
     &      - data_in(j,data_levels-1))/(r_at_data(j,data_levels)       &
     &      - r_at_data(j,data_levels-1))
          Else If (desired_r(j) == r_at_data(j,data_levels) ) Then
            data_out(j) = data_in(j,data_levels)
          End If

      ! If requested level is below bottom of model, do linear
      ! extrapolation using data on first and second levels.
          If ( desired_r(j)  <   r_at_data(j,1) ) Then
            data_out(j) = data_in(j,1) + (desired_r(j)                  &
     &      - r_at_data(j,1)) * (data_in(j,1)                           &
     &      - data_in(j,2))/(r_at_data(j,1) - r_at_data(j,2))
          Else If (desired_r(j) == r_at_data(j,1) ) Then
            data_out(j) = data_in(j,1)
          End If

      End Do

      else ! No linear extrapolation at top or bottom

        Do j = 1, data_points

!         Top : Set to top input data

          If ( desired_r(j)  >=  r_at_data(j,data_levels) ) Then
            data_out(j) = data_in(j,data_levels)
          Endif

!       Bottom : Set to bottom input data

          If ( desired_r(j)  <=  r_at_data(j,1) ) Then
            data_out(j) = data_in(j,1)
          Endif

        enddo

      endif

! ----------------------------------------------------------------------
! Section 2. Vertical interpolation.
! ----------------------------------------------------------------------

      Do j = 1, data_points

          If (level_below(j)  ==  0.or.                                 &
     &        level_below(j)  ==  data_levels) Then

             ! Output data has already been calculated - do nothing...
              data_out (j) = data_out (j)

          Else if ( abs( desired_r(j) - r_at_data(j, level_below(j)))   &
     &                         < Epsilon( desired_r(j))) Then

              ! Levels are essentially the same, copy in to out.
              data_out(j) = data_in(j,level_below(j))

          Else if (level_below(j)  ==  1 .or.                           &
     &        level_below(j)  ==  data_levels - 1                       &
     &        .or. interp_order  ==  1                                  &
     &        .or. interp_order  ==  2) Then
! linearly interpolate
            data_out (j) = ( (desired_r(j) -                            &
     &                          r_at_data(j,level_below(j)) )           &
     &                          * data_in (j,level_below(j)+1)          &
     &                         -(desired_r(j) -                         &
     &                           r_at_data(j,level_below(j)+1)) *       &
     &                           data_in (j,level_below(j)) ) /         &
     &                        ( r_at_data(j,level_below(j)+1) -         &
     &                          r_at_data(j,level_below(j)) )

          Else if (level_below(j)  ==  2 .or.                           &
     &             level_below(j)  ==  data_levels - 2                  &
     &             .or. interp_order  ==  3 ) Then

! cubicly interpolate

            r_here_minus = r_at_data(j,level_below(j)-1)
            r_here = r_at_data(j,level_below(j))
            r_here_plus = r_at_data(j,level_below(j)+1)
            r_here_plus2 = r_at_data(j,level_below(j)+2)

            data_out (j) = ( (desired_r(j) - r_here) *                  &
     &                           (desired_r(j) - r_here_plus )*         &
     &                           (desired_r(j) - r_here_plus2 ) ) /     &
     &                         ( (r_here_minus - r_here) *              &
     &                           (r_here_minus - r_here_plus )*         &
     &                           (r_here_minus - r_here_plus2 ) ) *     &
     &                         data_in (j,level_below(j)-1) +           &
     &                         ( (desired_r(j) - r_here_minus) *        &
     &                           (desired_r(j) - r_here_plus )*         &
     &                           (desired_r(j) - r_here_plus2 ) ) /     &
     &                         ( (r_here - r_here_minus) *              &
     &                           (r_here - r_here_plus )*               &
     &                           (r_here - r_here_plus2 ) ) *           &
     &                         data_in (j,level_below(j)) +             &
     &                         ( (desired_r(j) - r_here_minus) *        &
     &                           (desired_r(j) - r_here )*              &
     &                           (desired_r(j) - r_here_plus2 ) ) /     &
     &                         ( (r_here_plus - r_here_minus) *         &
     &                           (r_here_plus - r_here )*               &
     &                           (r_here_plus - r_here_plus2 ) ) *      &
     &                         data_in (j,level_below(j)+1) +           &
     &                         ( (desired_r(j) - r_here_minus) *        &
     &                           (desired_r(j) - r_here )*              &
     &                           (desired_r(j) - r_here_plus ) ) /      &
     &                         ( (r_here_plus2 - r_here_minus) *        &
     &                           (r_here_plus2 - r_here )*              &
     &                           (r_here_plus2 - r_here_plus ) ) *      &
     &                         data_in (j,level_below(j)+2)

          Else
! interpolate quinticly

            r_here_minus2 = r_at_data(j,level_below(j)-2)
            r_here_minus = r_at_data(j,level_below(j)-1)
            r_here = r_at_data(j,level_below(j))
            r_here_plus = r_at_data(j,level_below(j)+1)
            r_here_plus2 = r_at_data(j,level_below(j)+2)
            r_here_plus3 = r_at_data(j,level_below(j)+3)

            Data_out (j) = ((desired_r(j) - r_here_minus) *             &
     &                         (desired_r(j) - r_here )*                &
     &                         (desired_r(j) - r_here_plus )*           &
     &                         (desired_r(j) - r_here_plus2 )*          &
     &                         (desired_r(j) - r_here_plus3 ))/         &
     &                       ( (r_here_minus2 - r_here_minus) *         &
     &                         (r_here_minus2 - r_here )*               &
     &                         (r_here_minus2 - r_here_plus )*          &
     &                         (r_here_minus2 - r_here_plus2 )*         &
     &                         (r_here_minus2 - r_here_plus3 ) ) *      &
     &                         data_in (j,level_below(j)-2) +           &
     &                       ((desired_r(j) - r_here_minus2) *          &
     &                         (desired_r(j) - r_here )*                &
     &                         (desired_r(j) - r_here_plus )*           &
     &                         (desired_r(j) - r_here_plus2 )*          &
     &                         (desired_r(j) - r_here_plus3 ))/         &
     &                       ( (r_here_minus - r_here_minus2) *         &
     &                         (r_here_minus - r_here )*                &
     &                         (r_here_minus - r_here_plus )*           &
     &                         (r_here_minus - r_here_plus2 )*          &
     &                         (r_here_minus - r_here_plus3 ) ) *       &
     &                         data_in (j,level_below(j)-1) +           &
     &                       ((desired_r(j) - r_here_minus2) *          &
     &                         (desired_r(j) - r_here_minus )*          &
     &                         (desired_r(j) - r_here_plus )*           &
     &                         (desired_r(j) - r_here_plus2 )*          &
     &                         (desired_r(j) - r_here_plus3 ))/         &
     &                       ( (r_here - r_here_minus2) *               &
     &                         (r_here - r_here_minus )*                &
     &                         (r_here - r_here_plus )*                 &
     &                         (r_here - r_here_plus2 )*                &
     &                         (r_here - r_here_plus3 ) ) *             &
     &                         data_in (j,level_below(j)) +             &
     &                       ((desired_r(j) - r_here_minus2) *          &
     &                         (desired_r(j) - r_here_minus )*          &
     &                         (desired_r(j) - r_here )*                &
     &                         (desired_r(j) - r_here_plus2 )*          &
     &                         (desired_r(j) - r_here_plus3 ))/         &
     &                       ( (r_here_plus - r_here_minus2) *          &
     &                         (r_here_plus - r_here_minus )*           &
     &                         (r_here_plus - r_here )*                 &
     &                         (r_here_plus - r_here_plus2 )*           &
     &                         (r_here_plus - r_here_plus3 ) ) *        &
     &                         data_in (j,level_below(j)+1) +           &
     &                       ((desired_r(j) - r_here_minus2) *          &
     &                         (desired_r(j) - r_here_minus )*          &
     &                         (desired_r(j) - r_here )*                &
     &                         (desired_r(j) - r_here_plus )*           &
     &                         (desired_r(j) - r_here_plus3 ))/         &
     &                       ( (r_here_plus2 - r_here_minus2) *         &
     &                         (r_here_plus2 - r_here_minus )*          &
     &                         (r_here_plus2 - r_here )*                &
     &                         (r_here_plus2 - r_here_plus )*           &
     &                         (r_here_plus2 - r_here_plus3 ) ) *       &
     &                         data_in (j,level_below(j)+2) +           &
     &                       ((desired_r(j) - r_here_minus2) *          &
     &                         (desired_r(j) - r_here_minus )*          &
     &                         (desired_r(j) - r_here )*                &
     &                         (desired_r(j) - r_here_plus )*           &
     &                         (desired_r(j) - r_here_plus2 ))/         &
     &                       ( (r_here_plus3 - r_here_minus2) *         &
     &                         (r_here_plus3 - r_here_minus )*          &
     &                         (r_here_plus3 - r_here )*                &
     &                         (r_here_plus3 - r_here_plus )*           &
     &                         (r_here_plus3 - r_here_plus2 ) ) *       &
     &                         data_in (j,level_below(j)+3)

          End If

      End Do

! end of routine

      Return
      END SUBROUTINE vert_interp
#endif
