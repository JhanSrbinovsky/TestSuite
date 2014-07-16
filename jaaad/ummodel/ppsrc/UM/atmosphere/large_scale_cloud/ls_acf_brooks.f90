
!+    Brooks area cloud fraction parametrisation
!
!     *********************** COPYRIGHT *************************
!     Crown Copyright <year>, The Met. Office. All rights reserved.
!     *********************** COPYRIGHT *************************
!
!     Subroutine Interface:

      SUBROUTINE LS_ACF_Brooks (                                        &

! Parallel variables
     &  halo_i, halo_j, off_x, off_y                                    &

! model dimensions
     &, row_length, rows, model_levels, wet_model_levels                &

! in coordinate information
     &, r_theta_levels, delta_lambda, delta_phi                         &

! trig arrays
     &, FV_cos_theta_latitude                                           &

! in data fields
     &, bulk_cloud_fraction, cloud_fraction_liquid                      &
     &, cloud_fraction_frozen                                           &

! in logical control
     &, cumulus                                                         &

! out data fields
     &, area_cloud_fraction)

!     Description:
!       Calculates area_cloud_fraction from bulk_cloud_fraction
!
!     Method:
!       The calculation is  based on the parametrisation in
!       Brooks 2005 equations 2-3.
!       (Brooks et al, July 2005, JAS vol 62 pp 2248-2260)
!       The initial parametrisation uses the values in equations
!       4,5,7 and 8 of the paper for ice and liquid cloud without
!       wind shear.  For mixed phase clouds the maximum of the two
!       area_cloud_fractions resulting will be used.
!       Only area_cloud_fraction will be updated.
!       Grid box size is needed to be known.
!
!     Current Code Owner: LS cloud scheme owner
!
!     History:
!     Version  Date      Comment
!
!       6.4    19/09/05  Original code. (Amanda Kerr-Munslow, placed
!                                        into UM by D Wilson)
!
!     Code Description:
!       FORTRAN 77 with extensions recommended in the Met. Office
!       F77 Standard.
!
      IMPLICIT NONE

! Scalar arguments with INTENT(IN):
! Parallel setup variables

      Integer                                                           &
     &  halo_i                                                          &
                          ! Size of halo in i direction
     &, halo_j                                                          &
                          ! Size of halo in j direction
     &, off_x                                                           &
                          ! Size of small halo in i direction
     &, off_y             ! Size of small halo in j direction

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
                          ! number of points on a row
     &, rows                                                            &
                          ! number of rows in a theta field
     &, model_levels                                                    &
                          ! number of model levels
     &, wet_model_levels
                          ! number of model levels where moisture
                          ! variables are held

! Array Arguments with INTENT(IN)
! Co-ordinate arrays:
       Real                                                             &
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                   1-halo_j:rows+halo_j,0:model_levels)           &
                         ! height of theta levels (from centre of earth)
     &, FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
                         ! Finite volume cos(lat)
     &, delta_lambda                                                    &
                         ! EW (x) grid spacing in radians
     &, delta_phi        ! NS (y) grid spacing in radians

! Data arrays:
       Real                                                             &
     &  bulk_cloud_fraction (row_length, rows, wet_model_levels)        &
!       Cloud fraction at processed levels (decimal fraction).
     &, cloud_fraction_liquid (row_length, rows, wet_model_levels)      &
!       Liquid cloud fraction at processed levels (decimal fraction).
     &, cloud_fraction_frozen (row_length, rows, wet_model_levels)
!       Frozen cloud fraction at processed levels (decimal fraction).

       Logical, intent(in)::                                            &
     &  cumulus(row_length,rows)

! Arguments with INTENT(OUT):
! Data arrays:
       Real                                                             &
     &  area_cloud_fraction (row_length, rows, wet_model_levels)
!       Cloud fraction at processed levels (decimal fraction).

! Local Parameters:

! Parameters for liquid clouds, from Brooks 2005 equations 7 and 8
      Real power_law_gradient_liquid ! A
      Parameter (power_law_gradient_liquid = 0.1635)
      Real vert_fit_liquid  ! alpha
      Parameter (vert_fit_liquid = 0.6694)
      Real horiz_fit_liquid  ! beta
      Parameter (horiz_fit_liquid = -0.1882)

! Parameters for frozen clouds, from Brooks 2005 equations 4 and 5
      Real power_law_gradient_frozen ! A
      Parameter (power_law_gradient_frozen = 0.0880)
      Real vert_fit_frozen  ! alpha
      Parameter (vert_fit_frozen = 0.7679)
      Real horiz_fit_frozen  ! beta
      Parameter (horiz_fit_frozen = -0.2254)

! Local Scalars:
! Loop counters
       Integer                                                          &
     &  i, j, k

       Real                                                             &
     &  symmetric_adjustment_liquid                                     &
!    function f in eqn 7 in Brooks 2005
     &, symmetric_adjustment_frozen                                     &
!    function f in eqn 4 in Brooks 2005
     &, horiz_scale                                                     &
!    horizontal scale size of the grid box (m)
     &, vert_scale
!    vertical scale size of the grid box (m)


!  Local Arrays:
       Real                                                             &
     &  acf_liquid (row_length, rows, wet_model_levels)                 &
!    area cloud fraction based on liquid parameters
     &, acf_frozen (row_length, rows, wet_model_levels)
!    area cloud fraction based on frozen parameters

!-    End of header
! ----------------------------------------------------------------------

! ==Main Block==--------------------------------------------------------

! Initialise arrays and local variables to zero
        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
              area_cloud_fraction(i,j,k) = 0.0
              acf_liquid(i,j,k) = 0.0
              acf_frozen(i,j,k) = 0.0
            End Do
          End Do
        End Do
        horiz_scale = 0.0
        vert_scale = 0.0
        symmetric_adjustment_liquid = 0.0
        symmetric_adjustment_frozen = 0.0

        Do k = 1, wet_model_levels
          Do j = 1, rows
            Do i = 1, row_length
!
! Test if bulk_cloud_fraction is within bounds of possibility
!
            If ( bulk_cloud_fraction(i,j,k) <= 0.0 ) Then

              area_cloud_fraction(i,j,k) = 0.0

            Else If ( bulk_cloud_fraction(i,j,k) >= 1.0 ) Then

              area_cloud_fraction(i,j,k) = 1.0

            Else if (cumulus(i,j)) then
              ! This is a convective point so do not apply the
              ! area cloud representation
              area_cloud_fraction(i,j,k) =                              &
     &          bulk_cloud_fraction(i,j,k)

            Else

! Only calculate area_cloud_fraction if the bulk_cloud_fraction
! is between (not equal to) 0.0 and 1.0
!
! Calculate horizontal and vertical scales.
! The horizontal scale is taken as the square root of the
! area of the grid box.
! The vertical scale is taken as the difference in radius
! from the centre of the Earth between the upper and lower
! boundaries of the grid box.

              horiz_scale = SQRT (                                      &
     &                          r_theta_levels(i,j,k)                   &
     &                          * r_theta_levels(i,j,k)                 &
     &                          * delta_lambda * delta_phi              &
     &                          * FV_cos_theta_latitude(i,j) )
              If (k .eq. wet_model_levels) then
                ! Assume top layer thickness is the same as the
                ! thickness of the layer below
                vert_scale =  r_theta_levels(i,j,k)                     &
     &                       - r_theta_levels(i,j,k-1)
              Else 
                vert_scale =  r_theta_levels(i,j,k+1)                   &
     &                       - r_theta_levels(i,j,k)
              End if  ! k eq wet_model_levels

! Calculate the symmetric_adjustment (f).
! This parameter controls the extent to which the area cloud fraction
! is greater than the bulk cloud fraction.  If f = 0, they are equal.

              symmetric_adjustment_liquid =                             &
     &                   power_law_gradient_liquid                      &
     &                   * ( vert_scale ** vert_fit_liquid )            &
     &                   * ( horiz_scale ** horiz_fit_liquid )
              symmetric_adjustment_frozen =                             &
     &                   power_law_gradient_frozen                      &
     &                   * ( vert_scale ** vert_fit_frozen )            &
     &                   * ( horiz_scale ** horiz_fit_frozen )

! Calculate the area cloud fractions for liquid and frozen cloud
! Calculate the liquid and frozen fractions separately to
! allow for greatest flexibility in future choice of decisions
! regarding mixed phase cloud.

              acf_liquid(i,j,k) = 1./                                   &
     &            ( 1. + ( exp(-1.*symmetric_adjustment_liquid)         &
     &                     * ( 1./bulk_cloud_fraction(i,j,k) - 1.) ) )
              acf_frozen(i,j,k) = 1./                                   &
     &            ( 1. + ( exp(-1.*symmetric_adjustment_frozen)         &
     &                     * ( 1./bulk_cloud_fraction(i,j,k) - 1.) ) )

! Calculate the final area cloud fraction for each grid box
! Currently this is based on which there is more of, ice or liquid.

              If ( cloud_fraction_frozen(i,j,k) == 0.0 ) Then
                If ( cloud_fraction_liquid(i,j,k) == 0.0 ) Then

! If there is no liquid or frozen cloud, there should be no area cloud
                  area_cloud_fraction(i,j,k) = 0.0

                Else

! If there is no frozen cloud but there is liquid cloud,
! then the area cloud fraction is given by the liquid parametrisation
! 0 no cloud, 1 either, 2 liq, 3 ice'
                  area_cloud_fraction(i,j,k) = acf_liquid(i,j,k)
                End If

              Else ! cloud_fraction_frozen

                If ( cloud_fraction_liquid(i,j,k) == 0.0 ) Then

! If there is frozen cloud but there is no liquid cloud,
! then the area cloud fraction is given by the frozen parametrisation
                  area_cloud_fraction(i,j,k) = acf_frozen(i,j,k)

                Else

! If there is frozen cloud and there is liquid cloud,
! then the area cloud fraction is given by the maximum of the two
! parametrisations
                  area_cloud_fraction(i,j,k) =                          &
     &               MAX( acf_liquid(i,j,k),acf_frozen(i,j,k) )

                End If

              End If ! cloud_fraction_frozen

            End If ! bulk_cloud_fraction between 0.0 and 1.0

            End Do
          End Do
        End Do

      RETURN
      END SUBROUTINE LS_ACF_Brooks
! ======================================================================
