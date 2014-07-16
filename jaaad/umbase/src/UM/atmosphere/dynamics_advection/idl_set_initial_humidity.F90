#if defined(A12_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_set_initial_humidity

      Subroutine idl_set_initial_humidity(                              &
     &                      p_zero, recip_Kappa                         &
     &,                     Earth_radius, height_domain                 &
     &,                     row_length, rows                            &
     &,                     model_levels, wet_model_levels              &
     &,                     me, halo_i, halo_j                          &
     &,                     qprofile_number, L_dry, q1                  &
     &,                     max_num_profile_data, num_profile_data      &
     &,                     zprofile_data, qprofile_data                &
     &,                     zprofile_orog, idl_interp_option, hf        &
     &,                     r_theta_levels, eta_theta_levels            &
     &,                     q, theta, exner_theta_levels                &
     &,                     q_ref, theta_ref, exner_ref                 &
     &,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Input from idealised namelist.
!
!
!
! Original Programmer: Richard Forbes
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! ----     -------     -----
! 6.2      10/10/01    New deck created from code in IDL_INITIAL_DATA.
!                      R. Forbes
!
! Code Description:
!   Language: FORTRAN 77 + extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent (In)

      !  Physical constants
      Real, Intent(In) :: recip_kappa         ! 1/Kappa
      Real, Intent(In) :: p_zero              ! Reference pressure (hPa)
      Real, Intent(In) :: Earth_radius        ! Earth Radius (m)
      Real, Intent(In) :: height_domain       ! Height of top of domain

      ! Grid dimensions
      Integer, Intent(In) :: row_length       ! No. of points on a row
      Integer, Intent(In) :: rows             ! No. of rows (theta)
      Integer, Intent(In) :: model_levels     ! No. of model levels
      Integer, Intent(In) :: wet_model_levels ! No. of wet model levels
      Integer, Intent(In) :: halo_i           ! Size of halo in x
      Integer, Intent(In) :: halo_j           ! Size of halo in y

      ! Multi-processor
      Integer, Intent(In) :: me               ! My processor number

      ! Idealised options
      Logical, Intent(In) :: L_code_test      ! User switch
      Logical, Intent(In) :: L_dry            ! Dry model switch

      Integer, Intent(In) :: max_num_profile_data ! max no. profile data
      Integer, Intent(In) :: num_profile_data ! actual no. data values
      Integer, Intent(In) :: qprofile_number  ! Humidity profile option
      Integer, Intent(In) :: idl_interp_option ! Profile interp option

      Real,    Intent(In) :: q1               ! For constant RH
      Real,    Intent(In) :: zprofile_orog    ! Orog hgt of initial prof
      Real,    Intent(In) :: hf               ! Hgt for constant interp

      ! Idealised namelist data profile
      Real, Intent(In) :: zprofile_data(max_num_profile_data) ! Heights
      Real, Intent(In) :: qprofile_data(max_num_profile_data) ! Humidity

      ! Vertical co-ordinate information
      Real, Intent(In) :: r_theta_levels(1-halo_i:row_length+halo_i,    &
     &                       1-halo_j:rows+halo_j,0:model_levels)
      Real, Intent(In) :: eta_theta_levels(0:model_levels)

! Arguments with Intent (InOut)

      ! 3D prognostic field arrays
      Real, Intent(InOut) ::                                            &
     &  theta(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,          &
     &                                   model_levels)                  &
     &, exner_theta_levels(1-halo_i:row_length+halo_i,                  &
     &                     1-halo_j:rows+halo_j, model_levels)          &
     &, q(1-halo_i:row_length+halo_i,                                   &
     &          1-halo_j:rows+halo_j, wet_model_levels)

      ! Vertical reference profile
      Real, Intent(InOut) :: theta_ref(model_levels)
      Real, Intent(InOut) :: exner_ref(model_levels+1)
      Real, Intent(InOut) :: q_ref(wet_model_levels)

! Local variables

      ! Work arrays
      Real :: p_row(1-halo_i:row_length+halo_i)  ! pressure on a row
      Real :: t_row(1-halo_i:row_length+halo_i)  ! temperature on a row
      Real :: qs_row(1-halo_i:row_length+halo_i) ! qsat on a row
      Real :: p_ref(model_levels)      ! reference pressure profile
      Real :: t_ref(model_levels+1)    ! reference temperature profile
      Real :: qs_ref(wet_model_levels) ! reference saturation profile
      Real :: weight                   ! interpolation weight
      Real :: z_at_theta               ! height asl of theta level
      Real :: z_at_orog                ! height of orog
      Real :: hs                       ! Temporary variable
      Real :: eta_model                ! eta coord for model levels
      Real :: eta_profile(num_profile_data) ! eta coord for input prof

      Integer :: i, j, k, k2           ! Loop counters
      Integer :: length                ! array length

      ! Error reporting
      Character (Len=*),  Parameter :: RoutineName=                     &
     &                                 'idl_set_initial_humidity'
      Character (Len=256)           :: Cmessage
      Integer                       :: ErrorStatus

#include "qprofile.h"
!-------------------------------------------------------------------

      If (me == 0) Then
        Write (6,*) ' '
        Write (6,*) ' HUMIDITY PROFILE  '
      End If

!-------------------------------------------------------------------
!
!                      Dry. Humidity = 0
!
!-------------------------------------------------------------------

      If (L_dry .or. (qprofile_number == qp_dry)) Then

        ! set moisture fields to zero.
        If (me == 0) Then
          Write (6,*) '   Dry model run, q = 0.0'
        End If   !(me == 0)

        q(:,:,:) = 0.0
        q_ref(:) = 0.0

!-------------------------------------------------------------------
!
!                  Set constant relative humidity
!
!-------------------------------------------------------------------
      Else If (qprofile_number == qp_qsat) Then

        If (me == 0) Then
          Write (6,*) '   Setting constant relative humidity',          &
     &                ' ',q1,'% wrt water T>0degC, wrt ice T<0degC'
        End If

        ! set q field to constant relative humidity
        length = row_length + halo_i + halo_i
        Do k = 1, wet_model_levels
          Do j = 1-halo_j, rows+halo_j

            ! Calculate temperature and pressure for QSAT
            Do i = 1-halo_i, row_length+halo_i
              t_row(i)=theta(i,j,k) * exner_theta_levels(i,j,k)
              p_row(i)=p_zero*exner_theta_levels(i,j,k)**recip_Kappa
            End Do

! DEPENDS ON: qsat
            Call qsat( qs_row(1-halo_i), t_row(1-halo_i),               &
     &                  p_row(1-halo_i), length)

            Do i = 1-halo_i, row_length+halo_i
              q(i,j,k) = q1 * qs_row(i)
            End Do

          End Do
        End Do

        ! Set reference profile
        Do k = 1, wet_model_levels
          t_ref(k) = theta_ref(k) * exner_ref(k)
          p_ref(k) = p_zero*exner_ref(k)**recip_Kappa
        End Do
        length = wet_model_levels
! DEPENDS ON: qsat
        Call qsat( qs_ref, t_ref, p_ref, length)
        Do k = 1, wet_model_levels
          q_ref(k) = q1 * qs_ref(k)
        End Do

!-------------------------------------------------------------------
!
!            Set field from profile in idealised namelist
!
!-------------------------------------------------------------------
      Else If (qprofile_number == qp_namelist .or.                      &
     &         qprofile_number == qp_namelist_rh) Then

        If (me == 0) Then
          Write (6,*) '   Setting humidity from namelist profile'
        End If

        ! Check to make sure the namelist profile data extends
        ! to the top of the model.
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            If (r_theta_levels(i,j,model_levels) - Earth_radius         &
     &           >   zprofile_data(num_profile_data)) Then
              Write(Cmessage,*)                                         &
     &          'Idealised namelist vertical profile data'              &
     &          //'does not extend to the top of the model.'            &
     &          //'Please modify the namelist data.'
              ErrorStatus = 1
! DEPENDS ON: ereport
              Call Ereport( RoutineName, ErrorStatus, Cmessage )
            End If
          End Do
        End Do

        ! Interpolate q from namelist profile to model levels
        Do k = 1, model_levels
          Do k2 = 1, num_profile_data-1
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i

                z_at_theta = r_theta_levels(i,j,k) - Earth_radius

                If (z_at_theta  >   zprofile_data(k2) .and.             &
     &              z_at_theta  <=  zprofile_data(k2+1)) Then

                  weight = (z_at_theta - zprofile_data(k2))             &
     &                   /(zprofile_data(k2+1) - zprofile_data(k2))

                  q(i,j,k) = qprofile_data(k2) + weight*                &
     &                   (qprofile_data(k2+1) - qprofile_data(k2))
                End If

              End Do
            End Do
          End Do
        End Do


        !---------------------------------------------------------------
        ! Alternative interpolation of initial profile over orography
        !---------------------------------------------------------------
        !
        !  idl_interp_option = 1: constant on height levels
        !                         (default above, no need to modify)
        !  idl_interp_option = 2: hybrid height everywhere up to a
        !                         specified height hf (the same as modl)
        !                         levels are defined). hs=height_domain
        !  idl_interp_option = 3: as option 2 but only when orography is
        !                         less than input profile orography.
        !                         hs = zprofile_orog
        !
        ! Sets up eta coord for each level of the initial profile data
        ! and an eta coordinate for each model column, and interpolates
        ! in eta space if the model level height is less than "hf" and
        ! the model orography height is less than "hs"
        !---------------------------------------------------------------

        If (idl_interp_option == 2 .or. idl_interp_option == 3) Then

          If (idl_interp_option == 2) hs = height_domain
          If (idl_interp_option == 3) hs = zprofile_orog

          ! Set up eta coord for each level of the initial profile data
          eta_profile(1) = 0.0
          Do k2 = 2, num_profile_data
            eta_profile(k2) =  (zprofile_data(k2) - zprofile_orog)      &
     &                        /(hf - zprofile_orog)
          End Do

          ! Interpolate in eta space and overwrite q where appropriate
          Do k = 1, wet_model_levels
            Do k2 = 1, num_profile_data - 1
              Do j = 1-halo_j, rows+halo_j
                Do i = 1-halo_i, row_length+halo_i
                  z_at_theta = r_theta_levels(i,j,k) - Earth_radius
                  z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
                  eta_model = (z_at_theta - z_at_orog)/(hf - z_at_orog)

                  If ( (z_at_orog <= hs ) .and. (z_at_theta < hf) .and. &
     &                 (eta_model > eta_profile(k2))  .and.             &
     &                 (eta_model <= eta_profile(k2+1)) ) Then

                    weight = (eta_model - eta_profile(k2))              &
     &                       /(eta_profile(k2+1) - eta_profile(k2))
                    q(i,j,k) = qprofile_data(k2) + weight*              &
     &                        (qprofile_data(k2+1) - qprofile_data(k2))
                  End If
                End Do
              End Do
            End Do
          End Do

        End If ! on idl_interp_option = 2 or 3


        !---------------------------------------------
        ! Set up reference q profile (assumes no orography)
        !---------------------------------------------
        Do k = 1, model_levels
          Do k2 = 1, num_profile_data-1

                z_at_theta = eta_theta_levels(k)*height_domain

                If (z_at_theta  >   zprofile_data(k2) .and.             &
     &              z_at_theta  <=  zprofile_data(k2+1)) Then

                  weight = (z_at_theta - zprofile_data(k2))             &
     &                   /(zprofile_data(k2+1) - zprofile_data(k2))

                  q_ref(k) = qprofile_data(k2) + weight*                &
     &                      (qprofile_data(k2+1) - qprofile_data(k2))
                End If

          End Do
        End Do


        !---------------------------------------------
        ! If input data is relative humidity (%/100)
        ! then convert to q (kg/kg)
        !---------------------------------------------
        If (qprofile_number == qp_namelist_rh) Then

          ! set q field to constant relative humidity wrt water
          If (me == 0) Then
            Write (6,*) '   Relative humidity (wrt water) data '
          End If

          ! Calculate humidity in kg/kg
          ! q(kg/kg) = RH (held in variable q) * qsat

          length = row_length + halo_i + halo_i
          Do k = 1, wet_model_levels
            Do j = 1-halo_j, rows+halo_j

              ! Calculate temperature and pressure for QSAT
              Do i = 1-halo_i, row_length+halo_i
                t_row(i)=theta(i,j,k) * exner_theta_levels(i,j,k)
                p_row(i)=p_zero*exner_theta_levels(i,j,k)**recip_Kappa
              End Do

! DEPENDS ON: qsat_wat
              Call qsat_wat( qs_row(1-halo_i), t_row(1-halo_i),         &
     &                        p_row(1-halo_i), length )

              Do i = 1-halo_i, row_length+halo_i
                q(i,j,k) = q(i,j,k) * qs_row(i)
              End Do

            End Do
          End Do

          ! Set reference profile
          Do k = 1, wet_model_levels
            t_ref(k) = theta_ref(k) * exner_ref(k)
            p_ref(k) = p_zero*exner_ref(k)**recip_Kappa
          End Do
          length = wet_model_levels
! DEPENDS ON: qsat
          Call qsat( qs_ref, t_ref, p_ref, length)
          Do k = 1, wet_model_levels
            q_ref(k) = q_ref(k) * qs_ref(k)
          End Do

        End If

!-------------------------------------------------------------------
!
!            Leave q exactly as it is (e.g. from initial dump)
!
!-------------------------------------------------------------------
      Else If (qprofile_number == qp_dump) Then

        Write (6,*) '   Start dump 3D humidity field will be used'

!-------------------------------------------------------------------
!
!              qprofile_number not recognised !
!
!-------------------------------------------------------------------
      Else

        Write(Cmessage,*)                                               &
     &          'qprofile_number not recognised.'                       &
     &        //' Please modify the idealised namelist data.'
        ErrorStatus = 1
! DEPENDS ON: ereport
        Call Ereport( RoutineName, ErrorStatus, Cmessage )

      End If

      Return
      END SUBROUTINE idl_set_initial_humidity

      !  End subroutine IDL_set_initial_humidity
#endif
