
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Force_Lbc
!

      Subroutine IDL_Force_Lbc(                                         &
     &             R, g, Cp, kappa, epsilon, p_zero                     &
     &,            row_length, rows, off_x, off_y                       &
     &,            halo_i, halo_j, Earth_radius,  LENRIMA               &
     &,            timestep, timestep_number                            &
     &,            model_levels, wet_levels                             &
     &,            max_model_levels, max_num_force_times                &
     &,            u_lbc, v_lbc, theta_lbc                              &
     &,            q_lbc, u_adv_lbc, v_adv_lbc                          &
     &,            exner_rho_levels_LBC                                 &
     &,            r_theta_levels, r_rho_levels                         &
     &,            eta_theta_levels, eta_rho_levels                     &
     &,            height_domain, theta_surface                         &
     &,            pforce_option                                        &
     &,            tforce_option, qforce_option, uvforce_option         &
     &,            num_pforce_times                                     &
     &,            num_tforce_times, num_qforce_times                   &
     &,            num_uvforce_times                                    &
     &,            pforce_time_interval                                 &
     &,            tforce_time_interval, qforce_time_interval           &
     &,            uvforce_time_interval                                &
     &,            p_surface_data                                       &
     &,            tforce_data_modlev, qforce_data_modlev               &
     &,            uforce_data_modlev, vforce_data_modlev               &
     &,            newtonian_timescale )


      Implicit None

! Purpose: To apply idealised LBC forcing.
!
! Method: Resets LBC data to profile data timeseries specified in
!         idealised namelist.
!
! Original Programmer:
! Current code owner:
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 6.2    01/03/06   Original code.  Yongming Tang
!
! Code Description:
!   Language: FORTRAN 77 + common extensions
!   This code is written to UMDP3 programming standards.

! Parameters required for dimensioning some of the arguments
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end

! Variables with Intent (In)

      ! Physical constants
      Real, Intent(In) ::                                               &
     &  R                                                               &
     &, g                                                               &
     &, kappa                                                           &
     &, epsilon                                                         &
     &, Cp                                                              &
     &, p_zero                                                          &
     &, Earth_radius

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, off_x                                                           &
                         ! Size of small halo in i
     &, off_y                                                           &
                         ! Size of small halo in j.
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_levels                                                      &
     &, num_pforce_times                                                &
                            ! No. of times in P forcing data
     &, num_tforce_times                                                &
                            ! No. of times in T forcing data
     &, num_qforce_times                                                &
                            ! No. of times in Q forcing data
     &, num_uvforce_times                                               &
                            ! No. of times in UV forcing data
     &, pforce_option                                                   &
     &, tforce_option                                                   &
     &, qforce_option                                                   &
     &, uvforce_option                                                  &
     &, timestep_number                                                 &
                             ! Model timestep number in run&,
     &, max_model_levels                                                &
                             ! Max number of model levels
     &, max_num_force_times                                             &
                             ! Max number of times in forcing data
     &, lenrima(Nfld_max,NHalo_max)
                             ! IN : Size of single level of LBC

      Real, Intent(In) ::                                               &
     &  timestep                                                        &
                             ! Length of timestep in seconds
     &, pforce_time_interval                                            &
                             ! Forcing data time interval
     &, tforce_time_interval                                            &
                             ! Forcing data time interval
     &, qforce_time_interval                                            &
     &, uvforce_time_interval                                           &
     &, height_domain                                                   &
     &, theta_surface                                                   &
     &, newtonian_timescale  ! Relaxation timescale

      ! Forcing data interpolated to model levels
      Real, Intent(In) ::                                               &
     &  p_surface_data(max_num_force_times)                             &
     &, tforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, qforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, uforce_data_modlev(max_model_levels, max_num_force_times)       &
     &, vforce_data_modlev(max_model_levels, max_num_force_times)

! Variables with Intent (InOut)

      Real, Intent(InOut) ::                                            &
     &  u_lbc(lenrima(fld_type_u,halo_type_extended),                   &
     &            model_levels)                                         &
     &, v_lbc(lenrima(fld_type_v,halo_type_extended),                   &
     &            model_levels)                                         &
     &, theta_lbc(lenrima(fld_type_p,halo_type_extended),               &
     &            model_levels)                                         &
     &, q_lbc(lenrima(fld_type_p,halo_type_extended),                   &
     &        wet_levels)                                               &
     &, u_adv_lbc(lenrima(fld_type_u,halo_type_extended),               &
     &            model_levels)                                         &
     &, v_adv_lbc(lenrima(fld_type_v,halo_type_extended),               &
     &            model_levels)                                         &
     &, exner_rho_levels_lbc(lenrima(fld_type_p,halo_type_extended),    &
     &            model_levels+1)                                       &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, r_theta_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j, &
     &         0: model_levels)                                         &
     &, r_rho_levels(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,   &
     &         model_levels)

! Start C_IDL_FORCE_OPTIONS

! Description: Include file containing idealised forcing options
!
! Author : R. Forbes
!
! History:
! Version  Date      Comment
! -------  ----      -------
! 6.1      01/08/04  Original code     R. Forbes

      INTEGER, PARAMETER :: no_forcing      = 0
      INTEGER, PARAMETER :: force_increment = 1
      INTEGER, PARAMETER :: force_relax     = 2
      INTEGER, PARAMETER :: force_reset     = 3

! End C_IDL_FORCE_OPTIONS

! Local variables
      Integer i, k
      Integer Lbc_sizes_u, Lbc_sizes_v, Lbc_sizes_p

!-----------------------------------------------------------------------

      Lbc_sizes_p = lenrima(fld_type_p,halo_type_extended)
      Lbc_sizes_u = lenrima(fld_type_u,halo_type_extended)
      Lbc_sizes_v = lenrima(fld_type_v,halo_type_extended)

      !----------------------------------------------
      !              Temperature forcing
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_p, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_tforce_times, tforce_time_interval               &
     &,            tforce_data_modlev, theta_lbc )

      !----------------------------------------------
      !              Humidity forcing
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_p, wet_levels                              &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_qforce_times, qforce_time_interval               &
     &,            qforce_data_modlev, q_lbc )

!      !----------------------------------------------
!      !              Exner Pressure
!      !   (Reset data to specified forcing data)
!      !----------------------------------------------
!      ! Not yet included
!
!        Call IDL_LBC_Pressure(
!     &             g, Cp, kappa, epsilon,  p_zero
!     &,            row_length, rows, halo_i, halo_j
!     &,            Earth_radius, Lbc_sizes_p
!     &,            model_levels, wet_levels
!     &,            theta_surface
!     &,            r_theta_levels, r_rho_levels
!     &,            eta_theta_levels, eta_rho_levels
!     &,            theta_lbc, q_lbc
!     &,            timestep, timestep_number
!     &,            max_num_force_times
!     &,            num_pforce_times, pforce_time_interval
!     &,            p_surface_data
!     &,            exner_rho_levels_LBC )


      !----------------------------------------------
      !           Horizontal Wind Forcing
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_u, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            uforce_data_modlev, u_lbc )

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_v, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            vforce_data_modlev, v_lbc )

      !----------------------------------------------
      !             Horizontal Wind (_adv)
      !   (Reset data to specified forcing data)
      !----------------------------------------------

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_u, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            uforce_data_modlev, u_adv_lbc )

! DEPENDS ON: idl_lbc_reset
        Call IDL_LBC_Reset(                                             &
     &             Lbc_sizes_v, model_levels                            &
     &,            timestep, timestep_number                            &
     &,            max_model_levels, max_num_force_times                &
     &,            num_uvforce_times, uvforce_time_interval             &
     &,            vforce_data_modlev, v_adv_lbc )


! End of routine.

      Return
      END SUBROUTINE IDL_Force_Lbc

