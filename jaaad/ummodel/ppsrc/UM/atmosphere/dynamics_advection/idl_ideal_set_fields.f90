
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_ideal_set_fields

      Subroutine IDL_ideal_set_fields(                                  &
     &                      R, g, kappa, epsilon, Cp                    &
     &,                     p_zero, Earth_radius, Pi                    &
     &,                     model_domain, row_length, rows, n_rows      &
     &,                     model_levels, wet_model_levels              &
     &,                     off_x, off_y, halo_i, halo_j                &
     &,                     me, n_proc, at_extremity                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     theta_eh, rho_eh, exner_rho_levels_eh       &
     &,                     u, v, w, u_adv, v_adv, w_adv                &
     &,                     theta, q, qcl, qcf, rho                     &
     &,                     exner_rho_levels, exner_theta_levels        &
     &,                     p, p_theta_levels, p_star                   &
     &,                     L_include_halos)

! Purpose:
!      1. Copy large halo initial data from work arrays into data arrays
!      2. Set work arrays (qcl, qcf, w_adv) and w to zero
!      3. Initialise pressure arrays!
!
! Original Progammer: T. Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version    Date     Comment
! ----     -------     -------
! 5.3       10/10/01    This deck created.      Terry Davies
! 6.2       31/01/06    Pass in extended halo theta,rho,exner. R. Forbes
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      Implicit None

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
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                             ! Local number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, off_x                                                           &
                             ! Size of small halo in i
     &, off_y                ! Size of small halo in j.

      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, n_proc     ! Total number of processors

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_include_halos

! Include physical constants
      Real                                                              &
           ! physical constants
     &  R                                                               &
     &, g                                                               &
     &, Cp                                                              &
     &, epsilon                                                         &
     &, kappa                                                           &
     &, p_zero                                                          &
               ! reference pressure
     &, Earth_radius                                                    &
     &, Pi

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)

      ! Work arrays with extended halos (_eh)
      ! Needed so that external halo values in LAMS can be set correctly
      ! for the lateral boundary arrays.
      REAL, Intent (InOut) ::                                           &
     &  theta_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &           model_levels)                                          &
     &, rho_eh(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &           model_levels)                                          &
     &, exner_rho_levels_eh(1-halo_i:row_length+halo_i,                 &
     &           1-halo_j:rows+halo_j, model_levels+1)

! Primary Arrays
      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    model_levels)                                                 &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &    model_levels)                                                 &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    0:model_levels)                                               &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                             &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &        model_levels)                                             &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        0:model_levels)                                           &
     &, rho(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &      model_levels)                                               &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &    model_levels+1)                                               &
     &, p_star(row_length, rows)                                        &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        model_levels)                                             &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &    wet_model_levels)                                             &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      wet_model_levels)                                           &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &      wet_model_levels)                                           &
     &, exner_rho_levels(1-off_x:row_length+off_x,                      &
     &                   1-off_y:rows+off_y, model_levels+1)            &
     &, exner_theta_levels(1-off_x:row_length+off_x,1-off_y:rows+off_y, &
     &        model_levels)                                             &
     &, p_theta_levels(1-off_x:row_length+off_x,  1-off_y:rows+off_y,   &
     &        model_levels)

! local variables
      Integer                                                           &
     & i, j, k   ! loop counters

      Real                                                              &
     & recip_kappa

      recip_kappa = 1. / kappa

! First copy from big haloed work arrays into normal arrays
! u,v fields are in u_adv, v_adv
! DEPENDS ON: copy_field
      CALL COPY_FIELD(THETA_EH, THETA                                   &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(EXNER_RHO_LEVELS_EH, EXNER_RHO_LEVELS             &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels+1, model_levels+1, 1, model_levels+1 &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(RHO_EH, RHO                                       &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(U_ADV, U                                          &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_u, .true., .false., .true.)
! DEPENDS ON: copy_field
      CALL COPY_FIELD(V_ADV, V                                          &
     &,               row_length, row_length, n_rows, n_rows            &
     &,               model_levels, model_levels, 1, model_levels       &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_v, .true., .false., .true.)

! Set required arrays to zero, ( w_adv, qcl, qcf)
      Do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          w_adv(i,j,0) = 0.0
        End Do
      End Do
      Do k = 1, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            w_adv(i,j,k) = 0.0
            qcl(i,j,k) = 0.0
            qcf(i,j,k) = 0.0
          End Do
        End Do
      End Do   !  k = 1, model_levels

! Initialise remaining field w=0
! DEPENDS ON: copy_field
      CALL COPY_FIELD(W_ADV, W                                          &
     &,               row_length, row_length, rows, rows                &
     &,               model_levels+1, model_levels+1, 1, model_levels+1 &
     &,               halo_i, halo_j, off_x, off_y                      &
     &,               fld_type_p, .true., .false., .false.)
!
! Must now re-calculate the pressure-based variables, namely pressure
! on both rho and theta levels, exner on theta levels and p_star so that
! they are all consistent with the new LBC-updated values of exner on
! rho levels.
!
! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta(                                         &
     &                      r_theta_levels, r_rho_levels,               &
     &                      Exner_rho_levels,                           &
     &                      row_length, rows, model_levels,             &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      Exner_theta_levels, L_include_halos)

! Calculate pressure from Exner at rho levels.
      Do k = 1, model_levels+1
        Do j = 1, rows
          Do i = 1, row_length
            p(i,j,k) = p_zero * exner_rho_levels(i,j,k)**recip_kappa
          End Do
        End Do
      End Do
! Calculate p_theta_levels from Exner at theta levels.
      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            p_theta_levels(i,j,k) = p_zero *                            &
     &                      exner_theta_levels(i,j,k)**recip_kappa
          End Do
        End Do
      End Do
! Halos updated
! DEPENDS ON: swap_bounds
      call Swap_Bounds(p,                                               &
     &                   row_length, rows, model_levels+1,              &
     &                   off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
      call Swap_Bounds(p_theta_levels,                                  &
     &                   row_length, rows, model_levels,                &
     &                   off_x, off_y, fld_type_p, .false.)

! DEPENDS ON: calc_p_star
      Call Calc_P_star(                                                 &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   p, rho, g, row_length, rows, model_levels,     &
     &                   off_x, off_y, halo_i, halo_j,                  &
     &                   p_star)

      return
      END SUBROUTINE IDL_ideal_set_fields
