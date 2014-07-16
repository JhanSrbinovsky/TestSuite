
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SL_tracer2
!
      Subroutine SL_tracer2(                                            &
     &                           super_array_size,                      &
     &                           super_array, super_tracer_phys2,       &
     &                           eta_theta_levels,                      &
     &                           r_rho_levels, r_theta_levels,          &
     &                           rho_n, rho_np1,                        &
     &                           row_length, rows, model_levels,        &
     &                           delta_lambda, delta_phi,               &
     &                           glambda_p, phi_p,                      &
     &                           grecip_dlamp, recip_dphip,             &
     &                           lambda_p_rm, lambda_p_rp,              &
     &                           phi_p_rm, phi_p_rp,                    &
     &                           recip_lambda_p_m, recip_lambda_p_0,    &
     &                           recip_lambda_p_p, recip_lambda_p_p2,   &
     &                           recip_phi_p_m, recip_phi_p_0,          &
     &                           recip_phi_p_p, recip_phi_p_p2,         &
     &                           base_lambda, base_phi,                 &
     &                           recip_dlam, recip_dphi, max_look,      &
     &                           look_lam, look_phi,                    &
     &                           halo_lam, halo_phi,                    &
     &                           FV_cos_theta_latitude,                 &
     &                           wet_to_dry_n, wet_to_dry_np1,          &
! LAM bit
     &                           me, n_proc, n_procx, n_procy,          &
     &                           halo_i, halo_j, l_datastart,           &
     &                           g_i_pe, at_extremity,                  &
     &                           global_row_length,                     &
     &                           global_rows,                           &
     &                           gc_all_proc_group,                     &
     &                           proc_row_group,                        &
     &                           proc_col_group, off_x, off_y,          &
     &                           L_regular, L_sl_halo_reprod,           &
     &                           high_order_scheme_moist,               &
     &                           monotone_scheme_moist,                 &
     &                           model_domain, L_high_moist,            &
     &                           L_mono_moist, L_conserv_moist,         &
     &                           check_bottom_levels,                   &
     &                           interp_vertical_search_tol,            &
     &                           first_constant_r_rho_level,            &
     &                           depart_lambda, depart_phi,             &
     &                           depart_r,                              &
     &                           CO2, L_CO2_interactive,                &
     &                           murk, L_murk_advect,                   &
     &                           soot_new, soot_agd, soot_cld, L_soot,  &
     &                           bmass_new, bmass_agd, bmass_cld,       &
     &                           L_biomass,                             &
     &                           ocff_new, ocff_agd, ocff_cld, l_ocff,  &
     &                           DUST_DIV1,DUST_DIV2,DUST_DIV3,         &
     &                           DUST_DIV4,DUST_DIV5,DUST_DIV6,         &
     &                           L_DUST,                                &
     &                           so2, so4_aitken, so4_accu,             &
     &                           so4_diss, nh3, dms,                    &
     &                           L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms, &
     &                           tracers, tr_levels, tr_vars,           &
     &                           tracers_ukca, tr_ukca,                 &
     &                           tracer_phys1, tracer_phys2,            &
     &                           i_start, i_end, j_start, j_end,        &
     &                           L_USE_CARIOLLE, OZONE_TRACER,          &
     &                           Error_Code)

! Purpose:
!          Performs semi-Lagrangian advection of tracers
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
! Version   Date       Comment
! ----     -------     -------
! 5.2      15/11/00     original deck based on sl_thermo  Andy Malcolm
! 5.3      25/04/01     Call added to fill external halos for murk
!                       variable in the mes model.        S. Cusack
! 5.3      11/07/01     bug fix                          Andy Malcolm
! 5.3      24/10/01  don't recalculate departure point       A.Malcolm
! 5.4      03/04/02  correct conservation constraint by using
!                    rho_dry rather than rho                 A.Malcolm
! 5.5      05/02/03  Add code to advect 3 modes of biomass aerosol.
!                                                     P Davison
! 5.5      14/02/03  improve tracer conservation properties A.Malcolm
! 5.5      25/02/03  tidy up code, remove moisture in       A. Malcolm
! 5.5  12/02/03  Include code for mineral dust scheme.      S Woodward
! 6.1      13/05/04   changes to use super_array_size      Andy Malcolm
! 6.2      16/12/05   Set external halo values for LAM's   Andy Malcolm
!
!  6.2    25/12/05  Variable resolution changes         Yongming Tang
!  6.4    17/12/06  Major speedup changes               Andy Malcolm
!  6.4    16/02/07  Fix Tracers in LAM's                Andy Malcolm
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
     &, rows                                                            &
                     ! number of rows.
     &, model_levels                                                    &
                     ! Number of model levels.
     &, me                                                              &
                     ! My processor number
     &, n_proc                                                          &
                     ! Total number of processors
     &, n_procx                                                         &
                     ! Number of processors in longitude
     &, n_procy                                                         &
                     ! Number of processors in latitude
     &, halo_i                                                          &
                     ! Size of large halo in i.
     &, halo_j                                                          &
                     ! Size of large halo in j.
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, l_datastart(3)                                                  &
                       ! First gridpoints held by this processor.
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, proc_col_group                                                  &
                       ! Group id for processors on the same column
     &, max_look                                                        &
                             ! max size of look-up arrays for searches
     &, gc_all_proc_group                                               &
     &, global_row_length                                               &
                            ! global number of points on a row
     &, global_rows                                                     &
                            ! global number of rows
     &, g_i_pe(1-halo_i:global_row_length+halo_i)
                             ! processor on my processor-row
                             ! holding a given value in i direction

      Integer                                                           &
     &  tr_levels                                                       &
     &, tr_vars                                                         &
     &, tr_ukca

      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Logical                                                           &
     &  L_sl_halo_reprod ! if true then sl code bit repoducible with
                         ! any sensible halo size

      Integer                                                           &
     &  high_order_scheme_moist                                         &
                                 ! a code saying which high order
                           ! scheme to use for moist variables.
     &, monotone_scheme_moist                                           &
                              ! a code saying which monotone
                           ! scheme to use for moist variables.
     &, interp_vertical_search_tol                                      &
                                   ! used in interpolation code.
     &, check_bottom_levels ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Integer   i_start, i_end, j_start, j_end


      Logical                                                           &
     &  L_high_moist                                                    &
                       ! True, if high order interpolation required
                       !       for moist variables.
     &, L_mono_moist                                                    &
                       ! True, if interpolation required to be monotone
                       !       for moist variables.
     &, L_conserv_moist                                                 &
                        ! True, if interpolation to be monotone and
                       !       conservative for moist variables.
     &, L_regular
 
      Real                                                              &
     &  delta_lambda                                                    &
                      ! holds spacing between points in the i
                      ! direction for the input data field.
     &, delta_phi                                                       &
                      ! holds spacing between points in the j
                      ! direction for the input data field.
     &, recip_dlam                                                      &
     &, recip_dphi                                                      &
     &, base_lambda                                                     &
     &, base_phi

! look-up table halos
      Integer                                                           &
     &  halo_lam                                                        &
     &, halo_phi

!VarRes horizontal co-ordinate look-up table
      Integer                                                           &
     &  look_lam(1-halo_lam : max_look-halo_lam)                        &
     &, look_phi(1-halo_phi : max_look-halo_phi)

!VarRes horizontal co-ordinate information
      Real                                                              &
     &  glambda_p(1-halo_i : global_row_length+halo_i)                  &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, grecip_dlamp(1-halo_i : global_row_length+halo_i)               &
     &, recip_dphip( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, lambda_p_rm(1-halo_i : row_length+halo_i)                       &
     &, lambda_p_rp(1-halo_i : row_length+halo_i)                       &
     &, phi_p_rm   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, phi_p_rp   ( 1-halo_i : row_length + halo_i                     &
     &,              1-halo_j : rows+halo_j )                           &
     &, recip_lambda_p_m(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_0(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p(1-halo_i : row_length+halo_i)                  &
     &, recip_lambda_p_p2(1-halo_i : row_length+halo_i)                 &
     &, recip_phi_p_m ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_0 ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p ( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )                        &
     &, recip_phi_p_p2( 1-halo_i : row_length + halo_i                  &
     &,                 1-halo_j : rows+halo_j )


      Logical, Intent(In) ::                                            &
     &  L_CO2_interactive                                               &
     &, L_murk_advect                                                   &
     &, L_Soot                                                          &
     &, L_biomass                                                       &
     &, L_ocff                                                          &
     &, L_DUST                                                          &
     &, L_sulpc_so2, L_sulpc_nh3, l_sulpc_dms                           &
     &, L_USE_CARIOLLE
     

      Real, Intent(InOut) ::                                            &
     & CO2 (1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &,murk(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &,soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,soot_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         model_levels)                                            &
     &,bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &,bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &,bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &,ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          model_levels)                                           &
     &,ocff_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         model_levels)                                            &
     &, DUST_DIV1(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV2(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV3(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV4(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV5(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, DUST_DIV6(1-OFF_X:ROW_LENGTH+OFF_X, 1-OFF_Y:ROWS+OFF_Y,         &
     &           MODEL_LEVELS)                                          &
     &, so2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &           model_levels)                                          &
     &, so4_aitken(1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &           model_levels)                                          &
     &, so4_accu(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &, so4_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &           model_levels)                                          &
     &, nh3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &           model_levels)                                          &
     &, dms(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
     &           model_levels)                                          &
     &, tracers(1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          tr_levels,tr_vars)                                      &
     &, tracers_ukca(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &          tr_levels,tr_ukca)                                      &
! Add cariolle specific parameters for ozone tracer     
     &, OZONE_TRACER(1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &           model_levels)                                          

      Real                                                              &
     &  r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, eta_theta_levels(0:model_levels)

      Real                                                              &
     &  rho_n   (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, model_levels)                      &
     &, rho_np1 (1-off_x:row_length+off_x,                              &
     &           1-off_y:rows+off_y, model_levels)                      &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)

      Real                                                              &
                      ! Trig functions.
     &  FV_cos_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)

      Real                                                              &
     &  depart_lambda (row_length, rows, model_levels)                  &
                                                      ! Lambda
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_phi (row_length, rows, model_levels)                     &
                                                    ! Phi Co-ordinate of
                                                      ! co-ordinate of
                                                      ! departure point.
     &, depart_r (row_length, rows, model_levels)     ! Vertical
                                                      ! co-ordinate of
                                                      ! departure point.

! Arguments with Intent OUT. ie: Output variables.

      Integer                                                           &
     &  Error_Code     ! Non-zero on exit if error detected.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

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

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k, l                                                      &
                       ! Loop indices
     &, temp                                                            &
     &, count                                                           &
     &, tr_start

      Logical                                                           &
     &  L_vector       ! True, if data is a horizontal vector component,
                       ! False, then data is a scalar.

! arrays

      Real                                                              &
     &  work(row_length, rows, model_levels)                            &
     &, work_np1(row_length, rows, model_levels)                        &
     &, drk, drkp1                                                      &
     &, work4(row_length, rows)

      integer super_array_size, array_size_count, super_array_size_qpos

      Real :: super_array(1-halo_i:row_length+halo_i,                   &
                          1-halo_j:rows+halo_j,                         &
                          model_levels, super_array_size)     
      Real :: super_tracer_phys2(row_length,                            &
                                 rows,                                  &
                                 model_levels, super_array_size)         
      Real :: tracer_phys1(1-halo_i:row_length+halo_i,                  &
                           1-halo_j:rows+halo_j,                        &
                           tr_levels, tr_vars)   
      Real :: tracer_phys2(row_length, rows, tr_levels, tr_vars)   
      Real :: data_out_super(row_length, rows, model_levels, super_array_size)

      Integer                                                           &
     &  i_out (row_length, rows, model_levels)                          &
     &, j_out (row_length, rows, model_levels)

      Real                                                              &
     &  weight_lambda (row_length, rows, model_levels)                  &
     &, weight_phi    (row_length, rows, model_levels)

! ----------------------------------------------------------------------
!  Section 0.    Initialise array_size_count
! ----------------------------------------------------------------------

      array_size_count=0
!    qpos is not called for murk 
      super_array_size_qpos=super_array_size
      if (L_murk_advect)super_array_size_qpos=super_array_size_qpos-1

!qcon block start
! ----------------------------------------------------------------------
! Section 3.1  Set appropriate weighted rho*r*r*delta_r at data points
! and at time t and t+deltat.
! ----------------------------------------------------------------------

! Weights come from rewriting: Sum_k(qbarr * rho*r*r*dr) =
!                              Sum_k(q * weighted average of rho*r*r*dr)

! Store in work (for current timestep n) and work_np1 (for next
! timestep n+1) -  'work' array therefore reused

! It is assumed that rho_n holds r-squared scaled current value of rho
! and that rho_np1 holds r-squared scaled value of rho at next timestep

        k = 1
! Note it is assumed that q(0) = q(1) and so q(0) contribution has
! been absorbed into that of q(1), hence different form of drk.
! This is not essential part of alogorithm and could be changed.
        Do j = 1, rows
          Do i = 1, row_length
            drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
            drk   = r_theta_levels(i,j,k)   - r_theta_levels(i,j,k-1)
          work(i,j,k)     = wet_to_dry_n(i,j,k+1)*rho_n(i,j,k+1)*drkp1  &
     &                   +  wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
        work_np1(i,j,k) = wet_to_dry_np1(i,j,k+1)*rho_np1(i,j,k+1)*drkp1&
     &                  + wet_to_dry_np1(i,j,k)  *rho_np1(i,j,k)  *drk
          End Do
        End Do

        Do k = 2, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              drkp1 = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
              drk   = r_rho_levels(i,j,k)     - r_theta_levels(i,j,k-1)
          work(i,j,k)     = wet_to_dry_n(i,j,k+1)*rho_n(i,j,k+1)*drkp1  &
     &                    + wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
       work_np1(i,j,k) = wet_to_dry_np1(i,j,k+1)*rho_np1(i,j,k+1)*drkp1 &
     &                 + wet_to_dry_np1(i,j,k)  *rho_np1(i,j,k)  *drk
            End Do
          End Do
        End Do

        k = model_levels
        Do j = 1, rows
          Do i = 1, row_length
              drk   = r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
              work(i,j,k)   = wet_to_dry_n(i,j,k)  *rho_n(i,j,k)  *drk
            work_np1(i,j,k) = wet_to_dry_np1(i,j,k)*rho_np1(i,j,k)*drk
          End Do
        End Do

!qcon block ends

! ----------------------------------------------------------------------
! Section 1.0  do interpolation of tracers
! ----------------------------------------------------------------------

        L_vector = .false.

! DEPENDS ON: calc_index
          Call Calc_Index(                                              &
     &                    row_length, rows, model_levels,               &
     &                    delta_lambda, delta_phi,                      &
     &                    base_lambda, base_phi,                        &
     &                    glambda_p, phi_p, grecip_dlamp, recip_dphip,  &
     &                    recip_dlam, recip_dphi, max_look,             &
     &                    look_lam, look_phi, halo_lam, halo_phi,       &
     &                    L_regular, depart_lambda, depart_phi,         &
     &                    halo_i, halo_j,                               &
     &                    global_row_length,                            &
     &                    row_length, rows, l_datastart,                &
     &                    i_out, j_out,                                 &
     &                    weight_lambda, weight_phi)

! super array set up in first tr_set_phys

      if(super_array_size >= 1)then

! DEPENDS ON: interpolation_qcon_multi
         Call Interpolation_qcon_multi(                                 &
     &                     super_array,                                 &
     &                     eta_theta_levels(1),                         &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     work, work_np1,                              &
     &                     r_theta_levels(1-halo_i,1-halo_j,1),         &
     &                     super_array_size, check_bottom_levels,       &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, model_levels,              &
     &                     rows,                                        &
     &                     row_length, rows, model_levels,              &
     &                     glambda_p, phi_p,                            &
     &                     lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,&
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out, j_out,                                &
     &                     weight_lambda, weight_phi,                   &
     &                     high_order_scheme_moist,                     &
     &                     monotone_scheme_moist,                       &
     &                     FV_cos_theta_latitude, L_regular,            &
     &                     L_vector, model_domain, L_high_moist,        &
     &                     L_mono_moist, L_conserv_moist,               &
     &                     depart_r, depart_lambda, depart_phi,         &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j, global_row_length,           &
     &                     global_rows, row_length, rows,               &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 1,           &
     &                     0, 0, L_sl_halo_reprod,                      &
     &                     off_x, off_y,                                &
     &                     data_out_super, Error_Code )

      endif

!  add on atmos_physics2 increment to readvected tracer
      Do l=1,super_array_size
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              super_tracer_phys2(i,j,k,l) = data_out_super(i,j,k,l) +   &
     &                  super_tracer_phys2(i,j,k,l)
            End Do
          End Do
        End Do
      End Do

      if (super_array_size_qpos >0)then
! DEPENDS ON: q_pos_ctl
        call Q_Pos_Ctl(                                                 &
     &                   super_tracer_phys2, row_length, rows,          &
     &                   model_levels*super_array_size_qpos,            &
     &                   global_row_length, global_rows,                &
     &                   me, n_proc, 0,0,                               &
     &                   gc_all_proc_group,                             & 
     &                   model_domain,                                  &
     &                   halo_type_no_halo, .true. ,  0.0               &
     &                   )
      endif      ! super_array_size_qpos >0

! ----------------------------------------------------------------------
! section 2.0  set end of timestep tracers
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.1  carbon cycle.
! ----------------------------------------------------------------------
      array_size_count=0
      if(l_CO2_interactive)then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              co2(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      endif  ! l_CO2_interactive

! ----------------------------------------------------------------------
! Section 2.2  Soot cycle.
! ----------------------------------------------------------------------
      If (l_Soot) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              soot_new(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              soot_agd(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              soot_cld(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      endif  ! l_soot

! ----------------------------------------------------------------------
! Section 2.3  Biomass aerosol.
! ----------------------------------------------------------------------
      If (l_Biomass) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              bmass_new(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              bmass_agd(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              bmass_cld(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      endif  ! l_biomass

! ----------------------------------------------------------------------
! Section 2.4  sulphur cycle.
! ----------------------------------------------------------------------
      If (l_Sulpc_so2) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              so4_aitken(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              so4_accu(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              so4_diss(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
              so2(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+3)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +3

        if(L_sulpc_nh3)then
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = j_start,j_end
              Do i = i_start,i_end
                nh3(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        endif

        if(L_sulpc_dms)then
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = j_start,j_end
              Do i = i_start,i_end
                dms(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        endif

      endif  ! l_sulpc_so2

! ----------------------------------------------------------------------
! Section 2.5  mineral dust.
! ----------------------------------------------------------------------
      IF (L_DUST) THEN
        array_size_count=array_size_count +1
        DO K = 1, MODEL_LEVELS
          Do j = j_start,j_end
            Do i = i_start,i_end
              DUST_DIV1(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              DUST_DIV2(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              DUST_DIV3(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
              DUST_DIV4(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+3)
              DUST_DIV5(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+4)
              DUST_DIV6(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+5)
            END DO
          END DO
        END DO
        array_size_count=array_size_count +5
      ENDIF  ! L_DUST

! ----------------------------------------------------------------------
! New addition  Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
      If (L_OCFF) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              ocff_new(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
              ocff_agd(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+1)
              ocff_cld(i,j,k) = super_tracer_phys2(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      endif  ! l_ocff

! ----------------------------------------------------------------------
! Section 2.5.1  Cariolle ozone tracer.
! ----------------------------------------------------------------------
      IF (L_USE_CARIOLLE) THEN
        array_size_count=array_size_count +1
        DO K = 1, MODEL_LEVELS
          Do j = j_start,j_end
            Do i = i_start,i_end
       OZONE_TRACER(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
            END DO
          END DO
        END DO
      ENDIF  ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 2.6.a  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF (model_levels==tr_levels.and. tr_vars>0) THEN
        do count=1,tr_vars
          array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS
            Do j = j_start,j_end
              Do i = i_start,i_end
              tracers(i,j,k,count) = super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End Do
      End IF  ! tr_vars>0

! ----------------------------------------------------------------------
! Section 2.7  UKCA tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF (model_levels==tr_levels.and. tr_ukca>0) THEN
        do count=1,tr_ukca
          array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS
            Do j = j_start,j_end
              Do i = i_start,i_end
              tracers_ukca(i,j,k,count) =                               &
     &           super_tracer_phys2(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End Do
      End IF  ! tr_ukca>0

! ----------------------------------------------------------------------
! Section 2.99  Murk cycle.  This must be the last Full level field in
!                           the super_array
! ----------------------------------------------------------------------
      IF (L_Murk_advect) then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = j_start,j_end
            Do i = i_start,i_end
              murk(i,j,k) = super_tracer_phys2(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      End IF  ! L_Murk_advect

! ----------------------------------------------------------------------
! Section 3.0  free tracers.   (model_levels/=tr_levels)
! ---------------------------------------------------------------------
! ----------------------------------------------------------------------
! Section 3.1  Trajectory limit at bottom tracer_level required
! ----------------------------------------------------------------------

      If (tr_vars > 0 .and. tr_levels /= model_levels) Then

! calculate trajectory limit at departure point

! assume once one level has no data below bottom then no higher level
! has either.
         k = model_levels - tr_levels + 1
         temp=1
         Do while ( temp  >   0 )
           temp = 0

! DEPENDS ON: bi_linear_h
           Call bi_linear_h(                                            &
     &                      r_theta_levels(1-halo_i,1-halo_j,           &
     &                                     model_levels-tr_levels+1),   &
     &                      depart_lambda(1,1,k),                       &
     &                      depart_phi(1,1,k),                          &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows, 1,                        &
     &                      row_length, rows,                           &
     &                      i_out(1,1,k), j_out(1,1,k),                 &
     &                      weight_lambda(1,1,k), weight_phi(1,1,k),    &
     &                      model_domain, me, n_procx,                  &
     &                      halo_i, halo_j, l_datastart,                &
     &                      global_row_length, g_i_pe, at_extremity,    &
     &                      1, proc_row_group,                          &
     &                      L_sl_halo_reprod, L_regular,                &
     &                      work4 )

           Do j = 1, rows
             Do i = 1, row_length

               If (depart_r(i,j,k) <   work4(i,j) ) Then
! move trajectory up to lowest level of data.
                 depart_r(i,j,k) = work4(i,j)
                 temp = temp + 1
               End If

             End Do     ! i loop
           End Do       ! j loop
           k = k + 1
           call gc_isum (1,n_proc,error_code,temp)
         End Do         ! do while

! ----------------------------------------------------------------------
! Section 3.2  do advection of tracer fields
! ----------------------------------------------------------------------
         tr_start = model_levels - tr_levels + 1
! DEPENDS ON: interpolation_qcon_multi
         Call Interpolation_qcon_multi(                                 &
     &                     tracer_phys1,                                &
     &                     eta_theta_levels(tr_start),                  &
     &                     r_theta_levels(1-halo_i,1-halo_j,tr_start),  &
     &                     work(1,1,tr_start), work_np1(1,1,tr_start),  &
     &                     r_theta_levels(1-halo_i,1-halo_j,tr_start),  &
     &                     tr_vars, check_bottom_levels,                &
     &                     interp_vertical_search_tol,                  &
     &                     first_constant_r_rho_level,                  &
     &                     row_length, rows, tr_levels,                 &
     &                     rows,                                        &
     &                     row_length, rows, tr_levels,                 &
     &                     glambda_p, phi_p,                            &
     &                     lambda_p_rm, lambda_p_rp, phi_p_rm, phi_p_rp,&
     &                     recip_lambda_p_m, recip_lambda_p_0,          &
     &                     recip_lambda_p_p, recip_lambda_p_p2,         &
     &                     recip_phi_p_m, recip_phi_p_0,                &
     &                     recip_phi_p_p, recip_phi_p_p2,               &
     &                     i_out(1,1,tr_start), j_out(1,1,tr_start),    &
     &                     weight_lambda(1,1,tr_start),                 &
     &                     weight_phi(1,1,tr_start),                    &
     &                     high_order_scheme_moist,                     &
     &                     monotone_scheme_moist,                       &
     &                     FV_cos_theta_latitude, L_regular,            &
     &                     L_vector, model_domain, L_high_moist,        &
     &                     L_mono_moist, L_conserv_moist,               &
     &                     depart_r(1,1,tr_start),                      &
     &                     depart_lambda(1,1,tr_start),                 &
     &                     depart_phi(1,1,tr_start),                    &
     &                     me, n_proc, n_procx, n_procy,                &
     &                     halo_i, halo_j, global_row_length,           &
     &                     global_rows, row_length, rows,               &
     &                     l_datastart, at_extremity, g_i_pe,           &
     &                     proc_row_group, proc_col_group, 1,           &
     &                     off_x, off_y, L_sl_halo_reprod,              &
     &                     off_x, off_y,                                &
     &                     tracers, Error_Code)

!  set end of timestep tracer
        Do l=1,tr_vars
          Do k = 1, tr_levels
            Do j = j_start,j_end
              Do i = i_start,i_end
                tracers(i,j,k,l) = tracers(i,j,k,l) + tracer_phys2(i,j,k,l)
              End Do
            End Do
          End Do
        End Do

!   reset tracers to be non-negative
! DEPENDS ON: q_pos_ctl
        call Q_Pos_Ctl(                                                 &
     &                 tracers, row_length, rows,                       &
     &                 tr_levels*tr_vars,                               &
     &                 global_row_length, global_rows,                  &
     &                 me, n_proc, off_x, off_y,                        &
     &                 gc_all_proc_group,                               & 
     &                 model_domain,                                    &
     &                 halo_type_single, .true. ,  0.0                  &
     &                 )

      End If    ! tr_vars >0 and tr_levels /= model_levels

      return
      END SUBROUTINE SL_tracer2
