#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE CHECK_IDEALISE                                         &
     &  ( LEN_INTHD, LEN_REALHD, LEN_LEVDEPC1, LEN_LEVDEPC2,            &
     &  INTHD, REALHD, LEVDEPC )

      Implicit None

!  Subroutine  - check the idealised namelist.
!
! Description:
!   Read the namelist, assigning and freeing Fortran units as
!   required.
!
! Method:
!   The namelists is provided in one file
!
! Current Code Owner: T. Davies
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   12/05/01   Set up to read idealise namelist. T. Davies
!   5.4   29/08/02   Add variables for convect-equil    Carol Roadnight
!  5.4  29/07/02   More appropriate formatting for some prints
!                                                         T. Davies
!  5.4  28/08/02   Move L_free_slip to Run_dyn namelist     Andy Malcolm
!  5.4  28/08/02   set unit number to read file IDEALISE    Andy Malcolm
!  5.5  03/02/03   Include checks on u and v profiles. R.M.Forbes
!  5.5  20/02/03   Add 'Fmt=' keyword to Write statements   P.Dando
!  6.0   18/08/03  Add extra control variables    Terry Davies
!  6.1   01/08/04  Add idealised forcing variables  C.Halliwell/R.Forbes
!  6.2   01/03/06  Add idealised pressure forcing variables. Y.M.Tang
!  6.1   02/08/04  Fix minor bug                          Terry Davies
!  6.2   13/01/06  Add control logicals for idealised baroclinic wave
!                  and cyclone options                    Bob Beare
!  6.2   31/01/06  Add idealised bubble variables. R.Forbes
!  6.2   31/01/06  Add idealised perturbation variables. R.Forbes
!  6.2   05/01/06  Setup new variables for uv profile heights. R.Forbes
!  6.2   21/02/06  Add lbc forcing     Yongming Tang
!  6.2  31/01/06  Initialise IdlSurfFluxSea variables. R.Forbes
!  6.2   21/03/06  Included nstypes.h J Ridley.
!  6.2   14/02/06  Default setting of dynamics control variables
!                     moved to deck READLSA2             Terry Davies
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Arguments
      Integer :: Len_IntHd
      Integer :: Len_RealHd
      Integer :: Len_LevDepC1
      Integer :: Len_LevDepC2

      Integer :: IntHd   (Len_IntHd)
      Real    :: RealHd  (Len_RealHd)
      Real    :: LevDepC (Len_LevDepC1, Len_LevDepC2)

! Local Variables/Paramters
      Character (Len=*), Parameter :: RoutineName =                     &
     &                                'Check_Idealise'
      Integer, Parameter           :: max_levels = 160
      Integer, Parameter           :: jtheta     = 1
      Integer, Parameter           :: jrho       = 2

      Character (Len=80)  :: Cmessage
      Character (Len=80)  :: FileName
      Integer             :: ErrorStatus
      Integer             :: status
      Integer             :: model_levels
      Integer             :: i             ! looper
      Integer             :: j             ! looper
      Integer             :: nft
      Logical             :: l_exist

! Comdecks
#include "nstypes.h"
#include "c_mdi.h"
#include "cprintst.h"
#include "parvars.h"
#include "cmaxsize.h"
#include "cntlatm.h"
#include "cruntimc.h"
#include "problem.h"
#include "tprofile.h"
#include "qprofile.h"
#include "uvhoriz.h"


      Namelist /Idealise/ L_idealised_data,                             &
     &       L_initialise_data, surface_type, grid_number, grid_flat,   &
     &       first_theta_height, thin_theta_height, height_domain,      &
     &       tprofile_number, qprofile_number, uvprofile_number,        &
     &       theta_surface, p_surface, Brunt_Vaisala,                   &
     &       dtheta_dz1, height_dz1,                                    &
     &       t_horizfn_number, t_horizfn_data,                          &
     &       uv_horizfn_number, u_in, v_in, height_u_in, q1,            &
     &       u_ramp_start, u_ramp_end, ujet_lat, ujet_width, r_plane,   &
     &       h_o, grow_steps, Witch_power,                              &
     &       lambda_fraction, phi_fraction,                             &
     &       half_width_x, half_width_y, plat_size_x, plat_size_y,      &
     &       first_constant_r_rho_level_new, L_constant_dz,             &
     &       L_polar_wind_zero, L_rotating,                             &
     &       L_trivial_trigs, f_plane, ff_plane, L_vert_Coriolis,       &
     &       L_fixed_lbcs, L_pressure_balance,                          &
     &       L_wind_balance, L_rotate_winds,                            &
     &       IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                &
     &       L_spec_z0, roughlen_z0m, roughlen_z0h,                     &
     &       L_inviscid, L_deep, L_hydrostatic,                         &
     &       L_geostrophic, L_solver, L_trap_errors, L_uv_zero,         &
     &       big_layers, transit_layers, mod_layers, big_factor, mag,   &
     &       SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,  &
     &       SuHe_pole_equ_deltaT, SuHe_static_stab,                    &
     &       base_frictional_timescale, SuHe_sigma_cutoff,              &
     &       SuHe_relax, SuHe_fric,                                     &
     &       L_SH_williamson, Instability_diagnostics,                  &
     &       L_simple_friction, tropics_deg,                            &
     &       num_profile_data,                                          &
     &       zprofile_data, tprofile_data, qprofile_data,               &
     &       num_uvprofile_data,                                        &
     &       z_uvprofile_data, uprofile_data, vprofile_data,            &
     &       tforce_option, qforce_option, uvforce_option,              &
     &       num_tforce_levels, num_tforce_times,                       &
     &       tforce_time_interval,                                      &
     &       num_qforce_levels, num_qforce_times,                       &
     &       qforce_time_interval,                                      &
     &       num_uvforce_levels, num_uvforce_times,                     &
     &       uvforce_time_interval,newtonian_timescale,                 &
     &       z_tforce_data, tforce_data,                                &
     &       z_qforce_data, qforce_data,                                &
     &       z_uvforce_data, uforce_data, vforce_data,                  &
     &       tforce_data_modlev, qforce_data_modlev,                    &
     &       uforce_data_modlev, vforce_data_modlev,                    &
     &       pforce_option,                                             &
     &       num_pforce_times, pforce_time_interval,                    &
     &       p_surface_data, L_pforce,                                  &
     &       L_perturb_t, perturb_magnitude_t                           &
     &,      L_perturb_q, perturb_magnitude_q                           &
     &,      L_perturb_correlate_tq                                     &
     &,      L_perturb_correlate_vert                                   &
     &,      L_perturb_correlate_time                                   &
     &,      perturb_type, perturb_height                               &
     &,      L_code_test, L_force, cool_rate                            &
     &,      L_cyclone, L_baroclinic                                    &
     &,       L_fix_orog_hgt_lbc, orog_hgt_lbc, L_force_lbc             &
     &,       idl_interp_option, zprofile_orog, hf                      &
     &,      idl_bubble_option, idl_bubble_max                          &
     &,      idl_bubble_width,  idl_bubble_depth, idl_bubble_height     &
     &,      idl_bubble_xoffset,idl_bubble_yoffset                      &
     &,      L_idl_bubble_saturate                                      &
     &,      L_damp, L_geo_for,  L_bomex                                &
     &,      DMPTIM, HDMP, ZDMP                                         &
     &,      u_geo, v_geo

      model_levels = IntHd(8)

      nft=106      ! set unit number explicitly

! Header values used in subroutine
!     IntHd(13)    : Number of boundary layer levels
!     IntHd(24)    : First rho level at which height is constant

! Set model defaults
      L_idealised_data = .false.
      L_initialise_data = .false.
      Instability_diagnostics = 0
      L_simple_friction = .false.
       L_inviscid = .false.
       L_deep = .true.
       L_hydrostatic = .false.
       L_geostrophic = .false.
       L_solver = .true.
       L_trap_errors = .false.
       L_uv_zero = .false.
       tropics_deg = 30.0
      L_code_test = .false.

! Set IDEALISE defaults
      surface_type = 10
      grid_number = 10
      first_theta_height = 10.0
      thin_theta_height = 1.0
      grid_flat = 3
      first_constant_r_rho_level_new = -1
      tprofile_number = 10
      qprofile_number = 10
      uvprofile_number = 10
      lambda_fraction = 0.5
      phi_fraction = 0.5
      h_o = 0.0
      grow_steps = 0
      plat_size_x = 0.0
      plat_size_y = 0.0
      Witch_power = 1.5
! Half-widths for Witch of Agnesi.  Also use to mask real orography
      half_width_x = 2500000.0
      half_width_y = 2500000.0
      theta_surface = 280.0
      p_surface = 100000.0
      height_domain = 10000.
      Brunt_Vaisala =  0.01
      L_constant_dz = .true.
      L_trivial_trigs = .false.
      L_vert_Coriolis = .false.
      r_plane = -90.0
      f_plane = -90.0
      ff_plane = -90.0
      L_rotating = .true.
      L_fixed_lbcs = .false.
      L_force_lbc = .false.
      L_fix_orog_hgt_lbc = .false.
      orog_hgt_lbc = 0.0
      idl_interp_option = 1  ! Interpolation -> const on height levels
      zprofile_orog = 0.0    ! Sea level
      hf = 0.0
      L_wind_balance = .false.
      L_rotate_winds = .false.
      L_pressure_balance = .false.
      L_polar_wind_zero= .false.
      L_perturb = .false.
      perturb_factor = 1.0
      L_perturb_t         = .False.
      perturb_magnitude_t = 0.5      ! Kelvin
      L_perturb_q         = .False.
      perturb_magnitude_q = 0.5E-3   ! kg/kg
      L_perturb_correlate_tq   = .True.
      L_perturb_correlate_vert = .True.
      L_perturb_correlate_time = .True.
      perturb_type        = 1        ! random
      perturb_height      = 0.0
      L_force = .false.
      L_cyclone = .false.
      L_baroclinic = .false.
      L_damp = .false.
      L_geo_for = .false.
      L_bomex = .false.
      DMPTIM = 0.0
      HDMP = 0.0
      ZDMP = 0.0 
      u_geo = 0.0 
      v_geo = 0.0 

      cool_rate = 0.0
      q1 = 70.0  ! 70% relative humidity
      t_horizfn_number = 0
      Do i = 1, 10
        t_horizfn_data(i) = 0.0
      End Do
      uv_horizfn_number = 0
      Do i = 1, 4
        u_in(i) = 0.0
        v_in(i) = 0.0
      End Do
      u_ramp_start = 0.0
      u_ramp_end = -90.0
      ujet_lat = -90.0
      ujet_width = 0.0
      Do i = 1, 3
        height_u_in(i) = 0.0  ! dummy value - not used for constant u
        dtheta_dz1(i) = 0.01
        height_dz1(i) = 0.0
      End Do
      big_layers = 0
      transit_layers = 0
      mod_layers = 0
      big_factor = 1.0
      mag = 1.0
      num_profile_data   = 0
      Do i = 1, max_num_profile_data
        zprofile_data(i) = 0.0
        tprofile_data(i) = 0.0
        qprofile_data(i) = 0.0
      End Do
      num_uvprofile_data = 0
      Do i = 1, max_num_profile_data
        z_uvprofile_data(i) = 0.0
        uprofile_data(i)    = 0.0
        vprofile_data(i)    = 0.0
      End Do
      pforce_option = 0
      num_pforce_times = 1
      pforce_time_interval = 600.0
      L_pforce = .false.
      tforce_option = 0
      qforce_option = 0
      uvforce_option = 0
      num_tforce_levels = 1
      num_tforce_times = 1
      tforce_time_interval = 600.0
      num_qforce_levels = 1
      num_qforce_times = 1
      qforce_time_interval = 600.0
      num_uvforce_levels = 1
      num_uvforce_times = 1
      uvforce_time_interval = 600.0
      newtonian_timescale = 3600.0
      Do i = 1, max_num_profile_data
        z_tforce_data(i) = 0.0
        z_qforce_data(i) = 0.0
        z_uvforce_data(i) = 0.0
        Do j = 1, max_num_force_times
          tforce_data(i,j)=0.0
          qforce_data(i,j)=0.0
          uforce_data(i,j)=0.0
          vforce_data(i,j)=0.0
        End Do
      End Do
      Do j = 1, max_num_force_times
        Do i = 1, max_model_levels
          tforce_data_modlev(i,j)= 0.0
          qforce_data_modlev(i,j)= 0.0
          uforce_data_modlev(i,j)= 0.0
          vforce_data_modlev(i,j)= 0.0
        End Do
      End Do
      Do j = 1, max_num_force_times
        p_surface_data(j)=0.0
      End Do

      Do i=1,idl_max_num_bubbles
        idl_bubble_option(i) = 0      ! Default no bubbles
        idl_bubble_max(i)    = 1.0    ! 1 K
        idl_bubble_height(i) = 1000.  ! 1 km
        idl_bubble_xoffset(i)= 0.5    ! Centre of domain
        idl_bubble_yoffset(i)= 0.5    ! Centre of domain
        idl_bubble_width(i)  = 1000.  ! 1 km
        idl_bubble_depth(i)  = 1000.  ! 1 km
        L_idl_bubble_saturate(i) = .False.
      End Do
! Set dynamical core defaults
      L_SH_williamson = .false.
      SuHe_pole_equ_deltaT = 60.
      SuHe_static_stab = 10.
      SuHe_sigma_cutoff = 0.7
      base_frictional_timescale = 1.1574e-5
      SuHe_newtonian_timescale_ka = 2.893e-7
      SuHe_newtonian_timescale_ks = 2.893e-6
      SuHe_relax = 2
      SuHe_fric = 2

      IdlSurfFluxSeaOption = 0
      IdlSurfFluxSeaParams(:) = 0.0

! defaults for specification of roughness length
      L_spec_z0 = .false.
      roughlen_z0m = 0.0
      roughlen_z0h = 0.0

      if(Problem_number  /=  0) then

! Now find the namelist filename from Env Vars

      Call Fort_Get_Env( 'IDEALISE', 8, FileName, 80, ErrorStatus )

      If ( ErrorStatus /= 0 ) Then
        ErrorStatus = 10
        Cmessage =                                                      &
     &  'Unable to Obtain Idealise Filename from Environment'
! DEPENDS ON: ereport
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

      FileName = Trim( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together

      Inquire( file=FileName, exist=l_exist )

      If ( .Not. l_exist ) Then
        Write (6,*) 'Idealise file: ',FileName
        ErrorStatus = 20
        Cmessage = ' Idealise Namelist file does not exist!'
! DEPENDS ON: ereport
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

! Open the file containing Idealised model settings
      Open( Unit=nft, File=FileName, IOstat=status )

      If ( PrintStatus >= PrStatus_Oper ) Then
        Write (6,*) '-Idealise settings file: ',FileName
      End If

! Quick error check to make sure parameter max_levels isn't too small
      If ( max_levels < model_levels ) Then
        ErrorStatus = 10
        Cmessage = 'Internal paramter max_levels is too small!'
! DEPENDS ON: ereport
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

! Read Idealise Namelist
      Read( Unit = nft, Nml=Idealise )

!     Write out namelist variables for diagnostic
        Write ( Unit = 6, fmt=*) 'Values in IDEALISE Namelist.'
        Write ( Unit = 6, fmt=*) 'L_initialise_data ',L_initialise_data
        If( L_inviscid ) then
          Write ( Unit = 6, fmt=*) 'Inviscid lower boundary condition ' &
     &,  ' on theta L_inviscid = ',L_inviscid
        EndIf ! L_inviscid
        If( L_deep ) then
          Write ( Unit = 6, fmt=*) 'Deep atmosphere equations '         &
     &,  ' being used '
        else
          Write ( Unit = 6, fmt=*) 'SHALLOW atmosphere equations '      &
     &,  ' being used '
        EndIf ! L_deep
        If( L_hydrostatic ) then
          Write ( Unit = 6, fmt=*) 'Hydrostatic switch is ACTIVE '
        EndIf ! L_hydrostatic
        If( .not. L_solver ) then
          Write ( Unit = 6, fmt=*) 'Solver is not to be called '
        EndIf ! .not. L_solver
        If( L_geostrophic ) then
          Write ( Unit = 6, fmt=*) 'Geostrophic switch is ACTIVE '
        EndIf ! L_geostrophic
        If( L_trap_errors ) then
          Write ( Unit = 6, fmt=*) 'Error trapping  is ACTIVE '
        EndIf ! L_trap_errors
        If( L_uv_zero ) then
          Write ( Unit = 6, fmt=*) ' Setting u,v zero is ACTIVE '
        EndIf ! L_uv_zero
        Write ( Unit = 6, fmt=*) 'L_code_test ',L_code_test
        If( L_initialise_data )then
        Write ( Unit = 6, fmt=*) 'surface_type ',surface_type
        Write ( Unit = 6, fmt=*) 'grid_number ',grid_number
        Write ( Unit = 6, fmt=*) 'grid_flat ',grid_flat
        Write ( Unit = 6, fmt=*) 'height_domain ',height_domain
        Write ( Unit = 6, fmt=*)'first_theta_height ',first_theta_height
        Write ( Unit = 6, fmt=*) 'thin_theta_height ',thin_theta_height
        Write ( Unit = 6, fmt=*) 'first_constant_r_rho_level_new ',     &
     &                          first_constant_r_rho_level_new
        Write ( Unit = 6, fmt=*) 'big_layers ',big_layers
        Write ( Unit = 6, fmt=*) 'transit_layers ',transit_layers
        Write ( Unit = 6, fmt=*) 'mod_layers ',mod_layers
        Write ( Unit = 6, fmt=*) 'big_factor ',big_factor
        Write ( Unit = 6, fmt=*) 'mag ',mag
        Write ( Unit = 6, fmt=*) 'tprofile_number ',tprofile_number
        Write ( Unit = 6, fmt=*) 'qprofile_number ',qprofile_number
        Write ( Unit = 6, fmt=*) 'uvprofile_number ',uvprofile_number
        Write ( Unit = 6, fmt=*) 'theta_surface ',theta_surface
        Write ( Unit = 6, fmt=*) 'p_surface ',p_surface
        Write ( Unit = 6, fmt=*) 'Brunt_Vaisala ',Brunt_Vaisala
        Write ( Unit = 6, fmt=*) 'dtheta_dz1'
        Write ( Unit = 6, Fmt='(3F10.7)' )( dtheta_dz1(i),i=1,3 )
        Write ( Unit = 6, fmt=*) 'height_dz1'
        Write ( Unit = 6, Fmt='(3F10.3)' )( height_dz1(i),i=1,3 )
        Write ( Unit = 6, fmt=*) 't_horizfn_number ',t_horizfn_number
        If (t_horizfn_number  /=  0) Then
          Write ( Unit = 6, fmt=*) 't_horizfn_data'
          Write ( Unit = 6, Fmt='(10F10.7)' )(t_horizfn_data(i),i=1,10)
        End If
        Write ( Unit = 6, fmt=*) 'uv_horizfn_number ',uv_horizfn_number
        Write ( Unit = 6, fmt=*) 'height_u_in'
        Write ( Unit = 6, Fmt='(3F10.3)' )( height_u_in(i),i=1,3 )
        Write ( Unit = 6, fmt=*) 'u_in'
        Write ( Unit = 6, Fmt='(4F10.3)' )( u_in(i),i=1,4 )
        Write ( Unit = 6, fmt=*) 'v_in'
        Write ( Unit = 6, Fmt='(4F10.3)' )( v_in(i),i=1,4 )
        Write ( Unit = 6, fmt=*) 'q1 ',q1
        Write ( Unit = 6, fmt=*) 'orog_height h_o ',h_o
        Write ( Unit = 6, fmt=*) 'grow_steps ',grow_steps
        Write ( Unit = 6, fmt=*) 'Witch_power ',Witch_power
        Write ( Unit = 6, fmt=*) 'lambda_fraction ',lambda_fraction
        Write ( Unit = 6, fmt=*) 'phi_fraction ',phi_fraction
        Write ( Unit = 6, fmt=*) 'half_width_x ',half_width_x
        Write ( Unit = 6, fmt=*) 'half_width_y ',half_width_y
        Write ( Unit = 6, fmt=*) 'plat_size_x ',plat_size_x
        Write ( Unit = 6, fmt=*) 'plat_size_y ',plat_size_y
        Write ( Unit = 6, fmt=*) 'L_constant_dz ',L_constant_dz
          Write ( Unit = 6, fmt=*) 'u_ramp_start ',u_ramp_start
          Write ( Unit = 6, fmt=*) 'u_ramp_end ',u_ramp_end
        Write ( Unit = 6, fmt=*) 'ujet_lat ',ujet_lat
        Write ( Unit = 6, fmt=*) 'ujet_width ',ujet_width
        Write ( Unit = 6, fmt=*) 'r_plane ',r_plane
        Write ( Unit = 6, fmt=*) 'f_plane ',f_plane
        Write ( Unit = 6, fmt=*) 'ff_plane ',ff_plane
        Write ( Unit = 6, fmt=*) 'L_trivial_trigs ',L_trivial_trigs
        Write ( Unit = 6, fmt=*) 'L_rotating ',L_rotating
        Write ( Unit = 6, fmt=*) 'L_vert_Coriolis ',L_vert_Coriolis
        Write ( Unit = 6, fmt=*) 'L_fixed_lbcs ',L_fixed_lbcs
        Write (Unit = 6, fmt=*) 'L_force_lbc ', L_force_lbc
        Write (Unit=6, fmt=*) 'L_fix_orog_hgt_lbc ',L_fix_orog_hgt_lbc
        Write (Unit = 6, fmt=*) 'orog_hgt_lbc ', orog_hgt_lbc
        Write (Unit = 6, fmt=*) 'idl_interp_option ',idl_interp_option
        Write (Unit = 6, fmt=*) 'zprofile_orog ', zprofile_orog
        Write (Unit = 6, fmt=*) 'hf ', hf
        Write (Unit = 6, fmt=*) 'p_surface_data ', p_surface_data
        Write ( Unit = 6, fmt=*) 'L_wind_balance ',L_wind_balance
        Write ( Unit = 6, fmt=*) 'L_rotate_winds ',L_rotate_winds
        Write ( Unit = 6, fmt=*)'L_pressure_balance ',L_pressure_balance
        Write ( Unit = 6, fmt=*) 'L_polar_wind_zero ',L_polar_wind_zero
        Write ( Unit = 6, fmt=*) 'L_perturb ',L_perturb
        Write ( Unit = 6, fmt=*) 'tropics_deg ',tropics_deg

        If (L_perturb) Then
          Write ( Unit = 6, fmt=*) 'perturb_factor ',perturb_factor
          Write ( Unit = 6, fmt=*) 'perturb_height ',perturb_height
        EndIf
        If (tprofile_number == tp_namelist                              &
     &      .or. qprofile_number == qp_namelist                         &
     &      .or. qprofile_number == qp_namelist_rh) Then

          Write (Unit=6, Fmt=*) 'num_profile_data ',num_profile_data

          ! Check to make sure data points in profile is less than max
          If ((num_profile_data == 0) .or.                              &
     &        (num_profile_data > max_num_profile_data)) Then
            Write (6,*) 'max_num_profile_data ',max_num_profile_data
            Write(Cmessage,*)                                           &
     &        'Idealised namelist vertical profile data:'               &
     &        //'Zero or too many points. '                             &
     &        //'num_profile_data must be 0 <= max_num_profile_data'
            ErrorStatus = 1
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          Write (Unit=6, Fmt=*) 'zprofile_data'
          Write (Unit=6, Fmt='(F10.3)')                                 &
     &              (zprofile_data(i),i=1,num_profile_data)
        End If
        If (tprofile_number == tp_namelist) Then
          Write (Unit=6, Fmt=*) 'tprofile_data'
          Write (Unit=6, Fmt='(F10.3)')                                 &
     &              (tprofile_data(i),i=1,num_profile_data)
        End If
        If (qprofile_number == qp_namelist                              &
     &      .or. qprofile_number == qp_namelist_rh) Then
          Write (Unit=6, Fmt=*) 'qprofile_data'
          Write (Unit=6, Fmt='(F10.7)')                                 &
     &              (qprofile_data(i),i=1,num_profile_data)
        End If


        If (uvprofile_number == uv_vert_namelist) Then

          ! If num_uvprofile_data not set in namelist then assume
          ! data is on same levels as t,q data
          If (num_uvprofile_data == 0 .AND. num_profile_data /= 0)      &
     &     Then
            Write (Unit=6, Fmt=*) 'Assuming uv data is on the same'     &
     &                          //' height levels as tq data.'
            num_uvprofile_data = num_profile_data
            Do i = 1, max_num_profile_data
              z_uvprofile_data(i) = zprofile_data(i)
            End Do
          End If

          Write (Unit=6, Fmt=*) 'num_uvprofile_data ',                  &
     &                                   num_uvprofile_data

          ! Check to make sure no. points in profile less than max
          If ((num_uvprofile_data == 0) .or.                            &
     &        (num_uvprofile_data > max_num_profile_data)) Then
            Write (6,*) 'max_num_profile_data ',max_num_profile_data
            Write(Cmessage,*)                                           &
     &        'Idealised namelist vertical profile data:'               &
     &        //'Zero or too many points. '                             &
     &        //'num_uv_profile_data must be 0 <= max_num_profile_data'
            ErrorStatus = 1
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          ! Write out arrays
          Write (Unit=6, Fmt=*) 'z_uvprofile_data'
          Write (Unit=6, Fmt='(F10.3)')                                 &
     &              (z_uvprofile_data(i),i=1,num_uvprofile_data)
          Write (Unit=6, Fmt=*) 'uprofile_data'
          Write (Unit=6, Fmt='(F10.4)')                                 &
     &              (uprofile_data(i),i=1,num_uvprofile_data)
          Write (Unit=6, Fmt=*) 'vprofile_data'
          Write (Unit=6, Fmt='(F10.4)')                                 &
     &              (vprofile_data(i),i=1,num_uvprofile_data)

        End If
        Write ( Unit = 6, fmt=*) 'L_force =',L_force
        Write ( Unit = 6, fmt=*) 'L_cyclone =',L_cyclone
        Write ( Unit = 6, fmt=*) 'L_baroclinic =',L_baroclinic
        Write ( Unit = 6, fmt=*) 'L_damp =',L_damp
        Write ( Unit = 6, fmt=*) 'L_geo_for =',L_geo_for
        Write ( Unit = 6, fmt=*) 'L_bomex =',L_bomex
        Write ( Unit = 6, fmt=*) 'Damping layer settings'
        Write ( Unit = 6, fmt=*) 'DMPTIM =',DMPTIM
        Write ( Unit = 6, fmt=*) 'HDMP =',HDMP
        Write ( Unit = 6, fmt=*) 'ZDMP =',ZDMP
        Write ( Unit = 6, fmt=*) 'Geostrophic forcings'
        Write ( Unit = 6, fmt=*) 'u_geo =',u_geo
        Write ( Unit = 6, fmt=*) 'v_geo =',v_geo
        ! -------------------------------------------------
        ! Theta forcing
        ! -------------------------------------------------

        If (tforce_option  ==  1) Then
          Write (Unit = 6, fmt=*)                                       &
     &          'Forcing increments added to theta field'
        Elseif (tforce_option  ==  2) Then
          Write ( Unit = 6, fmt=*)                                      &
     &          'Relaxation forcing for theta field'
! Option not yet implemented
!          Elseif (tforce_option  ==  3) Then
!            Write (Unit = 6, fmt=*)
!     &            'Theta reset to specified forcing data'
        Else
          Write ( Unit = 6, fmt=*)                                      &
     &     'No forcing of theta'
        Endif

        ! -------------------------------------------------
        ! Humidity forcing
        ! -------------------------------------------------

        If (qforce_option  ==  1) Then
          Write ( Unit = 6, fmt=*)                                      &
     &          'Forcing increments added to q field'
        Elseif (qforce_option  ==  2) Then
          Write ( Unit = 6, fmt=*)                                      &
     &          'Relaxation forcing for q field'
! Option not yet implemented
!          Elseif (qforce_option  ==  3) Then
!            Write ( Unit = 6, fmt=*)
!     &            'q reset to specified forcing data'
        Else
          Write ( Unit = 6, fmt=*) 'No forcing of q'
        Endif

        ! -------------------------------------------------
        ! Horizontal wind forcing
        ! -------------------------------------------------

        If (uvforce_option  ==  1) Then
          Write ( Unit = 6, fmt=*)                                      &
     &          'Forcing increments added to u and v fields'
        Elseif (uvforce_option  ==  2) Then
          Write ( Unit = 6, fmt=*)                                      &
     &          'Relaxation forcing for u and v fields'
! Option not yet implemented
!          Elseif (uvforce_option  ==  3) Then
!            Write ( Unit = 6, fmt=*)
!     &            'u and v reset to specified forcing data'
        Else
          Write ( Unit = 6, fmt=*)                                      &
     &       'No forcing of winds'
        Endif

      ElseIf( Problem_number  ==  dynamical_core)then

        Write ( Unit = 6, fmt=*) 'L_SH_williamson ',L_SH_williamson
        Write ( Unit = 6, fmt=*) 'SuHe_newtonian_timescale_ka ',        &
     &                         SuHe_newtonian_timescale_ka
        Write ( Unit = 6, fmt=*) 'SuHe_newtonian_timescale_ks ',        &
     &                         SuHe_newtonian_timescale_ks
        Write ( Unit = 6, fmt=*) 'SuHe_pole_equ_deltaT ',               &
     &                         SuHe_pole_equ_deltaT
        Write ( Unit = 6, fmt=*) 'SuHe_static_stab ',                   &
     &                         SuHe_static_stab
        Write ( Unit = 6, fmt=*) 'base_frictional_timescale ',          &
     &                         base_frictional_timescale
        Write ( Unit = 6, fmt=*) 'SuHe_sigma_cutoff ',                  &
     &                         SuHe_sigma_cutoff
        Write ( Unit = 6, fmt=*) 'SuHe_relax ',SuHe_relax
        Write ( Unit = 6, fmt=*) 'SuHe_fric ',SuHe_fric
        EndIf ! L_initialise_data
        Write ( Unit = 6, fmt=*) 'Instability_diagnostics ',            &
     &                          Instability_diagnostics
        Write ( Unit = 6, fmt=*) 'L_simple_friction ',L_simple_friction


      Close( Unit=nft )

      endif

      Return
      END SUBROUTINE CHECK_IDEALISE
#endif
