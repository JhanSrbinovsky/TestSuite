#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine SETCONA ------------------------------------------------
!LL
!LL Purpose : to calculate additional constants derived from
!LL      information in the input dump and pass out into the model via
!LL      the argument lists in argcona and arglndm.
!LL
!LL Method:
!LL    0. Initialise vertical coordinates from eta values and surface
!LL       orography.
!LL    0.1 Initialise data arrays held in secondary dump space.
!LL    1. Trigonometric functions:
!LL        Convert from degrees to radians.
!LL        Calculate trig functions for this grid.
!LL        Update halos for trig functions.
!LL    2. Set up semi-lagrangian advection options.
!LL       Initialise land/soil points in index arrays.
!LL    3. Determine model level boundaries for L/M/H cloud diagnostics.
!LL    4. Add locally derived parameters.
!LL
!LL Language: FORTRAN 77 + common extensions also in Fortran 90.
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version 7.2, dated 5/2/98
!LL
!LL Documentation : Unified Model Documentation Paper No P0
!LL
!LLEND----------------------------------------------------------------
!
!*L arguments

      SUBROUTINE Setcona                                                &
! Input data (including some logical fields)
     &     (eta_theta_in,eta_rho_in,Smvcst,Land_sea_mask,Orog,          &
     &      grad_x, grad_y,                                             &
     &      rho,exner_rho_levels,                                       &
     &      orog_lbc,exner_lbc,                                         &
! Input size and control variables
     &      LENRIMA,LBC_SIZEA,LBC_STARTA,                               &
     &      RIMWIDTHA, RIMWEIGHTSA,                                     &
     &      global_row_length, global_rows,                             &
     &      MODEL_LEVELS,ROWS,N_ROWS,ROW_LENGTH,                        &
     &      LAND_FIELD,WET_LEVELS,boundary_layer_levels,                &
     &      first_constant_r_rho_level,cloud_levels,z_top_of_model,     &
     &       tr_levels, tr_vars, tr_ukca,                               &
     &      height_gen_method,                                          &
! other constants
     &      Delta_lambda_in,Delta_phi_in,Base_phi_in,Base_lambda_in,    &
     &      lat_rot_NP_in,long_rot_NP_in,                               &
! VarRes grid info in degrees
     &      Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in,               &
     &      RowDepCStart,                                               &
! Initialise variables in secondary D1:
     &      exner_theta_levels,p_theta_levels,p,p_star,                 &
! Output data
#include "argcona.h"
#include "arglndm.h"
     &     icode,Cmessage,isSTASH)  ! FLUME-STASH

      Use solinc_data, Only: L_orog, slope_angle
Use level_heights_Mod
Use trignometric_Mod
Use dyn_coriolis_Mod
Use dyn_var_res_Mod
Use diff_coeff_Mod
Use rad_mask_trop_Mod
Use rot_coeff_Mod
Use volcts_Mod
Use scvary_Mod

      IMPLICIT NONE

!     INCLUDED COMDECKS
#include "nstypes.h"
#include "parvars.h"
#include "rimtypes.h"
#include "cmaxsize.h"
#include "cruntimc.h"
#include "problem.h"
#include "cprintst.h"
#include "natforce.h"

      INTEGER                                                           &
     &       LENRIMA(Nfld_max,NHalo_max,Nrima_max),                     &
                                                     ! IN LBC data len
     &       LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max),                 &
                                                       ! IN LBC size
     &       LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max),                &
                                                         ! IN LBC start
     &       RIMWIDTHA(Nrima_max),                                      &
                                    ! IN RIM width
     &       global_row_length,                                         &
                                 ! IN total number of point in a row
     &       global_rows,                                               &
                           ! IN total number of rows in model
     &       MODEL_LEVELS,ROWS,N_ROWS,ROW_LENGTH,ICODE,                 &
     &       LAND_FIELD,                                                &
                           ! IN Number of land points in model from umui
     &       WET_LEVELS,                                                &
                           ! IN Number of wet levels in model, from umui
     &       first_constant_r_rho_level,                                &
                                         !(IN) 1st constant height level
     &       boundary_layer_levels,                                     &
                                    ! (IN) Num. of boundary layer levels
     &       height_gen_method,                                         &
                                  ! (IN) Method for generating heights
     &       CLOUD_LEVELS                                               &
                           ! IN No of cloudy levels in the model
     &      ,tr_levels, tr_vars, tr_ukca                                &
     &      ,yvolc, mvolc, ysol

      INTEGER RowDepCStart
                        ! IN Start of Row dependent constant
      CHARACTER*(80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='SETCONA')

      REAL                                                              &
            ! Input arguments (grid constants)
     &       Delta_lambda_in                                            &
                             ! EW (x) grid spacing in degrees
     &      ,Delta_phi_in                                               &
                             ! NS (y) grid spacing in degrees
     &      ,Base_phi_in                                                &
                             ! Latitude of first theta point in degrees
     &      ,Base_lambda_in                                             &
                             ! Longitude of first theta point in degs
     &      ,lat_rot_NP_in                                              & 
                           ! Real latitude of 'pseudo' N pole in degs
     &      ,long_rot_NP_in                                             & 
                           ! Real longitude of 'pseudo' N pole in degs
     &      ,z_top_of_model                                             &
                            ! (IN) Height of top of model in metres
     &,      RIMWEIGHTSA(RIMWIDTHA(rima_type_norm))  ! IN RIM weights

      REAL                                                              &
            ! Input VarRes grid info in degrees
     &   Lambda_p_in(global_row_length)                                 &
                                               !IN EW and NS VarRes grid
     &  ,Lambda_u_in(global_row_length)                                 &
                                               !IN EW and NS u,v grid
     &  ,Phi_p_in(global_rows)                                          &
                                               !location in degrees
     &  ,Phi_v_in(global_rows)                 !location in degrees

      REAL                                                              &
            ! Array arguments with intent(IN):
     &     Smvcst(land_field)                                           &
                                      ! IN Volumetric saturation point
     &    ,eta_theta_in(0:model_levels)                                 &
                                        ! IN Eta values for theta levs
     &    ,eta_rho_in(model_levels)                                     &
                                        ! IN Eta values for rho levels
     &    ,Orog(row_length, rows)                                       &
                                        ! IN Orography (on all points)
     &    ,grad_x(row_length*rows)                                      &
                                        ! IN Orographic X-gradient
     &    ,grad_y(row_length*rows)                                      &
                                        ! IN Orographic Y-gradient
     &    ,rho(1-offx:row_length+offx, 1-offy:rows+offy,                &
     &       model_levels)                                              &
                           ! density*(r**2): used in call to Calc_P_star
     &    ,exner_rho_levels(1-offx:row_length+offx,                     &
                                                     ! Exner at p levels
     &                      1-offy:rows+offy, model_levels+1)           &
     &,      orog_lbc(LENRIMA(fld_type_p,halo_type_extended,            &
     &                        rima_type_orog))                          &
                                                 ! Orography LBC
     &,      exner_lbc(LENRIMA(fld_type_p,halo_type_extended,           &
                                                              !Exner LBC
     &                         rima_type_norm),MODEL_LEVELS+1)

      REAL                                                              &
            ! Output args (secondary arrays calculated from dump)
     &  p_theta_levels(1-offx:row_length+offx,                          &
                                               ! press on theta levs Pa
     &                 1-offy:rows+offy, model_levels)                  &
     &, p(1-offx:row_length+offx, 1-offy:rows+offy, model_levels+1)     &
                                                              ! in Pa
     &, p_star(row_length, rows)                                        &
                                 ! surface pressure (Pa)
     &, exner_theta_levels(1-offx:row_length+offx,                      &
                                                   ! Exner at theta levs
     &                     1-offy:rows+offy, model_levels)

      LOGICAL                                                           &
     &       Land_sea_mask(row_length,rows)

      REAL                                                              &
     &  orog_halos(1-halo_i:ROW_LENGTH+halo_i,                          &
     &           1-halo_j:ROWS+halo_j)

      LOGICAL       isSTASH  ! FLUME-STASH

#include "typcona.h"
#include "typlndm.h"
#include "cphyscon.h"
#include "c_eta_pmsl.h"
#include "acparm.h"
#include "mppac.h"
#include "cntlatm.h"
#include "trophgt1.h"
#include "highos.h"


! Local variables
      REAL                                                              &
     &   Delta_lambda_wk                                                &
                             ! EW (x) grid spacing in degrees
     &  ,Delta_phi_wk                                                   &
                             ! NS (y) grid spacing in degrees
     &  ,Base_phi_wk                                                    &
                             ! Latitude of first theta point in degrees
     &  ,Base_lambda_wk                                                 &
                             ! Longitude of first theta point in degrees
     &  ,latitude_rot(row_length, rows)                                 &
                                         ! rot. latit. in degrees (LAM)
     &  ,longitude_rot(row_length, rows)                                &
                                         ! rot. longit. in degrees (LAM)
     &  ,true_latitude_b(row_length, rows)                              &
                                           ! true latitude on B-grid
     &  ,true_longitude_b(row_length,rows)                              &
                                           ! true longitude on B-grid
     &  ,bear_rot_NP(row_length, rows)                                  & 
                                         ! Bearing of 'pseudo' N pole
     &  ,temp1                                                          &
                   ! Used in calculation of Coriolis terms (LAM only)
     &  ,temp2                                                          &
                   ! Used in calculation of Coriolis terms (LAM only)
     &  ,constant                                                       &
                   ! 1/kappa used to calc pressure on model levels
     &   , f_plane_rad                                                  &
                        ! f_plane latitude in radians
     &   , ff_plane_rad                                                 &
                         ! f_plane latitude in radians
     &   , f1_temp                                                      &
                    ! f1 term for f-plane or ff-plane
     &   , f2_temp                                                      &
                    ! f2 term for f-plane or ff-plane
     &   , f3_temp                                                      &
                    ! f3 term for f-plane or ff-plane
     &  , scale                                                         &
     &  ,CLOUD_BOUND(NUM_CLOUD_TYPES+1)                                 &
                                        ! boundaries of cloud types
     &  ,r_ref_theta(model_levels)                                      &
                                   ! Local dynamic array
     &  ,r_ref_rho(model_levels)   ! Local dynamic array

      LOGICAL                                                           &
     &   landmask(1-offx:row_length+offx, 1-offy:rows+offy)             &
     &  ,L_error                                                        &
     &  ,L_varres

      INTEGER                                                           &
              ! Mostly loop counters, but j also used for interp points
     &   I,KK,                                                          &
     &   J,GI,GJ,                                                       &
     &   j0,j1,k,                                                       &
     &   LEVEL                                                          &
                         ! Used to set up cloud type boundaries
     &  ,first_row                                                      &
     &  ,last_row                                                       &
     &  ,info                                                           &
     &  , active_levels

      INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
      INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

!  External subroutines called :
      EXTERNAL                                                          &
     &   EQTOLL                                                         &
     &  ,W_Coeff                                                        &
     &  ,Ereport                                                        &
     &  ,Calc_Exner_at_theta                                            &
     &  ,Calc_P_from_Exner                                              &
     &  ,Calc_P_star                                                    &
     &  ,Swap_Bounds
      EXTERNAL pole_bearing,aspang

! Set flag for orography correction in module solinc_data
      L_orog = L_use_orog_corr .or. L_use_grad_corr
! Grid type in dump: L_varres = T / F for variable / unif grid
      L_varres = .FALSE.
      If (RowDepCStart >= 1) L_varres = .TRUE.
! ----------------------------------------------------------------------
! 0. Set-up vertical co-ordinate arrays
! ----------------------------------------------------------------------
! Set up orog_tmp - copy of orog but with wide halos
      DO j=1-halo_j,ROWS+halo_j
        DO i=1-halo_i,ROW_LENGTH+halo_i
          orog_halos(i,j)=0.0
        END DO
      END DO

      DO j=1,ROWS
        DO i=1,ROW_LENGTH
          orog_halos(i,j)=orog(i,j)
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(orog_halos,row_length,rows,1,                    &
     &                 halo_i,halo_j,fld_type_p,.FALSE.)

! ----------------------------------------------------------------------
! 0.1 Write out some useful information at start of run
! ----------------------------------------------------------------------

#if defined(REPROD)
      write(6,*) '  '
      write(6,*) '  ****  This run uses bit-reproducible code  ****'
#else
      write(6,*) '  '
      write(6,*) '  ****  This run uses fast code, therefore  ****'
      write(6,*) '  ****  different processor configurations  ****'
      write(6,*) '  ****  might not bit-reproduce             ****'
#endif

      write(6,*) '  '
      write(6,*) '  ****  PE configuration for this run  ****   '
      if ( nproc_x > 1 .and. nproc_y > 1 ) then
        write(6,*) '      ', nproc_x, ' processors East-West ',         &
     &                     nproc_y, ' processors North South'
      else if ( nproc_x > 1 .and. nproc_y == 1 ) then
        write(6,*) '      ', nproc_x, ' processors East-West ',         &
     &                     nproc_y, ' processor North South'
      else if ( nproc_x == 1 .and. nproc_y > 1 ) then
        write(6,*) '      ', nproc_x, ' processor East-West ',          &
     &                     nproc_y, ' processors North South'
      else ! nproc_x = 1 .and. nproc_y =1 ) then
        write(6,*) '       Single processor 1x1 '
      end if ! nproc_x > 1 .and. nproc_y > 1 )

      if ( offx >= halo_i .or. offy >= halo_j ) then
        if (.not.isSTASH ) then                 ! FLUME-STASH
          write(6,*) '  '
          write(6,*) '  ****  DANGER  ****   '
          write(6,*) '  *** SMALL halo >= LARGE halo  *****   '
          write(6,*) '  This could result in overwriting  '
          ICODE    = 20
          CMESSAGE ='SETCONA: Large halo size is too small'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        end if
      end if ! offx >= halo_i .or. offy >= halo_j

      If (L_mix_ratio) Then
        write(6,*) '  '
        write(6,*) '  ***  This run uses mixing ratios for      ***'
        write(6,*) '  ***  the moist variables in the dynamics  ***'
      Else
        write(6,*) '  '
        write(6,*) '  ***  This run uses specific humidities for  ***'
        write(6,*) '  ***  the moist variables in the dynamics    ***'
      End If !  L_mix_ratio

      If (model_domain  ==  mt_LAM) Then
      
        write(6,*) '  '
        Write ( Unit = 6, fmt=*) ' ***   LAM set-up ***' 
        
        Write ( Unit = 6, fmt=*) 'LBC frame size for solver, '          &
     &,                          ' n_rims_to_do = ', n_rims_to_do
        If( L_LBC_balance ) then
          Write ( Unit = 6, fmt=*) 'L_LBC_balance = ', L_LBC_balance    &
     &,          ' Impose vertically balanced Exner pressures'          &
     &,          ' and rho and set w=0 in lbcs '
        Else !  L_LBC_balance = .false.
          Write ( Unit = 6, fmt=*) 'L_LBC_balance = ', L_LBC_balance    &
     &,                            ' No balancing of lbcs '
        EndIf ! L_LBC_balance      
        If( L_lbc_new ) then
          Write ( Unit = 6, fmt=*) 'L_lbc_new = ', L_lbc_new            &
     &,                              'Use new lbc algorithm '
        Else !  L_lbc_new = .false.
          Write ( Unit = 6, fmt=*) 'L_lbc_new = ', L_lbc_new            &
     &,                              'Use old lbc algorithm '
        EndIf ! L_lbc_new
        If( L_int_uvw_lbc) then
          Write ( Unit = 6, fmt=*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
     &,          ' Advecting winds interpolated in lateral boundaries '
        Else !  L_int_uvw_lbc = .false.
          Write ( Unit = 6, fmt=*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
     &,                ' Extrapolated advecting winds from lbc file '
        EndIf ! L_int_uvw_lbc
        write(6,*) '  '

! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &    ROW_LENGTH,ROWS,halo_i,halo_j,1,fld_type_p,orog_halos,        &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_orog),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_orog),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_orog),   &
     &    halo_i,halo_j,orog_lbc,RIMWIDTHA(rima_type_orog),             &
     &    RIMWIDTHA(rima_type_orog),RIMWEIGHTSA,                        &
     &    at_extremity,                                                 &
     &    .FALSE.,.TRUE.)

! DEPENDS ON: set_lateral_boundaries
        CALL SET_LATERAL_BOUNDARIES(                                    &
     &    ROW_LENGTH,ROWS,Offx,Offy,MODEL_LEVELS+1,fld_type_p,          &
     &    exner_rho_levels,                                             &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i,halo_j,exner_lbc,RIMWIDTHA(rima_type_norm),            &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    at_extremity,                                                 &
     &    .FALSE.,.TRUE.)

      END IF ! IF (model_domain  ==  mt_LAM)

! Set reference height profile

!!! Allocate arrays for level_heights_mod module
      If (.not.allocated(eta_theta_levels)) Then
        Allocate (eta_theta_levels(0:model_levels))
      End If
      If (.not.allocated(eta_rho_levels)) Then
        Allocate (eta_rho_levels(model_levels))
      End If
      If (.not.allocated(r_theta_levels)) Then
        Allocate (r_theta_levels(1-halo_i:row_length+halo_i,           &
     &                           1-halo_j:rows+halo_j, 0:model_levels))
      End If
      If (.not.allocated(r_rho_levels)) Then
        Allocate (r_rho_levels(1-halo_i:row_length+halo_i,             &
     &                         1-halo_j:rows+halo_j, model_levels))
      End If

      eta_theta_levels(0) = eta_theta_in(0)
      Do k = 1, model_levels
        eta_theta_levels(k) = eta_theta_in(k)
        eta_rho_levels(k) = eta_rho_in(k)
        r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
      End Do
      Do k = 1, model_levels
        r_ref_rho(k) = eta_rho_levels(k) * z_top_of_model
      End Do
! set bottom level, ie: orography
      Do j = 1-halo_j, rows+halo_j
        Do i= 1-halo_i, row_length+halo_i
          r_theta_levels(i,j,0) = Orog_halos(i,j) + Earth_radius
        End Do
      End Do
! For constant levels set r to be a constant on the level
      Do k = first_constant_r_rho_level, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = Earth_radius + r_ref_theta(k)
            r_rho_levels(i,j,k) = Earth_radius + r_ref_rho(k)
          End Do
        End Do
      End Do

      Select Case( height_gen_method )
        Case( height_gen_original )
! The original version of height generation used in the SI dynamics
!
! For boundary layer levels set depth to be constant.
      Do k = 1, boundary_layer_levels
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = r_theta_levels(i,j,0)               &
     &                                 + r_ref_theta(k)
            r_rho_levels(i,j,k) = r_theta_levels(i,j,0)                 &
     &                               + r_ref_rho(k)
          End Do
        End Do
      End Do
! For intermediate levels use linear relaxation to constant value.
! set orographic heights.
      Do k = boundary_layer_levels+1,                                   &
     &       first_constant_r_rho_level-1
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
              r_rho_levels(i,j,k) =                                     &
     &          ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
     &            r_theta_levels(i,j,boundary_layer_levels) ) *         &
     &          ( eta_rho_levels(k) -                                   &
     &            eta_theta_levels(boundary_layer_levels) ) /           &
     &          ( eta_rho_levels(first_constant_r_rho_level) -          &
     &            eta_theta_levels(boundary_layer_levels) )             &
     &           +  r_theta_levels(i,j,boundary_layer_levels)
            r_theta_levels(i,j,k) =                                     &
     &          ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
     &            r_theta_levels(i,j,boundary_layer_levels) ) *         &
     &          ( eta_theta_levels(k) -                                 &
     &            eta_theta_levels(boundary_layer_levels) ) /           &
     &          ( eta_rho_levels(first_constant_r_rho_level) -          &
     &            eta_theta_levels(boundary_layer_levels) )             &
     &           +  r_theta_levels(i,j,boundary_layer_levels)
          End Do
        End Do
      End Do

        Case( height_gen_smooth )
! A smooth quadratic height generation
          Do k = 1, first_constant_r_rho_level-1
            Do j = 1-halo_j, rows+halo_j
              Do i= 1-halo_i, row_length+halo_i
              r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +&
     &         Earth_radius + Orog_halos(i,j) * (1.0 - eta_rho_levels(k)&
     &              /eta_rho_levels(first_constant_r_rho_level))**2
              r_theta_levels(i,j,k) = eta_theta_levels(k) *             &
     &             z_top_of_model + Earth_radius + Orog_halos(i,j) *    &
     &             (1.0 - eta_theta_levels(k) /                         &
     &              eta_rho_levels(first_constant_r_rho_level))**2
              End Do
            End Do
          End Do

        Case Default
          icode = 10
          Write (Cmessage,*) 'Unrecognised height generation method - ',&
     &                       'Dump needs to be reconfigured'
! DEPENDS ON: ereport
          Call Ereport( RoutineName, icode, Cmessage )
      End Select

! 0.1 Initialise secondary arrays.
! Exner at p (=rho) levels is obtained from dump.
! calculate p from exner_rho_levels
! [halos required for diagnostic calculations at T+0 in INITDIAG.]
      constant = 1./ kappa
      Do k = 1, model_levels+1
        Do j = 1-offy, rows+offy
          Do i = 1-offx, row_length+offx
            p(i,j,k)= (exner_rho_levels(i,j,k) ** constant)             &
     &                      * p_zero
          End Do
        End Do
      End Do

! calculate exner at theta_levels which is then used to get
! p at theta_levs
! DEPENDS ON: calc_exner_at_theta
      Call Calc_Exner_at_theta( r_theta_levels, r_rho_levels,           &
     &                exner_rho_levels,                                 &
     &                 row_length, rows, model_levels,                  &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 exner_theta_levels,.TRUE.)

! Calculate pressure from Exner at theta levels.

! DEPENDS ON: calc_p_from_exner
       call Calc_P_from_Exner(                                          &
     &                      p_theta_levels, kappa, p_zero,              &
     &                      row_length, rows, model_levels,             &
     &                      offx, offy,                                 &
     &                      exner_theta_levels,.TRUE.)

! calculate p_star using rho (from dump) and p on model levels
! DEPENDS ON: calc_p_star
      Call Calc_P_star (r_theta_levels, r_rho_levels, p, rho,           &
     &                  g, row_length, rows, model_levels,              &
     &                  offx, offy, halo_i, halo_j,                     &
     &                  p_star)

! ----------------------------------------------------------------------
! Initialise Constants and trigonometric functions.
! ----------------------------------------------------------------------

      If( L_rotating) Then
        OMEGA = 7.292116E-5 ! Angular speed of Earth's rotation
                            ! = 2*pi/siderial day (23h56m04s)
      Else
        OMEGA = 0.0         ! Zero angular speed of rotation for planet
      End If    ! L_rotating
      two_omega = 2. * OMEGA

! Convert grid-spacing and base values from degrees to radians
      If( L_regular ) then ! 
        base_lambda_wk = base_lambda_in
        base_phi_wk = base_phi_in
        delta_lambda_wk = delta_lambda_in
        delta_phi_wk = delta_phi_in
      Else
        base_lambda_wk = lambda_p_in(1)
        base_phi_wk = phi_p_in(1)
        delta_lambda_wk = (lambda_p_in(global_row_length) -             &
     &                    lambda_p_in(1))/(global_row_length - 1) 
        delta_phi_wk = (phi_p_in(global_rows) -                         &
     &                    phi_p_in(1))/(global_rows - 1) 
      End If 
      delta_lambda=delta_lambda_wk * Pi_over_180
      delta_phi=delta_phi_wk * Pi_over_180
      base_phi=base_phi_wk * Pi_over_180
      base_lambda=base_lambda_wk * Pi_over_180
      f_plane_rad = f_plane * Pi_over_180
      lat_rot_NP =lat_rot_NP_in  * Pi_over_180
      long_rot_NP=long_rot_NP_in * Pi_over_180

      If ( L_regular ) Then 

        Write ( Unit = 6, fmt=*) ' '
        Write ( Unit = 6, fmt=*) '*** The horizontal grid is regular'   &
     &                    , '  L_regular = ', L_regular, ' ***'
        Write ( Unit = 6, fmt=*) ' '

!   allocate small arrays to avoid out-of-bounds in atm_step 

        If (.not.allocated(glambda_p)) Then
          Allocate (glambda_p ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(glambda_u)) Then
          Allocate (glambda_u ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(gdlambda_p)) Then
          Allocate (gdlambda_p ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(gdlambda_u)) Then
          Allocate (gdlambda_u ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(grecip_dlamp)) Then
          Allocate (grecip_dlamp ( 1-halo_i : halo_i ))
        End If
        If (.not.allocated(grecip_dlamu)) Then
          Allocate (grecip_dlamu ( 1-halo_i : halo_i ))
        End If

      Else  !  Variable grid; set up parameters

        Write ( Unit = 6, fmt=*) ' '
        Write ( Unit = 6, fmt=*) '*** This run uses a variable '        &
     &           , ' horizontal grid L_regular = ', L_regular, ' ***'
        Write ( Unit = 6, fmt=*) ' '

! Allocate arrays for dyn_var_res_mod module
        If (.not.allocated(glambda_p)) Then
          Allocate (glambda_p ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(glambda_u)) Then
          Allocate (glambda_u ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(phi_p)) Then
          Allocate (phi_p ( 1-halo_i : row_length + halo_i,             &
     &                      1-halo_j : rows + halo_j ))
        End If
        If (.not.allocated(phi_v)) Then
          Allocate (phi_v ( 1-halo_i : row_length + halo_i,             &
     &                      1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(gdlambda_p)) Then
          Allocate (gdlambda_p ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(gdlambda_u)) Then
          Allocate (gdlambda_u ( 1-halo_i : global_row_length+halo_i ))
        End If
        If (.not.allocated(dphi_p)) Then
          Allocate (dphi_p ( 1-halo_i : row_length + halo_i,            &
     &                       1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(dphi_v)) Then
          Allocate (dphi_v ( 1-halo_i : row_length + halo_i,            &
     &                       1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(grecip_dlamp)) Then
          Allocate (grecip_dlamp ( 1-halo_i : global_row_length+halo_i))
        End If
        If (.not.allocated(grecip_dlamu)) Then
          Allocate (grecip_dlamu ( 1-halo_i : global_row_length+halo_i))
        End If
        If (.not.allocated(recip_dphip)) Then
          Allocate (recip_dphip ( 1-halo_i : row_length + halo_i,       &
     &                            1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_dphiv)) Then
          Allocate (recip_dphiv ( 1-halo_i : row_length + halo_i,       &
     &                            1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(wt_lambda_p)) Then
          Allocate (wt_lambda_p ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(wt_lambda_u)) Then
          Allocate (wt_lambda_u ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(wt_phi_p)) Then
          Allocate (wt_phi_p ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(wt_phi_v)) Then
          Allocate (wt_phi_v ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(lambda_p_rm)) Then
          Allocate (lambda_p_rm ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(lambda_p_rp)) Then
          Allocate (lambda_p_rp ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(lambda_u_rm)) Then
          Allocate (lambda_u_rm ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(lambda_u_rp)) Then
          Allocate (lambda_u_rp ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(phi_p_rm)) Then
          Allocate (phi_p_rm ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(phi_p_rp)) Then
          Allocate (phi_p_rp ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(phi_v_rm)) Then
          Allocate (phi_v_rm ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(phi_v_rp)) Then
          Allocate (phi_v_rp ( 1-halo_i : row_length + halo_i,          &
     &                         1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_lambda_p_m)) Then
          Allocate (recip_lambda_p_m ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_p_0)) Then
          Allocate (recip_lambda_p_0 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_p_p)) Then
          Allocate (recip_lambda_p_p ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_p_p2)) Then
          Allocate (recip_lambda_p_p2 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_m)) Then
          Allocate (recip_lambda_u_m ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_0)) Then
          Allocate (recip_lambda_u_0 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_p)) Then
          Allocate (recip_lambda_u_p ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_lambda_u_p2)) Then
          Allocate (recip_lambda_u_p2 ( 1-halo_i : row_length+halo_i ))
        End If
        If (.not.allocated(recip_phi_p_m)) Then
          Allocate (recip_phi_p_m ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_p_0)) Then
          Allocate (recip_phi_p_0 ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_p_p)) Then
          Allocate (recip_phi_p_p ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_p_p2)) Then
          Allocate (recip_phi_p_p2 ( 1-halo_i : row_length + halo_i,    &
     &                               1-halo_j : rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_m)) Then
          Allocate (recip_phi_v_m ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_0)) Then
          Allocate (recip_phi_v_0 ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_p)) Then
          Allocate (recip_phi_v_p ( 1-halo_i : row_length + halo_i,     &
     &                              1-halo_j : n_rows+halo_j ))
        End If
        If (.not.allocated(recip_phi_v_p2)) Then
          Allocate (recip_phi_v_p2 ( 1-halo_i : row_length + halo_i,    &
     &                               1-halo_j : n_rows+halo_j ))
        End If

! DEPENDS ON: set_var_grid
      Call Set_var_grid (                                               &
     &                   Lambda_p_in, Lambda_u_in,                      &
     &                   Phi_p_in, Phi_v_in,                            &
     &                   global_row_length, global_rows,                &
     &                   row_length, rows,  n_rows,                     &
     &                   halo_i, halo_j, L_varres,                      &
     &                   delta_lambda_wk, delta_phi_wk,                 &
     &                   base_lambda_wk, base_phi_wk,                   &
     &                   lam_var, phi_var,                              &
     &                   var_ratio, lam_ratio, phi_ratio,               &
     &                   lam_frac, phi_frac, Pi_over_180,               &
     &                   glambda_p, glambda_u, phi_p, phi_v,            &
     &                   gdlambda_p, gdlambda_u, dphi_p, dphi_v,        &
     &                   grecip_dlamp, grecip_dlamu,                    &
     &                   recip_dphip, recip_dphiv,                      &
     &                   wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v,  &
     &                   lambda_p_rm, lambda_p_rp,                      &
     &                   lambda_u_rm, lambda_u_rp,                      &
     &                   phi_p_rm, phi_p_rp, phi_v_rm, phi_v_rp,        &
     &                   lambda_p_end, lambda_u_end,                    &
     &                   phi_p_end, phi_v_end, dlambda_p_end,           &
     &                   dlambda_u_end, dphi_p_end, dphi_v_end,         &
     &                   recip_lambda_p_m, recip_lambda_p_0,            &
     &                   recip_lambda_p_p, recip_lambda_p_p2,           &
     &                   recip_lambda_u_m, recip_lambda_u_0,            &
     &                   recip_lambda_u_p, recip_lambda_u_p2,           &
     &                   recip_phi_p_m, recip_phi_p_0,                  &
     &                   recip_phi_p_p, recip_phi_p_p2,                 &
     &                   recip_phi_v_m, recip_phi_v_0,                  &
     &                   recip_phi_v_p, recip_phi_v_p2,                 &
     &                   max_look, recip_dlam, recip_dphi,              &
     &                   halo_lam, halo_phi, look_lam, look_phi,        &
     &                   model_domain, datastart )

      End If   !  L_regular

! 1. set trig fields and Coriolis components

!!! Allocate latitude arrays for trignometric_mod module
      If (.not.allocated(cos_theta_latitude)) Then
        Allocate (cos_theta_latitude (1-Offx:row_length+Offx,           &
     &                                1-Offy:rows+Offy))
      End If
      If (.not.allocated(sec_theta_latitude)) Then
        Allocate (sec_theta_latitude (1-Offx:row_length+Offx,           &
     &                                1-Offy:rows+Offy))
      End If
      If (.not.allocated(FV_cos_theta_latitude)) Then
        Allocate (FV_cos_theta_latitude (1-Offx:row_length+Offx,        &
     &                                   1-Offy:rows+Offy))
      End If
      If (.not.allocated(FV_sec_theta_latitude)) Then
        Allocate (FV_sec_theta_latitude (1-Offx:row_length+Offx,        &
     &                                   1-Offy:rows+Offy))
      End If
      If (.not.allocated(sin_theta_latitude)) Then
        Allocate (sin_theta_latitude (row_length, rows))
      End If
      If (.not.allocated(tan_theta_latitude)) Then
        Allocate (tan_theta_latitude (row_length, rows))
      End If
      If (.not.allocated(sin_v_latitude)) Then
        Allocate (sin_v_latitude (row_length, n_rows))
      End If
      If (.not.allocated(tan_v_latitude)) Then
        Allocate (tan_v_latitude (row_length, n_rows))
      End If
      If (.not.allocated(cos_v_latitude)) Then
        Allocate (cos_v_latitude (1-Offx:row_length+Offx,               &
     &                            1-Offy:n_rows+Offy))
      End If
      If (.not.allocated(sec_v_latitude)) Then
        Allocate (sec_v_latitude (1-Offx:row_length+Offx,               &
     &                            1-Offy:n_rows+Offy))
      End If

!!! Allocate longitude and 'true' arrays for trignometric_mod module
      If (.not.allocated(cos_theta_longitude)) Then
        Allocate (cos_theta_longitude (row_length, rows))
      End If
      If (.not.allocated(sin_theta_longitude)) Then
        Allocate (sin_theta_longitude (row_length, rows))
      End If
      If (.not.allocated(cos_u_longitude)) Then
        Allocate (cos_u_longitude (row_length, rows))
      End If
      If (.not.allocated(sin_u_longitude)) Then
        Allocate (sin_u_longitude (row_length, rows))
      End If
      If (.not.allocated(true_latitude)) Then
        Allocate (true_latitude (row_length, rows))
      End If
      If (.not.allocated(true_longitude)) Then
        Allocate (true_longitude (row_length, rows))
      End If

      If(model_domain  /=  mt_global) Then
!  For LAMs only, Trigs may be set to equatorial values
        If (L_trivial_trigs) then
        If(mype  ==  0) write(6,*)'WARNING: trig values set to trivial'
        Do j = 1, rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            cos_theta_latitude(i, j) = 1.
            sec_theta_latitude(i, j) = 1.
        sin_theta_latitude(i, j) = 0.0
!       sin_theta_latitude(i, j) = sin(f_plane_rad)
            tan_theta_latitude(i, j) = 0.
            FV_cos_theta_latitude(i, j) = 1.
            FV_sec_theta_latitude(i, j) = 1.
          End Do
        End Do

        Do j = 1, n_rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
        sin_v_latitude(i, j) = 0.0
!       sin_v_latitude(i, j) = sin(f_plane_rad)
            cos_v_latitude(i, j) = 1.
            tan_v_latitude(i, j) = 0.
            sec_v_latitude(i, j) = 1.
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            sin_theta_longitude(i, j) = sin ((gi-1)*delta_lambda)
            cos_theta_longitude(i, j) = cos ((gi-1)*delta_lambda)
            sin_u_longitude(i, j) = sin ((gi-.5)*delta_lambda)
            cos_u_longitude(i, j) = cos ((gi-.5)*delta_lambda)
          End Do
        End Do

        End If        ! L_trivial_trigs = .true.
      End If        !model_domain  /=  mt_global
      If (model_domain  ==  mt_global .or. .not.L_trivial_trigs) Then

! non-trivial trigonometric info.
        Do j = 1, rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            cos_theta_latitude(i, j) = cos(Base_phi+(gj-1)*delta_phi)
            sin_theta_latitude(i, j) = sin(Base_phi+(gj-1)*delta_phi)
            FV_cos_theta_latitude(i, j) = cos_theta_latitude(i, j)
          End Do
        End Do

        If (model_domain  ==  mt_global) Then

          j0 = 1
          j1 = rows
          If (at_extremity(PNorth)) Then
            Do i = 1, row_length
              cos_theta_latitude(i, rows) = 0.
              sin_theta_latitude(i, rows) = 1.
              FV_cos_theta_latitude(i, rows) = delta_phi/8.
              sec_theta_latitude(i, rows) = rmdi
              tan_theta_latitude(i, rows) = rmdi
            End Do

            j1 = rows - 1

          End If

          If (at_extremity(PSouth)) Then
            Do i = 1, row_length
              cos_theta_latitude(i, 1) = 0.
              sin_theta_latitude(i, 1) = -1.
              FV_cos_theta_latitude(i, 1) = delta_phi/8.
              sec_theta_latitude(i, 1) = rmdi
              tan_theta_latitude(i, 1) = rmdi
            End Do

            j0 = 2

          End If


          Do j = j0, j1
            Do i = 1, row_length
              sec_theta_latitude(i, j) = 1./cos_theta_latitude(i,j)
              tan_theta_latitude(i, j) = sin_theta_latitude(i,j)        &
     &                                   /cos_theta_latitude(i,j)
            End Do
          End Do


        Else
! Limited area model
          Do j = 1, rows
            Do i = 1, row_length
              sec_theta_latitude(i, j) = 1./cos_theta_latitude(i, j)
              tan_theta_latitude(i, j) = sin_theta_latitude(i, j)       &
     &                                   /cos_theta_latitude(i, j)
            End Do
          End Do
        End If ! model_domain

        Do j = 1, rows
          Do i = 1, row_length
            FV_sec_theta_latitude(i, j) = 1./FV_cos_theta_latitude(i,j)
          End Do
        End Do

        Do j = 1, n_rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            sin_v_latitude(i, j) = sin(Base_phi + (gj-.5)*delta_phi)
            cos_v_latitude(i, j) = cos(Base_phi + (gj-.5)*delta_phi)
            sec_v_latitude(i, j) = 1./cos_v_latitude(i, j)
            tan_v_latitude(i, j) = sin_v_latitude(i, j)                 &
     &                             /cos_v_latitude(i, j)
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            sin_theta_longitude(i, j) =                                 &
     &                      sin (Base_lambda+(gi-1)*delta_lambda)
            cos_theta_longitude(i, j) =                                 &
     &                      cos (Base_lambda+(gi-1)*delta_lambda)
            sin_u_longitude(i, j) =                                     &
     &                      sin (Base_lambda+(gi-.5)*delta_lambda)
            cos_u_longitude(i, j) =                                     &
     &                      cos (Base_lambda+(gi-.5)*delta_lambda)
          End Do
        End Do

      End If !model_domain  ==  mt_global .or. .not.L_trivial_trigs

!!! Allocate arrays for dyn_coriolis_mod module
      If (.not.allocated(f1_at_v)) Then
        Allocate (f1_at_v (row_length, 0:n_rows+1))
      End If
      If (.not.allocated(f2_at_u)) Then
        Allocate (f2_at_u (0:row_length+1, rows))
      End If
      If (.not.allocated(f3_at_u)) Then
      Allocate (f3_at_u (1-offx:row_length+offx, 1-offy:rows+offy))
      End If
      If (.not.allocated(f3_at_v)) Then
      Allocate (f3_at_v (1-offx:row_length+offx, 1-offy:n_rows+offy))
      End If

      If (model_domain  ==  mt_global) Then
! global model
        Do j = 1, n_rows
          Do i = 1, row_length
            f1_at_v(i,j) = 0.
            f3_at_v(i,j) = two_omega * sin_v_latitude(i,j)
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            f2_at_u(i,j) = two_omega * cos_theta_latitude(i,j)
            f3_at_u(i,j) = two_omega * sin_theta_latitude(i,j)
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            true_longitude(i,j) = (gi-1)*delta_lambda
          End Do
        End Do
! set polar values of true longitude to be all the same

        If (at_extremity(PNorth)) then
          Do i = 1, row_length
            true_longitude(i,rows) = 0.
          End Do
        End If

        If (at_extremity(PSouth)) then
          Do i = 1, row_length
            true_longitude(i,1) = 0.
          End Do
        End If

! get parameters for common latlonmax within mppac
! these need to be in degrees
#if defined(MPP)
      lat_n = Base_phi_wk + (datastart(2)+rows-2) * Delta_phi_wk
      lat_s = Base_phi_wk

      long_w       = Base_lambda_wk
      long_w_model = long_w
      long_e       = Base_lambda_wk +                                   &
     &              (datastart(1)+row_length-2) * Delta_lambda_wk
      long_e_model = long_e

      if(long_w >  180.0) long_w = long_w-360.0
      if(long_e >  180.0) long_e = long_e-360.0

#endif
      Else
! limited area model
        temp1 = cos( lat_rot_NP)
        temp2 = sin( lat_rot_NP)

      If( L_trivial_trigs .or. f_plane  >   -89.0) Then
        f_plane_rad = f_plane * Pi_over_180
        ff_plane_rad = ff_plane * Pi_over_180
        f1_temp = - two_omega * temp1 * sin(ff_plane_rad)
        f2_temp = two_omega * (cos(f_plane_rad) * temp2                 &
     &               - sin(f_plane_rad) * temp1 * cos(ff_plane_rad))
        f3_temp = two_omega * ( sin(f_plane_rad) * temp2                &
     &               + cos(f_plane_rad) * temp1 * cos(ff_plane_rad))
        If (L_vert_Coriolis) Then
          f1_temp = 0.0
          f2_temp = 0.0
          f3_temp = two_omega * ( sin(f_plane_rad) * temp2              &
     &               + cos(f_plane_rad) * temp1 * cos(ff_plane_rad))
        End If  !  L_vert_Coriolis

        Do j = 1, n_rows
          Do i = 1, row_length
            f1_at_v(i,j) = f1_temp
            f3_at_v(i,j) = f3_temp
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            f2_at_u(i,j) = f2_temp
            f3_at_u(i,j) = f3_temp
          End Do
        End Do

      Else        !  L_trivial_trigs = .false.
        Do j = 1, n_rows
          Do i = 1, row_length
            f1_at_v(i,j) = - two_omega * temp1 *                        &
     &                       sin_theta_longitude(i,j)
            f3_at_v(i,j) = two_omega * ( sin_v_latitude(i,j) * temp2    &
     &                                  +cos_v_latitude(i,j) * temp1    &
     &                                  *cos_theta_longitude(i,j) )
          End Do
        End Do
        Do j = 1, rows
          Do i = 1, row_length
            f2_at_u(i,j) = two_omega * ( cos_theta_latitude(i,j) * temp2&
     &                                  -sin_theta_latitude(i,j) * temp1&
     &                                  *cos_u_longitude(i,j) )
            f3_at_u(i,j) = two_omega * ( sin_theta_latitude(i,j) * temp2&
     &                                  +cos_theta_latitude(i,j) * temp1&
     &                                  *cos_u_longitude(i,j) )
          End Do
        End Do
      End If   !  L_trivial_trigs

! calculate true longitude in radians
        Do j = 1, rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            longitude_rot(i,j) = (base_lambda + (gi-1) * delta_lambda)  &
     &                         / Pi_over_180
            latitude_rot(i,j) = (base_phi + (gj-1) * delta_phi)         &
     &                         / Pi_over_180
          End Do
        End Do

! DEPENDS ON: eqtoll
        Call eqtoll(latitude_rot, longitude_rot                         &
     &,             true_latitude, true_longitude                       &
     &,             lat_rot_NP_in, long_rot_NP_in, rows*row_length)

        Do j = 1, rows
          Do i = 1, row_length
            true_longitude(i,j) = true_longitude(i,j) * Pi_over_180
            true_latitude(i,j) = true_latitude(i,j) * Pi_over_180
          End Do
        End Do

! get parameters for common latlonmax within mppac
! these need to be in degrees
#if defined(MPP)
      lat_n = latitude_rot(1,rows)
      lat_s = latitude_rot(1,1)

      long_w       = longitude_rot(1,1)
      long_w_model = long_w
      long_e       = longitude_rot(row_length,1)
      long_e_model = long_e

      if(long_w >  180.0) long_w = long_w-360.0
      if(long_e >  180.0) long_e = long_e-360.0
      if(long_w_model >= 360.0)long_w_model = long_w_model-360.
      if(long_e_model >= 360.0)long_e_model = long_e_model-360.

#endif
! calculate lat/longitude for points on equatorial grid for B grid
        Do j = 1, n_rows
          gj = datastart(2) + j - 1
          Do i = 1, row_length
            gi = datastart(1) + i - 1
            longitude_rot(i,j) = (base_lambda + (gi-.5) * delta_lambda) &
     &                         / Pi_over_180
            latitude_rot(i,j) = (base_phi + (gj-.5) * delta_phi)        &
     &                         / Pi_over_180
          End Do
        End Do

! DEPENDS ON: eqtoll
        Call eqtoll(latitude_rot, longitude_rot                         &
     &,             true_latitude_b, true_longitude_b                   &
     &,             lat_rot_NP_in, long_rot_NP_in                       &
     &,             n_rows*row_length)

! Calculate rotation coefficients for wind

!!! Allocate arrays for rot_coeff_mod module
      If (.not.allocated(rot_coeff1)) Then
        Allocate (rot_coeff1 ( row_length, n_rows ))
      End If
      If (.not.allocated(rot_coeff2)) Then
        Allocate (rot_coeff2 ( row_length, n_rows ))
      End If

! DEPENDS ON: w_coeff
        Call w_coeff(rot_coeff1, rot_coeff2,                            &
     &               true_longitude_b, longitude_rot,                   &
     &              lat_rot_NP_in, long_rot_NP_in,                      &
     &               n_rows*row_length )

      End If ! model_domain

! Call swap_bounds for those trig fields which require halos

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (FV_cos_theta_latitude, row_length, rows, 1,    &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(FV_cos_theta_latitude,row_length,rows,   &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (FV_sec_theta_latitude, row_length, rows, 1,    &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(FV_sec_theta_latitude,row_length,rows,   &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (cos_theta_latitude, row_length, rows, 1,       &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(cos_theta_latitude,row_length,rows,      &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (sec_theta_latitude, row_length, rows, 1,       &
     &                   Offx, Offy, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(sec_theta_latitude,row_length,rows,      &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (cos_v_latitude, row_length, n_rows, 1,         &
     &                   Offx, Offy, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(cos_v_latitude,row_length,n_rows,        &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (sec_v_latitude, row_length, n_rows, 1,         &
     &                   Offx, Offy, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(sec_v_latitude,row_length,n_rows,        &
     &                         1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f2_at_u, row_length, rows, 1,                  &
     &                   1, 0, fld_type_u, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f2_at_u,row_length,rows,1,1,0)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f1_at_v, row_length, n_rows, 1,                &
     &                   0, 1, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f1_at_v,row_length,n_rows,1,0,1)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f3_at_u, row_length, rows, 1,                  &
     &                   Offx, Offy, fld_type_u, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f3_at_u,row_length,rows,1,Offx,Offy)

! DEPENDS ON: swap_bounds
      call Swap_Bounds  (f3_at_v, row_length, n_rows, 1,                &
     &                   Offx, Offy, fld_type_v, .false.)
! DEPENDS ON: fill_external_halos
      CALL FILL_EXTERNAL_HALOS(f3_at_v,row_length,n_rows,1,Offx,Offy)

! ----------------------------------------------------------------------
! 1. Set polar filtering and diffusion
! ----------------------------------------------------------------------
!!! Allocate arrays for diff_coeff_mod module
      If (.not.allocated(diff_coeff_u)) Then
        Allocate (diff_coeff_u ( 1-Offx : row_length+Offx,              &
     &                           1-Offy : rows+Offy ))
      End If
      If (.not.allocated(diff_coeff_v)) Then
        Allocate (diff_coeff_v ( 1-Offx : row_length+Offx,              &
     &                           1-Offy : n_rows+Offy ))
      End If

      If ( L_pofil_new ) Then

! DEPENDS ON: setdiff
        Call Setdiff(                                                   &
! Input size and control variables
     &      global_rows, row_length, rows, n_rows, model_levels,        &
     &      model_domain, at_extremity, datastart,                      &
     &      offx, offy, mype, nproc_y,                                  &
     &      max_model_levels, max_121_rows, max_sweeps,                 &
! other constants
     &      delta_lambda, delta_phi,                                    &
     &      pi, Pi_over_180, polar_cap, scale_ratio,                    &
     &      ref_lat_deg, diff_coeff_ref,                                &
     &      cos_theta_latitude, sin_theta_latitude,                     &
     &      cos_v_latitude, sin_v_latitude,                             &
! Output data
     &       global_u_filter, global_v_filter,                          &
     &       u_sweeps, v_sweeps,                                        &
     &       u_begin, u_end, v_begin, v_end,                            &
     &       diff_coeff_u, diff_coeff_v, diff_coeff_phi,                &
     &       diffusion_coefficient_thermo, diffusion_order_thermo,      &
     &       diffusion_coefficient_wind, diffusion_order_wind,          &
     &       diff_order_thermo, diff_order_wind,                        &
     &       diff_timescale_thermo, diff_timescale_wind,                &
     &       L_sponge, sponge_power, sponge_ew, sponge_ns,              &
     &       sponge_wts_ew, sponge_wts_ns,                              &
     &       L_diffusion, L_cdiffusion, L_filter, L_diff_auto,          &
     &       L_pfcomb, L_pfexner, L_pftheta, L_pfuv, L_pfw, L_pfincs,   &
     &       L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs)

      Else ! original combi/diffusion setting

! DEPENDS ON: setdiff_old
        Call Setdiff_old(                                               &
! Input size and control variables
     &      global_rows, row_length, rows, n_rows,                      &
     &      model_domain, at_extremity, datastart,                      &
     &      offx, offy, mype, nproc_y, max_121_rows,                    &
! other constants
     &      delta_lambda, delta_phi,                                    &
     &      pi, Pi_over_180, polar_cap, scale_ratio, diff_coeff_ref,    &
     &      cos_theta_latitude, sin_theta_latitude,                     &
     &      cos_v_latitude, sin_v_latitude,                             &
! Output data
     &       global_u_filter, global_v_filter,                          &
     &       u_sweeps, v_sweeps,                                        &
     &       u_begin, u_end, v_begin, v_end,                            &
     &       diff_coeff_u, diff_coeff_v,                                &
     &       diff_coeff_thermo, diff_coeff_wind,                        &
     &       diff_order_thermo, diff_order_wind,                        &
     &       diff_timescale_thermo, diff_timescale_wind,                &
     &       L_filter, L_pfcomb, L_pftheta, L_pfuv, L_pfw, L_pfincs,    &
     &       L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs)

      End If ! L_pofil_new

! ----------------------------------------------------------------------
! 1.5 Set filtering control variable
! ----------------------------------------------------------------------

      If ( L_polar_filter .or. L_pfcomb ) Then
        L_filter = .true.
        if ( model_domain  /=  mt_global) then
          write(6,*)' ******   Polar filter NOT ALLOWED in LAMs ****** '
          write(6,*)' ******       RESET UMUI or NAMELIST       ****** '
          ICODE    = 1
          CMESSAGE ='SETCONA:Polar filter NOT ALLOWED in LAMs'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)
        else if ( L_pofil_new ) then
          write(6,*)' ******  Newest filtering code being used ****** '
        end if ! model_domain  /=  mt_global
      Else
        write(6,*) ' ******   Polar filter is OFF   ****** '
      End If ! L_polar_filter .or. L_pfcomb
      If ( L_diff_thermo .or. L_diff_wind .or. L_diff_w ) Then
        if ( L_pofil_new ) then
          write(6,*) ' ******  Newest diffusion code being used ****** '
        endif !  L_pofil_new
        L_filter = .true.
      End If ! L_diff_thermo .or. L_diff_wind .or. L_diff_w
      If ( L_diff_incs .or. L_pfincs) Then
        L_filter_incs = .true.
      End If ! L_diff_incs .or. L_pfincs
      If ( L_pfexner  ) Then
        write(6,*) ' ** Polar filtering of Exner pressure is active ** '
      End If ! L_pfexner
      If ( L_diff_exner  ) Then
        write(6,*) ' ** Diffusion of Exner pressure is active ** '
      End If ! L_diff_exner

      If( L_tardiff_q ) Then
        write(6,*) ' '
        write(6,*) '  ***   Targeted diffusion of q is ACTIVE  *** '
        write(6, 913) w_conv_limit, tardiffq_factor
        write(6, 914) tardiffq_test
        write(6, 915) tardiffq_start, tardiffq_end
        write(6, 916) tar_horizontal
      Else
        write(6,*) ' '
        write(6,*) '  ***   Targeted diffusion of q is NOT ACTIVE  *** '
      End If !   L_tardiff_q

!L ----------------------------------------------------------
!  Print information to tell user what is active
!L ----------------------------------------------------------
!  vertical diffusion coefficients
      If ( L_vdiff_uv ) Then
        vdiffuv_factor = 0.25 *                                         &
     &                     (1.0 - EXP( -1.0 / vdiffuv_timescale) )
       Write ( Unit = 6, fmt=*) 'A factor delta_z**2/timestep is '      &
     & ,' allowed for in the application of the vertical diffusion.'
       Write ( Unit = 6, fmt=*) 'Vertical diffusion of order 1 with'    &
     &        , ' vertical diffusion coefficient = ', vdiffuv_factor
       Write ( Unit = 6, fmt=*) 'Vertical diffusion timescale is '      &
     &                           , vdiffuv_timescale,' steps '
      End If ! L_vdiff_uv

      if ( top_filt_start > model_levels ) then
        Write ( Unit = 6, fmt=*) 'No additional upper-level diffusion'
      else  ! top_filt_start <= model_levels
        active_levels = top_filt_end - top_filt_start + 1
        if ( active_levels > max_updiff_levels) then
          top_filt_start  = top_filt_end - max_updiff_levels + 1
          write(6,*) ' Max of uppermost ', max_updiff_levels            &
     &          , ' levels allowed for upper-level'                     &
     &          , ' diffusion - start level reset to ', top_filt_start
          active_levels = max_updiff_levels
        end if ! active_levels > max_updiff_levels
        Write ( Unit = 6, fmt=*) 'Extra upper-level diffusion applied'
        scale = 1.0
        if ( L_upper_ramp ) then
          scale = up_diff_scale
          write(6, 911) top_diff, top_filt_end
          write(6, 912) scale, top_filt_start
        else  !
          write(6, 910) top_diff, top_filt_start, top_filt_end
        end if  !  L_upper_ramp
        kk = active_levels
        up_diff(active_levels) = top_diff
        do k = top_filt_end - 1, top_filt_start, - 1
          kk = kk - 1
          up_diff(kk) = scale * up_diff(kk + 1)
        enddo !  k = top_filt_end - 1, top_filt_start, - 1
        Do k = 1, active_levels
          kk = k + top_filt_start - 1
          write(6,*)'Level', kk,' Diffusion factor = ',up_diff(k)
        End do !  k = 1, active_levels
      end if ! top_filt_start > model_levels

      if (L_adjust_theta) then
        if( adjust_theta_start < 2 ) then
          adjust_theta_start = 2
          Write ( Unit = 6, fmt=*) 'Start level for convective '        &
     & , ' adjustment reset to ', adjust_theta_start
        end if ! adjust_theta_start < 2
        Write ( Unit = 6, fmt=*) 'Convective adjustment applied when'   &
     &    ,' lapse rate is negative above level ',  adjust_theta_start  &
     &    ,' up to level ', adjust_theta_end
      end if  ! L_adjust_theta

! ----------------------------------------------------------------------
! Section 1.5a . Print FORMATTING
! ----------------------------------------------------------------------

 910  FORMAT(' Diffusion coefficient = ',F5.2                           &
     &      ,' from level ', I3,' to level ', I3)
 911  FORMAT(' Diffusion coefficient = ',F5.2,' at level ', I3)
 912  FORMAT(' decreasing by ', F6.3, ' at each level down to ', I3)
 913  FORMAT('   Vertical velocity test value is  ',F5.2,' m/s and ',   &
     &                                 'diffusion factor = ',F6.3)
 914  FORMAT('   Start testing at level  ',I3)
 915  FORMAT('   Apply from level  ',I3,' to level ',I3)
 916  FORMAT('   Slope test applied up to level ',I4)


! ----------------------------------------------------------------------
! 2. Semi-Lagrangian scheme options.
! ----------------------------------------------------------------------


! j holds number of points required for interpolation scheme
      j = 2
      Do i = 1, 3
        If (high_order_scheme(i)  ==  quinticLagrange) Then
          j = 3
          If (halo_j  <   5) Then
            write(6,*)' error halo_j too small ',halo_j
          End If
        Else
          If (halo_j  <   4) Then
            write(6,*)' error halo_j too small ',halo_j
          End If
        End If
      End Do

! Calculate max cfl possible in LAM
! Y-direction also applies in global model though variable is not
! used in code.

      LAM_max_cfl(1) = halo_i - j
      LAM_max_cfl(2) = halo_j - j

! ----------------------------------------------------------------------
! Set up data.
! ----------------------------------------------------------------------

! Calculate land_index
          land_points = 0
          Do j =1, rows
            Do i = 1, row_length
              If (land_sea_mask(i,j)) then
                land_points = land_points + 1
                land_index(land_points) = (j-1)*row_length + i
              End If
            End Do
          End Do

! set-up code for hydrology
          soil_points = 0
          land_ice_points = 0
          Do i= 1, land_points
! Test on soil moisture concentration at saturation
            If (smvcst(i) >   0.0) Then       ! Soil points
              soil_points = soil_points+1
              soil_index(soil_points)=i
            Else If (smvcst(i)  ==  0.0) Then   ! Land-ice points
              land_ice_points = land_ice_points+1
              land_ice_index(land_ice_points)=i
            End If
          End Do


! calculate r_at_u points on rho levels and
! r_at_v points on rho levels.

!!! Allocate r_at_u/v arrays for level_heights_mod module
      If (.not.allocated(r_at_u)) Then
        Allocate (r_at_u (1-halo_i:row_length+halo_i,                   &
     &                    1-halo_j:rows+halo_j, model_levels))
      End If
      If (.not.allocated(r_at_v)) Then
        Allocate (r_at_v (1-halo_i:row_length+halo_i,                   &
     &                    1-halo_j:n_rows+halo_j, model_levels))
      End If

        Do k = 1, model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i - 1
              r_at_u (i,j,k) = .5 * (r_rho_levels(i,j,k) +              &
     &                                 r_rho_levels(i+1,j,k) )
            End Do
          End Do
        End Do

        Do k = 1, model_levels
          Do j = 1-halo_j, n_rows+halo_j - 1
            Do i = 1-halo_i, row_length+halo_i
              r_at_v (i,j,k) = .5 * (r_rho_levels(i,j,k) +              &
     &                               r_rho_levels(i,j+1,k) )
            End Do
          End Do
        End Do

! call swap_bounds to set extra points
! DEPENDS ON: swap_bounds
       call Swap_Bounds(r_at_u, row_length, rows, model_levels,         &
     &                  halo_i, halo_j, fld_type_u, .false.)

! DEPENDS ON: swap_bounds
       call Swap_Bounds(r_at_v, row_length, n_rows, model_levels,       &
     &                  halo_i, halo_j, fld_type_v, .false.)

      If ( L_diag_L2norms ) then
        Write ( Unit = 6, fmt=*) 'Printing of norms during timestep'    &
     &, ' activated'
      endIf ! L_diag_L2norms
      If ( L_diag_L2helm ) then
        Write ( Unit = 6, fmt=*) 'Printing of coefficient norms  '      &
     &, ' in solver'
      endIf ! L_diag_L2helm
!  Initialise diagnostic printing items
! Sampling interval must not be larger than print interval
      if ( diag_interval > print_step) then
        diag_interval = print_step
      end if !  diag_interval > print_step
! Set array sizes needed for chosen diagnostic printing
! These are needed to do sums and max/mins across processors
      rpemax = 0
      rpemin = 0
      ipesum = 0
      rpesum = 0
      if (L_print_theta1 ) then
        rpemin = rpemin + 1
        ipesum = ipesum + 1
        rpesum = rpesum + 3
      end if !  L_print_theta1
      if (L_print_lapse ) then
        rpemin = rpemin + (model_levels - 1)
        ipesum = ipesum + 2 * (model_levels - 1)
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_lapse
      if (L_print_div) then
        rpemax = rpemax + model_levels
        rpemin = rpemin + model_levels
        ipesum = ipesum + 2 * model_levels
        rpesum = rpesum + 4 * model_levels
      end if !  L_print_div
      if (L_print_w .or. L_print_wmax) then
        rpemax = rpemax + model_levels - 1
        ipesum = ipesum + 2 * (model_levels - 1)
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_w .or. L_print_wmax
      if (L_print_shear) then
        rpemax = rpemax + model_levels - 1
        ipesum = ipesum + model_levels - 1
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_shear
      if (L_print_max_wind) then
        rpemax = rpemax + model_levels
        ipesum = ipesum + model_levels
        rpesum = rpesum + 2 * model_levels
      end if !  L_print_max_wind
      if ( L_diag_wind ) then
        rpesum = rpesum + 2 * model_levels
      end if  ! L_diag_wind
        min_theta1_run = 1000.0
        time_theta1_min = 0
      Do k = 1, model_levels
         max_w_run(k) =  0.0
         time_w_max(k) = 0
         max_div_run(k) =  0.0
         time_div_max(k) = 0
         min_div_run(k) =  0.0
         time_div_min(k) = 0
        min_lapse_run(k) =  1.0e6
         time_lapse_min(k) = 0
        max_shear_run(k) = 0.0
        time_max_shear(k) = 0
        max_wind_run(k) = 0.0
        time_max_wind(k) = 0
        max_KE_run(k) = 0.0
        min_KE_run(k) = 1.0e30
        time_KE_max(k) = 0
        time_KE_min(k) = 0
        time_noise_max(k) = 0
        max_noise_run(k) = 0.0
      End Do  ! k = 1, model_levels
      max_KE_run(model_levels + 1) = 0.0
      min_KE_run(model_levels + 1) = 1.0e30
      time_KE_max(model_levels + 1) = 0
      time_KE_min(model_levels + 1) = 0

!L ----------------------------------------------------------
!L 3. Set up cloud type boundaries for low/medium/high cloud.
!L ----------------------------------------------------------
       IF (NUM_CLOUD_TYPES  >   3) THEN
         ICODE    = 1
         CMESSAGE = 'SETCONA: Parameter NUM_CLOUD_TYPES exceeds 3'
         WRITE(6,*) 'NUM_CLOUD_TYPES=',NUM_CLOUD_TYPES
! DEPENDS ON: ereport
         CALL Ereport(RoutineName,ICODE,Cmessage)

       END IF

!  Diagnostics of cloud amount for low/medium/high cloud are calculated
!  by finding the max cloud cover over a range of model levels. The
!  ranges are determined by heights [ 1->2:2->3:3->4 = L:M:H ] held in
!  array h_split by comparison with heights of model levels (above
!  surface) in r_ref_theta.

      DO KK = 1, NUM_CLOUD_TYPES + 1
        LEVEL = 1
!
!       r_ref_theta turns out to be A_theta(k) - Earth_radius at every
!       model level because z_top_of_model = z_rho_top is chosen here.
!
        DO WHILE ((r_ref_theta(LEVEL)  <=  h_split(KK)) .AND.           &
     &                                       (LEVEL  <=  model_levels))
          LEVEL = LEVEL + 1

        END DO

        IF (LEVEL  >   model_levels) THEN
          ICODE    = 1
          CMESSAGE ='SETCONA:Error in locating levels for cloud layers'
! DEPENDS ON: ereport
          CALL Ereport(RoutineName,ICODE,Cmessage)

        END IF
        CLOUD_BOUND(KK) = LEVEL

      END DO

      LOW_BOT_LEVEL  = CLOUD_BOUND(1)
      LOW_TOP_LEVEL  = CLOUD_BOUND(2) - 1
      MED_BOT_LEVEL  = CLOUD_BOUND(2)
      MED_TOP_LEVEL  = CLOUD_BOUND(3) - 1
      HIGH_BOT_LEVEL = CLOUD_BOUND(3)
      HIGH_TOP_LEVEL = CLOUD_BOUND(4) - 1

      IF (LOW_TOP_LEVEL  >   CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA: No of cloud levels less than Top of Low'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      IF (MED_TOP_LEVEL >  CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA:  No of cloud levels less than Top of Med'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      IF (HIGH_TOP_LEVEL  >   CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA: No of cloud levels less than Top of High'
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

!L ----------------------------------------------------------
!L 4. Set up locally-derived parameters.
!L ----------------------------------------------------------

!     The tropopause diagnosed for radiative purposes
!     divides theta-levels considered to lie in the stratosphere
!     from those considered to lie in the troposphere: the
!     tropopause is therefore taken as a rho-level. This level
!     is constrained to lie between heights of z_min_trop
!     and z_max_trop. The level is used in the setting of the
!     background aerosol climatology and of ozone profiles,
!     subject to the choice of appropriate options; additionally
!     it is used in the calculation of diagnostics defined at
!     the tropopause.
!
!     Start at the second rho-level because the first is omitted
!     from the grid seen in the physics.
      min_trop_level=2
      Do ; If ( (r_ref_rho(min_trop_level) >= z_min_trop) .OR.          &
     &          (min_trop_level == model_levels) ) Exit
        min_trop_level = min_trop_level+1
      End Do
!
      max_trop_level=min_trop_level
      Do ; If ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.           &
     &          (max_trop_level == model_levels) ) Exit
        max_trop_level = max_trop_level+1
      End Do
      max_trop_level = max_trop_level-1
!
! set up super tracer array size

          super_array_size=0

          if (l_CO2_interactive) super_array_size=super_array_size+1
          If (l_Soot) super_array_size=super_array_size+3
          If (l_Biomass) super_array_size=super_array_size+3
          If (l_Sulpc_so2)super_array_size=super_array_size+4
          if (L_sulpc_nh3)super_array_size=super_array_size+1
          if (L_sulpc_dms)super_array_size=super_array_size+1
          IF (L_DUST) super_array_size=super_array_size+6
          if (L_ocff) super_array_size=super_array_size+3
          IF (L_Murk_advect)super_array_size=super_array_size+1
          if (L_USE_CARIOLLE)super_array_size=super_array_size+1
          if (tr_levels==model_levels.and.tr_vars>0)                    &
     &            super_array_size=super_array_size+tr_vars
          if (tr_levels==model_levels.and.tr_ukca>0)                    &
     &            super_array_size=super_array_size+tr_ukca


          if(mype == 0)write(6,*)'super_array_size =  ',                &
     &                          super_array_size

! set up moisture array size

!        always do q,qcl,qcf
          moisture_array_size=3
          if (l_pc2) moisture_array_size=moisture_array_size+4
          if (l_mcr_qrain)moisture_array_size=moisture_array_size+1
          if (l_mcr_qcf2 )moisture_array_size=moisture_array_size+1
          if (l_mcr_qgraup)moisture_array_size=moisture_array_size+1
          if(mype == 0)write(6,*)'moisture_array_size =  ',             &
     &                          moisture_array_size

!      end of set up code

! ----------------------------------------------------------------
! 5. Calculate the coeffs for the interpolation of the radiation
!    quantities if spatial degradation of the radiation code is on.
!    Calculate the mask which ensures a chequer-board pattern over
!    the whole domain (and not just for a single PE).
! -----------------------------------------------------------------

! Create a land mask with a halo.
      IF (L_rad_deg) THEN

        landmask(1-offx:0, 1-offy: rows+offy) = .FALSE.
        landmask(row_length+1:row_length+offx, 1-offy:rows+offy)=.FALSE.
        landmask(1:row_length, 1-offy: 0) = .FALSE.
        landmask(1:row_length, rows+1: rows+offy) = .FALSE.

        DO j=1,ROWS
          DO i=1,ROW_LENGTH
            landmask(i,j)=Land_sea_mask(i,j)
          END DO
        END DO

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(landmask,row_length,rows,1,                    &
     &                   offx,offy,fld_type_p,.FALSE.)

        first_row=1
        last_row=rows
        IF (model_domain  ==  mt_global) THEN
          IF (at_extremity(PNorth)) THEN
            last_row=rows-1
          END IF
          IF (at_extremity(PSouth)) THEN
            first_row=2
          END IF
        END IF

!!! Allocate arrays for rad_mask_trop_mod module
      If (.not.allocated(es_space_interp)) Then
        Allocate (es_space_interp(4, row_length, rows))
      End If
      If (.not.allocated(rad_mask)) Then
        Allocate (rad_mask(row_length, rows))
      End If

! DEPENDS ON: coeffs_degrade
        CALL COEFFS_DEGRADE(ES_SPACE_INTERP, landmask,                  &
     &                      first_row, last_row, model_domain,          &
     &                      mt_lam, at_extremity(PNorth),               &
     &                      at_extremity(PSouth), at_extremity(PWest),  &
     &                      at_extremity(PEast),                        &
     &                      row_length*rows,row_length, rows,offx, offy)
! DEPENDS ON: rad_degrade_mask
        CALL RAD_DEGRADE_MASK(RAD_MASK, datastart,                      &
     &                      Ndim_max, first_row, last_row, offx,        &
     &                      row_length, rows)

      ENDIF  !  If L_rad_deg loop

        If ( thmono_height >= 0.5 * r_ref_theta(1) ) then
          If ( thmono_height < r_ref_theta(1) ) then
            k = 1
          else if ( thmono_height > r_ref_theta(model_levels) ) then
            k = model_levels
          else
            k = 1
            do  !  cycle to find nearest level to  thmono_height
              k = k + 1
              if ( thmono_height < r_ref_theta(k) ) then
                if ( r_ref_theta(k) - thmono_height >                   &
     &               thmono_height - r_ref_theta(k-1) ) k = k - 1
                EXIT
              end if ! thmono_height < r_ref_theta(k)
              CYCLE
            end do !  cycle to find nearest level to  thmono_height
          End If ! thmono_height < r_ref_theta(1)
          thmono_levels = k
          Write ( Unit = 6, fmt=*) '3D Monotone limiter will be applied'&
     &,                            ' to advection of theta'
          Write ( Unit = 6, fmt=*) 'thmono_height was set to '          &
     &,                             thmono_height,' in the UMUI'
          Write ( Unit = 6, fmt=*) 'Limiter will be applied up to level'&
     &,   thmono_levels,', the nearest level to',thmono_height,'metres'
        Else   ! thmono_height < 0.5 * r_ref_theta(1)
          thmono_levels = 0
        End If ! thmono_height >= 0.5 * r_ref_theta(1)
!L ----------------------------------------------------------
!L 5. Set up dynamical core parameters
!L ----------------------------------------------------------
!    Set frictional_timescale = 0 everywhere since this value
!    is passed into PE_HELMHOLTZ from ATMSTEP
!    The value used in simple friction is calculated internally
!           ( stored in friction_level)
!     as is SuHe_level_weight (temp1) in the temperature relaxation

      Do k = 1, model_levels
        frictional_timescale(k) = 0.0
      End Do

      If( problem_number  ==  dynamical_core) Then

! DEPENDS ON: idl_set_suhe_params
        Call IDL_Set_SuHe_params                                        &
     &                      (row_length, rows, model_levels             &
     &,                     kappa, Cp, p_zero, Earth_radius, pi, R, g   &
     &,                     offx, offy, halo_i, halo_j                  &
     &,                     mype, nproc, at_extremity                   &
     &,                     datastart, gc_all_proc_group                &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     SuHe_pole_equ_deltaT, SuHe_static_stab      &
     &,                     p, p_theta_levels                           &
     &,                     base_frictional_timescale                   &
     &,                     friction_level, SuHe_sigma_cutoff           &
     &,                     SuHe_level_weight )

      else if( problem_number  ==  idealised_problem) Then
!     Clear arrays
        Do k = 1, model_levels
          SuHe_level_weight(k) = 0.0
          friction_level(k) = 0.0
        End Do    ! k = 1, model_levels

      end if    ! problem_number  ==  dynamical_core
! ----------------------------------------------------------------------
! 6. GCR diagnostics initialisation.
! ----------------------------------------------------------------------
      if (GCR_diagnostics == 3 )then
        L_error = .false.
        Do i = 1, 3
          if(GCR_its_avg_step(i) == 0 )L_error = .true.
        EndDo  !  i = 1, 3
        if ( L_error ) then
          Write ( Unit = 6, fmt=*) 'WARNING GCR iteration counting at'  &
     &,         ' timestep 0 or interval of 0 NOT PERMITTED '
! Following values will output iteration info at 6, 12 and
!  1440 timesteps (30 days for a 30 minute timestep)
          GCR_its_avg_step(1) = 6
          GCR_its_avg_step(2) = 12
          GCR_its_avg_step(3) = 1440
          write(6,*) ' Iteration count diagnostics reset for timesteps '&
     &          , GCR_its_avg_step(1), GCR_its_avg_step(2)              &
     &                               , GCR_its_avg_step(3)
          write(6,*)  '  and at intervals of ', GCR_its_avg_step(3),    &
     &            ' timesteps thereafter.'
          write(6,*)  ' Change in UMUI if you want different values '
        else ! L_error .false.
          write(6,*)  ' '
          write(6,*)  'Iteration count diagnostics at timesteps '       &
     &          , GCR_its_avg_step(1), GCR_its_avg_step(2)              &
     &          , GCR_its_avg_step(3)
          write(6,*)  ' and at intervals of ', GCR_its_avg_step(3),     &
     &            ' timesteps thereafter.'
        endif ! L_error
        GCR_its_switch =  1
        GCR_sum_its = 0
        GCR_min_its = 1000
        GCR_max_its = 0
        GCR_max_time = 0
        GCR_min_time = 0
      endif ! GCR_diagnostics == 3

! ----------------------------------------------------------------------
! 7. Error trapping for advection choices
! ----------------------------------------------------------------------

      do i=1,3
        if (.not. L_mono(i) .and. .not. L_high(i)) then
          if(i == 1)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for theta'
          if(i == 2)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for moisture'
          if(i == 3)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for winds'
          Write ( Unit = 6, fmt=*)                                      &
     &      'Both L_high and L_mono switches set to false '

        endif
      end do

! ----------------------------------------------------------------------
! Orography slope angle and aspect for the radiation code.
! ----------------------------------------------------------------------

      If (L_orog) Then

         If (model_domain == mt_global) Then
            bear_rot_NP = 0.0
         Else
! DEPENDS ON: pole_bearing
            Call pole_bearing(row_length, rows,                         &
     &            lat_rot_NP, long_rot_NP, true_longitude,              &
     &            f3_at_u(1:row_length,1:rows), two_Omega,              &
     &            bear_rot_NP)
         Endif

         If (L_use_grad_corr) Then
! DEPENDS ON: aspang_ancil
            Call aspang_ancil(row_length, rows, land_points,            &
     &            land_sea_mask, grad_x, grad_y, bear_rot_NP)
         Else
! DEPENDS ON: aspang
            Call aspang(row_length,rows, Delta_lambda, Delta_phi,       &
     &            Earth_radius, Orog_halos(0:row_length+1,0:rows+1),    &
     &            cos_theta_latitude(1:row_length,1:rows),              &
     &            bear_rot_NP)
         Endif

         If (model_domain == mt_global) Then
            If (at_extremity(PNorth)) slope_angle(:,rows) = 0.0
            If (at_extremity(PSouth)) slope_angle(:,1) = 0.0
         Endif

      Endif


#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! ----------------------------------------------------------------------
! Boundary layer solver options
! ----------------------------------------------------------------------
      If ( L_us_blsol ) Then
        If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*)  '*** New stable and non-oscillatory boundary-layer&
     & solver is ACTIVE ***'
          write(6,*)  '     It will run with Puns=',Puns,'Pstb=',Pstb
        End If 
      Else
        If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*)'*** Standard boundary-layer solver is ACTIVE ***'
        End If 
      End If
#endif
  
! ----------------------------------------------------------------------
! read in volcanic forcing, VOLCTS
! ----------------------------------------------------------------------
      if (L_VOLCTS) then
      write(6,*) 'SETCONA: Reading in volcanic forcing from file: ',    &
     &            FILE_VOLCTS
      OPEN(UNIT=56,file=FILE_VOLCTS,status='old') 
      do i=1850,2300 
        do j=1,12 
          read(56,*) yvolc,mvolc,volcts(:,j,i) 
        enddo 
      enddo 
      CLOSE(56)
      endif

! ----------------------------------------------------------------------
! read in varying Solar constant, SCVARY
! ----------------------------------------------------------------------
!
! Year-to-year variation of the solar "constant" - from 
! Solanki & Krivova (2003). + solar cycle into 21stC  
!
      if (L_SCVARY) then
      write(6,*) 'SETCONA: Reading in solar forcing from file: ',    &
     &            FILE_SCVARY
      OPEN(UNIT=56,file=FILE_SCVARY,status='old')
      do i=1,601
        read(56,*) ysol,scvary(i)
      enddo
      CLOSE(56)
      endif

  
      RETURN
      END SUBROUTINE Setcona

#endif
