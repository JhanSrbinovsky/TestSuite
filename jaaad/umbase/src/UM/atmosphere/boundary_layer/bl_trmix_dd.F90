#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine BL_TRMIX_DD(                                           &
     &        row_length, rows, bl_levels, off_x, off_y, at_extremity   &
     &,       dtrdz_charney_grid,r_rho_levels,r_theta_levels            &
     &,       halo_i,halo_j                                             &
     &,       tr_levels, tr_vars, model_levels                          &
     &,       DIM_CS2                                                   &
     &,       alpha_cd, rhokh_mix                                       &
     &,       p_star, p, timestep                                       &
! Control logicals
     &,       L_MURK,L_MURK_ADVECT,L_BL_TRACER_MIX,L_DUST,L_CAM_DUST    &
     &,       L_SULPC_SO2,l_sulpc_nh3, l_sulpc_dms, l_soot, l_biomass   &
     &,       l_ocff, l_co2_interactive                                 &
! rml 1/7/13 logical to write co2 fluxes into passive tracers
     &,       l_co2_emits, L_CO2_TRACER, L_USE_CARIOLLE                 &
! Fields to mix
     &,       murk, free_tracers, OZONE_TRACER                          &
! Mineral Dust
     &,       DUST_DIV1,DUST_DIV2,DUST_DIV3                             &
     &,       DUST_DIV4,DUST_DIV5,DUST_DIV6                             &
! Sulphur cycle
     &,       so2, dms, so4_aitken, so4_accu, so4_diss, nh3             &
! Soot cycle
     &,       soot_new, soot_aged, soot_cld                             &
! Biomass aerosol
     &,       bmass_new, bmass_agd, bmass_cld                           &
! Fossil-fuel organic carbon aerosol
     &,       ocff_new, ocff_agd, ocff_cld                              &
! Carbon cycle
     &,       co2, co2_emits, co2flux, npp, resp_s                      &
     &,       co2_flux_tot, land_co2                                    &
! Tracer fluxes - kdcorbin, 05/10
     &,       tracer_flux1, tracer_flux2, tracer_flux3, tracer_flux4    &
     &,       tracer_flux5, tracer_flux6, tracer_flux7, tracer_flux8    &
     &,       tracer_flux9, tracer_flux10,tracer_flux11,tracer_flux12   &
     &,       tracer_flux13,tracer_flux14,tracer_flux15,tracer_flux16   &
     &,       tracer_flux17,tracer_flux18,tracer_flux19,tracer_flux20   &
! Emissions fields
     &,       DUST_FLUX                                                 &
     &,       so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em    &
     &,       ocff_hilem, ocff_em, drydep_str                           &
     &,       kent, we_lim, t_frac, zrzi                                &
     &,       kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                &
     &,       zh, zhsc, z_half                                          &
! IN for dry deposition of tracers
     &,       rho_aresist, aresist, resist_b                            &
     &,       R_B_DUST,TSTAR                                            &
     &,       land_points, land_index, ice_frac                         &
! IN variables for MOSES II
     &,       ntype, ntiles, tile_pts, tile_index, tile_frac            &
     &,       canopy, catch, snow_tile, gc                              &
     &,       aresist_tile, resist_b_tile, FLANDG                       &
     &,                                                                 &
! STASH related variables
#include "argsts.h"
     &        stashwork                                                 &
     &        )
!
! ---------------------------------------------------------------------
! Purpose: Interface level to tr_mix for tracer variables that
!          require boundary layer mixing and/or dry deposition
!          with MOSES II surface scheme.
!
!          Called by BL_TRACER_INTCTL.
!
! Current code owner: M Woodage
!
! History:
! Version   Date   Comment
! -------   ----   -------
!  5.3   20/10/01  New deck, based on vn5.2 deck BL_TRACER_MIX.
!                                                             M Woodage
!  5.4   05/09/02   Initialise RES_FACTOR array to zero to prevent
!                      unset values being used.              M Woodage
!  5.4   29/08/02  Extra variables passed through for changes
!                  to tr_mix.                   P. Davison.
!  5.4   10/07/02  Include mixing and dry deposition of 3 modes of soot
!                  and output diagnostics.                  P. Davison
!    5.5  26/02/03  Implement mixing for carbon cycle. C Jones
!  5.5   05/02/03  Include mixing and dry deposition of biomass smoke
!                  and output diagnostics.                  P. Davison
!  5.5   19/02/03  Remove redundent RML code. Adrian Lock
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!  6.2  15/02/06  Use optimised  trmix for dust. Stephanie Woodward

!  6.2   15/02/06  Fix to trmix arg list for dust     Stephanie Woodward
!
! 6.2      15/12/05  Correct dust diagnostics when substepping phys2.
!                                                        M. Diamantakis
!  6.2  01/03/06  pass CO2 fluxes back up for use in diagnostics
!                 routines.                             C.D. Jones
!
!
!  6.4  03/01/07  Use standard tr_mix, now it's been optimised, for
!                      dust.                        Stephanie Woodward
!
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      Implicit None

#include "c_dust_ndiv.h"

!
! Parallel setup variables
      Integer, Intent(In) ::                                            &
     &  off_x                                                           &
                              ! Size of small halo in i
     &, off_y                                                           &
                              ! Size of small halo in j.
     &, halo_i                                                          &
                              ! Size of halo in i direction
     &, halo_j                ! Size of halo in j direction

       Logical, Intent(In) ::                                           &
     &  at_extremity(4)   ! Indicates if this processor is at north,
                          ! south, east or west of the processor grid

! Model dimensions
      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, bl_levels                                                       &
     &, tr_levels                                                       &
                          ! number of free tracer levels
     &, tr_vars                                                         &
                          ! number of free tracer variables
     &, DIM_CS2         ! soil carbon dimension


! Model switches
      Logical, Intent(In) ::                                            &
     &  L_Murk                                                          &
                          ! Switch for (visibility) aerosol
     &, L_murk_advect                                                   &
                          ! Switch for advecting aerosol
     &, L_bl_tracer_mix                                                 &
                          ! Switch for BL mixing of free tracers
     &, L_DUST                                                          &
                          ! switch for mineral dust
     &, L_CAM_DUST                                                      &
                          !Old version of dust_uplift scheme used in CAM NWP models
     &, L_sulpc_so2                                                     &
                          ! Switch for Sulphur Cycle
     &, L_sulpc_nh3                                                     &
                          ! NH3 included in Sulphur Cycle
     &, L_sulpc_dms                                                     &
                          ! DMS included in Sulphur Cycle
     &, L_soot                                                          &
                          ! Switch for Soot Cycle
     &, L_biomass                                                       &
                          ! Switch for Biomass aerosol
     &, L_ocff                                                          &
                          ! Switch for Fossil-fuel OC aerosol
     &, L_co2_interactive                                               &
                          ! Switch for interactive CO2
     &, L_co2_emits                                                     &
                          ! Switch for anthro CO2 emissions
! rml 2/7/13
     &, L_CO2_TRACER                                                    &
                          ! Switch for putting co2 flux into passive tracer
     &, L_USE_CARIOLLE
                          ! Switch for cariolle ozone tracer scheme
! Model parameters
      Real, Intent(In) ::                                               &
     &  timestep                                                        &
     &, alpha_cd(bl_levels)                                             &
     &, rhokh_mix (row_length, rows, bl_levels)                         &
     &, we_lim(row_length,rows,3)                                       &
                                      !rho*entrainment rate implied by
                                      !placing of subsidence
     &, zrzi(row_length,rows,3)                                         &
                                      !(z-z_base)/(z_i-z_base)
     &, t_frac(row_length,rows,3)                                       &
                                      !a fraction of the timestep
     &, we_lim_dsc(row_length,rows,3)                                   &
                                      !rho*entrainment rate implied by
                                      !  placing of subsidence
     &, zrzi_dsc(row_length,rows,3)                                     &
                                      !(z-z_base)/(z_i-z_base)
     &, t_frac_dsc(row_length,rows,3)                                   &
                                      !a fraction of the timestep
     &, z_half(row_length,rows,BL_LEVELS)                               &
                                      !Z_HALF(*,K) is height of half
                                      !    level k-1/2.
     &, zhsc(row_length,rows)                                           &
                                      !Top of decoupled layer
     &, zh  (row_length,rows)                                           &
                                      !Top of surface mixed layer
     &, dtrdz_charney_grid(row_length,rows,BL_LEVELS)
                                      ! For tracer mixing

      Integer, Intent(In) ::                                            &
     &  kent(row_length,rows)                                           &
                                      !grid-level of SML inversion
     &, kent_dsc(row_length,rows)                                       &
                                      !grid-level of DSC inversion
     &, land_points                                                     &
                                      !No.of land points being processed
     &, land_index(land_points)                                         &
                                      !set from land_sea_mask
! For MOSES II
     &, ntype                                                           &
                                      !No. of tile types
     &, ntiles                                                          &
                                      !No.of land-surface tiles
     &, tile_pts(ntype)                                                 &
                                      !No.of tile points
     &, tile_index(land_points,ntype) !Index of tile points
!
! Required fields (input)
      Real, Intent(In) ::                                               &
     &  p_star(row_length, rows)                                        &
     &, p(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                   1-halo_j:rows+halo_j,0:bl_levels)              &
                                      ! Heights of theta levels
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, bl_levels)
                                      ! Heights of rho levels

! Emissions fields (input)
      REAL, INTENT(IN) ::                                               &
     &  DUST_FLUX ( ROW_LENGTH, ROWS, NDIV ) !dust emission flux

      Real, Intent(In) ::                                               &
     &  so2_hilem ( row_length, rows )                                  &
     &, so2_em    ( row_length, rows )                                  &
     &, nh3_em    ( row_length, rows )                                  &
     &, dms_em    ( row_length, rows )                                  &
     &, soot_hilem( row_length, rows )                                  &
     &, soot_em   ( row_length, rows )                                  &
     &, ocff_hilem( row_length, rows )                                  &
     &, ocff_em   ( row_length, rows )                                  &
     &, co2_emits ( row_length, rows )                                  &
                                        ! anthro CO2 emissions
                                        !  (Kg CO2 m-2 s-1)
     &, co2flux   ( row_length, rows )  ! ocean  CO2 emissions
                                        !  (Kg CO2 m-2 s-1)

! Tracer fluxes - kdcorbin, 05/10
! rml 12/9/11 allow tracer_flux1 to be InOut
! rml 1/7/13 change tracer_flux2 and tracer_flux3 to InOut
      Real, Intent(InOut) ::                                            &
    &     tracer_flux1(row_length, rows )                               &
    &,    tracer_flux2(row_length, rows )                               &
    &,    tracer_flux3(row_length, rows )
      Real, Intent(In) ::                                               &
    &                                      tracer_flux4(row_length, rows )  &
    &,    tracer_flux5(row_length, rows ), tracer_flux6(row_length, rows )  &
    &,    tracer_flux7(row_length, rows ), tracer_flux8(row_length, rows )  &
    &,    tracer_flux9(row_length, rows ), tracer_flux10(row_length,rows )  &
    &,    tracer_flux11(row_length, rows), tracer_flux12(row_length,rows )  &
    &,    tracer_flux13(row_length, rows), tracer_flux14(row_length,rows )  &
    &,    tracer_flux15(row_length, rows), tracer_flux16(row_length,rows )  &
    &,    tracer_flux17(row_length, rows), tracer_flux18(row_length,rows )  &
    &,    tracer_flux19(row_length, rows), tracer_flux20(row_length,rows ) 

! For dry deposition of tracers (input)
      Real                                                              &
     &  rho_aresist(row_length,rows)                                    &
                                       !RHOSTAR*CD_STD*VSHR
     &, aresist(row_length,rows)                                        &
                                       !1/(CD_STD*VSHR)
     &, resist_b(row_length,rows)                                       &
                                       !(1/CH-1/CD_STD)/VSHR
     &, R_B_DUST( ROW_LENGTH, ROWS,NDIV )                               &
                                          !surf layer res for dust
     &, TSTAR(ROW_LENGTH,ROWS)                                          &
                               ! temperature

     &, ice_frac(row_length,rows)                                       &
                                       !sea_ice fractn (ice/qrclim.ice)
! For MOSES II
     &, tile_frac(land_points,ntiles)                                   &
                                       !fractional coverage for each
!                                       surface tile
     &, canopy(land_points,ntiles)                                      &
                                       !Surface/canopy water (kg/m2)
     &, catch(land_points,ntiles)                                       &
                                       !Surface/canopy water capacity
!                                       of snow-free land tiles (kg/m2)
     &, snow_tile(land_points,ntiles)                                   &
                                       !Snow on tiles (kg/m2).
     &, gc(land_points,ntiles)                                          &
                                       !Stomatal conductance to evapn
!                                       for land tiles (m/s)
     &, aresist_tile(land_points,ntiles)                                &
!                                      ! 1/(CD_STD*VSHR) on land tiles
     &, resist_b_tile(land_points,ntiles)                               &
!                                    !(1/CH-1/CD_STD)/VSHR on land tiles
     &, FLANDG(row_length,rows)                                         &
                                     !land fraction of grid box
!                                       (0 if all sea, 1 if all land)
     &, npp(land_points)                                                &
                                       ! net primary productivity
                                       !  (Kg C m-2 s-1)
     &, resp_s(DIM_CS2)                                                 &
                                       ! soil respiration
                                       !  (Kg C m-2 s-1)
     &, co2_flux_tot(row_length,rows)                                   &
                                       ! total CO2 flux
                                       !  (Kg CO2 m-2 s-1)
     &, drydep_str(row_length, rows)                                    &

     &, land_co2(land_points)             ! terrestrial CO2 flux
                                       !  (Kg CO2 m-2 s-1)
! Tracer fields for mixing
      Real, Intent(InOut) ::                                            &
     &  murk        ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,free_tracers( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y,                         &
     &                tr_levels, tr_vars)

      REAL, INTENT(INOUT) ::                                            &
     &  DUST_DIV1   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV2   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV3   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV4   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV5   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )          &
     &, DUST_DIV6   ( 1 - OFF_X : ROW_LENGTH + OFF_X,                   &
     &                1 - OFF_Y : ROWS + OFF_Y, MODEL_LEVELS )
      Real, Intent(InOut) ::                                            &
     &  so2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,dms         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_aitken  ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_accu    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,so4_diss    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,nh3         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )

      Real, Intent(InOut) ::                                            &
     &  soot_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_aged   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,soot_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_new   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_agd   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,bmass_cld   ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_new    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_agd    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,ocff_cld    ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,co2         ( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )          &
     & ,OZONE_TRACER( 1 - off_x : row_length + off_x,                   &
     &                1 - off_y : rows + off_y, model_levels )  
     
      Real, Intent(InOut) ::                                            &
     &  stashwork(*)

! Argument comdeck declaration
#include "csubmodl.h"
#include "typsts.h"

! Local variables
      Integer, Parameter :: sect = 3               ! BL section
      Integer            :: item                   ! stash item code
      Integer            :: im_index               ! internal model
      Integer            :: n_tracer               ! looper
      Integer            :: ErrorStatus            ! error reporting
      INTEGER            :: IDIV            ! loop counter, dust divs
      INTEGER            :: LEV1            ! =1 no. levs for vgrav calc

      Real        :: zero_field(row_length, rows)  ! field of 0s
      Real        :: tr_flux( row_length, rows                          &
                                                   ! Fluxes returned
     &                      , bl_levels)           ! from mixing

      Character (Len=*), Parameter :: RoutineName='bl_tracer_mixing'
      Character (Len=80)           :: Cmessage
!
      Integer i, j, K, l                !loop counters
!
      Real                                                              &
     &  work_2d(row_length,rows)                                        &
                                        !2d work array for stash
     &, STR_RESIST_B(row_length,rows)                                   &
                                        !Rb for S Cycle dry dep
     &, STR_RESIST_S(row_length,rows)                                   &
                                        !Rs for S Cycle dry dep
     &, RES_FACTOR(row_length,rows)                                     &
                                        !Ra/(Ra+Rb+Rs) for dry dep
     &, RES_FACTOR_LAND(row_length,rows)                                &
                                        !Ra/(Ra+Rb+Rs) mean over land
     &, RESIST_S(row_length*rows)                                       &
                                        !stomatal resistance
     &, VSTOKES1(ROW_LENGTH,ROWS)                                       &
                                 !dust settling velocity, lowest layer
     &, DRYDEP_DUST(ROW_LENGTH,ROWS,NDIV)                               &
                                          ! dust deposition flux
     &, WORK1(ROW_LENGTH,ROWS,BL_LEVELS,NDIV)                           &
                                              ! workspace
     &, WORK2(ROW_LENGTH,ROWS)                                          &
     &, WORK3(ROW_LENGTH,ROWS)                                          &
     &, SNOW_F                          !calculated snow fraction

!kdcorbin, 05/10 - adding variable to store sources for tracers
     Real :: MYTRACER_FLUX(ROW_LENGTH,ROWS)

!kdcorbin, 04/10 - tracer flux information
INTEGER, parameter :: flux_vars=20
integer :: n_tr

!
! comdeck containing RESB and RESS parameters
#include "c_sulbdy.h"
#include "ccarbon.h"
#include "c_mdi.h"
#include "c_dustgen.h"
! comdeck containing RESB and RESS parameters for soot
#include "c_st_bdy.h"
! comdeck containing RESB and RESS parameters for biomass aerosol
#include "c_bm_bdy.h"
! comdeck containing RESB and RESS parameters for OCFF aerosol
#include "c_ocff_bdy.h"
!
! Externals
       External tr_mix
       External ereport
       External copydiag
       External copydiag_3d
       External sresfact
       EXTERNAL VGRAV

!
!---------------------------------------------------------------------
! Initial setup
!---------------------------------------------------------------------
      zero_field(:,:) = 0.0
      ErrorStatus = 0
      im_index = internal_model_index(atmos_im)
!
!  Initialise RES_FACTOR array to zero to prevent use of unset values
        Do j=1,rows
          Do i=1,row_length
            RES_FACTOR(i,j) = 0.0
          End Do
        End Do

!---------------------------------------------------------------------
! Mixing for Aerosol
!---------------------------------------------------------------------
      If ( l_murk_advect ) Then
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                         &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, murk( 1:row_length, 1:rows, 1:bl_levels )         &
     &      ,zero_field, zero_field, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )

        If (ErrorStatus /= 0) Then
          Cmessage = 'Problem with mixing aerosol'
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

        item = 129
        If ( sf(item,sect) ) Then      ! diagnostic flux required
! DEPENDS ON: copydiag_3d
          Call copydiag_3d( stashwork(si(item,sect,im_index)),          &
     &        tr_flux,                                                  &
     &        row_length,rows,bl_levels,0,0,0,0, at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus ,cmessage)
        End If

        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

      End If

!
!------------------------------------------------------------------
! Mixing for mineral dust, including dry deposition
!------------------------------------------------------------------
!
      IF (L_DUST) THEN

        LEV1 = 1

        DO K = 1, BL_LEVELS
          DO J = 1,ROWS
            DO I = 1,ROW_LENGTH
              WORK1(I,J,K,1)=DUST_DIV1(I,J,K)
              WORK1(I,J,K,2)=DUST_DIV2(I,J,K)
              WORK1(I,J,K,3)=DUST_DIV3(I,J,K)
              WORK1(I,J,K,4)=DUST_DIV4(I,J,K)
              WORK1(I,J,K,5)=DUST_DIV5(I,J,K)
              WORK1(I,J,K,6)=DUST_DIV6(I,J,K)
            ENDDO !ROW_LENGTH
          ENDDO !ROWS
        ENDDO !BL_LEVELS

        DO IDIV = 1,NDIV

! DEPENDS ON: vgrav
          CALL VGRAV(                                                   &
     &  ROW_LENGTH,ROWS,LEV1,DREP(IDIV),                                &
     &  RHOP,P_STAR,TSTAR,                                              &
     &  VSTOKES1,WORK2,WORK3                                            &
     &  )
         IF (L_CAM_DUST) THEN 
!           Using NWP version of the dust scheme, which moves
!           settling from level 1 to gravsett
            DO J = 1,ROWS
               DO I = 1,ROW_LENGTH
                  RES_FACTOR(I,J)=ARESIST(I,J)/                         & 
     &             (ARESIST(I,J)+R_B_DUST(I,J,IDIV))
               ENDDO !ROW_LENGTH
            ENDDO !ROWS
         ELSE
!           Using HadGEM version of the tracer scheme
            DO J = 1,ROWS
               DO I = 1,ROW_LENGTH
                  RES_FACTOR(I,J)=ARESIST(I,J) * ( VSTOKES1(I,J) + 1./  &
     &             (ARESIST(I,J)+R_B_DUST(I,J,IDIV)+                    &
     &              ARESIST(I,J)*R_B_DUST(I,J,IDIV)*VSTOKES1(I,J)) )
               ENDDO !ROW_LENGTH
            ENDDO !ROWS
         ENDIF

! DEPENDS ON: tr_mix
       CALL TR_MIX(                                                     &
     &       HALO_I, HALO_J, ROW_LENGTH, ROWS, BL_LEVELS                &
     &      ,OFF_X, OFF_Y, ALPHA_CD                                     &
     &      ,RHOKH_MIX(1,1,2), RHO_ARESIST                              &
     &      ,DTRDZ_CHARNEY_GRID, R_RHO_LEVELS, R_THETA_LEVELS           &
     &      ,TIMESTEP                                                   &
     &      ,TR_FLUX, WORK1( 1:ROW_LENGTH, 1:ROWS, 1:BL_LEVELS,IDIV )   &
     &      ,DUST_FLUX(1:ROW_LENGTH, 1:ROWS, IDIV), RES_FACTOR          &
     &      ,DRYDEP_STR                                                 &
     &      ,KENT, WE_LIM, T_FRAC, ZRZI                                 &
     &      ,KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                 &
     &      ,ZH ,ZHSC, Z_HALF                                           &
     &      ,ERRORSTATUS, .FALSE.                                       &
     &       )
!
         IF (ERRORSTATUS /=0) THEN
! DEPENDS ON: ereport
           CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
         ENDIF

         ITEM = 440 + IDIV      !dust dry dep flux, from lowest layer
           IF ( SF(ITEM,SECT) ) THEN
!
! Change sign of dry dep flux (otherwise negative)
             DO J=1,ROWS
               DO I=1,ROW_LENGTH
                 WORK_2D(I,J) = -DRYDEP_STR(I,J)
               END DO
             END DO
!
! DEPENDS ON: copydiag
             CALL COPYDIAG(STASHWORK(SI(ITEM,SECT,IM_INDEX)),WORK_2D,   &
     &        ROW_LENGTH,ROWS,0,0,0,0,AT_EXTREMITY,                     &
     &        ATMOS_IM,SECT,ITEM,                                       &
     &        ERRORSTATUS,CMESSAGE)
!
            IF (ERRORSTATUS /=0) THEN
! DEPENDS ON: ereport
              CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
            ENDIF
!
          ENDIF !stashflag
        ENDDO !NDIV

        DO K = 1, BL_LEVELS
          DO J = 1,ROWS
            DO I = 1,ROW_LENGTH
              DUST_DIV1(I,J,K)=WORK1(I,J,K,1)
              DUST_DIV2(I,J,K)=WORK1(I,J,K,2)
              DUST_DIV3(I,J,K)=WORK1(I,J,K,3)
              DUST_DIV4(I,J,K)=WORK1(I,J,K,4)
              DUST_DIV5(I,J,K)=WORK1(I,J,K,5)
              DUST_DIV6(I,J,K)=WORK1(I,J,K,6)
            ENDDO !ROW_LENGTH
          ENDDO !ROWS
        ENDDO !BL_LEVELS

      ENDIF !L_DUST
!------------------------------------------------------------------
! Mixing for Sulphur Cycle variables, including dry deposition
!------------------------------------------------------------------
!
      IF (l_sulpc_so2) THEN
!
!  Calculate resistance factors for species to be dry deposited.
!
! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.
!
!   For RESIST_B check to eliminate possible negative values
        Do j=1,rows
          Do i=1,row_length
            If (RESIST_B(i,j)  <   0.0)     Then
              RESIST_B(i,j) = 0.0
            End If
          End Do
        End Do
!
!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.
!
!   Initialise 'stomatal' resistance array to zero:
!
        Do K=1,row_length*rows
          RESIST_S(K)=0.0
        End Do
!
!
! CODE FOR SO2
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_SO2*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_SO2*RESIST_S(K)
!
!  Need to set STR_RESIST_S = R_SNOW for SO2 over sea ice
!  (0 is acceptable over sea).
                If (ice_frac(i,j)  >   0.0 )  Then
                  STR_RESIST_S(i,j)=R_SNOW
                End If
!
!  Calculate RES_FACTOR for SO2 (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for SO2 to calculate correct RES_FACTOR_LAND values
! if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &                 NTILES, TILE_INDEX, TILE_PTS, .TRUE.,            &
     &                 CANOPY, CATCH,                                   &
     &                 GC, TILE_FRAC,                                   &
     &                 SNOW_TILE,                                       &
     &                 ARESIST, ARESIST_TILE, RESIST_B_TILE,            &
     &                 RESB_SO2, RESS_SO2, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for SO2 combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
!
! Mix SO2   (Note emiss added in AERO_CTL via TRSRCE call)
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, so2( 1:row_length, 1:rows, 1:bl_levels )          &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 270                      !dry dep flux SO2
        IF ( sf(item,sect) ) THEN
!
! Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF

! CODE FOR NH3 (if present)
!
        If (l_sulpc_nh3) Then
!
!  Calculate RES_FACTOR for NH3 to allow dry deposition in same way
!  as for SO2 (including code for snow and ice)
!
          Do j=1,rows
            Do i=1,row_length
              K=(j-1)*row_length + i
!
              IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
                STR_RESIST_B(i,j) = RESB_NH3*RESIST_B(i,j)
                STR_RESIST_S(i,j) = RESS_NH3*RESIST_S(K)
!
!  Need to set STR_RESIST_S = R_SNOW for SO2 over sea ice
!  (0 is acceptable over sea).
                  If (ice_frac(i,j)  >   0.0 )  Then
                    STR_RESIST_S(i,j)=R_SNOW
                  End If
!
!  Calculate RES_FACTOR for NH3  (only correct for sea points)
!
                RES_FACTOR(i,j) = aresist(i,j) /                        &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
              END IF
!
            End Do
          End Do
!
! Call SRESFACT for NH3 to calculate correct RES_FACTOR_LAND values
! if land present:
          IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
            Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,   &
     &                 NTILES, TILE_INDEX, TILE_PTS, .TRUE.,            &
     &                 CANOPY, CATCH,                                   &
     &                 GC, TILE_FRAC,                                   &
     &                 SNOW_TILE,                                       &
     &                 ARESIST, ARESIST_TILE, RESIST_B_TILE,            &
     &                 RESB_NH3, RESS_NH3, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for NH3 combining sea and land values
            Do j=1,rows
              Do i=1,row_length
                RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +   &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
              End Do
            End Do
!
          END IF
!
! Mix NH3   (Note emiss added in AERO_CTL via TRSRCE call)
!
! DEPENDS ON: tr_mix
          Call tr_mix(                                                  &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, nh3( 1:row_length, 1:rows, 1:bl_levels )          &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
! Write diagnostics to STASH
!
          item = 300                      !dry dep flux NH3
          IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
            Do j=1,rows
              Do i=1,row_length
                work_2d (i,j) = -drydep_str(i,j)
              End Do
            End Do
!
! DEPENDS ON: copydiag
            Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,    &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
            If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
              Call Ereport( RoutineName, ErrorStatus, Cmessage )
            End If
!
          END IF
!
        End If                        ! End If nh3 condition

! CODE FOR SO4_AIT
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_SO4_AIT*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_SO4_AIT*RESIST_S(K)
!
!  Calculate RES_FACTOR for SO4_AIT  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
         End Do
!
! Call SRESFACT for SO4_AIT to calculate correct RES_FACTOR_LAND values
! if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_SO4_AIT, RESS_SO4_AIT, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for SO4_AIT combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix so4_aitken
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, so4_aitken( 1:row_length, 1:rows, 1:bl_levels )   &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 271                      !dry dep flux SO4_AIT
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF
!
! CODE FOR SO4_ACC
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_SO4_ACC*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_SO4_ACC*RESIST_S(K)
!
!  Calculate RES_FACTOR for SO4_ACC  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for SO4_ACC to calculate correct RES_FACTOR_LAND values
! if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_SO4_ACC, RESS_SO4_ACC, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for SO4_ACC combining sea and land values
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix so4_accu
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, so4_accu( 1:row_length, 1:rows, 1:bl_levels )     &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 272                      !dry dep flux SO4_ACC
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF
!
! CODE FOR SO4_DIS
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_SO4_DIS*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_SO4_DIS*RESIST_S(K)
!
!  Calculate RES_FACTOR for SO4_DIS  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for SO4_ACC to calculate correct RES_FACTOR_LAND values
! if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_SO4_DIS, RESS_SO4_DIS, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for SO4_DIS combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix so4_diss
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, so4_diss( 1:row_length, 1:rows, 1:bl_levels )     &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 273                      !dry dep flux SO4_DIS
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF
!
! CODE FOR DMS (if present)
!
        If (l_sulpc_dms) Then
!
! (Note no dry dep for DMS, and emiss added in AERO_CTL via TRSRCE call)
!
! Mix DMS
!
! DEPENDS ON: tr_mix
          Call tr_mix(                                                  &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                         &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, dms( 1:row_length, 1:rows, 1:bl_levels )          &
     &      ,zero_field, zero_field, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        End If                          ! End If dms condition
!
      END IF                            ! END IF l_sulpc_so2
!
!---------------------------------------------------------------------
! Mixing for Carbon cycle
!---------------------------------------------------------------------
      If ( l_co2_interactive ) Then

        DO I=1,ROW_LENGTH
          DO J=1,ROWS
            CO2_FLUX_TOT(I,J)=0.0

!  (i) CO2 emissions from ancillary file.
            IF(L_CO2_EMITS) THEN
              IF ( CO2_EMITS(I,J)  /=  RMDI ) THEN
                CO2_FLUX_TOT(I,J)=CO2_FLUX_TOT(I,J) + CO2_EMITS(I,J)
              ENDIF
            ENDIF

!  (ii) CO2 flux from ocean. (+ve implies air to sea)
            IF ( CO2FLUX(I,J)  /=  RMDI ) THEN
              CO2_FLUX_TOT(I,J)=CO2_FLUX_TOT(I,J) - CO2FLUX(I,J)
            ENDIF

          ENDDO
        ENDDO

!  (iii) CO2 flux from land processes. (+ve implies biosphere to atmos)
        DO L=1,LAND_POINTS
          J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I =  LAND_INDEX(L) - (J-1)*ROW_LENGTH
          LAND_CO2(L) = (RESP_S(L) - NPP(L)) * M_CO2 / M_CARBON
          CO2_FLUX_TOT(I,J)=CO2_FLUX_TOT(I,J) + LAND_CO2(L)
        ENDDO

! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                         &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, co2( 1:row_length, 1:rows, 1:bl_levels )          &
     &      ,co2_flux_tot, zero_field, drydep_str                       &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )

        If (ErrorStatus /= 0) Then
          Cmessage = 'Problem with CO2 mixing'
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

      End If   ! C-cycle
!
!
!------------------------------------------------------------------
! Mixing for Soot variables, including dry deposition
!------------------------------------------------------------------
!
! Note: Emissions of soot are dealt with in aero_ctl. Here we consider
! just boundary layer mixing and dry deposition of the three modes
! of soot.
!
      IF (L_SOOT)  THEN   ! If soot modelling is included
!
!  Calculate resistance factors for species to be dry deposited.
!
! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.
!
!   For RESIST_B check to eliminate possible negative values
        Do j=1,rows
          Do i=1,row_length
            If (RESIST_B(i,j)  <   0.0)     Then
              RESIST_B(i,j) = 0.0
            End If
          End Do
        End Do
!
!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.
!
!   Initialise 'stomatal' resistance array to zero:
!
        Do K=1,row_length*rows
          RESIST_S(K)=0.0
        End Do
!
! CODE FOR FRESH SOOT
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_FreshSoot*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_Soot*RESIST_S(K)
!
!  Calculate RES_FACTOR for fresh soot  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for fresh soot to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_FreshSoot, RESS_Soot, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for fresh soot combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix fresh soot
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, soot_new( 1:row_length, 1:rows, 1:bl_levels )     &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 301                      !dry dep flux fresh soot
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
! CODE FOR AGED SOOT
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_AgedSoot*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_Soot*RESIST_S(K)
!
!  Calculate RES_FACTOR for aged soot  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for aged soot to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_AgedSoot, RESS_Soot, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for aged soot combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix aged soot
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, soot_aged( 1:row_length, 1:rows, 1:bl_levels )    &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 302                      !dry dep flux aged soot
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
! CODE FOR CLOUD SOOT
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_SootInCloud*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_Soot*RESIST_S(K)
!
!  Calculate RES_FACTOR for cloud soot  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for cloud soot to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_SootInCloud, RESS_Soot, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for cloud soot combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix cloud soot
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, soot_cld( 1:row_length, 1:rows, 1:bl_levels )     &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 303                   !dry (occult) dep flux cloud soot
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
      END IF  ! L_SOOT Soot modelling included
!
! END of soot section
!
!
!------------------------------------------------------------------
! Mixing for Biomass aerosol, including dry deposition
!------------------------------------------------------------------
!
! Note: Emissions of biomass aerosol are dealt with in aero_ctl.
! Here we consider just boundary layer mixing and dry deposition.
!
      IF (L_BIOMASS)  THEN   ! If biomass aerosol is included
!
!  Calculate resistance factors for species to be dry deposited.
!
! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.
!
!   For RESIST_B check to eliminate possible negative values
        Do j=1,rows
          Do i=1,row_length
            If (RESIST_B(i,j)  <   0.0)     Then
              RESIST_B(i,j) = 0.0
            End If
          End Do
        End Do
!
!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.
!
!   Initialise 'stomatal' resistance array to zero:
!
        Do K=1,row_length*rows
          RESIST_S(K)=0.0
        End Do
!
! CODE FOR FRESH BIOMASS
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_FreshBmass*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_Bmass*RESIST_S(K)
!
!  Calculate RES_FACTOR for fresh biomass  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for fresh biomass to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_FreshBmass, RESS_Bmass, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for fresh biomass combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix fresh biomass smoke
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, bmass_new( 1:row_length, 1:rows, 1:bl_levels )    &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 396                      !dry dep flux fresh biomass
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
! CODE FOR AGED BIOMASS
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_AgedBmass*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_Bmass*RESIST_S(K)
!
!  Calculate RES_FACTOR for aged biomass  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for aged biomass to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_AgedBmass, RESS_Bmass, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for aged smoke combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix aged smoke
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, bmass_agd( 1:row_length, 1:rows, 1:bl_levels )    &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 397                      !dry dep flux aged biomass
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
! CODE FOR CLOUD BIOMASS
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j) <  1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_BmassInCloud*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_Bmass*RESIST_S(K)
!
!  Calculate RES_FACTOR for cloud smoke  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
     &           (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for cloud smoke to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS >  0) THEN
!
! DEPENDS ON: sresfact
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
     &               NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
     &               CANOPY, CATCH,                                     &
     &               GC, TILE_FRAC,                                     &
     &               SNOW_TILE,                                         &
     &               ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
     &               RESB_BmassInCloud, RESS_Bmass, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for cloud smoke combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
     &                         FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix cloud smoke
!
! DEPENDS ON: tr_mix
        Call tr_mix(                                                    &
     &       halo_i, halo_j, row_length, rows, bl_levels                &
     &      ,off_x, off_y, alpha_cd                                     &
     &      ,rhokh_mix(1,1,2), rho_aresist                              &
     &      ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
     &      ,timestep                                                   &
     &      ,tr_flux, bmass_cld( 1:row_length, 1:rows, 1:bl_levels )    &
     &      ,zero_field, RES_FACTOR, drydep_str                         &
     &      ,kent, we_lim, t_frac, zrzi                                 &
     &      ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
     &      ,zh ,zhsc, z_half                                           &
     &      ,ErrorStatus, .false.                                       &
     &       )
!
        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 398                   !dry (occult) dep flux cloud smoke
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
! DEPENDS ON: copydiag
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
     &        atmos_im,sect,item,                                       &
     &        ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
      END IF  ! L_BIOMASS Biomass smoke modelling included
!
! END of biomass section
!
!
!
!------------------------------------------------------------------
! Mixing for Fossil-fuel organic carbon, including dry deposition
!------------------------------------------------------------------
!
! Note: Emissions of ocff are dealt with in aero_ctl. Here we consider
! just boundary layer mixing and dry deposition of the three modes
! of ocff.
!
      IF (L_ocff)  THEN   ! If ocff modelling is included
!
!  Calculate resistance factors for species to be dry deposited.
!
! Note that for MOSES II with coastal tiling, at coastal points
! (i.e. those containing both land and sea):
! RESIST_B here contains  values approriate for SEA only (i.e. no
! meaning including land tiles has been done in SFEXCH8A), but
! ARESIST and RHO_ARESIST contain grid box mean values calculated
! in SFEXCH8A taking account of land/sea fractions.
!
!   For RESIST_B check to eliminate possible negative values
        Do j=1,rows
          Do i=1,row_length
            If (RESIST_B(i,j) .LT. 0.0)     Then
              RESIST_B(i,j) = 0.0
            End If
          End Do
        End Do
!
!   For STR_RESIST_S values depend on surface type (land, sea,
!    snow,ice) as well as tracer identity.
!
!   Initialise 'stomatal' resistance array to zero:
!
        Do K=1,row_length*rows
          RESIST_S(K)=0.0
        End Do
!
! CODE FOR FRESH OCFF
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j).LT.1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_Freshocff*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_ocff*RESIST_S(K)
!
!  Calculate RES_FACTOR for fresh ocff  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
                 (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for fresh ocff to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS.GT.0) THEN
!
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
                     NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
                     CANOPY, CATCH,                                     &
                     GC, TILE_FRAC,                                     &
                     SNOW_TILE,                                         &
                     ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
                     RESB_Freshocff, RESS_ocff, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for fresh ocff combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
                               FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix fresh ocff
!
        Call tr_mix(                                                    &
             halo_i, halo_j, row_length, rows, bl_levels                &
            ,off_x, off_y, alpha_cd                                     &
            ,rhokh_mix(1,1,2), rho_aresist                              &
            ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
            ,timestep                                                   &
            ,tr_flux, ocff_new( 1:row_length, 1:rows, 1:bl_levels )     &
            ,zero_field, RES_FACTOR, drydep_str                         &
            ,kent, we_lim, t_frac, zrzi                                 &
            ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
            ,zh ,zhsc, z_half                                           &
            ,ErrorStatus, .false.                                       &
             )
!
        If (ErrorStatus /=0) Then
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 407                        !dry dep flux fresh ocff
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
              row_length,rows,0,0,0,0,at_extremity,                     &
              atmos_im,sect,item,                                       &
              ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
! CODE FOR AGED OCFF
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j).LT.1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_Agedocff*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_ocff*RESIST_S(K)
!
!  Calculate RES_FACTOR for aged ocff  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
                 (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for aged ocff to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS.GT.0) THEN
!
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
                     NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
                     CANOPY, CATCH,                                     &
                     GC, TILE_FRAC,                                     &
                     SNOW_TILE,                                         &
                     ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
                     RESB_Agedocff, RESS_ocff, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for aged ocff combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
                               FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix aged OCFF
!
        Call tr_mix(                                                    &
             halo_i, halo_j, row_length, rows, bl_levels                &
            ,off_x, off_y, alpha_cd                                     &
            ,rhokh_mix(1,1,2), rho_aresist                              &
            ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
            ,timestep                                                   &
            ,tr_flux, ocff_agd( 1:row_length, 1:rows, 1:bl_levels )     &
            ,zero_field, RES_FACTOR, drydep_str                         &
            ,kent, we_lim, t_frac, zrzi                                 &
            ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
            ,zh ,zhsc, z_half                                           &
            ,ErrorStatus, .false.                                       &
             )
!
        If (ErrorStatus /=0) Then
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 408                  ! dry dep flux aged ocff
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
              row_length,rows,0,0,0,0,at_extremity,                     &
              atmos_im,sect,item,                                       &
              ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
! CODE FOR CLOUD OCFF
!
        Do j=1,rows
          Do i=1,row_length
            K=(j-1)*row_length + i
!
            IF (FLANDG(i,j).LT.1.0) THEN
!  For grid boxes containing sea:
              STR_RESIST_B(i,j) = RESB_ocffInCloud*RESIST_B(i,j)
              STR_RESIST_S(i,j) = RESS_ocff*RESIST_S(K)
!
!  Calculate RES_FACTOR for cloud ocff  (only correct for sea points)
!
              RES_FACTOR(i,j) = aresist(i,j) /                          &
                 (aresist(i,j)+STR_RESIST_B(i,j)+STR_RESIST_S(i,j))
!
            END IF
!
          End Do
        End Do
!
! Call SRESFACT for cloud ocff to calculate correct RES_FACTOR_LAND
! values if land present:
        IF (LAND_POINTS.GT.0) THEN
!
          Call SRESFACT (ROW_LENGTH, ROWS, LAND_POINTS, LAND_INDEX,     &
                     NTILES, TILE_INDEX, TILE_PTS, .FALSE.,             &
                     CANOPY, CATCH,                                     &
                     GC, TILE_FRAC,                                     &
                     SNOW_TILE,                                         &
                     ARESIST, ARESIST_TILE, RESIST_B_TILE,              &
                     RESB_ocffInCloud, RESS_ocff, RES_FACTOR_LAND)
!
! Recalculate RES_FACTOR for cloud ocff combining sea and land values
!
          Do j=1,rows
            Do i=1,row_length
              RES_FACTOR(i,j) = (1.0-FLANDG(i,j))*RES_FACTOR(i,j) +     &
                               FLANDG(i,j)*RES_FACTOR_LAND(i,j)
            End Do
          End Do
!
        END IF
!
! Mix cloud ocff
!
        Call tr_mix(                                                    &
             halo_i, halo_j, row_length, rows, bl_levels                &
            ,off_x, off_y, alpha_cd                                     &
            ,rhokh_mix(1,1,2), rho_aresist                              &
            ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
            ,timestep                                                   &
            ,tr_flux, ocff_cld( 1:row_length, 1:rows, 1:bl_levels )     &
            ,zero_field, RES_FACTOR, drydep_str                         &
            ,kent, we_lim, t_frac, zrzi                                 &
            ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
            ,zh ,zhsc, z_half                                           &
            ,ErrorStatus, .false.                                       &
             )
!
        If (ErrorStatus /=0) Then
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If
!
! Write diagnostics to STASH
!
        item = 409                 !MSR dry (occult) dep flux cloud ocff
        IF ( sf(item,sect) ) THEN
!
!  Change sign of dry dep flux (otherwise negative)
          Do j=1,rows
            Do i=1,row_length
              work_2d(i,j) = -drydep_str(i,j)
            End Do
          End Do
!
          Call copydiag(STASHwork(si(item,sect,im_index)),work_2d,      &
              row_length,rows,0,0,0,0,at_extremity,                     &
              atmos_im,sect,item,                                       &
              ErrorStatus,Cmessage)
!
          If (ErrorStatus /=0) Then
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
!
        END IF  ! SF(item,sect)
!
      END IF  ! L_ocff ocff modelling included
!
! END of fossil-fuel organic carbon section
!
!---------------------------------------------------------------------
! Mixing for Cariolle Ozone tracer
!---------------------------------------------------------------------
      If ( L_USE_CARIOLLE ) Then
! DEPENDS ON: tr_mix
        Call tr_mix(                                                   &
            halo_i, halo_j, row_length, rows, bl_levels                &
           ,off_x, off_y, alpha_cd                                     &
           ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                         &
           ,dtrdz_charney_grid, r_rho_levels, r_theta_levels           &
           ,timestep                                                   &
           ,tr_flux, OZONE_TRACER( 1:row_length, 1:rows, 1:bl_levels ) &
           ,zero_field, zero_field, drydep_str                         &
           ,kent, we_lim, t_frac, zrzi                                 &
           ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
           ,zh ,zhsc, z_half                                           &
           ,ErrorStatus, .false.                                       &
            )

        If (ErrorStatus /= 0) Then
          Cmessage = 'Problem with mixing aerosol'
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

        item = 480
        If ( sf(item,sect) ) Then      ! diagnostic flux required
! DEPENDS ON: copydiag_3d
          Call copydiag_3d( stashwork(si(item,sect,im_index)),         &
             tr_flux,                                                  &
             row_length,rows,bl_levels,0,0,0,0, at_extremity,          &
             stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
             stash_levels,num_stash_levels+1,                          &
             atmos_im,sect,item,                                       &
             ErrorStatus ,cmessage)
        End If

        If (ErrorStatus /=0) Then
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

      End If

!
!---------------------------------------------------------------------
! Mixing for Free Tracers
! Will mimic 4.5 here, but note that free tracer levels exist only
! at the top of the model so this may be incorrect if
! model_levels /= tr_levels!
!---------------------------------------------------------------------

      If ( l_bl_tracer_mix ) Then

        If (tr_levels /= model_levels) Then
          Write(Cmessage,*) 'Tracer mixing gives unexpected results',   &
     &                      'when model_levels /= tr_levels'
          ErrorStatus = -10
! DEPENDS ON: ereport
          Call Ereport( RoutineName, ErrorStatus, Cmessage )
        End If

       !-----------------------------
        !Adding in tracer fluxes to tracers - kdcorbin, 05/10
        !-----------------------------

        if (tr_vars .gt. flux_vars) then
           write(6,*)'ERROR IN BL_TRMIX_DD:  Too many tracers'
           write(6,*)'  Only using ',flux_vars,' fluxes'
           n_tr=flux_vars
        else
           n_tr=tr_vars
        endif
     
! rml 9/9/11 put land carbon fluxes into tracer 1
        if (L_CO2_TRACER.and.n_tr.ge.1) then
          DO L=1,LAND_POINTS
            J = (LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I =  LAND_INDEX(L) - (J-1)*ROW_LENGTH
            TRACER_FLUX1(I,J)= (RESP_S(L) - NPP(L)) * M_CO2 / M_CARBON
          ENDDO
        endif

! rml 1/7/13 put ocean carbon fluxes into tracer 2 (if running coupled model)
        if (L_CO2_TRACER.and.n_tr.ge.2) then
          DO I=1,ROW_LENGTH
            DO J=1,ROWS

!  (ii) CO2 flux from ocean. (+ve implies air to sea)
!! RML : WILL NEED TO CHECK SIGN CONVENTION
               IF ( CO2FLUX(I,J)  /=  RMDI ) THEN
                 TRACER_FLUX2(I,J)= -1. * CO2FLUX(I,J)
               ENDIF

            ENDDO
          ENDDO
        endif

! rml 1/7/13 if co2 (fossil) emissions are switched on, write into tracer 3
        if (L_CO2_TRACER.and.L_CO2_EMITS.and.n_tr.ge.3) then
          DO I=1,ROW_LENGTH
            DO J=1,ROWS
              IF ( CO2_EMITS(I,J)  /=  RMDI ) THEN
                TRACER_FLUX3(I,J) = CO2_EMITS(I,J)
              ENDIF
            ENDDO
          ENDDO
        endif

        Do n_tracer = 1, n_tr

            mytracer_flux = 0.0

          Select Case (n_tracer)

           Case (1)
                !Tracer flux 1
                mytracer_flux(:,:) = TRACER_FLUX1(:,:)

           Case (2)
                !Tracer flux 2
                mytracer_flux(:,:) = TRACER_FLUX2(:,:)

           Case (3)
                !Tracer flux 3
                mytracer_flux(:,:) = TRACER_FLUX3(:,:)
    
           Case (4)
                !Tracer flux 4
                mytracer_flux(:,:) = TRACER_FLUX4(:,:)
  
           Case (5)
                !Tracer flux 5
                mytracer_flux(:,:) = TRACER_FLUX5(:,:)

           Case (6)
                !Tracer flux 6
                mytracer_flux(:,:) = TRACER_FLUX6(:,:)

           Case (7)
                !Tracer flux 7
                mytracer_flux(:,:) = TRACER_FLUX7(:,:)
        
           Case (8)
                !Tracer flux 8
                mytracer_flux(:,:) = TRACER_FLUX8(:,:)

           Case (9)
                !Tracer flux 9
                mytracer_flux(:,:) = TRACER_FLUX9(:,:)

           Case (10)
                !Tracer flux 10
                mytracer_flux(:,:) = TRACER_FLUX10(:,:)

           Case (11)
                !Tracer flux 11
                mytracer_flux(:,:) = TRACER_FLUX11(:,:)

           Case (12)
                !Tracer flux 12
                mytracer_flux(:,:) = TRACER_FLUX12(:,:)

           Case (13)
                !Tracer flux 13
                mytracer_flux(:,:) = TRACER_FLUX13(:,:)

           Case (14)
                !Tracer flux 14
                mytracer_flux(:,:) = TRACER_FLUX14(:,:)

           Case (15)
                !Tracer flux 15
                mytracer_flux(:,:) = TRACER_FLUX15(:,:)
         
           Case (16)
                !Tracer flux 16
                mytracer_flux(:,:) = TRACER_FLUX16(:,:)

           Case (17)
                !Tracer flux 17
                mytracer_flux(:,:) = TRACER_FLUX17(:,:)

           Case (18)
                !Tracer flux 18
                mytracer_flux(:,:) = TRACER_FLUX18(:,:)

           Case (19)
                !Tracer flux 19
                mytracer_flux(:,:) = TRACER_FLUX19(:,:)

           Case (20)
                !Tracer flux 20
                mytracer_flux(:,:) = TRACER_FLUX20(:,:)

           End Select

        !kdcorbin, 05/10 - added tracer_flux to tr_mix call
! DEPENDS ON: tr_mix
          Call tr_mix(                                                  &
     &         halo_i, halo_j, row_length, rows, bl_levels              &
     &        ,off_x, off_y, alpha_cd                                   &
     &        ,rhokh_mix(1,1,2), rhokh_mix(1,1,1)                       &
     &        ,dtrdz_charney_grid, r_rho_levels, r_theta_levels         &
     &        ,timestep                                                 &
     &        ,tr_flux, free_tracers( 1:row_length, 1:rows,             &
     &                                1:bl_levels,  n_tracer )          &
     &        ,mytracer_flux, zero_field, drydep_str                       &
     &        ,kent, we_lim, t_frac, zrzi                               &
     &        ,kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc               &
     &        ,zh ,zhsc, z_half                                         &
     &        ,ErrorStatus, .false.                                     &
     &         )

          If (ErrorStatus /= 0) Then
            Cmessage = 'Problem with mixing tracers'
! DEPENDS ON: ereport
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

!kdcorbin, 05/10 - commenting out copy to tracer fluxes
!          item = 99 + n_tracer
!          If ( sf(item,sect) )  Then     ! diagnostic flux required
!! DEPENDS ON: copydiag_3d
!            Call copydiag_3d( stashwork(si(item,sect,im_index)),        &
!     &          tr_flux,                                                &
!     &          row_length,rows,bl_levels,0,0,0,0, at_extremity,        &
!     &          stlist(1,stindex(1,item,sect,im_index)),len_stlist,     &
!     &          stash_levels,num_stash_levels+1,                        &
!     &          atmos_im,sect,item,                                     &
!     &          ErrorStatus, cmessage)
!          End If

!          If (ErrorStatus /=0) Then
!! DEPENDS ON: ereport
!            Call Ereport( RoutineName, ErrorStatus, Cmessage )
!          End If

        End Do
      End If

      Return
      End Subroutine bl_trmix_dd
#endif
