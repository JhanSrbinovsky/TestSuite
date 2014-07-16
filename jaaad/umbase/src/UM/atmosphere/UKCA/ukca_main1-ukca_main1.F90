#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!  To gather required fields from D1, then call the UK Chemistry
!  and Aerosols submodel components.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Method:
! 1) Interface with atmosphere model:
!   Prognostic fields are selected from D1 using the item codes
!   defined in the routine UKCA_SetD1Defs, these include tracers
!   from section 34 (UKCA).
!   Diagnostic fields are selected from D1 using a tag of 98, these
!   need to be set in the STASH panel of the UMUI.
! 2) UKCA routines are called depending on the scheme selected:
!    - Woodward dust scheme
!    - Emissions routine
!    - photolysis routine
!    - Chemistry control routine
! 3) Updated tracer arrays are returned to the D1 array
!
! CONTAINED subroutines:  GETD1FLDS, PUTD1FLDS
!
!  Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                            Fiona O'Connor
!
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_MAIN1(                                           &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "arglndm.h"
#include "argppx.h"
#include "argptra.h"
         Idummy)

      USE ASAD_MOD,          ONLY: asad_mod_final, method
      USE UKCA_D1_DEFS
      USE UKCA_TRACER_VARS
      USE UKCA_CSPECIES
      USE UKCA_CONSTANTS
      USE UKCA_TROPOPAUSE
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE UKCA_CHEM1_DAT
      USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels
      USE trignometric_mod,  ONLY: true_longitude, sin_theta_longitude, &
                                   sin_theta_latitude,                  &
                                   cos_v_latitude,                      &
                                   FV_cos_theta_latitude
      USE dyn_coriolis_mod,  ONLY: f3_at_u
      IMPLICIT NONE

#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "typlndm.h"
#include "typptra.h"
#include "typbnd.h"
#include "typsts.h"
#include "typcona.h"
#include "nstypes.h"
! River routing
#include "typatcpl.h"
#include "chsunits.h"
#include "ccontrol.h"
#include "cruntimc.h"
#include "problem.h"
#include "surface.h"
#include "vgrid.h"
#include "o3intp.h"
#include "ctime.h"
#include "stparam.h"
#include "c_global.h"
#include "c_writd.h"
#include "ppxlook.h"
#include "cprintst.h"
#include "c_dust_ndiv.h"
#include "ctfilt.h"
! #include "model.h"
! #include "cstash.h"

      INTEGER, INTENT(INOUT)     :: Idummy       !

! ErrorStatus

      INTEGER                    :: errcode=0     ! Error flag (0 = OK)
      CHARACTER(72)              :: cmessage      ! Error return message

! Local scalars

      INTEGER    :: Nfields           ! fields counter
      INTEGER    :: Nukca_Fields      ! No. fields requested
      INTEGER    :: len               ! local dimension
      INTEGER    :: levs              ! number of levels
      INTEGER    :: StashCode         ! stash code
      INTEGER    :: section           ! stash section
      INTEGER    :: item              ! stash item
      INTEGER    :: addr              ! stash item
      INTEGER    :: field_typ         ! Field type
      INTEGER    :: halo_typ          ! halo type
      INTEGER    :: grid_typ          ! grid type
      INTEGER    :: tag               ! stash tag
      INTEGER    :: ptd1              ! D1 pointer
      INTEGER    :: I,I1,I2           ! loop variables
      INTEGER    :: J,J1,J2           ! loop variables
      INTEGER    :: K,L,N             ! loop variables
      INTEGER    :: m_atm_modl        ! sub model
      INTEGER    :: GET_FLD_TYPE      ! UM function
      INTEGER    :: timestep_number   ! no. of atmos timesteps since bas
      INTEGER    :: jradon_em         ! alternative D1 pointer for rn em
      INTEGER    :: nd_o3             ! size of um ozone array
      INTEGER    :: index             ! Indicies for controlling
      INTEGER    :: index2            ! albedo
      INTEGER    :: index3            ! calculation
      INTEGER    :: A_Steps_per_hr    ! Atmos steps per hour
      INTEGER    :: Radsteps_per_day  ! Radiation TS per day
      INTEGER    :: icnt              ! counter


      REAL       :: timestep                 ! atmosphere model TS
      REAL       :: secs_day=3600.0*24.0     ! Secs per day
      REAL       :: min_surf_albedo=0.1      ! fill-in value
      REAL       :: fx                       ! interpolation fraction
      REAL       :: eq_time                  ! Equation of time
      REAL       :: Sindec                   ! Solar declination
      REAL       :: secondssincemidnight     ! day time
      REAL       :: SCS                      ! solar constant scaling factor

      LOGICAL, SAVE :: first       = .TRUE.    ! true only on 1st call
      LOGICAL       :: L_moses_ii                 ! T if MOSES_II in use

! Local arrays

! MPP-related Arrays
! Table of number of points on a row
      INTEGER, DIMENSION(0:nproc-1) :: g_row_length
! Table number of rows in theta field
      INTEGER, DIMENSION(0:nproc-1) :: g_rows

! Tile index etc
      INTEGER, DIMENSION(land_points,ntype) :: tile_index !
      INTEGER, DIMENSION(ntype)             :: tile_pts   ! No of tile p

! arrays filled from D1, allocated dynamically with correct halo read
!  in from addressing array
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: conv_cloud_base
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: conv_cloud_top
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: kent
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: kent_dsc

      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: land_sea_mask

      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: all_tracers
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: trmol_post_atmstep
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: chem_diags
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: BE_fluxdiags
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: mode_diags
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: strat_fluxdiags
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: all_emissions
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: T_bias
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: Q_bias
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: aircraftems
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: SO2emiss_3D
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: tracer1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: chem_diag1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: BE_fluxdiag1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: strat_fluxdiag1
      REAL, DIMENSION(:,:),   ALLOCATABLE :: emission1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: theta
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_rho_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_theta_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: tnd
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: exner_theta_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: exner_rho_levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: pv_on_theta_mlevs
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: q
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rel_humid_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rho_r2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: qcl
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: qcf
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_ppn_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_liq_water
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cloud_ice_content
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_rain3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_snow3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ls_ppn3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_rain3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_snow3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_ppn3d
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: conv_cloud_amount
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: area_cloud_fraction
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: so4_aitken
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: so4_accum
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_wet_h2o2
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_wet_o3
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_dry_oh
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_drydep
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: delso2_wetdep
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: um_ozone   ! i.e. from D1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: um_ozone3d ! after O3_to_3D
      REAL, DIMENSION(:),     ALLOCATABLE :: um_ozone1d ! ip to O3_to_3D
      REAL, DIMENSION(:,:), ALLOCATABLE   :: pstar
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Tstar
      REAL, DIMENSION(:,:), ALLOCATABLE   :: u_10m      ! 10m wind U
      REAL, DIMENSION(:,:), ALLOCATABLE   :: v_10m      ! 10m wind V
      REAL, DIMENSION(:,:), ALLOCATABLE   :: U_scalar_10m
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Soil_clay
      REAL, DIMENSION(:,:), ALLOCATABLE   :: dust_mrel1
      REAL, DIMENSION(:,:), ALLOCATABLE   :: dust_mrel2
      REAL, DIMENSION(:,:), ALLOCATABLE   :: dust_mrel3
      REAL, DIMENSION(:,:), ALLOCATABLE   :: dust_mrel4
      REAL, DIMENSION(:,:), ALLOCATABLE   :: dust_mrel5
      REAL, DIMENSION(:,:), ALLOCATABLE   :: dust_mrel6
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: dust_ustar
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Soil_Layer_Moisture
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Tile_Frac
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Tstar_tile
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Snow_tile
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Frac_types
      REAL, DIMENSION(:,:), ALLOCATABLE   :: Rough_length
! surf_albedo is interpolated at every timestep from land_albedo_all
! which is calculated on radiation timesteps using the SW fluxes.
      REAL, DIMENSION(:,:), ALLOCATABLE   :: land_albedo
      REAL, DIMENSION(:,:,:), ALLOCATABLE, SAVE :: land_albedo_all  ! all r
      REAL, DIMENSION(:,:), ALLOCATABLE   :: surf_albedo  ! interpolated
      REAL, DIMENSION(:,:), ALLOCATABLE   :: net_surf_SW
      REAL, DIMENSION(:,:), ALLOCATABLE   :: tot_surf_SW
      REAL, DIMENSION(:), ALLOCATABLE     :: Fland
      REAL, DIMENSION(:,:), ALLOCATABLE   :: climoz2d !climatological O3
      REAL, DIMENSION(:,:), ALLOCATABLE   :: zbl      ! BL height
      REAL, DIMENSION(:,:), ALLOCATABLE   :: surf_hf
      REAL, DIMENSION(:,:), ALLOCATABLE   :: seaice_frac
      REAL, DIMENSION(:,:), ALLOCATABLE   :: conv_cloud_lwp
      REAL, DIMENSION(:,:), ALLOCATABLE   :: u_s     ! surf frict velocity
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: stcon       ! stomatal cond
      REAL, DIMENSION(:,:), ALLOCATABLE   :: laift_lp    ! LAI
      REAL, DIMENSION(:,:), ALLOCATABLE   :: canhtft_lp  ! canopy height
      REAL, DIMENSION(:,:), ALLOCATABLE   :: z0tile_lp
      REAL, DIMENSION(:,:), ALLOCATABLE   :: canwctile_lp ! canopy WC
      REAL, DIMENSION(:,:), ALLOCATABLE   :: ch4_wetl_emiss
      REAL, DIMENSION(:,:), ALLOCATABLE   :: theta_latitude ! gridbox lat
      REAL, DIMENSION(:,:), ALLOCATABLE   :: sinv_latitude  ! sin(boundary lat)
      REAL, DIMENSION(:,:), ALLOCATABLE   :: tropopause_height
      REAL, DIMENSION(:,:), ALLOCATABLE   :: cos_zenith_angle
      REAL, DIMENSION(:,:), ALLOCATABLE   :: ml_depth
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rhokh_mix
      REAL, DIMENSION(:,:), ALLOCATABLE   :: rho_aresist
      REAL, DIMENSION(:,:), ALLOCATABLE   :: aresist
      REAL, DIMENSION(:,:), ALLOCATABLE   :: resist_b
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: dtrdz_charney_grid
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: we_lim
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: t_frac
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: zrzi
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: we_lim_dsc
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: t_frac_dsc
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: zrzi_dsc
      REAL, DIMENSION(:,:), ALLOCATABLE   :: zhsc
      REAL, DIMENSION(:,:), ALLOCATABLE   :: rb_dust_div1
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rb_dust_ndivs
! solid angle of grid cell saved for mass calculation
      REAL, SAVE, ALLOCATABLE :: mass (:,:,:)
      REAL, SAVE, ALLOCATABLE :: solid_angle(:,:)
      REAL, SAVE, ALLOCATABLE :: z_half(:,:,:)
      REAL, SAVE, ALLOCATABLE :: volume(:,:,:)   ! gridbox volume
      REAL, SAVE, ALLOCATABLE :: area(:,:,:)     ! gridbox area

! Local dynamic arrays
      REAL, DIMENSION (:), ALLOCATABLE    :: STASHwork34

! Derived data arrays, allocatable if necessary
      REAL, DIMENSION(land_points)           :: Pstar_Land     ! surface
      REAL, DIMENSION(land_points)           :: Clay_Land      ! soil cl
      REAL, DIMENSION(land_points)           :: Rhostar_land   ! surface
      REAL, DIMENSION(land_points,ndiv)      :: Mrel_Land      ! relativ
      REAL, DIMENSION(land_points,ntiles)    :: ustar_tile     ! ustar o
      REAL, DIMENSION(row_length,rows)       :: land_fraction  !
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: t_theta_levels ! temp on
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: Thick_bl_levels    ! val
!                                               thickness in metres of t
!                                               distance from boundaries
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: p_layer_boundaries
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: dj     ! photolysis rates
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: so4_sa ! aerosol surface area
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: T_chem
      REAL, DIMENSION(:,:,:), ALLOCATABLE    :: Q_chem

      REAL :: delta_latitude
      REAL :: delta_longitude

! budget variables for stratospheric chemistry
      INTEGER :: n_opl = 0
      INTEGER, PARAMETER :: nphotf = 2
      INTEGER, PARAMETER :: nbimol = 5
      INTEGER, PARAMETER :: nabund = 2
      REAL,ALLOCATABLE :: semiss (:,:,:)
      REAL,ALLOCATABLE :: phot(:,:,:,:)
      REAL,ALLOCATABLE :: bimol(:,:,:,:)
      REAL,ALLOCATABLE :: termol(:,:,:)
      REAL,ALLOCATABLE :: hetero(:,:,:)
      REAL,ALLOCATABLE :: drydep(:,:,:)
      REAL,ALLOCATABLE :: wetdep(:,:,:)
      REAL,ALLOCATABLE :: noxlight(:,:,:)
      REAL,ALLOCATABLE :: abundance(:,:,:,:)
      REAL, ALLOCATABLE :: ozonebud(:,:,:,:)

      INTEGER, SAVE :: fastj_levels ! number of levels on which to do FAST-J

      REAL, SAVE :: lambda_aitken, lambda_accum ! parameters for computation
                                                ! of surface area density
      INTEGER, SAVE :: interval     ! interval in timesteps between calls to
                                    ! chemistry
      LOGICAL :: do_chemistry


! Data needed to prescribe lower boundary conditions for selected species
! needed in stratospheric chemistry
      INTEGER, PARAMETER :: n_boundary_vals = 29

      CHARACTER*10, SAVE :: lbc_spec(n_boundary_vals) =                 &
         (/'N2O       ','CF2Cl2    ','CFCl3     ','MeBr      ',         &
           'HCl       ','HOCl      ','Clx       ','OClO      ',         &
           'ClONO2    ','HBr       ','HOBr      ','BrCl      ',         &
           'Brx       ','BrONO2    ','CF2ClCFCl2','CHF2Cl    ',         &
           'CCl4      ','MeCl      ','MeCCl3    ','CF2ClBr   ',         &
           'CF3Br     ','H2        ','COS       ','Cl        ',         &
           'ClO       ','Cl2O2     ','Br        ','BrO       ',         &
           'CH4       ' /)

! Mixing ratios for longlived tracers (from WMO, 2002)
! Short-lived stratospheric tracers are forced to 0 at the surface.
      REAL, PARAMETER :: my_n2o_mmr = 4.7842E-7
      REAL, PARAMETER :: f12_mmr    = 2.2145E-9
      REAL, PARAMETER :: f11_mmr    = 1.3294E-9
      REAL, PARAMETER :: mebr_mmr   = 39.1206E-12
      REAL, PARAMETER :: hcl_mmr    = 1.5125E-10
      REAL, PARAMETER :: f113_mmr   = 5.3072E-10
      REAL, PARAMETER :: h22_mmr    = 4.1801E-10
      REAL, PARAMETER :: ccl4_mmr   = 5.3158E-10
      REAL, PARAMETER :: mecl_mmr   = 9.4133E-10
      REAL, PARAMETER :: meccl3_mmr = 2.3041E-10
      REAL, PARAMETER :: h1211_mmr  = 2.83858e-11
      REAL, PARAMETER :: h1301_mmr  = 1.59722e-11
      REAL, PARAMETER :: h2_mmr     =  3.4520E-8
      REAL, PARAMETER :: cos_mmr    = 1.0356E-9
      REAL, PARAMETER :: ch4_mmr    = 9.8530E-7

!  MMRs for named species
      REAL, SAVE :: lbc_mmr(n_boundary_vals) =                          &
        (/my_n2o_mmr, f12_mmr , f11_mmr  , mebr_mmr ,                   &
             hcl_mmr,       0.,      0.  ,        0.,                   &
                  0.,       0.,      0.  ,        0.,                   &
                  0.,       0., f113_mmr , h22_mmr  ,                   &
            ccl4_mmr, mecl_mmr,meccl3_mmr, h1211_mmr,                   &
           h1301_mmr, h2_mmr  ,cos_mmr   ,       0.0,                   &
           0.0      , 0.0     ,0.0       ,       0.0,                   &
           ch4_mmr /)

!- End of header

! Common Blocks

        INTERFACE UKCA_FASTJ
! DEPENDS ON: UKCA_FASTJ
          SUBROUTINE UKCA_FASTJ(                                       &
#include "arglndm.h"
            dj,                                                        &
            p_rho_levels, p0,                                          &
            p_theta_levels,                                            &
            t, tstar,                                                  &
            so4_aitken, so4_accum,                                     &
            q, qcl, qcf,                                               &
            conv_cloud_lwp, conv_cloud_top, conv_cloud_base,           &
            conv_cloud_amount,                                         &
            surf_alb,                                                  &
            nd_o3, um_ozone3d,                                         &
            land_frac_ctile,                                           &
            z_top_of_model,                                            &
            l_moses_ii)

            USE FASTJ_DATA,  ONLY: fastj_set_limits,                   &
                                fastj_allocate_memory,                 &
                                Blocking,kpcx,jpcl,                    &
                                jjpnl,jppj,nsl,sa,p,rz_3d,             &
                                tau,month,iday,nslon,nslat,            &
                                sza,sza_2d,SZAFAC,SZAFAC_2d,u0,        &
                                od,odw,odi,sulphur,fj_ozone,           &
                                t_fastj=>t

            USE FASTJ_MIE,   ONLY: NWB3,NWB2


#include "cmaxsize.h"
#include "nstypes.h"
#include "parvars.h"
#include "typsize.h"
#include "typcona.h"
#include "typlndm.h"
#include "chsunits.h"
#include "csubmodl.h"
#include "ccontrol.h"
#include "aercmp3a.h"
#include "aerprm3a.h"
#include "cphyscon.h"
#include "ctime.h"
!     DECLARE THE REDUCED SW SPECTRAL FILE AND ITS HOLDING COMMON BLOCK.
#include "mxsize3a.h"
#include "swspdl3a.h"
#include "swspcm3a.h"
#include "c_v_m.h"

            INTEGER, INTENT(in) :: nd_o3
            INTEGER, DIMENSION(row_length, rows),                      &
                     INTENT(in) :: conv_cloud_top
            INTEGER, DIMENSION(row_length, rows),                      &
                     INTENT(in) :: conv_cloud_base

            LOGICAL, INTENT(in) :: l_moses_ii

            REAL, INTENT(in) :: z_top_of_model
            REAL, DIMENSION(row_length, rows, model_levels, jppj),     &
                  INTENT(inout) :: dj
            REAL, DIMENSION(row_length, rows, model_levels+1),         &
                  INTENT(in) :: p_rho_levels
            REAL, DIMENSION(row_length, rows),                         &
                  INTENT(in) :: p0
            REAL, DIMENSION(row_length, rows, model_levels),           &
                  INTENT(in) :: p_theta_levels
            REAL, DIMENSION(row_length, rows, model_levels),           &
                  INTENT(in) :: t
            REAL, DIMENSION(row_length, rows),                         &
                  INTENT(in) :: tstar
            REAL, DIMENSION(row_length, rows, tr_levels),              &
                  INTENT(in) :: so4_aitken
            REAL, DIMENSION(row_length, rows, tr_levels),              &
                  INTENT(in) :: so4_accum
            REAL, DIMENSION(row_length, rows, wet_levels),             &
                  INTENT(in) :: q
            REAL, DIMENSION(row_length, rows, wet_levels),             &
                  INTENT(in) :: qcl
            REAL, DIMENSION(row_length, rows, wet_levels),             &
                  INTENT(in) :: qcf
            REAL, DIMENSION(row_length, rows),                         &
                  INTENT(in) :: conv_cloud_lwp
            REAL, DIMENSION(row_length, rows, wet_levels),             &
                  INTENT(in) :: conv_cloud_amount
            REAL, DIMENSION(row_length, rows),                         &
                  INTENT(in) :: surf_alb
            REAL, DIMENSION(row_length,rows,ozone_levels),             &
                  INTENT(in) :: um_ozone3d
            REAL, DIMENSION(land_points),                              &
                  INTENT(in) :: land_frac_ctile

          END SUBROUTINE UKCA_FASTJ
        END INTERFACE UKCA_FASTJ

! Initialisation

! timestep information

      m_atm_modl = SUBMODEL_FOR_SM(a_im)   ! submodel code
      timestep = SECS_PER_STEPim(atmos_im) ! timestep in seconds
      timestep_number = STEPim(atmos_im)   ! no. of steps since basis
      A_Steps_per_hr = 3600*STEPS_PER_PERIODim(a_im)/                  &
                            SECS_PER_PERIODim(a_im)
      Radsteps_per_day=A_Steps_per_hr*24/A_SW_Radstep
      call_chem_freq = 1            ! no of chemical steps per timestep


! Set mpp arrays from parvars.h information

      DO i= 0,nproc-1
        g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
        g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
      ENDDO ! i processors

! Allocate diagnostic space for STASH
      ALLOCATE (STASHwork34(STASH_maxlen(34,A_im)))

      IF (first) THEN
        WRITE(6,*) 'First call to ukca_main'
        WRITE(6,*) 'timestep = ',timestep
        WRITE(6,*) 'timestep number = ',timestep_number
        WRITE(6,*) 'call chem freq = ',CALL_CHEM_FREQ

! Define D1 definitions, emissions and diagnostics required
! DEPENDS ON: ukca_setd1defs
        CALL UKCA_SETD1DEFS(row_length,rows,n_rows, model_levels,      &
                 bl_levels,tr_levels,wet_levels,land_points,sm_levels, &
                 ntiles,tr_ukca)

! Initialise chemical definition arrays
        CALL UKCA_CHEM1_INIT()

! Define UKCA-mode setup based on value of i_mode_setup from UMUI
! below are called from modules UKCA_SETUP_INDICES and UKCA_MODE_SETUP
        IF (L_ukca_mode) THEN
          CALL UKCA_MODE_IMSCAVCOFF
          IF(i_mode_setup == 1) THEN
           CALL UKCA_SV1_SUSS_4MODE
           CALL UKCA_MODE_ALLCP_4MODE
          ELSE IF(i_mode_setup == 2) THEN
            CALL UKCA_ORGV1_SUSSBCOC_5MODE
            CALL UKCA_MODE_SUSSBCOC_5MODE
          ENDIF
        ENDIF
      ENDIF ! first

! Loop through all the objects in D1 to find address and characteristics
!  of the requested items
      IF (first) THEN
        WRITE(6,*) 'UKCA Prognostics and Diagnostics from D1:'
        WRITE(6,*) 'section,item,levels,length,address,halo_type,',    &
                 'grid_type,field_type'
        DO i=1,no_obj_D1(m_atm_modl)
          section = d1_addr(D1_section,i,m_atm_modl)
          item = d1_addr(D1_item,i,m_atm_modl)
          levs = d1_addr(d1_no_levels,i,m_atm_modl)
          len = d1_addr(d1_length,i,m_atm_modl)
          addr    = d1_addr(d1_address,i,m_atm_modl)
          halo_typ = d1_addr(d1_halo_type,i,m_atm_modl)
          grid_typ = d1_addr(d1_grid_type,i,m_atm_modl)
! DEPENDS ON: get_fld_type
          field_typ = get_fld_type(grid_typ)
          DO J=1,Nukca_D1items
            IF (UkcaD1Codes(J)%section == section .AND.                &
                UkcaD1Codes(J)%item == item .AND.                      &
                d1_addr(d1_object_type,i,m_atm_modl) == prognostic     &
                .AND. UkcaD1codes(J)%Prognostic) THEN
              UkcaD1Codes(J)%n_levels=levs
              UkcaD1Codes(J)%address=addr
              UkcaD1Codes(J)%length=len
              UkcaD1Codes(J)%halo_type=halo_typ
              UkcaD1Codes(J)%grid_type=grid_typ
              UkcaD1Codes(J)%field_type=field_typ
              WRITE(6,*) 'P',J, section,item,levs,len,addr,            &
                             halo_typ,grid_typ,field_typ
            ENDIF
          ENDDO
        ENDDO

!       Diagnostics loop through all stashlist items, check
!       for tag value, read adress etc into UkcaD1Codes.

        DO L=1,totitems
          tag=STLIST(st_macrotag,L)-(1000*(STLIST(st_macrotag,L)/1000))
          IF (tag == 98) THEN
            ptd1      = STLIST(st_D1pos,L)
            section   = d1_addr(d1_section,  ptd1,m_atm_modl)
            item      = d1_addr(d1_item,     ptd1,m_atm_modl)
            levs      = d1_addr(d1_no_levels,ptd1,m_atm_modl)
            len       = d1_addr(d1_length,   ptd1,m_atm_modl)
            addr      = d1_addr(d1_address,  ptd1,m_atm_modl)
            halo_typ  = d1_addr(d1_halo_type,ptd1,m_atm_modl)
            grid_typ  = d1_addr(d1_grid_type,ptd1,m_atm_modl)
! DEPENDS ON: get_fld_type
            field_typ = get_fld_type(grid_typ)
            stashcode = section*1000 + item
            DO J=1,Nukca_D1items
              IF (UkcaD1Codes(J)%section == section .AND.              &
                  UkcaD1Codes(J)%item == item .AND.                    &
                  .NOT. UkcaD1Codes(J)%Prognostic) THEN
                UkcaD1Codes(J)%n_levels  = levs
                UkcaD1Codes(J)%address   = addr
                UkcaD1Codes(J)%length    = len
                UkcaD1Codes(J)%halo_type = halo_typ
                UkcaD1Codes(J)%grid_type = grid_typ
                UkcaD1Codes(J)%field_type= field_typ
                WRITE(6,*) 'D',section,item,levs,len,addr,             &
                            halo_typ,grid_typ,field_typ
              ENDIF
            ENDDO
          ENDIF  ! tag
        ENDDO  ! L

! set constants needed in sulphur chemistry
      IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
         (L_ukca_stratcfc))  THEN
          lambda_aitken = 3./rho_so4 / chi /                            &
                         rad_ait * EXP(-2.5 * (ALOG(sigma)**2))         &
                         * 0.01 ! revert from m^-1 to cm^2/cm^3
          lambda_accum  = 3./rho_so4 / chi /                            &
                         rad_acc * EXP(-2.5 * (ALOG(sigma)**2))         &
                         * 0.01 ! revert from m^-2 to cm^2/cm^3
      ENDIF

! set up timestep counting. Interval depends on solver: IMPACT and Rosenbrock
! are run every dynamical timestep. Newton-Raphson is run every 2 / 3
! timesteps for a 30/20 minutes dynamical timestep.

        IF (method == 3) THEN
          IF (timestep == 1200.) THEN
            interval = 3
          ELSE
            interval = 2
          END IF
        ELSE
          interval = 1
        END IF

      ENDIF  ! first

! decide whether to do chemistry
      do_chemistry = (MOD(timestep_number, interval) == 0)


! Check if all items selected have been identified in D1

      DO i=1,Nukca_D1items
        IF (UkcaD1codes(i)%address == IMDI .AND.                       &
            UkcaD1codes(i)%required) THEN
          cmessage=' Item address not found in D1 array'
          WRITE(6,*) cmessage
! DEPENDS ON: ereport
          CALL EREPORT('UKCA_MAIN',                                    &
          UkcaD1Codes(i)%section*1000+UkcaD1Codes(I)%item,cmessage)
        ENDIF
      ENDDO

!     Copy fields from D1 array into named item, allocate arrays

      DO i=1,Nukca_D1items
        IF (UkcaD1Codes(i)%required) THEN
          CALL GETD1FLDS(i)
        ENDIF
      ENDDO

!     Transform tracers to ensure elemental conservation
      IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
         (L_ukca_stratcfc))  THEN
! DEPENDS ON: ukca_transform_halogen
        CALL ukca_transform_halogen(tr_ukca,rows,row_length,            &
                             model_levels, offx, offy, all_tracers,     &
                            halo_i, halo_j, q, .TRUE.,timestep_number)
      ENDIF

!     Calculate any derived variables, SO2 fluxes are required in chemistry_ctl

      i1 = lbound(theta,dim=1)
      i2 = ubound(theta,dim=1)
      j1 = lbound(theta,dim=2)
      j2 = ubound(theta,dim=2)
      ALLOCATE(t_theta_levels(I1:I2,J1:J2,model_levels))
      ALLOCATE(rel_humid_frac(I1:I2,J1:J2,wet_levels))
      ALLOCATE(Thick_bl_levels(1:row_length,1:rows,bl_levels))
      ALLOCATE(tile_frac(land_points,ntiles))              ! see below
      ALLOCATE(p_layer_boundaries(row_length,rows,0:model_levels))
      ALLOCATE(surf_albedo(row_length,rows))
      ALLOCATE(tnd(row_length,rows,model_levels))
      ALLOCATE(ls_ppn3d(row_length,rows,model_levels))
      ALLOCATE(conv_ppn3d(row_length,rows,model_levels))
      ALLOCATE(so4_sa(row_length, rows, model_levels)) ! sulphate area density
      ALLOCATE(delSO2_wet_h2o2(row_length,rows,model_levels))
      ALLOCATE(delSO2_wet_o3(row_length,rows,model_levels))
      ALLOCATE(delSO2_dry_oh(row_length,rows,model_levels))
      ALLOCATE(delSO2_drydep(row_length,rows,model_levels))
      ALLOCATE(delSO2_wetdep(row_length,rows,model_levels))

! Initialise fluxes to zero as may not be filled by chemistry
      delSO2_wet_h2o2(:,:,:)=0.0
      delSO2_wet_o3(:,:,:)=0.0
      delSO2_dry_oh(:,:,:)=0.0
      delSO2_drydep(:,:,:)=0.0
      delSO2_wetdep(:,:,:)=0.0

! Required in call to ukca_chemistry_ctl, but may be unallocated
      IF (.NOT. ALLOCATED(BE_fluxdiags)) THEN
        n_be_fluxdiags=1
        ALLOCATE(BE_fluxdiags(row_length,rows,model_levels,             &
                 n_be_fluxdiags))
        BE_fluxdiags=0.0
      ENDIF

! Required in call to UKCA_MODE, but may be unallocated
      IF (.NOT. ALLOCATED(mode_diags)) THEN
        n_mode_diags=1
        ALLOCATE(mode_diags(row_length,rows,model_levels,               &
                 n_mode_diags))
        mode_diags=0.0
      ENDIF

! Required in call to UKCA_EMISSION_CTL, but may be unallocated
      IF (.NOT. ALLOCATED(SO2emiss_3D)) THEN
        ALLOCATE(SO2emiss_3D(row_length,rows,model_levels))
      ENDIF


      IF (first) THEN

        ALLOCATE(theta_latitude(row_length, rows))
        ALLOCATE(sinv_latitude(row_length, 0:rows))
        ALLOCATE(solid_angle(row_length, rows))
        ALLOCATE(volume(row_length,rows,model_levels))
        ALLOCATE(area(row_length,rows,model_levels))
        ALLOCATE(mass(row_length, rows, model_levels))
        ALLOCATE(p_tropopause(row_length,rows))
        ALLOCATE(tropopause_level(row_length,rows))
        ALLOCATE(theta_trop(row_length,rows))
        ALLOCATE(pv_trop(row_length,rows))
        ALLOCATE(L_troposphere(row_length,rows,model_levels))

        theta_latitude  = asin(sin_theta_latitude)
        sinv_latitude(1:row_length,0:rows) =                            & 
     &      sin(acos(cos_v_latitude(1:row_length,0:rows))) 
                 
        solid_angle = delta_lambda *                                    & 
             (sinv_latitude(:,1:rows) - sinv_latitude(:,0:rows-1)) 

        DO k=1,model_levels-1
          volume(1:row_length,1:rows,k) = solid_angle / 3.0 *           &
             (r_theta_levels(1:row_length,1:rows,k)**3 -                &
              r_theta_levels(1:row_length,1:rows,k-1)**3)
        ENDDO

        DO j = 1, rows 
          DO i = 1, row_length 
            volume(i,j,model_levels) = solid_angle(i,j) / 3.0 *        & 
              (r_theta_levels(i,j,model_levels)**3 -                   & 
               r_rho_levels(i,j,model_levels)**3) 
          END DO 
        END DO 

        area(1:row_length,1:rows,model_levels) = 2*pi*                  &
             r_theta_levels(1:row_length,1:rows,model_levels)**2 *      &
             (sinv_latitude(1:row_length,1:rows) -                      &
              sinv_latitude(1:row_length,0:rows-1))/                    &
              global_row_length

        DEALLOCATE(theta_latitude)

        WRITE(6,*) 'Radsteps_per_day: ',Radsteps_per_day,              &
                    a_sw_radstep,a_steps_per_hr              !debug
        ALLOCATE(land_albedo_all(row_length,rows,Radsteps_per_day))
        land_albedo_all=min_surf_albedo

       ALLOCATE(z_half(row_length,rows,bl_levels))
       DO k = 1,bl_levels
         DO j = 1, rows
           DO i= 1, row_length
              z_half(i,j,k) = r_rho_levels(i,j,k) - r_theta_levels(i,j,0)
            ENDDO
          ENDDO
        ENDDO

      ENDIF

      index  = MOD(timestep_number,a_sw_radstep)
      index2 = 1+MOD(timestep_number,a_sw_radstep*Radsteps_per_day)/   &
               a_sw_radstep
      index3 = 1+MOD(index2,8)

      IF (L_SW_Radiate) THEN  ! Update land_albedo_all
        WHERE (tot_surf_sw > 0.0 .and. tot_surf_sw < 1.5e3)
          land_albedo=(tot_surf_sw-net_surf_sw)/tot_surf_sw
        ELSEWHERE
          land_albedo=min_surf_albedo
        ENDWHERE
        land_albedo_all(:,:,index2)=land_albedo(:,:)
        WRITE(6,*) 'land_albedo: ',land_albedo_all(9,9,:)   !debug
      ENDIF

!     Interpolate surf_albedo to timestep value

      fx=(REAL(index)+0.5)/a_sw_radstep
      surf_albedo(:,:)=fx*land_albedo_all(:,:,index2)+                 &
                       (1.0-fx)*land_albedo_all(:,:,index3)
!      WRITE(6,*) 'surf_albedo: ',minval(surf_albedo),                 &
!     &            maxval(surf_albedo)   !debug
!      WRITE(6,*) 'surf_albedo: ',minloc(surf_albedo),                 &
!     &     minval(surf_albedo),maxloc(surf_albedo),maxval(surf_albedo)
      DEALLOCATE(net_surf_sw)
      DEALLOCATE(tot_surf_sw)

      t_theta_levels=exner_theta_levels*theta
      Thick_bl_levels(1:row_length,1:rows,1) = 2.0*                    &
            (r_theta_levels(1:row_length,1:rows,1) -                   &
             r_theta_levels(1:row_length,1:rows,0))
      DO k=2,bl_levels
        Thick_bl_levels(1:row_length,1:rows,k) =                       &
              r_rho_levels(1:row_length,1:rows,k+1) -                  &
              r_rho_levels(1:row_length,1:rows,k)
      ENDDO
      p_layer_boundaries(1:row_length,1:rows,0)                =       &
              pstar(1:row_length,1:rows)
      p_layer_boundaries(1:row_length,1:rows,model_levels)     =       &
              p_theta_levels(1:row_length,1:rows,model_levels)
      p_layer_boundaries(1:row_length,1:rows,1:model_levels-1) =       &
              p_rho_levels(1:row_length,1:rows,2:model_levels)


! Mass:
! We make the hydrostatic assumption in the diagnostic mass calculation
! here so that mass = -b*(solid_angle/(3*g))*(r_top^3 - r_bottom^3)
! where   b =                                                           &
!    &    (p_layer_boundaries(:,:,k) - p_layer_boundaries(:,:,k-1))/    &
!    &    (r_theta_levels(:,:,k)     - r_theta_levels(:,:,k-1))         &

        IF (L_ukca_strat .OR. L_ukca_stratcfc .OR.                      &
             L_ukca_strattrop .AND. first)  THEN
         DO k=1,model_levels
           DO j=1,rows
            DO i=1,row_length
              mass(i,j,k) = (-solid_angle(i,j) / (3.*g)) *              &
            ((p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1))/ &
             (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)) ) *      &
             (r_theta_levels(i,j,k)**3. - r_theta_levels(i,j,k-1)**3.)
            END DO
          END DO
         END DO
        ENDIF      ! L_ukca_strat .... first

      DO k=1,model_levels
        DO j=1,rows
         DO i=1,row_length
           tnd(i,j,k) = p_theta_levels(i,j,k)                         &
                       /(zboltz*t_theta_levels(i,j,k))  ! molec/m3
            conv_ppn3d(i,j,k) = conv_rain3d(i,j,k) + conv_snow3d(i,j,k)
            ls_ppn3d(i,j,k)   = ls_rain3d  (i,j,k) + ls_snow3d  (i,j,k)

            conv_ppn3d(i,j,k) = MAX(conv_ppn3d(i,j,k), 0.0)
            ls_ppn3d  (i,j,k) = MAX(ls_ppn3d  (i,j,k), 0.0)
          ENDDO
        ENDDO
      ENDDO

!     Diagnostics for UKCA_MODE and aerosol chemistry

      IF (L_ukca_mode .OR. L_ukca_aerchem) THEN
        i1=lbound(rel_humid_frac,dim=1)
        j1=lbound(rel_humid_frac,dim=2)
        DO k = 1, wet_levels
! DEPENDS ON: qsat
          CALL QSAT(rel_humid_frac(i1,j1,k),t_theta_levels(i1,j1,k),   &
                p_theta_levels(i1,j1,k),SIZE(t_theta_levels(:,:,k)))

! Convert to relative humidity fraction
          rel_humid_frac(1:row_length,1:rows,k) =                      &
                q(1:row_length,1:rows,k)/                              &
                rel_humid_frac(1:row_length,1:rows,k)
        ENDDO
        WHERE (rel_humid_frac < 0.0)
          rel_humid_frac=0.0           ! remove negatives
        ENDWHERE
        WHERE (rel_humid_frac >= 1.0)
          rel_humid_frac=0.999         ! remove values >= 1.0
        ENDWHERE
      ENDIF

      IF (L_ukca_mode) THEN
        ALLOCATE(U_scalar_10m(row_length,rows))
! grid staggering, at N. Pole set final row=1st of penultimate row
        IF (at_extremity(PNorth) ) THEN
          DO i=1,rows-1
            U_scalar_10m(:,i)= SQRT(u_10m(1:row_length,i)**2 +          &
                                    v_10m(1:row_length,i)**2)
          ENDDO
          U_scalar_10m(:,rows)= U_scalar_10m(1,rows-1)
! at S. Pole set first row=1st of 2nd row
        ELSEIF (at_extremity(PSouth) ) THEN
          DO i=2,rows
            U_scalar_10m(:,i)= SQRT(u_10m(1:row_length,i)**2 +          &
                                    v_10m(1:row_length,i)**2)
          ENDDO
          U_scalar_10m(:,1)= U_scalar_10m(1,2)
        ELSE   ! non-polar points
          U_scalar_10m(:,:)= SQRT(u_10m(1:row_length,1:rows)**2 +       &
                                  v_10m(1:row_length,1:rows)**2)
        ENDIF

      ENDIF      ! L_ukca_mode

! Set up land fraction (required for MODE and DRY DEPN)
      land_fraction(:,:)=0.0
      DO l = 1, land_points
        j = (land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length
        land_fraction(i,j) = fland(l)
      END DO

      IF (first) THEN
        WRITE(6,*) ' ==================================='
        WRITE(6,*) '  MAX and MIN of UKCA INPUTS from D1'
        WRITE(6,*) ' ==================================='
        WRITE(6,*) 'Soil_Layer_Moisture: ',mype,                       &
           MAXVAL(soil_layer_moisture),MINVAL(soil_layer_moisture)
        DO j=1,n_use_tracers
         WRITE(6,*) 'troptracer: ',mype,j,                            &
           MAXVAL(all_tracers(:,:,:,j)),MINVAL(all_tracers(:,:,:,j))
        ENDDO
        DO j=1,n_chem_emissions
          WRITE(6,*) 'emission: ',mype,j,                              &
           MAXVAL(all_emissions(:,:,j)),MINVAL(all_emissions(:,:,j))
        ENDDO
       DO j=1,n_chem_diags
          WRITE(6,*) 'chemdiag1: ',mype, j,                            &
           MAXVAL(chem_diags(:,:,:,j)),MINVAL(chem_diags(:,:,:,j))
        ENDDO
        IF (L_ukca_BEflux) THEN
          DO j=1,n_BE_fluxdiags
            WRITE(6,*) 'BE_fluxdiag: ',mype, j,                        &
            MAXVAL(BE_fluxdiags(:,:,:,j)),MINVAL(BE_fluxdiags(:,:,:,j))
          ENDDO
       ENDIF
        IF (n_strat_fluxdiags > 0) THEN
          DO j=1,n_strat_fluxdiags
            WRITE(6,*) 'strat_fluxdiag: ',mype, j,                     &
                        MAXVAL(strat_fluxdiags(:,:,:,j)),              &
                        MINVAL(strat_fluxdiags(:,:,:,j))
          ENDDO
        ENDIF

        IF (L_ukca_dust) THEN
          WRITE(6,*) 'frac_types: ',mype,                              &
              MAXVAL(frac_types),MINVAL(frac_types)
          WRITE(6,*) 'tstar_tile: ',mype,                              &
              MAXVAL(tstar_tile),MINVAL(tstar_tile)
          WRITE(6,*) 'snow_tile: ',mype,MAXVAL(snow_tile),             &
              MINVAL(snow_tile)
          WRITE(6,*) 'soil_clay: ',mype,                               &
              MAXVAL(soil_clay),MINVAL(soil_clay)
          WRITE(6,*) 'dust_mrel1: ',mype,                              &
              MAXVAL(dust_mrel1),MINVAL(dust_mrel1)
          WRITE(6,*) 'dust_mrel2: ',mype,                              &
              MAXVAL(dust_mrel2),MINVAL(dust_mrel2)
          WRITE(6,*) 'dust_mrel3: ',mype,                              &
              MAXVAL(dust_mrel3),MINVAL(dust_mrel3)
          WRITE(6,*) 'dust_mrel4: ',mype,                              &
              MAXVAL(dust_mrel4),MINVAL(dust_mrel4)
          WRITE(6,*) 'dust_mrel5: ',mype,                              &
              MAXVAL(dust_mrel5),MINVAL(dust_mrel5)
          WRITE(6,*) 'dust_mrel6: ',mype,                              &
              MAXVAL(dust_mrel6),MINVAL(dust_mrel6)
          WRITE(6,*) 'dust_ustar: ',mype,                              &
              MAXVAL(dust_ustar),MINVAL(dust_ustar)
        ENDIF

        WRITE(6,*) 'Rough_length: ',mype,                              &
              MAXVAL(Rough_length),MINVAL(Rough_length)
        WRITE(6,*) 'Thick_bl_levels (level 1): ',mype,                 &
              MAXVAL(Thick_bl_levels(:,:,1)),                          &
              MINVAL(Thick_bl_levels(:,:,1))
        WRITE(6,*) 'conv_cloud_base: ',mype,                           &
              MAXVAL(conv_cloud_base(:,:)),                            &
              MINVAL(conv_cloud_base(:,:))
        WRITE(6,*) 'conv_cloud_top: ',mype,                            &
              MAXVAL(conv_cloud_top(:,:)),                             &
              MINVAL(conv_cloud_top(:,:))
        WRITE(6,*) 'land_sea_mask(1,1): ',land_sea_mask(1,1)
        WRITE(6,*) 'f3_at_u: ',mype,                                   &
              MAXVAL(f3_at_u(:,:)),                                    &
              MINVAL(f3_at_u(:,:))
        WRITE(6,*) 'FV_cos_theta_latitude: ',mype,                     &
              MAXVAL(FV_cos_theta_latitude(:,:)),                      &
              MINVAL(FV_cos_theta_latitude(:,:))

        DO k=1,model_levels
          WRITE(6,*) 'LEVEL: ',k
          WRITE(6,*) 'p_rho_levels:   ',mype,k,                        &
                      MAXVAL(p_rho_levels(1:row_length,1:rows,k)),     &
                      MINVAL(p_rho_levels(1:row_length,1:rows,k))
          WRITE(6,*) 'p_theta_levels: ',mype,k,                        &
                      MAXVAL(p_theta_levels(1:row_length,1:rows,k)),   &
                      MINVAL(p_theta_levels(1:row_length,1:rows,k))
          WRITE(6,*) 't_theta_levels: ',mype,k,                        &
                      MAXVAL(t_theta_levels(1:row_length,1:rows,k)),   &
                      MINVAL(t_theta_levels(1:row_length,1:rows,k))
          WRITE(6,*) 'rho_r2:         ',mype,k,                        &
                      MAXVAL(rho_r2(1:row_length,1:rows,k)),           &
                      MINVAL(rho_r2(1:row_length,1:rows,k))
          WRITE(6,*) 'ls_ppn_frac:    ',mype,k,                        &
                      MAXVAL(ls_ppn_frac(1:row_length,1:rows,k)),      &
                      MINVAL(ls_ppn_frac(1:row_length,1:rows,k))
          WRITE(6,*) 'ls_rain3d:      ',mype,k,                        &
                      MAXVAL(ls_rain3d(1:row_length,1:rows,k)),        &
                      MINVAL(ls_rain3d(1:row_length,1:rows,k))
          WRITE(6,*) 'ls_snow3d:      ',mype,k,                        &
                      MAXVAL(ls_snow3d(1:row_length,1:rows,k)),        &
                      MINVAL(ls_snow3d(1:row_length,1:rows,k))
         WRITE(6,*) 'conv_rain3d:      ',mype,k,                       &
                      MAXVAL(conv_rain3d(1:row_length,1:rows,k)),      &
                      MINVAL(conv_rain3d(1:row_length,1:rows,k))
         WRITE(6,*) 'conv_snow3d:      ',mype,k,                       &
                      MAXVAL(conv_snow3d(1:row_length,1:rows,k)),      &
                      MINVAL(conv_snow3d(1:row_length,1:rows,k))
         WRITE(6,*) 'aircraftems:      ',mype,k,                       &
                      MAXVAL(aircraftems(1:row_length,1:rows,k)),      &
                      MINVAL(aircraftems(1:row_length,1:rows,k))
          IF (L_ukca_aerchem) THEN
            WRITE(6,*) 'Natural SO2 emissions: ',mype,k,               &
                      MAXVAL(SO2emiss_3D(1:row_length,1:rows,k)),      &
                      MINVAL(SO2emiss_3D(1:row_length,1:rows,k))
          ENDIF
         WRITE(6,*) 'PV_on_theta_mlevs:      ',mype,k,                 &
                     MAXVAL(PV_on_theta_mlevs(1:row_length,1:rows,k)), &
                     MINVAL(PV_on_theta_mlevs(1:row_length,1:rows,k))
       ENDDO

        DO k=1,wet_levels
          WRITE(6,*) 'WET LEVEL: ',k
          WRITE(6,*) 'q:              ',mype,k,                        &
                      MAXVAL(q(1:row_length,1:rows,k)),                &
                      MINVAL(q(1:row_length,1:rows,k))
          WRITE(6,*) 'qcl:            ',mype,k,                        &
                      MAXVAL(qcl(1:row_length,1:rows,k)),              &
                      MINVAL(qcl(1:row_length,1:rows,k))
          WRITE(6,*) 'qcf:            ',mype,k,                        &
                      MAXVAL(qcf(1:row_length,1:rows,k)),              &
                      MINVAL(qcf(1:row_length,1:rows,k))
          IF (L_ukca_mode .OR. L_ukca_aerchem)                         &
            WRITE(6,*) 'rel_humid_frac: ',mype,k,                      &
                MAXVAL(rel_humid_frac(1:row_length,1:rows,k)),         &
                MINVAL(rel_humid_frac(1:row_length,1:rows,k))
        ENDDO
        DO k=1,model_levels+1
          WRITE(6,*) 'EXNER_rho_levels: ',mype,k,                      &
                      MAXVAL(exner_rho_levels(1:row_length,1:rows,k)), &
                      MINVAL(exner_rho_levels(1:row_length,1:rows,k))
        ENDDO
       DO k=1,bl_levels
         WRITE(6,*) 'BL LEVEL: ',k
         WRITE(6,*) 'rhokh_mix:     ',mype,k,                         &
                      MAXVAL(rhokh_mix(1:row_length,1:rows,k)),        &
                      MINVAL(rhokh_mix(1:row_length,1:rows,k))
          WRITE(6,*) 'dtrdz_charney_grid:     ',mype,k,                &
                MAXVAL(dtrdz_charney_grid(1:row_length,1:rows,k)),     &
                MINVAL(dtrdz_charney_grid(1:row_length,1:rows,k))
         WRITE(6,*) 'z_half:     ',mype,k,                             &
                MAXVAL(z_half(1:row_length,1:rows,k)),                 &
                MINVAL(z_half(1:row_length,1:rows,k))
        ENDDO
        IF (L_ukca_dust) THEN
         DO k=1,ndiv
           WRITE(6,*) 'DUST DIV: ',k
           WRITE(6,*) 'rb_dust_ndivs: ', mype,k,                      &
                      MAXVAL(rb_dust_ndivs(1:row_length,1:rows,k)),    &
                      MINVAL(rb_dust_ndivs(1:row_length,1:rows,k))
          ENDDO
        ENDIF
       WRITE(6,*) 'rho_aresist: ',MAXVAL(rho_aresist),                &
                                   MINVAL(rho_aresist)
        WRITE(6,*) 'aresist:  ', mype,MAXVAL(aresist), MINVAL(aresist)
        WRITE(6,*) 'resist_b: ', mype,MAXVAL(resist_b),MINVAL(resist_b)
        WRITE(6,*) 'kent:     ', mype,MAXVAL(kent),    MINVAL(kent)
       WRITE(6,*) 'kent_dsc: ', mype,MAXVAL(kent_dsc),MINVAL(kent_dsc)
       WRITE(6,*) 'zhsc:     ', mype,MAXVAL(zhsc),    MINVAL(zhsc)
        WRITE(6,*) 'ml_depth: ', mype,MAXVAL(ml_depth),MINVAL(ml_depth)
        DO k=1,3
         WRITE(6,*) 'we_lim: ', mype,k,                               &
                      MAXVAL(we_lim(1:row_length,1:rows,k)),           &
                      MINVAL(we_lim(1:row_length,1:rows,k))
          WRITE(6,*) 't_frac: ', mype,k,                               &
                      MAXVAL(t_frac(1:row_length,1:rows,k)),           &
                      MINVAL(t_frac(1:row_length,1:rows,k))
          WRITE(6,*) 'zrzi: ', mype,k,                                 &
                      MAXVAL(zrzi(1:row_length,1:rows,k)),             &
                      MINVAL(zrzi(1:row_length,1:rows,k))
          WRITE(6,*) 'we_lim_dsc: ', mype,k,                           &
                      MAXVAL(we_lim_dsc(1:row_length,1:rows,k)),       &
                      MINVAL(we_lim_dsc(1:row_length,1:rows,k))
          WRITE(6,*) 't_frac_dsc: ', mype,k,                           &
                      MAXVAL(t_frac_dsc(1:row_length,1:rows,k)),       &
                      MINVAL(t_frac_dsc(1:row_length,1:rows,k))
          WRITE(6,*) 'zrzi_dsc: ', mype,k,                             &
                      MAXVAL(zrzi_dsc(1:row_length,1:rows,k)),         &
                      MINVAL(zrzi_dsc(1:row_length,1:rows,k))
        ENDDO
        WRITE(6,*) 'pstar: ',    mype,MAXVAL(pstar),   MINVAL(pstar)
        WRITE(6,*) 'Tstar: ',    mype,MAXVAL(Tstar),   MINVAL(Tstar)
        WRITE(6,*) 'fland: ',    mype,MAXVAL(fland),   MINVAL(fland)
        WRITE(6,*) 'u_s  : ',    mype,MAXVAL(u_s  ),   MINVAL(u_s  )

!       Warning statements about possible mismatch between interactive
!       methane emissions and emissions ancillary

        IF (L_ukca_qch4inter) THEN
          WRITE(6,*) '******************WARNING*********************'
          WRITE(6,*) 'WARNING: INTERACTIVE CH4 EMISSIONS FROM WETLANDS'
          WRITE(6,*) 'ARE SWITCHED ON. Make sure that the surface'
          WRITE(6,*) 'methane emissions ancillary DOES NOT CONTAIN a'
          WRITE(6,*) 'contribution from wetlands and rice paddies.'
          WRITE(6,*) '******************WARNING*********************'
        ELSE
          WRITE(6,*) '******************WARNING*********************'
          WRITE(6,*) 'WARNING: INTERACTIVE CH4 EMISSIONS FROM WETLANDS'
          WRITE(6,*) 'ARE SWITCHED OFF. Make sure that the surface'
          WRITE(6,*) 'methane emissions ancillary CONTAINS a '
          WRITE(6,*) 'contribution from wetlands and rice paddies.'
          WRITE(6,*) '******************WARNING*********************'
       ENDIF

      ENDIF   ! End of IF (first) statement

!     Set up tile info

! DEPENDS ON: tilepts
      CALL TILEPTS(land_points,frac_types,tile_pts,tile_index)

!     Set tile fractions to 1 if aggregate tiles are used (NTILES=1).
!     Otherwise, set tile fractions to surface type fractions.

      tile_frac=0.0
      IF (ntiles == 1) THEN
        DO l=1,land_points
          tile_frac(l,1) = 1.
        ENDDO
      ELSE
        DO n=1,ntiles
          DO j=1,tile_pts(n)
            l = tile_index(j,n)
            tile_frac(l,n) = frac_types(l,n)
          ENDDO
        ENDDO
      ENDIF

      IF (L_ukca_dust) THEN

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA DUST MODEL',5)

!       Put fields into land arrays

        DO l = 1,land_points
          j = (land_index(l)-1)/row_length + 1
          i = land_index(l) - (j-1)*row_length
          pstar_land(l)   = pstar(i,j)
          rhostar_land(l) = pstar_land(l)/(tstar_tile(l,soil)*r)
          clay_land(l)    = soil_clay(i,j)
          DO n=1,ntiles
            ustar_tile(l,n) = dust_ustar(i,j,n)
            IF (ustar_tile(L,N) > 5.0) THEN
              WRITE(6,*) 'ustar_tile: ',ustar_tile(l,n),l,n
              WRITE(6,*) 'dust_ustar: ',dust_ustar(i,j,soil),i,j,soil
            ENDIF
          ENDDO
        ENDDO           ! land_Points

        DO l = 1,land_points
          j = (land_index(l)-1)/row_length + 1
          i = land_index(l) - (j-1)*row_length
          mrel_land(l,1) = dust_mrel1(i,j)
          mrel_land(l,2) = dust_mrel2(i,j)
          mrel_land(l,3) = dust_mrel3(i,j)
          mrel_land(l,4) = dust_mrel4(i,j)
          mrel_land(l,5) = dust_mrel5(i,j)
          mrel_land(l,6) = dust_mrel6(i,j)
        ENDDO            ! land_points

! DEPENDS ON: ukca_dust_ctl
        CALL UKCA_DUST_CTL(timestep,                                   &
          land_points, p_rho_levels, p_theta_levels, pstar, tstar,     &
          t_theta_levels(1:row_length,1:rows,:),                       &
          r_theta_levels, r_rho_levels,                                &
          rho_r2(1:row_length,1:rows,:), q(1:row_length,1:rows,:),     &
          qcl(1:row_length,1:rows,:), qcf(1:row_length,1:rows,:),      &
          ls_rain3d(1:row_length,1:rows,:),                            &
          ls_snow3d(1:row_length,1:rows,:),                            &
          all_tracers,                                                 &
          tile_pts,tile_index,tile_frac,fland,                         &
          pstar_land,tstar_tile,rhostar_land,                          &
          soil_layer_moisture,snow_tile,                               &
          ustar_tile,mrel_land,clay_land,                              &
          z_half, alpha_cd(1:bl_levels), ml_depth,                     &
          rhokh_mix, rho_aresist, aresist, resist_b,                   &
          dtrdz_charney_grid, kent, we_lim(:,:,1:3),                   &
          t_frac(:,:,1:3), zrzi(:,:,1:3), kent_dsc,                    &
          we_lim_dsc(:,:,1:3), t_frac_dsc(:,:,1:3),                    &
          zrzi_dsc(:,:,1:3), zhsc, rb_dust_ndivs,                      &
#include "arglndm.h"
#include "argsts.h"
        STASHwork34                                                    &
          )

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER('UKCA DUST MODEL',6)
      ENDIF    ! L_UKCA_dust

      IF (L_ukca_chem) THEN
! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA CHEMISTRY MODEL',5)

!       Initialise ASAD chemistry

        IF (FIRST) THEN
          ALLOCATE(c_species(n_chem_tracers+n_aero_tracers))
! DEPENDS ON: ukca_iniasad
          CALL UKCA_INIASAD(theta_field_size)
          CALL UKCA_CALC_CSPECIES
          IF (ANY(c_species < 0.001)) THEN
! check that ratio of molec. wts. array has been intialised.
            cmessage='c_species array contains zero value'
! DEPENDS ON: ereport
            CALL EREPORT('UKCA_MAIN',1,cmessage)
          ENDIF
        ENDIF

!       Calculate tropopause pressure using a combined
!       theta and PV surface

        CALL UKCA_CALC_TROPOPAUSE(row_length, rows, model_levels,      &
             sin_theta_latitude(1:row_length,1:rows),                  &
             theta(1:row_length,1:rows,1:model_levels),                &
             pv_on_theta_mlevs(1:row_length,1:rows,1:model_levels),    &
             p_layer_boundaries(1:row_length,1:rows,0:model_levels),   &
             p_theta_levels(1:row_length,1:rows,1:model_levels),       &
             p_tropopause(1:row_length,1:rows),                        &
             tropopause_level(1:row_length,1:rows))
                        
!       Calculate the difference in tracers per timestep in the
!       troposphere (moles/s) due to transport using the
!       trmol_post_atmstep array from this timestep and the
!       trmol_post_chem array from the previous timestep.
!       Do this only for those stratospheric flux diagnostics
!       which are switched on.

        IF (first .AND. n_strat_fluxdiags > 0) THEN

         ALLOCATE(trmol_post_chem(row_length,rows,model_levels,       &
                                    n_chem_tracers))  ! moles
          trmol_post_chem = 0.0
                
          DO l=1,n_strat_fluxdiags
            DO k=1,model_levels
             DO j=1,rows
                    DO i=1,row_length
                  strat_fluxdiags(i,j,k,l) = 0.0
               ENDDO
             ENDDO
           ENDDO
          ENDDO

       ELSE IF ((.NOT. first) .AND. n_strat_fluxdiags > 0) THEN

              ALLOCATE(trmol_post_atmstep(row_length,rows,model_levels,     &
                                    n_chem_tracers))  ! moles
         trmol_post_atmstep = 0.0

         DO l=1,n_chem_tracers
           DO k=1,model_levels
               DO j=1,rows
              DO i=1,row_length
                trmol_post_atmstep(i,j,k,l) = all_tracers(i,j,k,l)    &
                                            *tnd(i,j,k)*volume(i,j,k)  &
                                            /(c_species(l)*avogadro)   ! moles
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        icnt = 0
        DO l=1,n_chem_tracers
          IF (UkcaD1codes(istrat_first+l-1)%required) THEN
            icnt = icnt + 1

             DO k=1,model_levels
               DO j=1,rows
                 DO i=1,row_length
                IF (L_troposphere(i,j,k)) THEN              ! troposphere
                  strat_fluxdiags(i,j,k,icnt) =                     &
                      (trmol_post_atmstep(i,j,k,l)-                    &
                       trmol_post_chem(i,j,k,l))/timestep      ! moles/sec
                   ELSE                                        ! stratosphere
                     strat_fluxdiags(i,j,k,icnt) = 0.0
                   ENDIF
                 ENDDO
               ENDDO
            ENDDO
          ENDIF
        ENDDO

        DEALLOCATE(trmol_post_atmstep)  ! moles

        ENDIF      ! End of IF first and n_strat_fluxdiags statement

!       Do pre-processing for photolysis scheme

! The following (FAST-J calculation) only needs to be done when
! chemistry is called:

        IF (do_chemistry) THEN

          ALLOCATE(um_ozone3d(row_length,rows,ozone_levels))
          DO i=1,row_length
            um_ozone3d(i,:,:)=um_ozone(1,:,:)
          ENDDO

!       Convert 2D or 3D ozone field from UM to a 1D field

          nd_o3=size(um_ozone)
          ALLOCATE(um_ozone1D(nd_o3))
          um_ozone1D = reshape(um_ozone,(/nd_o3/))

!       Set the MOSES II flag if 8A boundary layer selected

          IF(h_sect(3)=="08A".OR. h_sect(3)=="08B") THEN
            l_moses_ii = .true.
          ELSE
            l_moses_ii = .false.
          ENDIF

          ALLOCATE(dj(row_length,rows,model_levels,jppj))
          dj=0.0

          IF (L_ukca_fastj) THEN
! define number of levels on which to do FAST-J
            IF (model_levels == 38) THEN
              fastj_levels = model_levels
            ELSEIF (model_levels == 60) THEN
              fastj_levels = model_levels - 5
            ELSE
              cmessage='Fast-J levels are not defined'
! DEPENDS ON: ereport
              CALL ereport('UKCA_MAIN1',1,cmessage)
            END IF

! DEPENDS ON: ukca_fastj
            CALL UKCA_FASTJ(                                           &
#include "arglndm.h"
            dj(1:row_length,1:rows,1:model_levels,1:jppj),             &
            p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
            pstar(1:row_length,1:rows),                                &
            p_theta_levels(1:row_length,1:rows,1:model_levels),        &
            t_theta_levels(1:row_length,1:rows,1:model_levels),        &
            tstar(1:row_length,1:rows),                                &
            so4_aitken(1:row_length,1:rows,1:tr_levels),               &
            so4_accum(1:row_length,1:rows,1:tr_levels),                &
            q(1:row_length,1:rows,1:wet_levels),                       &
            qcl(1:row_length,1:rows,1:wet_levels),                     &
            qcf(1:row_length,1:rows,1:wet_levels),                     &
            conv_cloud_lwp(1:row_length,1:rows),                       &
            conv_cloud_top(1:row_length,1:rows),                       &
            conv_cloud_base(1:row_length,1:rows),                      &
            conv_cloud_amount(1:row_length,1:rows,1:wet_levels),       &
            surf_albedo(1:row_length,1:rows),                          &
            nd_o3, um_ozone3d(1:row_length,1:rows,1:ozone_levels),     &
            fland,                                                     &
            a_realhd(rh_z_top_theta),                                  &
            l_moses_ii)

          ENDIF  ! End of IF (L_ukca_fastj) statement

!       Clear workspace

          DEALLOCATE(um_ozone1D)
         DEALLOCATE(surf_albedo)

        END IF   ! do_chemistry

!       Allocate array for wetland ch4 emissions when
!       L_ukca_qch4inter is false. Set emissions to zero.

        IF (.NOT.(L_ukca_qch4inter)) THEN
         ALLOCATE(ch4_wetl_emiss(1:row_length,1:rows))
         ch4_wetl_emiss(1:row_length,1:rows) = 0.0
        ENDIF

! Emission part. Before calling emissions,
! set up array of lower-boundary mixing ratios for stratospheric species.
! Note that some mixing ratios (CH4, N2O, F11, F12, F113, H22) are defined
! in cruntimc.h through the UMUI, others are not. If for some species
! actual emissions (not prescribed lower boundary conditions) are required,
! remove these species from lbc_spec, e.g, for CH4.
! Mass mixing ratios are best estimates, taken from the WMO Scientific
! assessment of ozone depletion (2002), for 2000 conditions.
! H2 corresponds to 0.5 ppmv.
! Select lower boundary MMRs. Either from UMUI or predefined constants
! Note that for bromine compounds (MeBr, H1211, H1301) are increased
! by 24% to account for non-included species and HBr is set to 0.
!
        IF (L_ukca_strat .OR. L_ukca_strattrop) THEN
          IF (L_ukca_useumuivals) THEN
            lbc_mmr( 1: 3) = (/ n2ommr, c12mmr     , c11mmr /)
            lbc_mmr(15:16) = (/c113mmr,  hcfc22mmr          /)
            lbc_mmr(29   ) = ch4mmr
          ENDIF

          IF ((first) .AND.(.NOT. L_UKCA_stratcfc)) THEN
! Add effects of unincluded tracers onto MeBr, F11, and F12
! MeBr: add effects of CF2ClBr and CF3Br:
            lbc_mmr(4) = lbc_mmr(4) + c_mebr *                          &
                     (h1211_mmr / c_cf2clbr + h1301_mmr / c_cf3br)

! F12: add effects of F113 and F22 onto F12
            lbc_mmr(2) = lbc_mmr(2) + 0.5*c_cf2cl2*                     &
                (3.*lbc_mmr(15) / c_cf2clcfcl2 + lbc_mmr(16) / c_chf2cl)

! F11: add effects of CCl4, MeCl, MeCCl3, and CF2ClBr
            lbc_mmr(3) = lbc_mmr(3) + c_cfcl3 / 3. *                    &
                (4.*ccl4_mmr / c_ccl4 + mecl_mmr / c_mecl               &
                +3.*meccl3_mmr / c_meccl3 + h1211_mmr / c_cf2clbr)
          ENDIF
        ENDIF    ! L_ukca_strat etc

! Stratospheric budget terms, but always needed in calls to routines
        ALLOCATE( semiss (row_length, rows, model_levels))
        ALLOCATE( phot   (row_length, rows, model_levels, nphotf))
        ALLOCATE( bimol  (row_length, rows, model_levels, nbimol))
        ALLOCATE( termol (row_length, rows, model_levels))
        ALLOCATE( hetero (row_length, rows, model_levels))
        ALLOCATE( drydep (row_length, rows, model_levels))
        ALLOCATE( wetdep (row_length, rows, model_levels))
        ALLOCATE(noxlight(row_length, rows, model_levels))
        ALLOCATE(abundance(row_length, rows, model_levels, nabund))

        IF (l_ukca_o3budget) THEN
          n_opl = 20
        ELSE    ! In call, so needs to be allocated
          n_opl = 1
        ENDIF
        ALLOCATE(ozonebud(row_length, rows, model_levels, n_opl))

! DEPENDS ON: ukca_emission_ctl
        CALL UKCA_EMISSION_CTL(                                        &
           n_chem_tracers+n_aero_tracers, n_use_emissions,             &
           timestep,                                                   &
           em_chem_spec(1:n_use_emissions),                            &
           f3_at_u(1:row_length,1:rows),                               &
           r_rho_levels(1:row_length,1:rows,1:model_levels),           &
           r_theta_levels(1:row_length,1:rows,0:model_levels),         &
           sin_theta_latitude(1:row_length,1:rows),                    &
           FV_cos_theta_latitude(1:row_length,1:rows),                 &
           true_longitude, delta_lambda, delta_phi,                    &
           i_year, i_month, i_day,                                     &
           tropopause_height(1:row_length,1:rows),                     &
           land_sea_mask(1:row_length,1:rows),                         &
           conv_cloud_base(1:row_length,1:rows),                       &
           conv_cloud_top(1:row_length,1:rows),                        &
           theta(1:row_length,1:rows,1:model_levels),                  &
           q(1:row_length,1:rows,1:wet_levels),                        &
           qcl(1:row_length,1:rows,1:wet_levels),                      &
           qcf(1:row_length,1:rows,1:wet_levels),                      &
           exner_rho_levels(1:row_length,1:rows,1:model_levels+1),     &
           rho_r2(1:row_length,1:rows,1:model_levels),                 &
           p_layer_boundaries(1:row_length,1:rows,0:model_levels),     &
           p_theta_levels(1:row_length,1:rows,1:model_levels),         &
           all_emissions(1:row_length,1:rows,1:n_chem_emissions),      &
           aircraftems(1:row_length,1:rows,1:model_levels),            &
           SO2emiss_3D(1:row_length,1:rows,1:model_levels),            &
           z_half, alpha_cd(1:bl_levels), ml_depth,                    &
           rhokh_mix, rho_aresist, aresist, resist_b,                  &
           dtrdz_charney_grid, kent, we_lim(:,:,1:3),                  &
           t_frac(:,:,1:3), zrzi(:,:,1:3), kent_dsc,                   &
           we_lim_dsc(:,:,1:3), t_frac_dsc(:,:,1:3),                   &
           zrzi_dsc(:,:,1:3), zhsc, rb_dust_ndivs,                     &
           ch4_wetl_emiss(1:row_length,1:rows),                        &
           all_tracers(1:row_length,1:rows,:,1:n_chem_tracers+         &
                       n_aero_tracers),                                &
           semiss, noxlight,                                           &
           n_boundary_vals, lbc_spec, lbc_mmr,                         &
           mass)

! Do chemistry calculation here (at chemistry timesteps)
        IF (do_chemistry) THEN

!       Apply correction to Temperature and/or Specific Humidity Bias
          ALLOCATE(T_chem(row_length,rows,model_levels))
          ALLOCATE(Q_chem(row_length,rows,model_levels))

         IF (L_ukca_Tbias) THEN
           T_chem(:,:,:) = t_theta_levels(1:row_length,1:rows,:)       &
                          - T_bias(1:row_length,1:rows,:)
          ELSE
            T_chem(:,:,:) = t_theta_levels(1:row_length,1:rows,:)
          ENDIF

          IF (L_ukca_Qbias) THEN
            Q_chem(:,:,:) = q(1:row_length,1:rows,:)                    &
                          - Q_bias(1:row_length,1:rows,:)
            WHERE (Q_chem < 0.0)
              Q_chem = 0.0
           ENDWHERE
         ELSE
           Q_chem(:,:,:) = q(1:row_length,1:rows,:)
         ENDIF

!       Calculate solar zenith angle

          ALLOCATE(cos_zenith_angle(row_length, rows))

! DEPENDS ON: solpos
          CALL SOLPOS (PREVIOUS_TIME(7), PREVIOUS_TIME(1),              &
               LCAL360, L_SEC_VAR, L_EqT, Eq_Time, Sindec, SCS)

          secondssincemidnight = REAL(Previous_time(4)*3600             &
                               +      Previous_time(5)*60               &
                               +      Previous_time(6))

! DEPENDS ON: ukca_solang
          CALL UKCA_SOLANG(Sindec, secondssincemidnight,                &
                      timestep, Eq_Time,                                &
                      sin_theta_latitude, true_longitude,               &
                      theta_field_size, cos_zenith_angle )

! DEPENDS ON: ukca_chemistry_ctl
         CALL UKCA_CHEMISTRY_CTL(I_month, I_day_number, I_hour,        &
             I_minute - int(timestep)/60,                              &
             REAL(interval)*timestep,                                  &
             model_levels,                                             &
             n_chem_tracers+n_aero_tracers, n_chem_diags,              &
             n_BE_fluxdiags,                                           &
             f3_at_u(1:row_length,1:rows)/two_omega,                   &
             FV_cos_theta_latitude(1:row_length,1:rows),               &
             true_longitude,                                           &
             p_theta_levels(1:row_length,1:rows,1:model_levels),       &
             T_chem(1:row_length,1:rows,:),                            &
             Q_chem(1:row_length,1:rows,:),                            &
             qcf(1:row_length,1:rows,:),                               &
             qcl(1:row_length, 1:rows, :),                             &
             rel_humid_frac(1:row_length, 1:rows, :),                  &
             p_layer_boundaries(1:row_length,1:rows,0:model_levels),   &
             r_theta_levels(1:row_length, 1:rows, :),                  &
             cos_zenith_angle,                                         &
             all_tracers(1:row_length,1:rows,:,                        &
                         1:n_chem_tracers+n_aero_tracers),             &
             chem_diags(1:row_length,1:rows,:,:),                      &
             BE_fluxdiags(1:row_length,1:rows,:,:),                    &
             Tstar,                                                    &
             Thick_bl_levels,                                          &
             Rough_length,                                             &
             u_s,                                                      &
             ls_ppn3d, conv_ppn3d,                                     &
             cloud_frac(1:row_length, 1:rows, :),                      &
             dj(:,:,:,:),                                              &
             volume(:,:,:),                                            &
             mass(:,:,:),                                              &
! Extra variables for new dry dep scheme
             land_points, land_index,                                  &
             tile_pts, tile_index, tile_frac,                          &
             zbl, surf_hf, seaice_frac, stcon,                         &
             soil_layer_moisture(:,1), fland,                          &
             laift_lp, canhtft_lp,                                     &
             z0tile_lp, tstar_tile, canwctile_lp,                      &
             nbimol, nphotf, n_opl, nabund,                            &
             abundance, phot, bimol, termol, hetero,                   &
             drydep, wetdep,                                           &
             ozonebud,                                                 &
             pv_on_theta_mlevs,                                        &
             theta(1:row_length,1:rows,:),                             &
             um_ozone3d,                                               &
             delso2_wet_h2o2,                                          &
             delso2_wet_o3,                                            &
             delso2_dry_oh,                                            &
             delso2_drydep,                                            &
             delso2_wetdep,                                            &
             so4_sa                                                    &
          )

          DEALLOCATE(cos_zenith_angle)
          DEALLOCATE(T_chem)
          DEALLOCATE(Q_chem)

          IF (L_ukca_budget2) THEN
           IF (n_chem_diags > 5)                                        &
               chem_diags(1:row_length,1:rows,:, 5: 6) = phot
           IF (n_chem_diags > 10)                                       &
               chem_diags(1:row_length,1:rows,:, 7:11) = bimol
           IF (n_chem_diags > 11)                                       &
               chem_diags(1:row_length,1:rows,:,12   ) = termol
           IF (n_chem_diags > 12)                                       &
               chem_diags(1:row_length,1:rows,:,13   ) = hetero
           IF (n_chem_diags > 13)                                       &
               chem_diags(1:row_length,1:rows,:,14   ) = drydep
           IF (n_chem_diags > 14)                                       &
               chem_diags(1:row_length,1:rows,:,15   ) = wetdep
           IF (n_chem_diags > 16)                                       &
               chem_diags(1:row_length,1:rows,:,16:17) = abundance
           IF ((L_ukca_o3budget) .AND. (n_chem_diags >= n_opl))         &
! update special ozone budget
           chem_diags(1:row_length,1:rows,:,1:n_opl)  = ozonebud

             IF (n_chem_diags > 17) THEN 
               chem_diags(1:row_length,1:rows,:,18) = so4_sa
             ENDIF
           END IF

         END IF ! do_chemistry

! Keep information on emissions
         IF (L_ukca_budget2) THEN
           IF (n_chem_diags > 2)                                        &
             chem_diags(1:row_length,1:rows,:, 3   ) = semiss
           IF (n_chem_diags > 3)                                        &
             chem_diags(1:row_length,1:rows,:, 4   ) = noxlight
         END IF

        IF (n_strat_fluxdiags > 0) THEN
         DO l=1,n_chem_tracers
           DO k=1,model_levels
             DO j=1,rows
               DO i=1,row_length
                 trmol_post_chem(i,j,k,l) = all_tracers(i,j,k,l)      &
                                           *tnd(i,j,k)*volume(i,j,k)   &
                                           /(c_species(l)*avogadro)    ! moles
             ENDDO
             ENDDO
           ENDDO
         ENDDO
        ENDIF   ! End of IF n_strat_fluxdiags > 0 statement


! Transform halogen/nitrogen/hydrogen species back
      IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
       (L_ukca_stratcfc)) THEN
! DEPENDS ON: ukca_transform_halogen
        CALL ukca_transform_halogen(tr_ukca,rows,row_length,            &
                             model_levels, offx, offy, all_tracers,     &
                             halo_i, halo_j, q, .FALSE., 0)
      ENDIF

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA CHEMISTRY MODEL',6)
      ENDIF    ! End of IF (L_ukca_chem) statement

      IF(L_UKCA_MODE) THEN

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA AEROSOL MODEL  ',5)

! DEPENDS ON: ukca_aero_ctl
        CALL UKCA_AERO_CTL(i_month, i_day_number, i_hour,               &
             I_minute - int(timestep)/60,                               &
             timestep,                                                  &
             model_levels, rows, row_length,                            &
             wet_levels,                                                &
             global_row_length,global_rows,                             &
             n_aero_tracers-1,                                          &
             n_mode_tracers,                                            &
             area,                                                      &
             p_theta_levels(1:row_length,1:rows,1:model_levels),        &
             t_theta_levels(1:row_length,1:rows,:),                     &
             q(1:row_length,1:rows,:),                                  &
             rel_humid_frac(1:row_length,1:rows,:),                     &
             p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
             all_tracers(1:row_length,1:rows,:,                         &
                  n_chem_tracers+1:                                     &
                  n_chem_tracers+n_aero_tracers-1),                     &
             all_tracers(1:row_length,1:rows,:,                         &
                  n_chem_tracers+n_aero_tracers+1:                      &
                  n_chem_tracers+n_aero_tracers+                        &
                  n_mode_tracers),                                      &
             bl_levels,                                                 &
             Tstar(1:row_length,1:rows),                                &
             seaice_frac(1:row_length,1:rows),                          &
             Rough_length(1:row_length,1:rows),                         &
             u_s,                                                       &
             U_scalar_10m,                                              &
             ls_rain3d(1:row_length,1:rows,:),                          &
             conv_rain3d,                                               &
             land_fraction,                                             &
             theta_field_size,                                          &
             delso2_wet_h2o2,                                           &
             delso2_wet_o3,                                             &
             mode_diags(1:row_length,1:rows,:,:),                       &
             n_use_emissions,                                           &
             em_chem_spec(1:n_use_emissions),                           &
             all_emissions(1:row_length,1:rows,                         &
                           1:n_chem_emissions),                         &
             SO2emiss_3D(1:row_length,1:rows,                           &
                         1:model_levels),                               &
             cloud_frac(1:row_length, 1:rows, :) )

! DEPENDS ON: timer
        IF (LTIMER) CALL TIMER('UKCA AEROSOL MODEL  ',6)
      ENDIF   ! L_UKCA_MODE

!     Return fields to D1

      CALL PUTD1FLDS

! In case water vapour feedback is activated, call Q_POS_CNTL to remove points
! with low water vapour and set polar rows to uniform constants

      IF (L_ukca_h2o_feedback)                                          &
! DEPENDS ON: q_pos_ctl
        CALL Q_Pos_Ctl(  D1(jq(1)), row_length, rows, wet_levels,       &
                         global_row_length, global_rows,                &
                         mype, nproc, halo_i, halo_j,                   &
!                        gc_all_proc_group, model_domain,               &
                                            model_domain,               &
                         halo_type_extended, l_q_pos_local, qlimit)

! Remove negatives from all chemical tracers and reset polar rows. Note: This
! call replaces update my_tracer_qpos4.mf77 that would have placed this call
! into atmstep2 subroutine. The call is necessary here to avoid a crash
! if radiative feedback is selected.

! DEPENDS ON: q_pos_ctl
      CALL Q_Pos_Ctl(                                                   &
              D1(jtr_ukca(1,1)), row_length, rows, tr_ukca*tr_levels,   &
                         global_row_length, global_rows,                &
                         mype, nproc, offx, offy,                       &
!                        gc_all_proc_group,model_domain,                &
                                           model_domain,                &
                         halo_type_single, .true. ,  0.0 )

! Add output fields to STASHwork arrays.

!      DO I = A_max_ukcavars+1,ndiag
!        item = I
!        IF (sf(item,ukca_sect)) THEN
! DEPENDS ON: copydiag
!          CALL COPYDIAG (STASHwork34(si(item,ukca_sect,IM_INDEX)),
!     &      XXXXXXXXX(1:row_length,1:rows,IDIV),
!     &      row_length,rows,0,0,0,0, at_extremity,
!     &      atmos_im,sect,item,
!     &      errorstatus,cmessage)
!          IF (errorstatus.GT.0) THEN
!            cmessage=": ERROR IN COPYDIAG"
! DEPENDS ON: ereport
!            CALL EREPORT( 'UKCA_MAIN1', errorstatus, cmessage )
!          ENDIF
!        ENDIF
!      ENDDO

!     Deallocate arrays

!      CALL ASAD_MOD_FINAL

      DEALLOCATE(tnd)
      DEALLOCATE(um_ozone3d)
      DEALLOCATE(conv_cloud_base)
      DEALLOCATE(conv_cloud_top)
      DEALLOCATE(land_sea_mask)
      DEALLOCATE(aircraftems)
      DEALLOCATE(theta)
      DEALLOCATE(p_rho_levels)
      DEALLOCATE(p_theta_levels)
      DEALLOCATE(exner_theta_levels)
      DEALLOCATE(exner_rho_levels)
      DEALLOCATE(q)
      DEALLOCATE(rho_r2)
      DEALLOCATE(qcl)
      DEALLOCATE(qcf)
      DEALLOCATE(ls_ppn_frac)
      DEALLOCATE(cloud_liq_water)
      DEALLOCATE(cloud_ice_content)
      DEALLOCATE(ls_rain3d)
      DEALLOCATE(ls_snow3d)
      DEALLOCATE(ls_ppn3d)
      DEALLOCATE(conv_rain3d)
      DEALLOCATE(conv_snow3d)
      DEALLOCATE(conv_ppn3d)
      DEALLOCATE(conv_cloud_amount)
      DEALLOCATE(area_cloud_fraction)
      DEALLOCATE(so4_aitken)
      DEALLOCATE(so4_accum)
      DEALLOCATE(dj)
      DEALLOCATE(rel_humid_frac)
      IF (l_ukca_mode) THEN
        DEALLOCATE(u_10m)
        DEALLOCATE(v_10m)
        DEALLOCATE(U_scalar_10m)
      ENDIF
      DEALLOCATE(pstar)
      DEALLOCATE(tstar)
      IF (l_ukca_dust) THEN
        DEALLOCATE(Soil_clay)
        DEALLOCATE(dust_mrel1)
        DEALLOCATE(dust_mrel2)
        DEALLOCATE(dust_mrel3)
        DEALLOCATE(dust_mrel4)
        DEALLOCATE(dust_mrel5)
        DEALLOCATE(dust_mrel6)
        DEALLOCATE(dust_ustar)
      ENDIF
      DEALLOCATE(Soil_Layer_Moisture)
      DEALLOCATE(Tile_Frac)
      DEALLOCATE(Tstar_tile)
      DEALLOCATE(Snow_tile)
      DEALLOCATE(Frac_types)
      DEALLOCATE(Rough_length)
      DEALLOCATE(fland)
      DEALLOCATE(zbl)
      DEALLOCATE(surf_hf)
      DEALLOCATE(seaice_frac)
      DEALLOCATE(conv_cloud_lwp)
      DEALLOCATE(u_s)
      DEALLOCATE(stcon)
      DEALLOCATE(laift_lp)
      DEALLOCATE(canhtft_lp)
      DEALLOCATE(z0tile_lp)
      DEALLOCATE(canwctile_lp)
      DEALLOCATE(Thick_bl_levels)
      DEALLOCATE(STASHwork34)
      DEALLOCATE(p_layer_boundaries)
      DEALLOCATE(t_theta_levels)
      DEALLOCATE(ch4_wetl_emiss)
      DEALLOCATE(delSO2_wet_h2o2)
      DEALLOCATE(delSO2_wet_o3)
      DEALLOCATE(delSO2_dry_oh)
      DEALLOCATE(delSO2_drydep)
      DEALLOCATE(delSO2_wetdep)
      DEALLOCATE(SO2emiss_3D)
      DEALLOCATE( semiss )
      DEALLOCATE( phot   )
      DEALLOCATE( bimol  )
      DEALLOCATE( termol )
      DEALLOCATE( hetero )
      DEALLOCATE( drydep )
      DEALLOCATE( wetdep )
      DEALLOCATE(noxlight)
      DEALLOCATE(abundance)
      DEALLOCATE(ozonebud)
      DEALLOCATE(so4_sa)
      DEALLOCATE(tropopause_height)
      DEALLOCATE(mode_diags)

      first=.false.

      CONTAINS
! ######################################################################
      SUBROUTINE GETD1FLDS(N)

! Interface block to allow generic calls
      INTERFACE UKCA_EXTRACT_D1_DATA

! DEPENDS ON: UKCA_EXTRACT_D1_DATA1D
        SUBROUTINE UKCA_EXTRACT_D1_DATA1D(                             &
#include "argd1.h"
        first,N,X)
          USE UKCA_D1_DEFS
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          REAL, DIMENSION(:), INTENT(INOUT) :: X  ! extracted array
          REAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_DATA1D

! DEPENDS ON: UKCA_EXTRACT_D1_DATA2D
        SUBROUTINE UKCA_EXTRACT_D1_DATA2D(                             &
#include "argd1.h"
        first,N,X)
          USE UKCA_D1_DEFS
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          REAL, DIMENSION(:,:), INTENT(INOUT) :: X  ! extracted array
          REAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_DATA2D

! DEPENDS ON: UKCA_EXTRACT_D1_INTEGER_DATA2D
        SUBROUTINE UKCA_EXTRACT_D1_INTEGER_DATA2D(                     &
#include "argd1.h"
        first,N,X)
          USE UKCA_D1_DEFS
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          INTEGER, DIMENSION(:,:), INTENT(INOUT) :: X  ! extracted array
          INTEGER, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_INTEGER_DATA2D

! DEPENDS ON: UKCA_EXTRACT_D1_LOGICAL_DATA2D
        SUBROUTINE UKCA_EXTRACT_D1_LOGICAL_DATA2D(                     &
#include "argd1.h"
        first,N,X)
          USE UKCA_D1_DEFS
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          LOGICAL, DIMENSION(:,:), INTENT(INOUT) :: X  ! extracted array
          LOGICAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_LOGICAL_DATA2D

! DEPENDS ON: UKCA_EXTRACT_D1_DATA3D
        SUBROUTINE UKCA_EXTRACT_D1_DATA3D(                             &
#include "argd1.h"
        first,N,X)
          USE UKCA_D1_DEFS
#include "parvars.h"
#include "typsize.h"
#include "typd1.h"
          INTEGER, INTENT(IN)               :: N  ! id of array
          LOGICAL, INTENT(IN)               :: first
          REAL, DIMENSION(:,:,:), INTENT(INOUT) :: X  ! extracted array
          REAL, DIMENSION(:), ALLOCATABLE   :: data1
          CHARACTER(LEN=72)                 :: cmessage
          INTEGER                           :: ERRCODE
        END SUBROUTINE UKCA_EXTRACT_D1_DATA3D
      END INTERFACE

      INTERFACE
! DEPENDS ON: UKCA_SET_ARRAY_BOUNDS
        SUBROUTINE UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          USE UKCA_D1_DEFS
          INTEGER, INTENT(IN)  :: N            ! id of array
          INTEGER, INTENT(OUT) :: I1,I2,J1,J2  ! array bounds
        END SUBROUTINE UKCA_SET_ARRAY_BOUNDS
      END INTERFACE


      INTEGER, INTENT(IN)     :: N             ! id of array

      INTEGER, SAVE           :: next_tr       ! count tracers
      INTEGER, SAVE           :: next_em       ! count emissions
      INTEGER, SAVE           :: next_cd       ! count chem diags
      INTEGER                 :: I1,I2,J1,J2   ! array bounds


!     Items in range 1-150 fill all_tracers array in sequence

      IF (UkcaD1codes(N)%section == UKCA_sect .AND.                    &
          UkcaD1codes(N)%item <= n_all_tracers) THEN
        IF (.NOT. ALLOCATED(tracer1)) THEN
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tracer1(I1:I2,J1:J2,                                &
                              ukcaD1codes(N)%len_dim3))
          ALLOCATE(all_tracers(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3,n_use_tracers))
          next_tr = 1
        ENDIF

        IF (.NOT. ALLOCATED(tr_index)) ALLOCATE(tr_index(n_use_tracers))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
#include "argd1.h"
           first,N,tracer1)
        all_tracers(:,:,:,next_tr) = tracer1(:,:,:)
        tr_index(next_tr)          = UkcaD1codes(N)%item
        next_tr                    = next_tr+1

!     Items in range 150+23 - 150+49 fill all_emissions array in sequence

      ELSEIF ((UkcaD1codes(N)%section == UKCA_ems_sect .AND.           &
          UkcaD1codes(N)%item >= n_emiss_first .AND.                   &
          UkcaD1codes(N)%item <= n_emiss_last ) .OR.                   &
         (UkcaD1codes(N)%section == 0 .AND.                            &
         (UkcaD1codes(N)%item == 58  .OR.                              &
          UkcaD1codes(N)%item == 126 .OR.                              &
          UkcaD1codes(N)%item == 127 .OR.                              &
          UkcaD1codes(N)%item == 128 .OR.                              &
          UkcaD1codes(N)%item == 129 .OR.                              &
          UkcaD1codes(N)%item == 130 .OR.                              &
          UkcaD1codes(N)%item == 131)) .OR.                            &
         (UkcaD1codes(N)%section == 17 .AND.                           &
          UkcaD1codes(N)%item == 205 )) THEN
        IF (.NOT. ALLOCATED(emission1)) THEN
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(emission1(I1:I2,J1:J2))
          ALLOCATE(all_emissions(I1:I2,J1:J2,n_chem_emissions))
          next_em = 1
        ENDIF
        IF (.NOT. ALLOCATED(em_index))                                 &
                ALLOCATE(em_index(n_chem_emissions))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
#include "argd1.h"
           first,N,emission1)
        all_emissions(:,:,next_em) = emission1(:,:)
        em_index(next_em)          = UkcaD1codes(N)%item
        next_em                    = next_em+1

! Add 3-D SO2 emissions
      ELSEIF (UkcaD1codes(N)%section == 0 .AND.                        &
              UkcaD1codes(N)%item == 121) THEN
        CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        ALLOCATE(SO2emiss_3D(I1:I2,J1:J2,ukcaD1codes(N)%len_dim3))
        CALL UKCA_EXTRACT_D1_DATA(                                     &
#include "argd1.h"
           first,N,SO2emiss_3D)

! Add 3-D Aircraft emissions
      ELSEIF (UkcaD1codes(N)%section == 0 .AND.                        &
              UkcaD1codes(N)%item == 340) THEN
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(aircraftems(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,aircraftems)

!     Prognostics: fill appropriate array

      ELSEIF (UkcaD1codes(N)%section == 0 .AND.                        &
              UkcaD1codes(N)%prognostic) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(4)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(theta(I1:I2,J1:J2,                                  &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,theta)
        CASE(9)
          ALLOCATE(soil_layer_moisture(ukcaD1codes(N)%len_dim1,        &
                                       ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,soil_layer_moisture)
        CASE(10)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(q(I1:I2,J1:J2,                                      &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,q)

!         Fill advected chemical H2O tracer if required, n.b.
!         halo size of q is larger than the other tracers

          IF (l_UKCA_advh2o) THEN
            all_tracers(:,:,1:wet_levels,next_tr) =                    &
              q(1-offx:row_length+offx,1-offy:rows+offy,:)
            next_tr = next_tr + 1
          ENDIF
        CASE(12)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(qcf(I1:I2,J1:J2,                                    &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,qcf)
        CASE(14)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_base(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,conv_cloud_base)
        CASE(15)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_top(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,conv_cloud_top)
        CASE(16)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_lwp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,conv_cloud_lwp)
        CASE(24)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tstar(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,tstar)
        CASE(25)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zbl(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,zbl)
        CASE(26)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(Rough_length(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,Rough_length)
        CASE(30)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(land_sea_mask(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,land_sea_mask)
        CASE(31)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(seaice_frac(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,seaice_frac)
        CASE(60)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(um_ozone(I1:I2,J1:J2,                               &
                              ukcaD1codes(N)%len_dim3))
          um_ozone = RESHAPE(D1(UkcaD1Codes(N)%address:                &
             UkcaD1Codes(N)%address+UkcaD1Codes(N)%length-1),          &
             (/SIZE(um_ozone,DIM=1),SIZE(um_ozone,DIM=2),              &
             SIZE(um_ozone,DIM=3)/))
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
        CASE(103)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(so4_aitken(I1:I2,J1:J2,                             &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,so4_aitken)
        CASE(104)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(so4_accum(I1:I2,J1:J2,                              &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,so4_accum)
        CASE(211)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_cloud_amount(I1:I2,J1:J2,                      &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,conv_cloud_amount)
        CASE(216)
! DEPENDS ON: ukca_set_array_bounds
          ALLOCATE(frac_types(ukcaD1codes(N)%len_dim1,                 &
                              ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,frac_types)
        CASE(217)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(laift_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,laift_lp)
        CASE(218)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(canhtft_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,canhtft_lp)
        CASE(229)
! DEPENDS ON: ukca_set_array_bounds
           CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
           ALLOCATE(canwctile_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,canwctile_lp)
        CASE(233)
          ALLOCATE(tstar_tile(ukcaD1codes(N)%len_dim1,                 &
                              ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,tstar_tile)
        CASE(234)
! DEPENDS ON: ukca_set_array_bounds
           CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
           ALLOCATE(z0tile_lp(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,z0tile_lp)
        CASE(240)
          ALLOCATE(snow_tile(ukcaD1codes(N)%len_dim1,                  &
                             ukcaD1codes(N)%len_dim2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,snow_tile)
        CASE(253)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(rho_r2(I1:I2,J1:J2,                                 &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,rho_r2)
        CASE(254)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(qcl(I1:I2,J1:J2,                                    &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,qcl)
        CASE(255)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(exner_rho_levels(I1:I2,J1:J2,                       &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,exner_rho_levels)
        CASE(265)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(area_cloud_fraction(I1:I2,J1:J2,                    &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,area_cloud_fraction)
        CASE(266)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_frac(I1:I2,J1:J2,ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,cloud_frac)
        CASE(338)                        ! (HadGAM1-ERA) T bias
         CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
         ALLOCATE(T_bias(I1:I2,J1:J2,                                 &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
          first,N,T_bias)
        CASE(339)                        ! (HadGAM1-ERA) Sp. humidity bias
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(Q_bias(I1:I2,J1:J2,                                 &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,Q_bias)
        CASE(340)                                      ! aircraft emiss
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(aircraftems(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,aircraftems)
        CASE(418)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(soil_clay(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,soil_clay)
        CASE(421)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_mrel1(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_mrel1)
        CASE(422)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_mrel2(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_mrel2)
        CASE(423)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_mrel3(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_mrel3)
        CASE(424)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_mrel4(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_mrel4)
        CASE(425)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_mrel5(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_mrel5)
        CASE(426)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_mrel6(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_mrel6)
        CASE(505)
          ALLOCATE(fland(ukcaD1codes(N)%len_dim1))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,fland)
        CASE(510)
          IF (L_SW_Radiate) THEN   ! Only on radiation TS
! DEPENDS ON: ukca_set_array_bounds
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(land_albedo(I1:I2,J1:J2))
            CALL UKCA_EXTRACT_D1_DATA(                                 &
#include "argd1.h"
            first,N,land_albedo)
          ENDIF

        CASE DEFAULT
          cmessage='Item not found in prognostic case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',ukcaD1codes(N)%item,cmessage)
        END SELECT

! Diagnostics (section 0): fill appropriate array
      ELSEIF (UkcaD1codes(N)%section == 0 .AND.                        &
              .NOT. UkcaD1codes(N)%prognostic) THEN
        SELECT CASE(UkcaD1codes(N)%item)

        CASE(406)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(exner_theta_levels(I1:I2,J1:J2,                     &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,exner_theta_levels)
        CASE(407)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(p_rho_levels(I1:I2,J1:J2,                           &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,p_rho_levels)
        CASE(408)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(p_theta_levels(I1:I2,J1:J2,                         &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,p_theta_levels)
        CASE(409)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(pstar(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,pstar)
        CASE DEFAULT
          cmessage='N not found in diagnostic(0) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (section 1): fill appropriate array
      ELSEIF (UkcaD1codes(N)%section == 1) THEN
        SELECT CASE(UkcaD1codes(N)%item)

        CASE(201)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(net_surf_SW(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,net_surf_SW)
        CASE(235)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tot_surf_SW(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,tot_surf_SW)
        CASE DEFAULT
          cmessage='N not found in diagnostic(1) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

!     Diagnostics (section 3): fill appropriate array

      ELSEIF (UkcaD1codes(N)%section == 3) THEN
        SELECT CASE(UkcaD1codes(N)%item)
       CASE(25)
! DEPENDS ON: ukca_set_array_bounds
         CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
         ALLOCATE(ml_depth(I1:I2,J1:J2))
         CALL UKCA_EXTRACT_D1_DATA(                                    &
#include "argd1.h"
           first,N,ml_depth)
        CASE(60)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(rhokh_mix(I1:I2,J1:J2,                              &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,rhokh_mix)
        CASE(61)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(rho_aresist(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,rho_aresist)
        CASE(62)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(aresist(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,aresist)
        CASE(63)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(resist_b(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,resist_b)
        CASE(64)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dtrdz_charney_grid(I1:I2,J1:J2,                     &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dtrdz_charney_grid)
        CASE(65)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(kent(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,kent)
        CASE(66)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(we_lim(I1:I2,J1:J2,                                 &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,we_lim)
        CASE(67)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(t_frac(I1:I2,J1:J2,                                 &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,t_frac)
        CASE(68)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zrzi(I1:I2,J1:J2,                                   &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,zrzi)
        CASE(69)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(kent_dsc(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,kent_dsc)
        CASE(70)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(we_lim_dsc(I1:I2,J1:J2,                             &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,we_lim_dsc)
        CASE(71)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(t_frac_dsc(I1:I2,J1:J2,                             &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,t_frac_dsc)
        CASE(72)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zrzi_dsc(I1:I2,J1:J2,                               &
                             ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,zrzi_dsc)
        CASE(73)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(zhsc(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,zhsc)
        CASE(74:79)  ! r_b_dust for dust divisions
          IF (.NOT. ALLOCATED(rb_dust_div1)) THEN
! DEPENDS ON: ukca_set_array_bounds
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(rb_dust_div1(I1:I2,J1:J2))
            ALLOCATE(rb_dust_ndivs(I1:I2,J1:J2,ndiv))
           next_cd = 1
         ENDIF
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,rb_dust_div1)
           rb_dust_ndivs(:,:,next_cd) = rb_dust_div1(:,:)
           next_cd                    = next_cd + 1
        CASE(209)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(u_10m(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,u_10m)
        CASE(210)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(v_10m(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,v_10m)
        CASE(217)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(surf_hf(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,surf_hf)
        CASE(430)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(dust_ustar(I1:I2,J1:J2,UkcaD1Codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,dust_ustar)
        CASE(462)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(stcon(I1:I2,J1:J2,Ukcad1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,stcon)
        CASE(465)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(u_s(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,u_s)
        CASE DEFAULT
          cmessage='N not found in diagnostic(3) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (section 4): fill appropriate array
      ELSEIF (UkcaD1codes(N)%section == 4) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(205)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_liq_water(I1:I2,J1:J2,                        &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,cloud_liq_water)
        CASE(206)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(cloud_ice_content(I1:I2,J1:J2,                      &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,cloud_ice_content)
        CASE(222)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ls_rain3d(I1:I2,J1:J2,                              &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,ls_rain3d)
        CASE(223)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ls_snow3d(I1:I2,J1:J2,                              &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,ls_snow3d)
        CASE(227)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ls_ppn_frac(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,ls_ppn_frac)
        CASE DEFAULT
          cmessage='N not found in diagnostic(4) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',1,cmessage)
        END SELECT

! Diagnostics (section 5): fill appropriate array

      ELSEIF (UkcaD1codes(N)%section == 5) THEN
        SELECT CASE(UkcaD1codes(N)%item)
       CASE(227)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_rain3d(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,conv_rain3d)
        CASE(228)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(conv_snow3d(I1:I2,J1:J2,                            &
                              ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,conv_snow3d)
        CASE DEFAULT
          cmessage='N not found in diagnostic(5) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',1,cmessage)
        END SELECT

! Diagnostics (section 8): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 8) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(242)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(ch4_wetl_emiss(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,ch4_wetl_emiss)
        CASE DEFAULT
          cmessage='N not found in diagnostic(8) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',1,cmessage)
        END SELECT

! Diagnostics (section 15): fill appropriate array

      ELSE IF (UkcaD1codes(N)%section == 15) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(218)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(PV_on_theta_mlevs(I1:I2,J1:J2,                      &
                   ukcaD1codes(N)%len_dim3))
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,PV_on_theta_mlevs)
        CASE DEFAULT
          cmessage='N not found in diagnostic(8) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',1,cmessage)
       END SELECT

! Diagnostics (Section 30): fill appropriate array
      ELSEIF (UkcaD1codes(N)%section == 30) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(453)
! DEPENDS ON: ukca_set_array_bounds
          CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
          ALLOCATE(tropopause_height(I1:I2,J1:J2))
          CALL UKCA_EXTRACT_D1_DATA(                                    &
#include "argd1.h"
           first,N,tropopause_height)
        CASE DEFAULT
          cmessage='N not found in diagnostic(30) case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

! Diagnostics (Section 33): fill appropriate array

      ELSEIF (UkcaD1codes(N)%section == UKCA_sect) THEN
        SELECT CASE(UkcaD1codes(N)%item)
        CASE(251:300)                ! chemical diagnostics
          IF (.NOT. ALLOCATED(chem_diag1)) THEN
! DEPENDS ON: ukca_set_array_bounds
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(chem_diag1(I1:I2,J1:J2,                           &
                     ukcaD1codes(N)%len_dim3))
            ALLOCATE(chem_diags(I1:I2,J1:J2,                           &
                     ukcaD1codes(N)%len_dim3,n_chem_diags))
            next_cd = 1
          ENDIF
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,chem_diag1)
           chem_diags(:,:,:,next_cd) = chem_diag1(:,:,:)
           next_cd                   = next_cd + 1
        CASE(301:470)                ! BE flux diagnostics
          IF (.NOT. ALLOCATED(BE_fluxdiag1)) THEN
! DEPENDS ON: ukca_set_array_bounds
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(BE_fluxdiag1(I1:I2,J1:J2,                         &
                     ukcaD1codes(N)%len_dim3))
            ALLOCATE(BE_fluxdiags(I1:I2,J1:J2,                         &
                     ukcaD1codes(N)%len_dim3,n_BE_fluxdiags))
            next_cd = 1
          ENDIF
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,BE_fluxdiag1)
           BE_fluxdiags(:,:,:,next_cd) = BE_fluxdiag1(:,:,:)
           next_cd                   = next_cd + 1
        CASE(471:512)                ! stratospheric flux diagnostics
         IF (.NOT. ALLOCATED(strat_fluxdiag1)) THEN
           CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
           ALLOCATE(strat_fluxdiag1(I1:I2,J1:J2,                      &
                     ukcaD1codes(N)%len_dim3))
            ALLOCATE(strat_fluxdiags(I1:I2,J1:J2,                      &
                     ukcaD1codes(N)%len_dim3,n_strat_fluxdiags))
            next_cd = 1
          ENDIF
          CALL UKCA_EXTRACT_D1_DATA(                                   &
#include "argd1.h"
           first,N,strat_fluxdiag1)
           strat_fluxdiags(:,:,:,next_cd) = strat_fluxdiag1(:,:,:)
           next_cd                        = next_cd + 1
        CASE DEFAULT
          cmessage='N not found in ukca_section case statement'
! DEPENDS ON: ereport
          CALL EREPORT('GETD1FLDS',N,cmessage)
        END SELECT

      ELSE
        cmessage='N not located in IF statement'
! DEPENDS ON: ereport
        CALL EREPORT('GETD1FLDS',N,cmessage)
      ENDIF

      END SUBROUTINE GETD1FLDS
! ######################################################################
       SUBROUTINE PUTD1FLDS
! Put tracer and chemical diagnostic fields back into D1

       REAL, DIMENSION(:), ALLOCATABLE :: data1
       REAL, DIMENSION(:), ALLOCATABLE :: data2
       REAL, DIMENSION(:), ALLOCATABLE :: data3

! DEPENDS ON: ukca_set_array_bounds
       CALL UKCA_SET_ARRAY_BOUNDS(1,I1,I2,J1,J2)
       ALLOCATE(data1(ukcaD1codes(1)%length))

! DEPENDS ON: ukca_set_array_bounds
       CALL UKCA_SET_ARRAY_BOUNDS(idiag_first,I1,I2,J1,J2)
       ALLOCATE(data2(ukcaD1codes(idiag_first)%length))

       DO N=1,n_use_tracers
         IF ((ukcaD1codes(N)%item .NE. 10) .OR.                        &
             (ukcaD1codes(N)%section .NE. 0)) THEN
           tracer1(:,:,:)=all_tracers(:,:,:,N)
           data1=RESHAPE(tracer1,(/SIZE(data1)/))

!        Update D1
           D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+           &
              ukcaD1codes(N)%length-1)=data1
         ENDIF

       ENDDO   ! N tracer loop

       DO N=idiag_first,idiag_last
         chem_diag1(:,:,:)=chem_diags(:,:,:,N-idiag_first+1)
         data2=RESHAPE(chem_diag1,(/SIZE(data2)/))

!        Update D1
         D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+             &
            ukcaD1codes(N)%length-1)=data2

       ENDDO   ! N chemical diagnostics loop

       IF (L_ukca_h2o_feedback) THEN
! Copy water vapour back into D1 array.
         N = size(q)
         ALLOCATE(data3(N))
         data3=RESHAPE(q,(/N/))

! Update D1
         DO N=n_use_tracers+n_use_emissions+1,                         &
              n_use_tracers+n_use_emissions+n_in_progs
           IF ((ukcaD1codes(N)%item  == 10) .AND.                      &
               (ukcaD1codes(N)%section == 0 ))                         &
             D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+         &
               ukcaD1codes(N)%length-1)=data3
         END DO

         DEALLOCATE(data3)
       ENDIF

       IF (L_ukca_BEflux) THEN

! DEPENDS ON: ukca_set_array_bounds
         CALL UKCA_SET_ARRAY_BOUNDS(iflux_first,I1,I2,J1,J2)
         ALLOCATE(data3(ukcaD1codes(iflux_first)%length))

         DO N=iflux_first,iflux_last
           BE_fluxdiag1(:,:,:)=BE_fluxdiags(:,:,:,N-iflux_first+1)
           data3=RESHAPE(BE_fluxdiag1,(/SIZE(data3)/))

!          Update D1
           D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address+           &
            ukcaD1codes(N)%length-1)=data3

         ENDDO   ! N BE flux diagnostics loop

        DEALLOCATE(data3)
        DEALLOCATE(BE_fluxdiag1)
        DEALLOCATE(BE_fluxdiags)

       ENDIF

       IF (L_ukca_stratflux) THEN

         icnt = 0
        DO N=istrat_first,istrat_last
          IF (UkcaD1codes(N)%required ) THEN
            CALL UKCA_SET_ARRAY_BOUNDS(N,I1,I2,J1,J2)
            ALLOCATE(data3(ukcaD1codes(N)%length))
            icnt = icnt+1
            strat_fluxdiag1(:,:,:)=strat_fluxdiags(:,:,:,icnt)
            data3=RESHAPE(strat_fluxdiag1,(/SIZE(data3)/))

!            Update D1 if diagnostic is required
             D1(ukcaD1codes(N)%address:ukcaD1codes(N)%address +        &
                ukcaD1codes(N)%length-1)=data3
             DEALLOCATE(data3)
           ENDIF
         ENDDO   ! N stratospheric flux diagnostics loop

        IF (ALLOCATED(strat_fluxdiag1)) DEALLOCATE(strat_fluxdiag1)
        IF (ALLOCATED(strat_fluxdiags)) DEALLOCATE(strat_fluxdiags)

       ENDIF     ! End of IF L_ukca_stratflux statement

       DEALLOCATE(data1)
       DEALLOCATE(data2)
       DEALLOCATE(tracer1)
       DEALLOCATE(emission1)
       DEALLOCATE(chem_diag1)
       DEALLOCATE(all_tracers)
       DEALLOCATE(all_emissions)
       DEALLOCATE(chem_diags)

       END SUBROUTINE PUTD1FLDS
! ######################################################################

       END SUBROUTINE UKCA_MAIN1
#endif
