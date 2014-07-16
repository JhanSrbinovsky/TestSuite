#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  To introduce a new deck (REDIST_STOCHEM, contains STOCH_SHARE
!  and ARCHDUMP)
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.4    28/01/99  Created.  C.E. Johnson
!  5.1    20/03/00  Module definition removed. C.E. Johnson
!  5.2    25/04/00  NRUN from dumps via logical StochStartDump added.
!                   C.E. Johnson.
!  5.2    23/11/00  Convective diagnostics added, 3D clouds.
!                   C.E. Johnson
!  5.2    20/07/01  New Dynamics version
!  5.3    27/11/01  Reinstated TROP; CB,CT arrays removed.
!  5.5    26/04/04  Added logicals and arrays to feed STOCHEM ch4 and
!                   ozone concentrations to UM radiation scheme.
!                   C.E. Johnson
!  6.1    20/08/04  Introduce STOCHEM logicals.  C.E. Johnson
!  6.2    01/03/06  Added new variables for dry deposition and
!                   biogenic VOC emission scheme. M.G. Sanderson.
!

! ######################################################################
!+ To redistribute arrays from D1 onto the STOCHEM processor grid
!+ and then call STOCHEM.
!
! Subroutine Interface:
      SUBROUTINE REDIST_STOCHEM(                                        &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "argcona.h"
     &   l_use_stochem_ch4, l_use_stochem_o3,                           &
     &   i_dummy)

      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_CHM
      USE IN_STOCHEM_OUT
      USE IN_STOCHEM_INTF
      USE LEVEL_HEIGHTS_MOD
      IMPLICIT NONE
!
! Description:
!   1) At the start of a CRUN, read in previous dump of data and
!      redistribute this to the STOCHEM processor decomposition as
!      the first set of data.
!   2) To gather fields from D1 required to run STOCHEM on one PE and
!      redistribute onto the STOCHEM processor decomposition as the
!      second set of data.
!   3) To write a data dump as required, 5 steps before a UM dump.
!   4) To call STOCHEM.
!
! Method:
!   Prognostic fields are selected from D1 using the item codes
!   defined in the array StochProgCodes.
!   Diagnostic fields are selected by tagging those required in
!   the UMUI.  There is a check on the total number of fields
!   processed - see NstochFields.  If the wrong fields are
!   tagged, they will not be identified by STOCH_SHARE and an
!   Error message will result.  If a field is not tagged, then
!   the wrong number of fields will be processed which will also
!   generate an error.
!   Each field is gathered on PE 0 using GENERAL_GATHER_FIELD,
!   then redistributed and renamed using STOCH_SHARE.
!   Note that fields on the U grid are padded by a row of zeros
!   before redistribution.
!
!
! Code Description:
!   Language: FORTRAN 90 + common extensions.
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):

#include "typsize.h"
#include "typd1.h"
#include "typduma.h"
#include "csubmodl.h"
#include "typsts.h"
#include "typcona.h"
!#include <model/model.h>        ! for logical 3dozone

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN) :: i_dummy        ! not used
      LOGICAL, INTENT(INOUT) :: l_use_stochem_ch4
      LOGICAL, INTENT(INOUT) :: l_use_stochem_o3

!   ErrorStatus

      INTEGER       :: icode=0      ! Error flag (0 = OK)
      CHARACTER(72) :: cmessage     ! Error return message

! Local parameters:

! Number of prognostic and diagnostic stashcodes used
      INTEGER, PARAMETER :: nstochprogs=13
      INTEGER, PARAMETER :: nstochdiags=27
      INTEGER, PARAMETER :: nstochinprogs=2

! Local scalars:
      INTEGER :: um_month          ! Month for STOCHEM
      INTEGER :: um_year           ! Year for STOCHEM
      INTEGER :: a_step            ! UM step number
      INTEGER :: d_step            ! UM dump frequency
      INTEGER :: nfields           ! fields counter
      INTEGER :: nstochfields      ! No. fields requested
      INTEGER :: length            ! local dimension
      INTEGER :: xlen              ! local dimension
      INTEGER :: address           ! address in D1
      INTEGER :: levels            ! number of levels
      INTEGER :: stashcode         ! stash code
      INTEGER :: section           ! stash section
      INTEGER :: item              ! stash item
      INTEGER :: ftype             ! Field type
      INTEGER :: halo_type         ! halo type
      INTEGER :: grid_type         ! grid type
      INTEGER :: tag               ! stash tag
      INTEGER :: ptd1              ! pointer to D1
      INTEGER :: i                 ! loop variable
      INTEGER :: j                 ! loop variable
      INTEGER :: k
      INTEGER :: m                 ! sub model
      INTEGER :: ierr              ! return code
      INTEGER :: GET_FLD_TYPE      ! UM function

      REAL    :: um_day            ! Day for STOCHEM

      CHARACTER(14) :: filename    ! name of dump file

      LOGICAL, SAVE :: first=.true.              ! true only on first ca
      LOGICAL, SAVE :: stochdumprequest=.false.  !
      LOGICAL, SAVE :: archdumprequest=.false.   ! true to archive a dum

      TYPE CODE
       INTEGER :: icode       ! stash code
       INTEGER :: field_type  ! field grid type
       INTEGER :: n_levels    ! number of levels
      ENDTYPE CODE

      TYPE(CODE), DIMENSION(nstochprogs)   :: stochprogcodes
      TYPE(CODE), DIMENSION(nstochdiags)   :: stochdiagcodes
      TYPE(CODE), DIMENSION(nstochinprogs) :: stochinprogcodes

! Met arrays
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: cc_base_level
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: cc_top_level

      REAL, DIMENSION(:,:,:), ALLOCATABLE :: t      ! temperature
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: q      ! spec humidity
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: clw    ! cloud liq water
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: up     ! conv updraft massf
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: upe    ! conv updraft entra
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: upd    ! conv updraft detra
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: cc     ! conv cloud cover
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: dc     ! dynm cloud cover
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p_th   ! pressure on theta
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: o3um   ! UM ozone levels
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: p      ! pressure on theta
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: lnp    ! log pressure theta
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: w      ! vertical vel (m/s)
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: u      ! northward vel (m/s
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: v      ! eastward vel (m/s)
      REAL, DIMENSION(:,:),   ALLOCATABLE :: p0     ! surface pressure
      REAL, DIMENSION(:,:),   ALLOCATABLE :: t0     ! surface temperatur
      REAL, DIMENSION(:,:),   ALLOCATABLE :: trop   ! tropopause pressur
      REAL, DIMENSION(:,:),   ALLOCATABLE :: tropz  ! tropopause height
      REAL, DIMENSION(:,:),   ALLOCATABLE :: bl     ! boundary layer hei
      REAL, DIMENSION(:,:),   ALLOCATABLE :: sice   ! sea ice fraction
      REAL, DIMENSION(:,:),   ALLOCATABLE :: snowam ! snow amount
      REAL, DIMENSION(:,:),   ALLOCATABLE :: acp    ! conv precip
      REAL, DIMENSION(:,:),   ALLOCATABLE :: adp    ! dynm precip
      REAL, DIMENSION(:,:),   ALLOCATABLE :: tauu   ! surface stress u
      REAL, DIMENSION(:,:),   ALLOCATABLE :: tauv   ! surface stress v
      REAL, DIMENSION(:,:),   ALLOCATABLE :: hf     ! heat flux
      REAL, DIMENSION(:,:),   ALLOCATABLE :: acr    ! conv rain
      REAL, DIMENSION(:,:),   ALLOCATABLE :: acs    ! conv snow
      REAL, DIMENSION(:,:),   ALLOCATABLE :: adr    ! dynm rain
      REAL, DIMENSION(:,:),   ALLOCATABLE :: ads    ! dynm snow
      REAL, DIMENSION(:,:),   ALLOCATABLE :: orog   ! orography height
      REAL, DIMENSION(:,:),   ALLOCATABLE :: land_fraction ! land fracti
      REAL, DIMENSION(:,:),   ALLOCATABLE :: totpar ! Total surface PAR
      REAL, DIMENSION(:,:),   ALLOCATABLE :: dirpar ! Direct component o
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: soilmc ! Soil moisture cont
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: t0tile ! Surface temperatur
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: gsf    ! Surface tile fract
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: lai_ft ! LAI on vegetated t
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: canht  ! Canopy ht on veg t
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: canwc  ! Canopy water on ve
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: gc     ! Stomatal Cond on t
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: z0     ! Roughness length o
      REAL, DIMENSION(:,:),   ALLOCATABLE :: sicetemp ! Sea ice temperat
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ra     ! aerodynamic resist
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: rq     ! quasi-laminar resi
      REAL, DIMENSION(:,:),   ALLOCATABLE :: so4_vd ! SO4 aerosol dry de

! 3D arrays for feedback of STOCHEM CH4 and O3 to UM via the D1 array
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: o3mmr
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: ch4mmr

! Temporary fields for retrieving data from D1
      INTEGER, DIMENSION(:), ALLOCATABLE :: ifield
      REAL,    DIMENSION(:), ALLOCATABLE :: rfield
      REAL,    DIMENSION(:), ALLOCATABLE :: zfield

!- End of header

! Common Blocks

#include "decomptp.h"
#include "decompdb.h"
#include "ctime.h"
#include "stparam.h"
#include "cntlgen.h"
#include "c_mdi.h"
#include "c_soilh.h"

      IF (first) THEN
! DEPENDS ON: inigri
        CALL INIGRI(eta_theta_levels,eta_rho_levels,                    &
                                                     ! to set grids
     &    a_realhd(rh_tot_mass_init))
! DEPENDS ON: inibounds2
        CALL INIBOUNDS2  ! Sets up bounds array etc in module
      END IF

! set logicals for coupling - now passed in through CNTLATM namelist
!     l_use_stochem_ch4 = l_coupled_ch4
!     l_use_stochem_o3 = l_coupled_o3

! Met arrays
      ALLOCATE(cc_base_level(nlonpe,nlatpe))
      ALLOCATE(cc_top_level(nlonpe,nlatpe))
      ALLOCATE(    t(nlonpe,nlatpe,nmetlev))
      ALLOCATE(    q(nlonpe,nlatpe,nmetlev))
      ALLOCATE(  clw(nlonpe,nlatpe,nmetlev))
      ALLOCATE(   up(nlonpe,nlatpe,nmetlev))
      ALLOCATE(  upe(nlonpe,nlatpe,nmetlev))
      ALLOCATE(  upd(nlonpe,nlatpe,nmetlev))
      ALLOCATE(   cc(nlonpe,nlatpe,nmetlev))
      ALLOCATE(   dc(nlonpe,nlatpe,nmetlev))
      ALLOCATE( p_th(nlonpe,nlatpe,nmetlev))
      ALLOCATE( o3um(nlonpe,nlatpe,nmetlev))
      ALLOCATE(    p(nlonpe,nlatpe,0:nmetlev))
      ALLOCATE(  lnp(nlonpe,nlatpe,0:nmetlev))
      ALLOCATE(    w(nlonpe,nlatpe,nmetlev+1))
      ALLOCATE(    u(nlonpe,nlatpe,0:nmetlev))
      ALLOCATE(    v(nlonpe,nlatpe,0:nmetlev))
      ALLOCATE(   p0(nlonpe,nlatpe))
      ALLOCATE(   t0(nlonpe,nlatpe))
      ALLOCATE( trop(nlonpe,nlatpe))
      ALLOCATE(tropz(nlonpe,nlatpe))
      ALLOCATE(   bl(nlonpe,nlatpe))
      ALLOCATE(sice(nlonpe,nlatpe))
      ALLOCATE(snowam(nlonpe,nlatpe))
      ALLOCATE(  acp(nlonpe,nlatpe))
      ALLOCATE(  adp(nlonpe,nlatpe))
      ALLOCATE( tauu(nlonpe,nlatpe))
      ALLOCATE( tauv(nlonpe,nlatpe))
      ALLOCATE(   hf(nlonpe,nlatpe))
      ALLOCATE(  acr(nlonpe,nlatpe))
      ALLOCATE(  acs(nlonpe,nlatpe))
      ALLOCATE(  adr(nlonpe,nlatpe))
      ALLOCATE(  ads(nlonpe,nlatpe))
      ALLOCATE( orog(nlonpe,nlatpe))
      ALLOCATE(land_fraction(nlonpe,nlatpe))
      ALLOCATE(totpar(nlonpe,nlatpe))
      ALLOCATE(dirpar(nlonpe,nlatpe))
      ALLOCATE(t0tile(nlonpe,nlatpe,ntype))
      ALLOCATE(gsf(nlonpe,nlatpe,ntype))
      ALLOCATE(lai_ft(nlonpe,nlatpe,npft))
      ALLOCATE(canht(nlonpe,nlatpe,npft))
      ALLOCATE(canwc(nlonpe,nlatpe,ntype))
      ALLOCATE(gc(nlonpe,nlatpe,npft))
      ALLOCATE(z0(nlonpe,nlatpe,ntype))
      ALLOCATE(sicetemp(nlonpe,nlatpe))
      ALLOCATE(soilmc(nlonpe,nlatpe,sm_levels))
      ALLOCATE(ra(nlonpe,nlatpe,ntype))
      ALLOCATE(rq(nlonpe,nlatpe,nc))
      ALLOCATE(so4_vd(nlonpe,nlatpe))

! STOCHEM arrays
      ALLOCATE(ch4mmr(nlonpe,nlatpe,nmetlev))
      ALLOCATE(o3mmr(nlonpe,nlatpe,nmetlev))

! Specify the stash codes of prognostics, field type and
! whether they are 3D fields.
      stochprogcodes(1)=code(256,fld_type_u,nmetlev)     ! u
      stochprogcodes(2)=code(257,fld_type_v,nmetlev)     ! v
      stochprogcodes(3)=code(10,fld_type_p,nmetlev)      ! q
      stochprogcodes(4)=code(14,fld_type_p,1)            ! cc_base_level
      stochprogcodes(5)=code(15,fld_type_p,1)            ! cc_top_level
      stochprogcodes(6)=code(23,fld_type_p,1)            ! snowam
      stochprogcodes(7)=code(24,fld_type_p,1)            ! t0
      stochprogcodes(8)=code(25,fld_type_p,1)            ! bl
      stochprogcodes(9)=code(31,fld_type_p,1)            ! sice
      stochprogcodes(10)=code(49,fld_type_p,1)           ! sicetemp
      stochprogcodes(11)=code(60,fld_type_p,nmetlev)     ! o3um (may be
      stochprogcodes(12)=code(258,fld_type_p,nmetlev+1)  ! w
      stochprogcodes(13)=code(265,fld_type_p,nmetlev)    ! dc

! Specify the stash codes of diagnostics, field type and
! whether they are 3D fields
! n.b. 205,408 & 409 are in section 0 but are diagnostics.
      stochdiagcodes(1)=code(408,fld_type_p,nmetlev)     ! P_TH
      stochdiagcodes(2)=code(409,fld_type_p,1)           ! P0
      stochdiagcodes(3)=code(1290,fld_type_p,1)          ! TOTPAR
      stochdiagcodes(4)=code(1291,fld_type_p,1)          ! DIRPAR
      stochdiagcodes(5)=code(3217,fld_type_p,1)          ! HF
      stochdiagcodes(6)=code(3219,fld_type_p,1)          ! TAUU
      stochdiagcodes(7)=code(3220,fld_type_p,1)          ! TAUV
      stochdiagcodes(8)=code(3316,fld_type_p,ntype)      ! T0TILE
      stochdiagcodes(9)=code(3317,fld_type_p,ntype)      ! GSF
      stochdiagcodes(10)=code(3318,fld_type_p,npft)      ! LAI_FT
      stochdiagcodes(11)=code(3319,fld_type_p,npft)      ! CANHT
      stochdiagcodes(12)=code(3321,fld_type_p,ntype)     ! CANWC
      stochdiagcodes(13)=code(3324,fld_type_p,ntype)     ! Z0
      stochdiagcodes(14)=code(3395,fld_type_p,1)         ! LAND_FRACTION
      stochdiagcodes(15)=code(3462,fld_type_p,npft)      ! GC
      stochdiagcodes(16)=code(4203,fld_type_p,1)         ! ADR
      stochdiagcodes(17)=code(4204,fld_type_p,1)         ! ADS
      stochdiagcodes(18)=code(4205,fld_type_p,nmetlev)   ! CLW
      stochdiagcodes(19)=code(5205,fld_type_p,1)         ! ACR
      stochdiagcodes(20)=code(5206,fld_type_p,1)         ! ACS
      stochdiagcodes(21)=code(5209,fld_type_p,nmetlev)   ! T
      stochdiagcodes(22)=code(5212,fld_type_p,nmetlev)   ! CC
      stochdiagcodes(23)=code(5250,fld_type_p,nmetlev-1) ! UP
      stochdiagcodes(24)=code(5253,fld_type_p,nmetlev-1) ! UPD
      stochdiagcodes(25)=code(8223,fld_type_p,sm_levels) ! SOILMC
      stochdiagcodes(26)=code(30451,fld_type_p,1)        ! TROP
      stochdiagcodes(27)=code(30453,fld_type_p,1)        ! TROPZ

! These are the zonal mean fields that are put into D1
      stochinprogcodes(1)=code(99,fld_type_p,nmetlev)    ! CH4
      stochinprogcodes(2)=code(100,fld_type_p,nmetlev)   ! O3

! Calculate the total number of fields requested
      nstochfields = SUM(stochprogcodes%n_levels) +                     &
     &  SUM(stochdiagcodes%n_levels)
      WRITE(6,*) 'NstochFields Requested: ',nstochfields

! The diagnostics required by STOCHEM are specified in the UMUI
! and are sent to secondary storage with a Tag of 99.

!      STOCHEM inputs are (STASH Section, Item, Stochem_Name):

!  Y   0   2   U      U compnt of wind after timestep
!  Y   0   3   V      V compnt of wind after timestep
!  N   0   9   SM     Soil Moisture in a layer
!  Y   0  10   Q      Specific humidity after timestep
!  Y   0  14   CC_base_level   CC base level no.
!  Y   0  15   CC_top_level    CC top level no.
!  N   0  20   ST     Deep soil temperature after timestep
!  Y   0  23   SNOWAM Snow amount after timestep
!  Y   0  24   T0     Surface temperature after timestep
!  Y   0  25   BL     Boundary layer depth after timestep
!  Y   0  31   SICE   Sea ice fraction after timestep
!  Y   0  33   OROG   Orography
!  Y   0  49   SICETEMP   Sea ice temperature
!  Y   0  60   O3UM   UM Ozone
!  Y   0 150   W      W compnt of wind after TS
!  Y   0 205   FLC    Fractional land cover
!  Y   0 265   DC     Area cloud fraction in layers
!  N   0 407   P_RH   P at Rho levels after TS
!  Y   0 408   P_TH   P at Theta levels after TS
!  Y   0 409   P0     Surface Pressure after TS
!  Y   1 290   TOTPAR Total PAR flux at surface (W/m2)
!  Y   1 291   DIRPAR Direct component of total PAR flux (W/m2)
!  Y   3 217   HF     Surface and boundary layer heat fluxes
!  Y   3 219   TAUU   X-comp of surf & bl wind stress N/m2
!  Y   3 220   TAUV   Y-comp of surf & bl wind stress N/m2
!  Y   3 316   T0TILE Surface temperature on tiles (K)
!  Y   3 317   GSF    Surface tile fractions
!  Y   3 318   LAI_FT Leaf area indices on vegetated tiles
!  Y   3 319   CANHT  Canopy height on vegetated tiles
!  Y   3 321   CANWC  Canopy water content on vegetated tiles
!  Y   3 324   Z0     Rougness length on tiles (m)
!  Y   3 395   LAND_FRACTION Land fraction (0-1)
!  Y   3 462   GC     Stomatal conductance (m s-1)
!  Y   4 203   ADR    Large scale rainfall rate kg/m2/s
!  Y   4 204   ADS    Large scale snowfall rate kg/m2/s
!  Y   4 205   CLW    Cloud liquid water after LS precip
!  Y   5 205   ACR    Convective rainfall rate
!  Y   5 206   ACS    Convective snowfall rate
!  Y   5 209   T      Temperature after convection
!  Y   5 212   CC     Conv cloud amount on each model level
!  Y   5 250   UP     Updraught mass flux (Pa/s)
!  Y   5 254   UPD    Updraught detrainment rate (Pa/s)
!  Y   8 223   SOILMC Soil moisture content in a layer (mm)
!  Y  30 451   TROP   Pressure at tropopause level
!  Y  30 453   TROPZ  Tropopause height

! Inibounds2 calculates the processor map for STOCHEM and the
! limits for each PE, depending on the number of PEs.

      m = submodel_for_sm(a_im)
      a_step = stepim(a_im)
      d_step = dumpfreqim(a_im)

! Dumps are required before a UM dump.
      IF (MOD(a_step,d_step) == 0) StochDumpRequest=.true.

! Dumps are archived when:
      IF (archdump_offsetim(a_im) /= 0 .AND.                            &
     &  archdump_freqim(a_im) /= 0) THEN
        IF (a_step/d_step  >=  archdump_offsetim(a_im) .AND.            &
     &    MOD(a_step/d_step,archdump_freqim(a_im)) == 0)                &
     &    archdumprequest=.true.
!       WRITE(6,*) 'ArchDumpRequest: ',ArchDumpRequest
!       WRITE(6,*) 'D_STEP: ',d_step
!       WRITE(6,*) 'ARCHDUMP_OFFSET: ',archdump_offsetim(a_im)
!       WRITE(6,*) 'ARCHDUMP_FREQim: ',archdump_freqim(a_im)
      END IF

      um_year = i_year
      um_month = i_month
! STOCHEM takes fractional day
      um_day = REAL(i_day) + REAL(i_hour)/24.0 +                        &
     &  REAL(i_minute)/1440.0 + REAL(i_second)/daysec
! DEPENDS ON: set_daym
      CALL SET_DAYM(daym,um_year)
! Stochem is catching up the UM time so subtract astep from time etc.
      um_day = um_day-secs_per_stepim(a_im)/daysec
      IF (um_day < 1.0) THEN
        um_month = um_month -1
        IF (um_month < 1) THEN
          um_year = um_year -1
          um_month = 12
        END IF
        um_day = um_day + REAL(daym(um_month))
      END IF

! Set fields counter. This variable sums up the number of fields
! read in from the D1 array, and must match nstochfields
      nfields = 0

! Set first level of U and V to zero
      u(:,:,0) = 0.0
      v(:,:,0) = 0.0

! Prognostics loop
! Loop through all the objects in D1
      DO i=1,no_obj_D1(m)
        item = D1_ADDR(D1_item,i,m)
        IF (ANY(StochProgCodes%icode == item) .AND.                     &
     &    d1_addr(d1_object_type,i,m) == prognostic) THEN
          levels = d1_addr(d1_no_levels,i,m)
          length = d1_addr(d1_length,i,m)
          address = d1_addr(d1_address,i,m)
          halo_type = d1_addr(d1_halo_type,i,m)
          grid_type = D1_ADDR(d1_grid_type,i,m)
! DEPENDS ON: get_fld_type
          ftype = get_fld_type(grid_type)

!         WRITE(6,*) 'Prognostic:'
!         WRITE(6,*) ' item:        ',d1_addr(d1_item,i,m)
!         WRITE(6,*) ' object_type: ',d1_addr(d1_object_type,i,m)
!         WRITE(6,*) ' grid_type:   ',d1_addr(d1_grid_type,i,m)
!         WRITE(6,*) ' halo_type:   ',d1_addr(d1_halo_type,i,m)
!         WRITE(6,*) ' field_type:  ',ftype
!         WRITE(6,*) ' no_levels:   ',d1_addr(d1_no_levels,i,m)
!         WRITE(6,*) ' length:      ',d1_addr(d1_length,i,m)

          nfields = nfields + levels

          SELECT CASE(item)
          CASE(14,15) ! cc_base_level, cc_top_level
            ALLOCATE(rfield(length))
            ALLOCATE(ifield(length))
            ifield = id1(address:address+length-1)
            rfield = REAL(ifield)
            DEALLOCATE(ifield)
          CASE(60) ! Zonal mean ozone
            xlen = lasize(1,ftype,halo_type)
            ALLOCATE(rfield(length*xlen))
            ALLOCATE(zfield(length))
            zfield = d1(address:address+length-1)
! Spread out Zfield over all longitudes
            DO j=1,length
              rfield((j-1)*xlen+1:j*xlen) = zfield(j)
            END DO
            DEALLOCATE(zfield)
            length = length * xlen      ! redefine for rfield
          CASE DEFAULT
            ALLOCATE(rfield(length))
            rfield = d1(address:address+length-1)
          END SELECT
          CALL STOCH_SHARE(rfield,length,levels,item,ftype,halo_type)
          DEALLOCATE(rfield)
        END IF        ! item and object type match
      END DO          ! I no_obj_D1(M)
      orog = rmdi
      orog(1:lasize(1,fld_type_p,halo_type_extended),                   &
     &  1:lasize(2,fld_type_p,halo_type_extended)) =                    &
     &  r_theta_levels(:,:,0)-earth_radius

      WRITE(6,*) 'U   : ',mype,u(1,1,1), MAXVAL(u),MINVAL(u)
      WRITE(6,*) 'V   : ',mype,v(1,1,1), MAXVAL(v),MINVAL(v)
      WRITE(6,*) 'Q   : ',mype,q(1,1,1), MAXVAL(q),MINVAL(q)
      WRITE(6,*) 'CC_base_level: ',mype,cc_base_level(1,1),             &
     &  MAXVAL(cc_base_level),MINVAL(cc_base_level)
      WRITE(6,*) 'CC_top_level: ',mype,cc_top_level(1,1),               &
     &  MAXVAL(cc_top_level),MINVAL(cc_top_level)
      WRITE(6,*) 'T0  : ',mype,t0(1,1),  MAXVAL(t0),MINVAL(t0)
      WRITE(6,*) 'BL  : ',mype,bl(1,1),  MAXVAL(bl),MINVAL(bl)
      WRITE(6,*) 'O3UM: ',mype,o3um(1,1,1),MAXVAL(o3um),MINVAL(o3um)
      WRITE(6,*) 'SNOWAM : ',mype,snowam(1,1),  MAXVAL(snowam),         &
     &  MINVAL(snowam)
      WRITE(6,*) 'SICE: ',mype,sice(1,1),MAXVAL(sice),MINVAL(sice)
      WRITE(6,*) 'SICETEMP: ',mype,sicetemp(1,1),MAXVAL(sicetemp),      &
     &  MINVAL(sicetemp)
      WRITE(6,*) 'W   : ',mype,w(1,1,1), MAXVAL(w),MINVAL(w)
      WRITE(6,*) 'DC  : ',mype,dc(1,1,1),MAXVAL(dc),MINVAL(dc)
      WRITE(6,*) 'OROG: ',mype,orog(1,1),MAXVAL(orog),MINVAL(orog)

! Diagnostics loop through all stashlist items
      DO j=1,totitems
        tag = STLIST(st_macrotag,j)-(1000*(STLIST(st_macrotag,j)/1000))
! Check if tag is set.
        IF (tag == 99) THEN
          ptd1 = STLIST(st_D1pos,j)
          section=  d1_addr(d1_section,  ptd1,m)
          item=     d1_addr(d1_item,     ptd1,m)
          levels=   d1_addr(d1_no_levels,ptd1,m)
          length=   d1_addr(d1_length,   ptd1,m)
          address=  d1_addr(d1_address,  ptd1,m)
          halo_type=d1_addr(d1_halo_type,ptd1,m)
          grid_type=d1_addr(d1_grid_type,ptd1,m)
! DEPENDS ON: get_fld_type
          ftype = get_fld_type(grid_type)
          stashcode = section*1000 + item

!         WRITE(6,*) 'Diagnostic: Stashcode = ',stashcode
!         WRITE(6,*) ' item:        ',d1_addr(d1_item,ptd1,m)
!         WRITE(6,*) ' object_type: ',d1_addr(d1_object_type,ptd1,m)
!         WRITE(6,*) ' grid_type:   ',d1_addr(d1_grid_type,ptd1,m)
!         WRITE(6,*) ' field_type:  ',ftype
!         WRITE(6,*) ' no_levels:   ',d1_addr(d1_no_levels,ptd1,m)
!         WRITE(6,*) ' length:      ',d1_addr(d1_length,ptd1,m)

          nfields = nfields + levels
          CALL STOCH_SHARE(d1(address:address+length-1),                &
     &      length,levels,stashcode,ftype,halo_type)

        END IF      ! Tag selection
      END DO        ! J Totitems in stashlist

! Set last level of convective diagnostics to zero
      up(:,:,nmetlev)  = 0.0
      upd(:,:,nmetlev) = 0.0
      WRITE(6,*) 'P_TH: ',mype,p_th(1,1,1),MAXVAL(p_th),MINVAL(p_th)
      WRITE(6,*) 'P0  : ',mype,p0(1,1),MAXVAL(p0),MINVAL(p0)
      WRITE(6,*) 'T   : ',mype,t(1,1,1),MAXVAL(t),MINVAL(t)
      WRITE(6,*) 'CC  : ',mype,cc(1,1,1),MAXVAL(cc),MINVAL(cc)
      WRITE(6,*) 'UP  : ',mype,up(1,1,1),MAXVAL(up),MINVAL(up)
      WRITE(6,*) 'UPD : ',mype,upd(1,1,1),MAXVAL(upd),MINVAL(upd)
      WRITE(6,*) 'ACR : ',mype,acr(1,1),MAXVAL(acr),MINVAL(acr)
      WRITE(6,*) 'ACS : ',mype,acs(1,1),MAXVAL(acs),MINVAL(acs)
      WRITE(6,*) 'ADR : ',mype,adr(1,1),MAXVAL(adr),MINVAL(adr)
      WRITE(6,*) 'ADS : ',mype,ads(1,1),MAXVAL(ads),MINVAL(ads)
      WRITE(6,*) 'TAUU: ',mype,tauu(1,1),MAXVAL(tauu),MINVAL(tauu)
      WRITE(6,*) 'TAUV: ',mype,tauv(1,1),MAXVAL(tauv),MINVAL(tauv)
      WRITE(6,*) 'HF  : ',mype,hf(1,1),MAXVAL(hf),MINVAL(hf)
      WRITE(6,*) 'TROP: ',mype,trop(1,1),MAXVAL(trop),MINVAL(trop)
      WRITE(6,*) 'TROPZ: ',mype,tropz(1,1),MAXVAL(tropz),MINVAL(tropz)
      WRITE(6,*) 'TOTPAR : ',mype,totpar(1,1),MAXVAL(totpar),           &
     &  MINVAL(totpar)
      WRITE(6,*) 'DIRPAR : ',mype,dirpar(1,1),MAXVAL(dirpar),           &
     &  MINVAL(dirpar)
      WRITE(6,*) 'T0TILE : ',mype,t0tile(1,1,1),MAXVAL(t0tile),         &
     &  MINVAL(t0tile)
      WRITE(6,*) 'GSF : ',mype,gsf(1,1,1),MAXVAL(gsf),MINVAL(gsf)
      WRITE(6,*) 'LAI_FT : ',mype,lai_ft(1,1,1),MAXVAL(lai_ft),         &
     &  MINVAL(lai_ft)
      WRITE(6,*) 'CANHT : ',mype,canht(1,1,1),MAXVAL(canht),            &
     &  MINVAL(canht)
      WRITE(6,*) 'CANWC : ',mype,canwc(1,1,1),MAXVAL(canwc),            &
     &  MINVAL(canwc)
      WRITE(6,*) 'SOILMC : ',mype,soilmc(1,1,1),MAXVAL(soilmc),         &
     &  MINVAL(soilmc)
      WRITE(6,*) 'GC : ',mype,gc(1,1,1),MAXVAL(gc),MINVAL(gc)
      WRITE(6,*) 'Z0  : ',mype,z0(1,1,1),MAXVAL(z0),MINVAL(z0)
      WRITE(6,*) 'Land_fraction: ',mype,land_fraction(1,1),             &
     &  MAXVAL(land_fraction),MINVAL(land_fraction)

! Check that the number of fields processed is correct.
      IF (nfields /= nstochfields) THEN
        cmessage = 'REDIST_STOCHEM: Wrong Number of Fields '
        WRITE(6,*) 'Wrong Number of Fields in REDIST_STOCH'
        WRITE(6,*) 'REDIST_STOCH is expecting ',nstochfields,' fields'
        WRITE(6,*) 'but has found ',nfields,' fields.'
! DEPENDS ON: ereport
        CALL EREPORT('Redist_Stochem',1,cmessage)
      END IF

! DEPENDS ON: cpustats
      CALL CPUSTATS('PROCESS_MET ',1)
! DEPENDS ON: process_met
      CALL PROCESS_MET(p0,p_th,t0,tauu,tauv,hf,adr,ads,bl,              &
     &  land_fraction,cc,acr,acs,sice,snowam,orog,gsf,lai_ft,canht,z0,  &
     &  gc,o3um,t0tile,sicetemp,soilmc,p,lnp,acp,adp,ra,rq,so4_vd)
! DEPENDS ON: cpustats
      CALL CPUSTATS('PROCESS_MET ',2)

!     WRITE(6,*) 'mype=',mype
!     WRITE(6,*) 'halo_i=',halo_i,'  halo_j',halo_j
!     WRITE(6,*) 'nmetlong=',nmetlong,'  nmetlat=',nmetlat
!     WRITE(6,*) 'nlonpe=',nlonpe,'  nlatpe=',nlatpe
!     WRITE(6,*) 'nlong=',nlong,'  mnlat=',mnlat
!     WRITE(6,*) 'nlnpe=',nlnpe,'  nlpe=',nlpe
!     WRITE(6,*) 'row_length=',row_length,'  rows=',rows
!     WRITE(6,*) 'lobound=',lobound,'  lnbound=',lnbound

! Run STOCHEM

! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER('STOCHEM',5)
! DEPENDS ON: stochem
      CALL STOCHEM(u,v,q,snowam,t0,bl,sice,orog,w,dc,p,lnp,o3um,        &
     &             clw,acr,t,cc,                                        &
     &             acp,adp,tropz,up,upd,                                &
     &             cc_base_level,cc_top_level,                          &
     &             land_fraction,                                       &
     &             totpar, dirpar, t0tile, gsf, lai_ft,                 &
     &             soilmc(:,:,1),                                       &
     &             canwc(:,:,1:npft), gc, ra, rq, so4_vd,               &
     &             um_year,um_month,um_day,                             &
     &             StochDumpRequest,                                    &
     &             a_step,a_fixhd(12),Secs_per_stepim(a_im),            &
     &             a_realhd(rh_z_top_theta),a_inthd(ih_1_c_rho_level),  &
     &             o3mmr,ch4mmr)
! DEPENDS ON: timer
      IF (LTIMER) CALL TIMER('STOCHEM',6)

! Write O3mmr and CH4mmr to user prognostic space
      IF (l_use_stochem_ch4 .OR. l_use_stochem_o3) THEN
        DO i=1,no_obj_d1(m)
          item = d1_addr(d1_item,i,m)
          IF (ANY(stochinprogcodes%icode == item) .AND.                 &
     &      d1_addr(d1_object_type,i,m) == prognostic) THEN
            levels = d1_addr(d1_no_levels,i,m)
            length = d1_addr(d1_length,i,m)
            address = d1_addr(d1_address,i,m)
            halo_type = d1_addr(d1_halo_type,i,m)
            grid_type = d1_addr(d1_grid_type,i,m)
! DEPENDS ON: get_fld_type
            ftype = get_fld_type(grid_type)
            WRITE(6,*) item,address,length,levels
            ALLOCATE(rfield(length))
            CALL STOCH_UNSHARE(rfield,length,levels,item,ftype,         &
     &        halo_type)
            d1(address:address+length-1) = rfield
            DEALLOCATE(rfield)
          END IF
        END DO
      END IF

 ! Archive and delete dumps as required.
      IF (mype == 0) CALL ARCHDUMP

      first = .false.
      StochDumpRequest = .false.
      ArchDumpRequest = .false.

      DEALLOCATE(cc_base_level)
      DEALLOCATE(cc_top_level)
      DEALLOCATE(land_fraction)
      DEALLOCATE(t)
      DEALLOCATE(q)
      DEALLOCATE(clw)
      DEALLOCATE(up)
      DEALLOCATE(upe)
      DEALLOCATE(upd)
      DEALLOCATE(cc)
      DEALLOCATE(dc)
      DEALLOCATE(p_th)
      DEALLOCATE(p)
      DEALLOCATE(u)
      DEALLOCATE(v)
      DEALLOCATE(w)
      DEALLOCATE(p0)
      DEALLOCATE(t0)
      DEALLOCATE(trop)
      DEALLOCATE(tropz)
      DEALLOCATE(bl)
      DEALLOCATE(sice)
      DEALLOCATE(snowam)
      DEALLOCATE(acp)
      DEALLOCATE(adp)
      DEALLOCATE(tauu)
      DEALLOCATE(tauv)
      DEALLOCATE(hf)
      DEALLOCATE(acr)
      DEALLOCATE(acs)
      DEALLOCATE(adr)
      DEALLOCATE(ads)
      DEALLOCATE(orog)
      DEALLOCATE(ch4mmr)
      DEALLOCATE(o3mmr)
      DEALLOCATE(totpar)
      DEALLOCATE(dirpar)
      DEALLOCATE(t0tile)
      DEALLOCATE(gsf)
      DEALLOCATE(lai_ft)
      DEALLOCATE(canht)
      DEALLOCATE(canwc)
      DEALLOCATE(z0)
      DEALLOCATE(sicetemp)
      DEALLOCATE(soilmc)
      DEALLOCATE(gc)
      DEALLOCATE(ra)
      DEALLOCATE(rq)
      DEALLOCATE(so4_vd)

      CONTAINS
! ######################################################################
      SUBROUTINE ARCHDUMP
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Archive and delete dump files from
!-                         the chemical model.
!-
!-   Inputs  : ArchDumpRequest
!-   Outputs :
!-   Controls:
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  5.0    25/06/99  Created.  C.E. Johnson
!  5.5    06/08/04  Now calls UM_FORT_FLUSH for compatability.
!                   C.E. Johnson
!  6.2    01/03/06  Corrected call to UM_FORT_FLUSH. M.G. Sanderson
!
!-
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------

      CHARACTER (80), SAVE :: pipename

      CALL GET_FILE(8,pipename,80,icode)        ! Get name of pipe
      OPEN(8,FILE=pipename)

! Archive current dumps when requested
      IF (archdumprequest .AND. stochdumprequest) THEN
! DEPENDS ON: findname
        CALL FINDNAME('d','a','c',0,0,filename)
        WRITE(6,*) 'Archiving: ',filename
        WRITE(8,620) filename
#if defined(T3E)
! DEPENDS ON: um_fort_flush
        CALL UM_FORT_FLUSH(8,icode)
#else
        CLOSE(8)
        OPEN(8,FILE=pipename)
#endif
      END IF

! Delete dumps 70 days previous to current
      IF (stochdumprequest .AND. a_step > 69*AINT(daysec/               &
     &  secs_per_stepim(atmos_im))) THEN
! DEPENDS ON: findname
        CALL FINDNAME('d','a','c',0,-70,filename)
        WRITE(6,*) 'Deleting: ',filename
        WRITE(8,610) filename
#if defined(T3E)
! DEPENDS ON: um_fort_flush
        CALL UM_FORT_FLUSH(8,icode)
#else
        CLOSE(8)
        OPEN(8,FILE=PIPENAME)
#endif
      END IF

  999 CONTINUE
  610 FORMAT('%%% ',A14,' DELETE')
  620 FORMAT('%%% ',A14,' ARCHIVE DUMP')

      END SUBROUTINE ARCHDUMP
! ######################################################################
      SUBROUTINE STOCH_SHARE(data1,length,levels,stashcode,ftype,       &
     &  halo_type)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Share out 1 block of met data to
!-                         the STOCHEM processor grid.
!-                         CONTAINED in Redist_Stochem.
!-
!-   Inputs  : StashCode,level=k,FirstShare(T/F)
!-             bigfield,smallfield,bounds
!-   Outputs : met data: U, V, ...
!-   Controls:
!-
!
! Current Owner of Code: C.E. Johnson
!
! History:
! Version   Date                    Comment
!  4.5    14/12/98  Created.  D.S. Stevenson
!  5.0    27/11/01  Removed CB, CT arrays. W.J. Collins.
!  6.2    20/04/05  SWAP_BOUNDS called for all halo types.
!                                          R. Johanni
!  6.2    01/03/06  Extra variables added to CASE statement.
!                                          M.G. Sanderson
!
!-
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
#include "c_mdi.h"
      INTEGER, INTENT(IN) :: stashcode     ! Stash code of field
      INTEGER, INTENT(IN) :: length        ! length of data1 array
      INTEGER, INTENT(IN) :: levels        ! no. of vertical levels
      INTEGER, INTENT(IN) :: ftype         ! field type
      INTEGER, INTENT(IN) :: halo_type     ! halo type

      REAL, DIMENSION(length), INTENT(IN) :: data1

      INTEGER :: ix0
      INTEGER :: ix1
      INTEGER :: iy0
      INTEGER :: iy1

! Note: smallfield is bigger than bigfield!
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev+1) :: smallfield
      REAL, DIMENSION(:,:,:), ALLOCATABLE      :: bigfield

      ALLOCATE(bigfield(lasize(1,ftype,halo_type),                      &
     &  lasize(2,ftype,halo_type),levels))
      bigfield = RESHAPE(data1,(/lasize(1,ftype,halo_type),             &
     &  lasize(2,ftype,halo_type),levels/))
      ix0 = 1 + halo_i - halosize(1,halo_type)
      ix1 = ix0 + lasize(1,ftype,halo_type) - 1
      iy0 = 1 + halo_j - halosize(2,halo_type)
      iy1 = iy0 + lasize(2,ftype,halo_type) - 1
      smallfield = rmdi
      smallfield(ix0:ix1,iy0:iy1,1:levels) = bigfield

! Need to set up haloes. Need to call SWAP_BOUNDS even if
! halo type is extended.
! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(                                                 &
     &  smallfield(1:2*halo_i+row_length,1:2*halo_j+rows,1:levels),     &
     &  row_length,rows,levels,halo_i,halo_j,ftype,.FALSE.)

      SELECT CASE(stashcode)
      CASE(256)   ! U
        u(:,:,1:nmetlev)=smallfield(:,:,1:levels)
      CASE(257)   ! V
        v(:,:,1:nmetlev)=smallfield(:,:,1:levels)
      CASE(10)    ! Q
        q               =smallfield(:,:,1:levels)
      CASE(14)    ! CC_base_Level
        cc_base_level   =nint(smallfield(:,:,1))
      CASE(15)    ! CC_Top_Level
        cc_top_level    =nint(smallfield(:,:,1))
      CASE(23)    ! SNOWAM
        snowam          =smallfield(:,:,1)
      CASE(25)    ! BL
        bl              =smallfield(:,:,1)
      CASE(24)    ! T0
        t0              =smallfield(:,:,1)
      CASE(31)    ! SICE
        sice            =smallfield(:,:,1)
      CASE(33)    ! OROG
        orog            =smallfield(:,:,1)
      CASE(49)    ! SICETEMP
        sicetemp        =smallfield(:,:,1)
      CASE(60)    ! O3UM
        o3um(:,:,1:nmetlev)=smallfield(:,:,1:levels)
      CASE(258)   ! w
        w               =smallfield(:,:,1:levels)
      CASE(265)   ! DC
        dc              =smallfield(:,:,1:levels)
      CASE(408)   ! P_TH
        p_th            =smallfield(:,:,1:levels)
      CASE(409)   ! P0
        p0              =smallfield(:,:,1)
      CASE(1290)  ! TOTPAR
        totpar          =smallfield(:,:,1)
      CASE(1291)  ! DIRPAR
        dirpar          =smallfield(:,:,1)
      CASE(3217)  ! HF
        hf              =smallfield(:,:,1)
      CASE(3219)  ! Tau_U
        tauu            =smallfield(:,:,1)
      CASE(3220)  ! Tau_V
        tauv            =smallfield(:,:,1)
      CASE(3316)  ! T0TILE
        t0tile          =smallfield(:,:,1:ntype)
      CASE(3317)  ! GSF
        gsf             =smallfield(:,:,1:ntype)
      CASE(3318)  ! LAI_FT
        lai_ft          =smallfield(:,:,1:npft)
      CASE(3319)  ! CANHT
        canht           =smallfield(:,:,1:npft)
      CASE(3321)  ! CANWC
        canwc           =smallfield(:,:,1:ntype)
      CASE(3324)  ! Roughness length on tiles
        z0              =smallfield(:,:,1:ntype)
      CASE(3395)  ! Land fraction
        land_fraction   =smallfield(:,:,1)
      CASE(3462)  ! Stomatal conductance on tiles
        gc              =smallfield(:,:,1:npft)
      CASE(4203)  ! ADR
        adr             =smallfield(:,:,1)
      CASE(4204)  ! ADS
        ads             =smallfield(:,:,1)
      CASE(4205)  ! CLW
        clw             =smallfield(:,:,1:levels)
      CASE(5205)  ! ACR
        acr             =smallfield(:,:,1)
      CASE(5206)  ! ACS
        acs             =smallfield(:,:,1)
      CASE(5209)  ! T
        t               =smallfield(:,:,1:levels)
      CASE(5212)  ! CC
        cc              =smallfield(:,:,1:levels)
      CASE(5250)  ! UP
        up(:,:,1:levels)=smallfield(:,:,1:levels)
      CASE(5253)  ! UPD
        upd(:,:,1:levels)=smallfield(:,:,1:levels)
      CASE(8223)  ! SOILMC
        soilmc(:,:,1:sm_levels)=smallfield(:,:,1:sm_levels)
      CASE(30451) ! TROP (in Pa)
        trop            =smallfield(:,:,1)
      CASE(30453) ! TROP (in m)
        tropz           =smallfield(:,:,1)
      CASE DEFAULT
        cmessage = 'Failure in case select of stashcode'
        WRITE(6,*) cmessage,'STASH: ',Stashcode
! DEPENDS ON: ereport
        CALL EREPORT('STOCH_SHARE',1,cmessage)
      END SELECT

      END SUBROUTINE STOCH_SHARE
! ######################################################################
      SUBROUTINE STOCH_UNSHARE(data1,length,levels,stashcode,ftype,     &
     &  halo_type)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Share out STOCHEM data to
!-                         the met grid.
!-                         Hard coded for zonal means
!-                         CONTAINED in Redist_Stochem.
!-
!-   Inputs  : length,levels,StashCode,ftype,halo_type
!-
!-   Outputs : met data: data1
!-   Controls:
!-
! Current Owner of Code: W.J. Collins
!
! History:
! Version   Date                    Comment
!  5.5    14/04/03  Created. W.J. Collins
!  5.5    29/10/03  Created. W.J. Collins
!  6.1    20/10/04  No change.
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
#include "c_mdi.h"

      INTEGER, INTENT(IN) :: stashcode     ! Stash code of field
      INTEGER, INTENT(IN) :: length        ! length of data1 array
      INTEGER, INTENT(IN) :: levels        ! no. of vertical levels
      INTEGER, INTENT(IN) :: ftype         ! field type
      INTEGER, INTENT(IN) :: halo_type     ! halo type

      REAL, DIMENSION(length), INTENT(OUT) :: data1

      INTEGER :: ix0
      INTEGER :: ix1
      INTEGER :: iy0
      INTEGER :: iy1

! Note: smallfield is bigger than bigfield
      REAL, DIMENSION(nlonpe,nlatpe,nmetlev) :: smallfield
      REAL, DIMENSION(lasize(1,ftype,halo_type),                        &
     &                lasize(2,ftype,halo_type),levels) :: bigfield

      SELECT CASE(StashCode)
      CASE(99)       ! ch4mmr
        smallfield = ch4mmr
      CASE(100)      ! o3mmr
        smallfield = o3mmr
      CASE DEFAULT
        cmessage = 'Failure in case select of stashcode'
        WRITE(6,*) cmessage,'STASH: ',stashcode
! DEPENDS ON: ereport
        CALL EREPORT('STOCH_UNSHARE',1,cmessage)
      END SELECT
      iy0 = 1 + halo_j - halosize(2,halo_type)
      iy1 = iy0 + lasize(2,ftype,halo_type) - 1
      ix0 = 1 + halo_i - halosize(1,halo_type)
      ix1 = ix0 + lasize(1,ftype,halo_type) - 1
      bigfield = smallfield(ix0:ix1,iy0:iy1,1:levels)
      data1 = RESHAPE(bigfield,(/lasize(1,ftype,halo_type)*             &
     &                lasize(2,ftype,halo_type)*levels/))

      END SUBROUTINE STOCH_UNSHARE

      END SUBROUTINE REDIST_STOCHEM
#endif
