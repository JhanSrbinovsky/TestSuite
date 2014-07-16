#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SET_ATM_POINTERS ---------------------------------------
!LL
!LL  Set pointers for primary atmosphere fields
!LL Initialisation routine for CRAY YMP
!LL
!LL MC, CW      <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL  3.1    9/02/93  : added comdeck CHSUNITS to define NUNITS for
!LL                    comdeck CCONTROL.
!LL 3.1     11/03/93  Set JTRACER(1,1) to have a sensible address
!L                    even if there are no tracers to remove bounds
!L                    checking problems in later routines. R. Rawlins
!LL 3.2    27/03/93 Dynamic allocation of main data arrays. R. Rawlins
!LL 3.4    20/01/94 Changes to allow for non-consecutive tracers.
!LL                 M.Carter
!LL  3.3  22/11/93  Add aerosol ancillary fields.  R T H Barnes.
!LL  3.3   26/10/93  M. Carter. Part of an extensive mod that:
!LL                  1.Removes the limit on primary STASH item numbers.
!LL                  2.Removes the assumption that (section,item)
!LL                    defines the sub-model.
!LL                  3.Thus allows for user-prognostics.
!LL                   Remove a_max_variables.
!LL                   PPINDEX now read from UI in INITCTL.
!LL  3.4  05/09/94  Add murk & user ancillary fields.  RTHBarnes.
!LL  3.4   18/05/94  J.Thomson add pointers for slab temperature
!LL                            and u and v ice velocities.
!LL  3.5   19/05/95  Some pointers for level dependent constants
!LL                  array removed. Sub_model change. D. Robinson
!LL  4.0   06/09/95  Set up pointers correctly for SLAB fields.
!LL                  D. Robinson
!LL
!LL  4.0 26/7/94 R.E. Essery Extra prognostic sea-ice temperature.
!LL  4.1 10/1/96      Extra prognostics frozen and unfrozen soil
!LL                   moisture fractions, and canopy conductance
!LL                   plus 2 extra vegetation fields J.Smith
!LL  4.1 30/04/96     Add pointers for 6 new variables and 6 new
!LL                   ancillary fields for Sulphur Cycle   M Woodage
!LL  4.3  18/3/97     Add pointers for HadCM2 sulphate loading patterns
!LL                                                   William Ingram
!LL  4.2  16/08/96    Added MPP PARVARS comdeck, and defined rowdepc
!LL                   (filterable wave numbers array) to be globally
!LL                   sized.                                P.Burton
!LL  4.3  06/03/97    Dimension multi_level land fields by LAND_FIELD
!LL                   for MPP jobs. D. Robinson.
!LL  4.4  05/09/97    Add pointer for net energy flux prognostic
!LL                   S.D.Mullerworth
!LL  4.4  04/08/97  Generalise JQCF pointer for mixed phase
!LL                 precipitation scheme.  RTHBarnes.
!LL  4.4  05/08/97    Add pointer for convective cloud amount on
!LL                   model levels (3D CCA) if L_3D_CCA. J.M.Gregory
!LL  4.4  10/09/97    Added pointers for snow grain size and snow soot
!LL                   content used in prognostic snow albedo scheme
!LL                   R. Essery
!LL  4.4  16/09/97    Add call to NSTYPES and pointers for new
!LL                   vegetation and land surface prognostics. R.A.Betts
!LL  4.5   1/07/98    Add pointers for ocean CO2 flux and surface
!LL                     CO2 emissions. C.D.Jones
!LL  4.5  04/03/98   Add pointers for NH3 prognostic and NH3 surface
!LL                    emiss for S Cycle                 M Woodage
!LL                  Also add pointers for 3 soot prognostic variables
!LL                    and 2 soot emiss                    M Woodage
!LL  4.5  15/07/98  Add pointers for new 3D CO2 array. C.D.Jones
!    4.5  22/10/98    Add pointers for extra multi-layer user
!                     ancillary fields
!                     Author D.M. Goddard
!    4.5  29/04/98    Pointer to total soil moisture content to point
!                     to non prognostic space in MOSES.
!                     Author D.M. Goddard
!LL  4.5  19/01/98    Replace JVEG_FLDS and JSOIL_FLDS with
!LL                   individual pointers. D. Robinson
!LL  5.0  13/07/99    Extensive changes for C-P C grid upgrade
!LL                   D Goddard
!LL  5.1  22/02/00  Add PARVARS for TYPSIZE                 P.Burton
!LL  5.1  28/02/00    Add pointers for ND LBCs.
!LL                   Old JRIM pointers are kept in until all
!LL                   references in code are deleted at 5.2
!LL                                                         P.Burton
!  5.1  28/04/00   Extra pressure level above top theta level
!                  p_top_theta_level removed                A.Malcolm
!LL  5.2  13/09/00    Remove JSOOT, JSO4, JH2SO4 pointers and
!LL                   change sizes for tracers. P.Selwood
!LL
!LL  5.1  29/02/00    Add JNET_MFLUX R A Stratton.
!LL  5.2  18/10/00    Removed JUICE and JVICE.            K.D.Williams
!    5.2  31/08/00    Add JOROG_LBC pointer
!                     Remove levels from LBC pointers
!                     Change section/item for LBC pointers
!                     Remove JRIM and JRIM_TENDENCY
!                                                   P.Burton
!    5.2  25/09/00    Clear out dead legacy code (3D RHcrit not needed
!                     in dump).                           A.C.Bushell
!
!LL  5.2  15/11/00  Additional initialisation of pointers for
!LL                 MOSES 2 and dynamic vegetation         M. Best
!LL  5.2  09/03/01  Introduce extra pointers in level dependent
!                   constants for height definitions. R Rawlins
!LL  5.3  04/10/01  Added new Slab model prognostics and re-order above
!LL                 where Sect_No is set to 31.            K.Williams
!    5.3  20/07/01  Replace hardwires with access of SI index when
!                   calculating pointers for items 406 to 412
!                   S.D.Mullerworth
!LL  5.3  19/06/01   Add JTPPSOZONE Dave Tan
!    5.4  02/09/02  New pointer JSNODEP_SEA added        K.Williams
!    5.4  28/08/02  Set up pointers for canopy snow and snow
!                   beneath canopy.  R. Essery
!    5.4  27/07/02  Increase cloud fraction halo sizes for PC2
!                                                          D. Wilson
!    5.5  05/11/02  Large-scale hydrology scheme: Add pointers for
!                   gamma function, water table depth,
!                   surface saturated fraction and wetland fraction.
!                                                        N. Gedney
!    5.5  05/02/03  Add pointers for biomass aerosol scheme   P Davison
!    5.5  13/02/03  New pointers for multicatagory sea ice. J.Ridley
!
!    5.5  03/02/03  Add pointers for 3 moisture arrays
!                   (qcf2,qrain,qgraup) and lbc,lbc_tend     R.M.Forbes
!    5.5  12/02/03  Include code for mineral dust scheme. S Woodward
!
!    5.5  26/02/03  Add pointers for river routing. P.Selwood.
!  6.0  14/08/03  NEC SX-6 Fix - avoid undefined value?
!                 R Barnes & J-C Rioual.
!    6.0  30/07/03  Add pointers for cloud fraction lbcs. Damian Wilson
!  6.0    12/08/03 removes JRIV_GRIDBOX. C.Bunton.
!    6.2  11/11/05  Set pointers for section 34 UKCA tracers. R Barnes.
!    6.1  28/06/04  Set pointers for section 33 free tracers. R Barnes.
!  6.1    20/08/03  Add pointers for STOCHEM fields.  C. Johnson
!    6.1  07/04/04  Add pointers for DMS concentration.      A. Jones
!    6.2  13/07/05  Correct free tracer LBC pointers. R Barnes.
!    6.2  19/08/05  Fix out of bounds reference. P.Selwood.
!    6.2  15/08/05  Free format fixes. P.Selwood
! 6.2  02/03/06 Added pointer direct PAR flux. M.G. Sanderson
!    6.2  01/10/04  Add pointers for murk aerosol lbc.    R.M.Forbes
!    6.2  01/03/06  Add pointers for RothC pools and fluxes. C.D. Jones
!   6.2   21/2/06  Re-route outflow from inland basins to soil moisture
!                  P. Falloon

! 6.2  25/11/05 Functionality for improved time stepping, radiative
!               forcing and radiance code added for versions 3C
!               and 3Z of radiation code             (J.-C. Thelen)
!    6.2  24/02/06  Add pointers for ocean DMS flux             J.Gunson
!LL Programming Standard: Unified Model DP NO. 3, Version 3
!LL
!LL  Logical task: P0
!LL
!LL  System Components: C21 (Atmosphere part)
!LL
!LL  Purpose:   Sets integer pointers to atmospheric
!LL             variables from STASHIN addresses.
!LL
!LL  External documentation: UMDP NO. C4 Version NO. 4
!LL
!LLEND-------------------------------------------------------------

      SUBROUTINE SET_ATM_POINTERS(                                      &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
     &                  ICODE,CMESSAGE)

      USE CSENARIO_MOD
      IMPLICIT NONE

!L
!*L Arguments
!L
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "nstypes.h"
#include "typduma.h"
#include "typsts.h"
#include "typptra.h"
      INTEGER                                                           &
     &    ICODE                  ! OUT: Error return code
!
      CHARACTER*80                                                      &
     &    CMESSAGE               ! OUT: Error return message
#include "chsunits.h"
#include "ccontrol.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "ctracera.h"
#include "cmaxsize.h"
#include "cruntimc.h"

! local variables

      INTEGER                                                           &
     &        IVAR,                                                     &
                                 ! Loop counts
     &        JVAR,                                                     &
                                 ! Loop counts
     &        IFLD,                                                     &
     &        LEV                                                       &
     &       ,im_ident                                                  &
                            !  Internal Model Identifier
     &       ,im_index                                                  &
                            !  Internal Model Index in Stash arrays
     &       ,Sect_No       !  Stash section number

! Tracer fluxes stash section number - kdcorbin, 06/10
     INTEGER :: tflux_sect

!     Set to atmosphere internal model
      im_ident  = atmos_im
      im_index  = internal_model_index(im_ident)
      Sect_No   = 0


! Set pointers for atmospheric primary variables from STASH :-
      JU(1)                = SI(  2,Sect_No,im_index)
      JV(1)                = SI(  3,Sect_No,im_index)
      JTHETA(1)            = SI(  4,Sect_No,im_index)
      JQ(1)                = SI( 10,Sect_No,im_index)
      JQCF(1)              = SI( 12,Sect_No,im_index)
      JTSTAR               = SI( 24,Sect_No,im_index)
      JLAND                = SI( 30,Sect_No,im_index)
      JOROG                = SI( 33,Sect_No,im_index)
      JW(0)                = SI(150,Sect_No,im_index)
      JRHO(1)              = SI(253,Sect_No,im_index)
      JQCL(1)              = SI(254,Sect_No,im_index)
      JEXNER_RHO_LEVELS(1) = SI(255,Sect_No,im_index)
      JQCF2(1)             = SI(271,Sect_No,im_index)
      JQRAIN(1)            = SI(272,Sect_No,im_index)
      JQGRAUP(1)           = SI(273,Sect_No,im_index)

!LL  5.3     06/01  Introduce extra pointers for coastal tiling.
!                                                     Nic Gedney.
      JFRAC_LAND     = SI(505,Sect_No,im_index)
      JTSTAR_LAND    = SI(506,Sect_No,im_index)
      JTSTAR_SEA     = SI(507,Sect_No,im_index)
      JTSTAR_SICE    = SI(508,Sect_No,im_index)
! Set pointers for seaice and land albedos
      JSICE_ALB      = SI(509,Sect_No,im_index)
      JLAND_ALB      = SI(510,Sect_No,im_index)
!
!LL  5.3     06/01  Introduce extra pointers for large-scale hydrology.
!                                                     Nic Gedney.
      JTI_MEAN       = SI(274,Sect_No,im_index)
      JTI_SIG        = SI(275,Sect_No,im_index)
      JFEXP          = SI(276,Sect_No,im_index)
      JGAMMA_INT     = SI(277,Sect_No,im_index)
      JWATER_TABLE   = SI(278,Sect_No,im_index)
      JFSFC_SAT      = SI(279,Sect_No,im_index)
      JF_WETLAND     = SI(280,Sect_No,im_index)
      JSTHZW         = SI(281,Sect_No,im_index)
      JA_FSAT        = SI(282,Sect_No,im_index)
      JC_FSAT        = SI(283,Sect_No,im_index)
      JA_FWET        = SI(284,Sect_No,im_index)
      JC_FWET        = SI(285,Sect_No,im_index)
      DO LEV=1,MODEL_LEVELS
        JW(LEV)=JW(LEV-1)+theta_off_size
      END DO
      DO LEV=2,MODEL_LEVELS
        JU(LEV)    =JU(LEV-1)     + u_off_size
        JV(LEV)    =JV(LEV-1)     + v_off_size
        JRHO(LEV)  =JRHO(LEV-1)   + theta_off_size
        JTHETA(LEV)=JTHETA(LEV-1) + theta_off_size
      END DO
      DO LEV=2,MODEL_LEVELS+1
        JEXNER_RHO_LEVELS(LEV)=JEXNER_RHO_LEVELS(LEV-1)+theta_off_size
      END DO
      DO LEV=2,WET_LEVELS
        JQ(LEV)    =JQ(LEV-1)     + theta_halo_size
        JQCL(LEV)  =JQCL(LEV-1)   + theta_halo_size
        JQCF(LEV)  =JQCF(LEV-1)   + theta_halo_size
        JQCF2(LEV)   =JQCF2(LEV-1)   + theta_halo_size
        JQRAIN(LEV)  =JQRAIN(LEV-1)  + theta_halo_size
        JQGRAUP(LEV) =JQGRAUP(LEV-1) + theta_halo_size
      END DO

      ! Pointers required to save coupling fields as
      ! prognostics when employing OASIS as the coupler.
      ! In non-coupled models these pointers will
      ! simply end up with a value of 1.
      JC_SOLAR = SI(171,Sect_No,im_index)
      JC_BLUE =  SI(172,Sect_No,im_index)
      JC_DOWN =  SI(173,Sect_No,im_index)
      JC_LONGWAVE =  SI(174,Sect_No,im_index)
      JC_TAUX =  SI(176,Sect_No,im_index)
      JC_TAUY =  SI(177,Sect_No,im_index)
      JC_WINDMIX =  SI(178,Sect_No,im_index)
      JC_SENSIBLE =  SI(179,Sect_No,im_index)
      JC_SUBLIM =  SI(180,Sect_No,im_index)
      JC_EVAP =  SI(181,Sect_No,im_index)
      JC_BOTMELTN =  SI(184,Sect_No,im_index)
      JC_TOPMELTN =  SI(185,Sect_No,im_index)
      JC_LSRAIN =  SI(186,Sect_No,im_index)
      JC_LSSNOW =  SI(187,Sect_No,im_index)
      JC_CVRAIN =  SI(188,Sect_No,im_index)
      JC_CVSNOW =  SI(189,Sect_No,im_index)
      JC_RIVEROUT =  SI(192,Sect_No,im_index)

#if defined(ACCESS)    
! gol124: auscom coupling
! (the same as JPSTAR)
      JC_PRESS = SI(409,Sect_No,im_index)
#endif

! Set pointers for optional atmospheric primary variables.
!  To be enabled at a later version than 5.0
!     IF(bit compatability required)
        JZH           = SI( 25,Sect_No,im_index)
        JU_ADV(1)     = SI(256,Sect_No,im_index)
        JV_ADV(1)     = SI(257,Sect_No,im_index)
        JW_ADV(0)     = SI(258,Sect_No,im_index)
        JNTML         = SI(259,Sect_No,im_index)
        JNBDSC        = SI(260,Sect_No,im_index)
        JNTDSC        = SI(261,Sect_No,im_index)
        JCUMULUS      = SI(262,Sect_No,im_index)
        JT1_SD        = SI(263,Sect_No,im_index)
        JQ1_SD        = SI(264,Sect_No,im_index)
        JCF_AREA(1)   = SI(265,Sect_No,im_index)
        JCF_BULK(1)   = SI(266,Sect_No,im_index)
        JCF_LIQUID(1) = SI(267,Sect_No,im_index)
        JCF_FROZEN(1) = SI(268,Sect_No,im_index)

        IF (L_3D_CCA) THEN
          JCCA(1) = SI(211,Sect_No,im_index)
          DO LEV=2,N_CCA_LEV
            JCCA(LEV)=JCCA(LEV-1)+THETA_FIELD_SIZE
          ENDDO
        ELSE
          JCCA(1)      = SI( 13,Sect_No,im_index)
        ENDIF

        JCCB           = SI( 14,Sect_No,im_index)
        JCCT           = SI( 15,Sect_No,im_index)
        JCCLWP         = SI( 16,Sect_No,im_index)
        JLCBASE        = SI( 21,Sect_No,im_index)
        JCANOPY_WATER  = SI( 22,Sect_No,im_index)
        JCCW_RAD(1)    = SI(212,Sect_No,im_index)

        DO LEV=2, WET_LEVELS
            JCCW_RAD(LEV)=JCCW_RAD(LEV-1)+THETA_FIELD_SIZE
        ENDDO
!     Endif !IF(bit compatability required)

      DO LEV=1,MODEL_LEVELS
        JW_ADV(LEV)=JW_ADV(LEV-1)+theta_halo_size
      END DO
      DO LEV=2,MODEL_LEVELS
        JU_ADV(LEV)=JU_ADV(LEV-1)+u_halo_size
        JV_ADV(LEV)=JV_ADV(LEV-1)+v_halo_size
      END DO
      DO LEV=2,WET_LEVELS
        JCF_AREA(LEV)  =JCF_AREA(LEV-1)+THETA_FIELD_SIZE
        JCF_BULK(LEV)  =JCF_BULK(LEV-1)+THETA_halo_SIZE
        JCF_LIQUID(LEV)=JCF_LIQUID(LEV-1)+THETA_halo_SIZE
        JCF_FROZEN(LEV)=JCF_FROZEN(LEV-1)+THETA_halo_SIZE
      END DO

!  Set pointers for secondary fields in D1
      JEXNER_THETA_LEVELS(1)= SI(406,Sect_No,im_index)
      JP(1)                 = SI(407,Sect_No,im_index)
      JP_THETA_LEVELS(1)    = SI(408,Sect_No,im_index)
      JPSTAR                = SI(409,Sect_No,im_index)
      JSW_INCS(0)           = SI(410,Sect_No,im_index)
      JLW_INCS(0)           = SI(411,Sect_No,im_index)

      DO LEV=1,MODEL_LEVELS+1
        JSW_INCS(LEV)=JSW_INCS(LEV-1)+THETA_FIELD_SIZE
      END DO
      DO LEV=1,MODEL_LEVELS
        JLW_INCS(LEV)=JLW_INCS(LEV-1)+THETA_FIELD_SIZE
      END DO

! Direct PAR flux for STOCHEM
      JDIRPAR               = SI(460,Sect_no,im_index)

      DO LEV=2,MODEL_LEVELS
        JEXNER_THETA_LEVELS(LEV)=JEXNER_THETA_LEVELS(LEV-1)+            &
     &                             theta_off_size
        JP_THETA_LEVELS(LEV)    =JP_THETA_LEVELS(LEV-1)+theta_off_size
      END DO
      DO LEV=2,MODEL_LEVELS+1
        JP(LEV)                 =JP(LEV-1)+theta_off_size
      END DO

!  Set pointers for ancillary fields in D1 from STASH
!     Soil fields
      JSMCL(1)            = SI(  9,Sect_No,im_index)
      J_DEEP_SOIL_TEMP(1) = SI( 20,Sect_No,im_index)
      JVOL_SMC_WILT       = SI( 40,Sect_No,im_index)
      JVOL_SMC_CRIT       = SI( 41,Sect_No,im_index)
      JVOL_SMC_SAT        = SI( 43,Sect_No,im_index)
      JSAT_SOIL_COND      = SI( 44,Sect_No,im_index)
      JTHERM_CAP          = SI( 46,Sect_No,im_index)
      JTHERM_COND         = SI( 47,Sect_No,im_index)
      JSAT_SOILW_SUCTION  = SI( 48,Sect_No,im_index)
      JCLAPP_HORN         = SI(207,Sect_No,im_index)
      JSTHU(1)            = SI(214,Sect_No,im_index)
      JSTHF(1)            = SI(215,Sect_No,im_index)

      DO LEV=2,ST_LEVELS
        J_DEEP_SOIL_TEMP(LEV)=J_DEEP_SOIL_TEMP(LEV-1)+LAND_FIELD
      ENDDO

      DO LEV=2,SM_LEVELS
        JSMCL(LEV)=JSMCL(LEV-1)+LAND_FIELD
        JSTHU(LEV)=JSTHU(LEV-1)+LAND_FIELD
        JSTHF(LEV)=JSTHF(LEV-1)+LAND_FIELD
      END DO

!     Vegetation Fields
      JZ0          = SI( 26,Sect_No,im_index) ! roughness length
      JVEG_FRAC    = SI( 50,Sect_No,im_index) ! vegetation fraction
      JROOT_DEPTH  = SI( 51,Sect_No,im_index) ! root depth
      JSFA         = SI( 52,Sect_No,im_index) ! snow free albedo
      JMDSA        = SI( 53,Sect_No,im_index) ! deep snow albedo
      JSURF_RESIST = SI( 54,Sect_No,im_index) ! surface resistance
      JSURF_CAP    = SI( 55,Sect_No,im_index) ! surface capacity
      JINFILT      = SI( 56,Sect_No,im_index) ! infiltration factor
      JLAI         = SI(208,Sect_No,im_index) ! leaf area index
      JCANHT       = SI(209,Sect_No,im_index) ! canopy height
      JGS          = SI(213,Sect_No,im_index) ! stomatal conductance

! Orography fields
      JOROG_SIL      = SI(17,Sect_No,im_index)   ! Silhouette area
      JOROG_HO2      = SI(18,Sect_No,im_index)   ! Peak to trough ht.
      JOROG_SD       = SI(34,Sect_No,im_index)
      JOROG_GRAD_X   = SI( 5,Sect_No,im_index)
      JOROG_GRAD_Y   = SI( 6,Sect_No,im_index)
      JOROG_GRAD_XX  = SI(35,Sect_No,im_index)
      JOROG_GRAD_XY  = SI(36,Sect_No,im_index)
      JOROG_GRAD_YY  = SI(37,Sect_No,im_index)

! Sea/Sea Ice fields
      JU_SEA         = SI( 28,Sect_No,im_index)
      JV_SEA         = SI( 29,Sect_No,im_index)
      JICE_FRACTION  = SI( 31,Sect_No,im_index)
      JICE_THICKNESS = SI( 32,Sect_No,im_index)
      JTI            = SI( 49,Sect_No,im_index)
      JICE_FRACT_CAT = SI(413,Sect_No,im_index)
      JICE_THICK_CAT = SI(414,Sect_No,im_index)
      JTI_CAT        = SI(415,Sect_No,im_index)
      JU_0_P         = SI(269,Sect_No,im_index)
      JV_0_P         = SI(270,Sect_No,im_index)

! Snow fields
      JSNODEP        = SI( 23,Sect_No,im_index) ! Snow depth over land
      JSNODEP_SEA    = SI( 95,Sect_No,im_index) ! Snow depth on sea ice
      JSNODEP_SEA_CAT= SI(416,Sect_No,im_index) ! Snow depth on ice cats
      JSNSOOT        = SI(221,Sect_No,im_index) ! Snow soot content
      JCATCH_SNOW    = SI(241,Sect_No,im_index)
      JSNOW_GRND     = SI(242,Sect_No,im_index)

! Ozone
      JOZONE(1)     = SI(60,Sect_No,im_index)
! Check for zonal ozone and calculate pointers accordingly
      LEXPAND_OZONE=.FALSE.
      IF (A_LOOKUP(LBNPT,PPINDEX(60,im_index)) == 1) THEN
        LEXPAND_OZONE = .TRUE.
      ENDIF
      DO LEV=2,OZONE_LEVELS
        IF(LEXPAND_OZONE) THEN
          JOZONE(LEV)=JOZONE(LEV-1)+ROWS
        ELSE
          JOZONE(LEV)=JOZONE(LEV-1)+THETA_FIELD_SIZE
        END IF
      END DO

!! Tropopause-based Ozone
      IF (tpps_ozone_levels >  0) THEN
        JTPPSOZONE(1)     = SI(341,Sect_No,im_index)

        !Check for zonal tpps_ozone and calculate pointers accordingly
        LEXPAND_TPPS_OZONE=.FALSE.
        IF (A_LOOKUP(LBNPT,PPINDEX(341,im_index)) == 1) THEN
          LEXPAND_TPPS_OZONE = .TRUE.
        ENDIF
      END IF

      DO LEV=2,TPPS_OZONE_LEVELS
        IF(LEXPAND_TPPS_OZONE) THEN
          JTPPSOZONE(LEV)=JTPPSOZONE(LEV-1)+ROWS
        ELSE
          JTPPSOZONE(LEV)=JTPPSOZONE(LEV-1)+THETA_FIELD_SIZE
        END IF
      END DO

! Add prognostic ozone tracer and cariolle parameters to section 0
! 
        JOZONE_TRACER = SI(480,sect_no,im_index)
        JO3_PROD_LOSS = SI(481,sect_no,im_index)
        JO3_P_L_VMR   = SI(482,sect_no,im_index)
        JO3_VMR       = SI(483,sect_no,im_index)
        JO3_P_L_TEMP  = SI(484,sect_no,im_index)
        JO3_TEMP      = SI(485,sect_no,im_index)
        JO3_P_L_COLO3 = SI(486,sect_no,im_index)
        JO3_COLO3     = SI(487,sect_no,im_index)

        DO  LEV = 2, model_levels
          JOZONE_TRACER(LEV) = JOZONE_TRACER(LEV-1) + THETA_OFF_SIZE
          JO3_PROD_LOSS(LEV) = JO3_PROD_LOSS(LEV-1) + ROWS
          JO3_P_L_VMR(LEV)   = JO3_P_L_VMR(LEV-1) + ROWS
          JO3_VMR(LEV)       = JO3_VMR(LEV-1) + ROWS
          JO3_P_L_TEMP(LEV)  = JO3_P_L_TEMP(LEV-1) + ROWS 
          JO3_TEMP(LEV)      = JO3_TEMP(LEV-1) + ROWS
          JO3_P_L_COLO3(LEV) = JO3_P_L_COLO3(LEV-1) + ROWS
          JO3_COLO3(LEV)     = JO3_COLO3(LEV-1) + ROWS

        END DO


! STOCHEM fields
! Add stochem fields
      JCH4_STOCH(1)=SI(99,Sect_No,im_index)
      JO3_STOCH(1)=SI(100,Sect_No,im_index)

      DO LEV=2,MODEL_LEVELS
        JCH4_STOCH(LEV)=JCH4_STOCH(LEV-1)+THETA_FIELD_SIZE
        JO3_STOCH(LEV)=JO3_STOCH(LEV-1)+THETA_FIELD_SIZE
      ENDDO

! Add sources and aerosol ancillaries

      JMURK_SOURCE(1) = SI(57,Sect_No,im_index) ! Murk source
      JSO2_EM   = SI(58,Sect_No,im_index)   ! Sulphur dioxide emiss.
      JDMS_EM   = SI(59,Sect_No,im_index)   ! Dimethyl sulphide emiss.
      JMURK(1)  = SI(90,Sect_No,im_index)   ! Murk concentration

! Add for Sulphur Cycle
      JSO2(1)       =SI(101,Sect_No,im_index) !Sulphur dioxide gas
      JDMS(1)       =SI(102,Sect_No,im_index) !Dimethyl sulphide gas
      JSO4_AITKEN(1)=SI(103,Sect_No,im_index) !Aitken mode SO4 aerosol
      JSO4_ACCU(1)  =SI(104,Sect_No,im_index) !Accumulation mode SO4 aer
      JSO4_DISS(1)  =SI(105,Sect_No,im_index) !Dissolved SO4 aerosol
      JH2O2(1)      =SI(106,Sect_No,im_index) !Hydrogen peroxide mmr
      JNH3(1)       =SI(107,Sect_No,im_index)  !Ammonia gas
      JSOOT_NEW(1)  =SI(108,Sect_No,im_index)  !Fresh soot
      JSOOT_AGD(1)  =SI(109,Sect_No,im_index)  !Aged soot
      JSOOT_CLD(1)  =SI(110,Sect_No,im_index)  !Soot in cloud
      JBMASS_NEW(1) =SI(111,Sect_No,im_index)  !Fresh biomass smoke
      JBMASS_AGD(1) =SI(112,Sect_No,im_index)  !Aged biomass smoke
      JBMASS_CLD(1) =SI(113,Sect_No,im_index)  !Cloud biomass smoke
      JOCFF_NEW(1)  =SI(114,Sect_No,im_index)  !Fresh ocff
      JOCFF_AGD(1)  =SI(115,Sect_No,im_index)  !Aged socff
      JOCFF_CLD(1)  =SI(116,Sect_No,im_index)  !Ocff in cloud
      JSO2_NATEM(1) =SI(121,Sect_No,im_index)  !Natural SO2 emissions
      JOH(1)        =SI(122,Sect_No,im_index)  !OH 3_D ancillary
      JHO2(1)       =SI(123,Sect_No,im_index)  !HO2 3_D ancillary
      JH2O2_LIMIT(1)=SI(124,Sect_No,im_index)  !H2O2 LIMIT 3_D ancillary
      JO3_CHEM(1)   =SI(125,Sect_No,im_index)  !O3 for chemistry 3_D anc
      JSO2_HILEM    =SI(126,Sect_No,im_index)  !High level SO2 emissions
      JNH3_EM       =SI(127,Sect_No,im_index)  !Ammonia surface emiss
      JSOOT_EM      =SI(128,Sect_No,im_index)  !Fresh soot surf emiss
      JSOOT_HILEM   =SI(129,Sect_No,im_index)  !Fresh soot high emiss
      JBMASS_EM     =SI(130,Sect_No,im_index)  !Fresh bmass surf emiss
      JBMASS_HILEM  =SI(131,Sect_No,im_index)  !Elevated bmass emiss
      JDMS_CONC     =SI(132,Sect_No,im_index)  !DMS conc in seawater
      JDMS_OFLUX    =SI(133,Sect_No,im_index)  !DMS flux from ocean
      JOCFF_EM      =SI(134,Sect_No,im_index)  !Fresh OCFF surf emiss
      JOCFF_HILEM   =SI(135,Sect_No,im_index)  !Fresh OCFF high emiss

! Aerosol climatologies
      JARCLBIOG_BG  =SI(351,Sect_No,im_index)  ! Biogenic aerosol climatology
      JARCLBIOM_FR  =SI(352,Sect_No,im_index)  ! Biomass burning (fresh) aerosol clim
      JARCLBIOM_AG  =SI(353,Sect_No,im_index)  ! Biomass burning (aged) aerosol clim
      JARCLBIOM_IC  =SI(354,Sect_No,im_index)  ! Biomass burning (in-cloud) aerosol clim
      JARCLBLCK_FR  =SI(355,Sect_No,im_index)  ! Black carbon (fresh) aerosol clim
      JARCLBLCK_AG  =SI(356,Sect_No,im_index)  ! Black carbon (aged) aerosol clim
      JARCLSSLT_FI  =SI(357,Sect_No,im_index)  ! Sea salt (film mode) aerosol clim 
      JARCLSSLT_JT  =SI(358,Sect_No,im_index)  ! Sea salt (jet mode) aerosol clim
      JARCLSULP_AC  =SI(359,Sect_No,im_index)  ! Sulphate (accumulation mode) aero clim
      JARCLSULP_AK  =SI(360,Sect_No,im_index)  ! Sulphate (Aitken mode) aerosol clim 
      JARCLSULP_DI  =SI(361,Sect_No,im_index)  ! Sulphate (dissolved) aerosol clim
      JARCLDUST_B1  =SI(362,Sect_No,im_index)  ! Dust (bin 1) aerosol climatology 
      JARCLDUST_B2  =SI(363,Sect_No,im_index)  ! Dust (bin 2) aerosol climatology 
      JARCLDUST_B3  =SI(364,Sect_No,im_index)  ! Dust (bin 3) aerosol climatology 
      JARCLDUST_B4  =SI(365,Sect_No,im_index)  ! Dust (bin 4) aerosol climatology 
      JARCLDUST_B5  =SI(366,Sect_No,im_index)  ! Dust (bin 5) aerosol climatology 
      JARCLDUST_B6  =SI(367,Sect_No,im_index)  ! Dust (bin 6) aerosol climatology 
      JARCLOCFF_FR  =SI(368,Sect_No,im_index)  ! Org carbon fossil fuel (fresh) aero clim
      JARCLOCFF_AG  =SI(369,Sect_No,im_index)  ! Org carbon fossil fuel (aged) aero clim
      JARCLOCFF_IC  =SI(370,Sect_No,im_index)  ! Org carbon fossil fuel (in-cloud) aero clim
      JARCLDLTA_DL  =SI(371,Sect_No,im_index)  ! Delta aerosol climatology

      DO LEV=2,MODEL_LEVELS
        JARCLBIOG_BG(LEV) = JARCLBIOG_BG(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_FR(LEV) = JARCLBIOM_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_AG(LEV) = JARCLBIOM_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLBIOM_IC(LEV) = JARCLBIOM_IC(LEV-1)+THETA_FIELD_SIZE
        JARCLBLCK_FR(LEV) = JARCLBLCK_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLBLCK_AG(LEV) = JARCLBLCK_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLSSLT_FI(LEV) = JARCLSSLT_FI(LEV-1)+THETA_FIELD_SIZE
        JARCLSSLT_JT(LEV) = JARCLSSLT_JT(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_AC(LEV) = JARCLSULP_AC(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_AK(LEV) = JARCLSULP_AK(LEV-1)+THETA_FIELD_SIZE
        JARCLSULP_DI(LEV) = JARCLSULP_DI(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B1(LEV) = JARCLDUST_B1(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B2(LEV) = JARCLDUST_B2(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B3(LEV) = JARCLDUST_B3(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B4(LEV) = JARCLDUST_B4(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B5(LEV) = JARCLDUST_B5(LEV-1)+THETA_FIELD_SIZE
        JARCLDUST_B6(LEV) = JARCLDUST_B6(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_FR(LEV) = JARCLOCFF_FR(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_AG(LEV) = JARCLOCFF_AG(LEV-1)+THETA_FIELD_SIZE
        JARCLOCFF_IC(LEV) = JARCLOCFF_IC(LEV-1)+THETA_FIELD_SIZE
        JARCLDLTA_DL(LEV) = JARCLDLTA_DL(LEV-1)+THETA_FIELD_SIZE
      END DO

! Mineral dust scheme

      JSOIL_CLAY   =SI(418,Sect_No,im_index)  ! soil clay fraction
      JSOIL_SILT   =SI(419,Sect_No,im_index)  ! soil silt fraction
      JSOIL_SAND   =SI(420,Sect_No,im_index)  ! soil sand fraction

      JDUST_MREL1=SI(421,Sect_No,im_index) !relative soil mass in div1
      JDUST_MREL2=SI(422,Sect_No,im_index) !relative soil mass in div2
      JDUST_MREL3=SI(423,Sect_No,im_index) !relative soil mass in div3
      JDUST_MREL4=SI(424,Sect_No,im_index) !relative soil mass in div4
      JDUST_MREL5=SI(425,Sect_No,im_index) !relative soil mass in div5
      JDUST_MREL6=SI(426,Sect_No,im_index) !relative soil mass in div6


      JDUST_DIV1(1)=SI(431,Sect_No,im_index)  ! dust mmr, division 1
      JDUST_DIV2(1)=SI(432,Sect_No,im_index)  ! dust mmr, division 2
      JDUST_DIV3(1)=SI(433,Sect_No,im_index)  ! dust mmr, division 3
      JDUST_DIV4(1)=SI(434,Sect_No,im_index)  ! dust mmr, division 4
      JDUST_DIV5(1)=SI(435,Sect_No,im_index)  ! dust mmr, division 5
      JDUST_DIV6(1)=SI(436,Sect_No,im_index)  ! dust mmr, division 6

      DO LEV = 2,MODEL_LEVELS
       JDUST_DIV1(LEV)=JDUST_DIV1(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV2(LEV)=JDUST_DIV2(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV3(LEV)=JDUST_DIV3(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV4(LEV)=JDUST_DIV4(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV5(LEV)=JDUST_DIV5(LEV-1)+THETA_OFF_SIZE
       JDUST_DIV6(LEV)=JDUST_DIV6(LEV-1)+THETA_OFF_SIZE
      ENDDO

! HadCM2 sulphate loading patterns
      JHadCM2_SO4(1)=SI(160,Sect_No,im_index)
      DO LEV=2, NSULPAT
        JHadCM2_SO4(LEV)=JHadCM2_SO4(LEV-1)+THETA_FIELD_SIZE
      ENDDO

! Add for Carbon cycle
      J_CO2FLUX = SI(250,Sect_No,im_index)
      J_CO2_EMITS  = SI(251,Sect_No,im_index)
      JCO2(1)      = SI(252,Sect_No,im_index)
      DO LEV=2,MODEL_LEVELS
        JMURK_SOURCE(LEV) = JMURK_SOURCE(LEV-1)+THETA_FIELD_SIZE
        JMURK(LEV) = JMURK(LEV-1)+THETA_OFF_SIZE

! For Sulphur Cycle variables
        JSO2(LEV)=JSO2(LEV-1)+THETA_OFF_SIZE
        JDMS(LEV)=JDMS(LEV-1)+THETA_OFF_SIZE
        JSO4_AITKEN(LEV)=JSO4_AITKEN(LEV-1)+THETA_OFF_SIZE
        JSO4_ACCU(LEV)=JSO4_ACCU(LEV-1)+THETA_OFF_SIZE
        JSO4_DISS(LEV)=JSO4_DISS(LEV-1)+THETA_OFF_SIZE
        JH2O2(LEV)=JH2O2(LEV-1)+THETA_OFF_SIZE
        JSO2_NATEM(LEV)=JSO2_NATEM(LEV-1)+THETA_FIELD_SIZE
        JOH(LEV) = JOH(LEV-1)+THETA_FIELD_SIZE
        JHO2(LEV) = JHO2(LEV-1)+THETA_FIELD_SIZE
        JNH3(LEV)      = JNH3(LEV-1)+THETA_OFF_SIZE
        JSOOT_NEW(LEV) = JSOOT_NEW(LEV-1)+THETA_OFF_SIZE
        JSOOT_AGD(LEV) = JSOOT_AGD(LEV-1)+THETA_OFF_SIZE
        JSOOT_CLD(LEV) = JSOOT_CLD(LEV-1)+THETA_OFF_SIZE
        JBMASS_NEW(LEV) = JBMASS_NEW(LEV-1)+THETA_OFF_SIZE
        JBMASS_AGD(LEV) = JBMASS_AGD(LEV-1)+THETA_OFF_SIZE
        JBMASS_CLD(LEV) = JBMASS_CLD(LEV-1)+THETA_OFF_SIZE
        JOCFF_NEW(LEV) = JOCFF_NEW(LEV-1)+THETA_OFF_SIZE
        JOCFF_AGD(LEV) = JOCFF_AGD(LEV-1)+THETA_OFF_SIZE
        JOCFF_CLD(LEV) = JOCFF_CLD(LEV-1)+THETA_OFF_SIZE
        JH2O2_LIMIT(LEV)=JH2O2_LIMIT(LEV-1)+THETA_FIELD_SIZE
        JO3_CHEM(LEV)=JO3_CHEM(LEV-1)+THETA_FIELD_SIZE

! For Carbon Cycle variables
        JCO2(LEV)=JCO2(LEV-1)+THETA_OFF_SIZE
      END DO

! Tracer fluxes - kdcorbin, 06/10
      tflux_sect = 3
      JTRACER_FLUX1 = SI(100,tflux_sect,im_index)
      JTRACER_FLUX2 = SI(101,tflux_sect,im_index)
      JTRACER_FLUX3 = SI(102,tflux_sect,im_index)
      JTRACER_FLUX4 = SI(103,tflux_sect,im_index)
      JTRACER_FLUX5 = SI(104,tflux_sect,im_index)
      JTRACER_FLUX6 = SI(105,tflux_sect,im_index)
      JTRACER_FLUX7 = SI(106,tflux_sect,im_index)
      JTRACER_FLUX8 = SI(107,tflux_sect,im_index)
      JTRACER_FLUX9 = SI(108,tflux_sect,im_index)
      JTRACER_FLUX10= SI(109,tflux_sect,im_index)
      JTRACER_FLUX11= SI(110,tflux_sect,im_index)
      JTRACER_FLUX12= SI(111,tflux_sect,im_index)
      JTRACER_FLUX13= SI(112,tflux_sect,im_index)
      JTRACER_FLUX14= SI(113,tflux_sect,im_index)
      JTRACER_FLUX15= SI(114,tflux_sect,im_index)
      JTRACER_FLUX16= SI(115,tflux_sect,im_index)
      JTRACER_FLUX17= SI(116,tflux_sect,im_index)
      JTRACER_FLUX18= SI(117,tflux_sect,im_index)
      JTRACER_FLUX19= SI(118,tflux_sect,im_index)
      JTRACER_FLUX20= SI(119,tflux_sect,im_index)

!jhan(re-assign atm_pointers
! Lest Sept 2012 - changed items from 300s (user ancils) to 800s
#     ifdef CABLE_SOIL_LAYERS
      !--- allow for 6 layers
      JTSOIL_TILE(1)   = SI(801,Sect_No,im_index)
      JTSOIL_TILE(2)   = SI(802,Sect_No,im_index)
      JTSOIL_TILE(3)   = SI(803,Sect_No,im_index)
      JTSOIL_TILE(4)   = SI(804,Sect_No,im_index)
      JTSOIL_TILE(5)   = SI(805,Sect_No,im_index)
      JTSOIL_TILE(6)   = SI(806,Sect_No,im_index)
      JSMCL_TILE(1)    = SI(807,Sect_No,im_index)
      JSMCL_TILE(2)    = SI(808,Sect_No,im_index)
      JSMCL_TILE(3)    = SI(809,Sect_No,im_index)
      JSMCL_TILE(4)    = SI(810,Sect_No,im_index)
      JSMCL_TILE(5)    = SI(811,Sect_No,im_index)
      JSMCL_TILE(6)    = SI(812,Sect_No,im_index)
      JSTHF_TILE(1)    = SI(813,Sect_No,im_index)
      JSTHF_TILE(2)    = SI(814,Sect_No,im_index)
      JSTHF_TILE(3)    = SI(815,Sect_No,im_index)
      JSTHF_TILE(4)    = SI(816,Sect_No,im_index)
      JSTHF_TILE(5)    = SI(817,Sect_No,im_index)
      JSTHF_TILE(6)    = SI(818,Sect_No,im_index)
      JSNOW_DEPTH3L(1) = SI(819,Sect_No,im_index)
      JSNOW_DEPTH3L(2) = SI(820,Sect_No,im_index)
      JSNOW_DEPTH3L(3) = SI(821,Sect_No,im_index)
      JSNOW_MASS3L(1)  = SI(822,Sect_No,im_index)
      JSNOW_MASS3L(2)  = SI(823,Sect_No,im_index)
      JSNOW_MASS3L(3)  = SI(824,Sect_No,im_index)
      JSNOW_TMP3L(1)   = SI(825,Sect_No,im_index)
      JSNOW_TMP3L(2)   = SI(826,Sect_No,im_index)
      JSNOW_TMP3L(3)   = SI(827,Sect_No,im_index)
      JSNOW_RHO3L(1)   = SI(828,Sect_No,im_index)
      JSNOW_RHO3L(2)   = SI(829,Sect_No,im_index)
      JSNOW_RHO3L(3)   = SI(830,Sect_No,im_index)
      JSNOW_RHO1L      = SI(831,Sect_No,im_index)
      JSNOW_AGE        = SI(832,Sect_No,im_index)
      JSNOW_FLG3l      = SI(833,Sect_No,im_index)
#     else
      !--- original 4 layer structure
      JTSOIL_TILE(1)   = SI(801,Sect_No,im_index)
      JTSOIL_TILE(2)   = SI(802,Sect_No,im_index)
      JTSOIL_TILE(3)   = SI(803,Sect_No,im_index)
      JTSOIL_TILE(4)   = SI(804,Sect_No,im_index)
      JSMCL_TILE(1)    = SI(805,Sect_No,im_index)
      JSMCL_TILE(2)    = SI(806,Sect_No,im_index)
      JSMCL_TILE(3)    = SI(807,Sect_No,im_index)
      JSMCL_TILE(4)    = SI(808,Sect_No,im_index)
      JSTHF_TILE(1)    = SI(809,Sect_No,im_index)
      JSTHF_TILE(2)    = SI(810,Sect_No,im_index)
      JSTHF_TILE(3)    = SI(811,Sect_No,im_index)
      JSTHF_TILE(4)    = SI(812,Sect_No,im_index)
      JSNOW_DEPTH3L(1) = SI(813,Sect_No,im_index)
      JSNOW_DEPTH3L(2) = SI(814,Sect_No,im_index)
      JSNOW_DEPTH3L(3) = SI(815,Sect_No,im_index)
      JSNOW_MASS3L(1)  = SI(816,Sect_No,im_index)
      JSNOW_MASS3L(2)  = SI(817,Sect_No,im_index)
      JSNOW_MASS3L(3)  = SI(818,Sect_No,im_index)
      JSNOW_TMP3L(1)   = SI(819,Sect_No,im_index)
      JSNOW_TMP3L(2)   = SI(820,Sect_No,im_index)
      JSNOW_TMP3L(3)   = SI(821,Sect_No,im_index)
      JSNOW_RHO3L(1)   = SI(822,Sect_No,im_index)
      JSNOW_RHO3L(2)   = SI(823,Sect_No,im_index)
      JSNOW_RHO3L(3)   = SI(824,Sect_No,im_index)
      JSNOW_RHO1L      = SI(825,Sect_No,im_index)
      JSNOW_AGE        = SI(826,Sect_No,im_index)
      JSNOW_FLG3l      = SI(827,Sect_No,im_index)
#     endif

! Lestevens Sept 2012 - START CASA
! sect_no equals 0
      JCPOOL_TILE(1)  =SI(851,Sect_No,im_index)  !
      JCPOOL_TILE(2)  =SI(852,Sect_No,im_index)  !
      JCPOOL_TILE(3)  =SI(853,Sect_No,im_index)  !
      JCPOOL_TILE(4)  =SI(854,Sect_No,im_index)  !
      JCPOOL_TILE(5)  =SI(855,Sect_No,im_index)  !
      JCPOOL_TILE(6)  =SI(856,Sect_No,im_index)  !
      JCPOOL_TILE(7)  =SI(857,Sect_No,im_index)  !
      JCPOOL_TILE(8)  =SI(858,Sect_No,im_index)  !
      JCPOOL_TILE(9)  =SI(859,Sect_No,im_index)  !
      JCPOOL_TILE(10) =SI(860,Sect_No,im_index)  !
      JNPOOL_TILE(1)  =SI(861,Sect_No,im_index)  !
      JNPOOL_TILE(2)  =SI(862,Sect_No,im_index)  !
      JNPOOL_TILE(3)  =SI(863,Sect_No,im_index)  !
      JNPOOL_TILE(4)  =SI(864,Sect_No,im_index)  !
      JNPOOL_TILE(5)  =SI(865,Sect_No,im_index)  !
      JNPOOL_TILE(6)  =SI(866,Sect_No,im_index)  !
      JNPOOL_TILE(7)  =SI(867,Sect_No,im_index)  !
      JNPOOL_TILE(8)  =SI(868,Sect_No,im_index)  !
      JNPOOL_TILE(9)  =SI(869,Sect_No,im_index)  !
      JNPOOL_TILE(10) =SI(870,Sect_No,im_index)  !
      JPPOOL_TILE(1)  =SI(871,Sect_No,im_index)  !
      JPPOOL_TILE(2)  =SI(872,Sect_No,im_index)  !
      JPPOOL_TILE(3)  =SI(873,Sect_No,im_index)  !
      JPPOOL_TILE(4)  =SI(874,Sect_No,im_index)  !
      JPPOOL_TILE(5)  =SI(875,Sect_No,im_index)  !
      JPPOOL_TILE(6)  =SI(876,Sect_No,im_index)  !
      JPPOOL_TILE(7)  =SI(877,Sect_No,im_index)  !
      JPPOOL_TILE(8)  =SI(878,Sect_No,im_index)  !
      JPPOOL_TILE(9)  =SI(879,Sect_No,im_index)  !
      JPPOOL_TILE(10) =SI(880,Sect_No,im_index)  !
      JPPOOL_TILE(11) =SI(881,Sect_No,im_index)  !
      JPPOOL_TILE(12) =SI(882,Sect_No,im_index)  !
      JSOIL_ORDER     =SI(883,Sect_No,im_index)  !
      JNIDEP          =SI(884,Sect_No,im_index)  !
      JNIFIX          =SI(885,Sect_No,im_index)  !
      JPDUST          =SI(887,Sect_No,im_index)  !
      JPWEA           =SI(888,Sect_No,im_index)  !
      JGLAI           =SI(893,Sect_No,im_index)  !
      JPHENPHASE      =SI(894,Sect_No,im_index)  !

!     JNFERT          =SI(886,Sect_No,im_index)  !
!     JPFERT          =SI(889,Sect_No,im_index)  !
!     JCBAL_SUM       =SI(890,Sect_No,im_index)  !
!     JNBAL_SUM       =SI(891,Sect_No,im_index)  !
!     JPBAL_SUM       =SI(892,Sect_No,im_index)  !
! Lestevens Sept 2012 - END CASA

      JCABLE_LAI(1)   = SI(895,Sect_No,im_index)   

! Add user ancillaries

      JUSER_ANC1  = SI(301,Sect_No,im_index)
      JUSER_ANC2  = SI(302,Sect_No,im_index)
      JUSER_ANC3  = SI(303,Sect_No,im_index)
      JUSER_ANC4  = SI(304,Sect_No,im_index)
      JUSER_ANC5  = SI(305,Sect_No,im_index)
      JUSER_ANC6  = SI(306,Sect_No,im_index)
      JUSER_ANC7  = SI(307,Sect_No,im_index)
      JUSER_ANC8  = SI(308,Sect_No,im_index)
      JUSER_ANC9  = SI(309,Sect_No,im_index)
      JUSER_ANC10 = SI(310,Sect_No,im_index)
      JUSER_ANC11 = SI(311,Sect_No,im_index)
      JUSER_ANC12 = SI(312,Sect_No,im_index)
      JUSER_ANC13 = SI(313,Sect_No,im_index)
      JUSER_ANC14 = SI(314,Sect_No,im_index)
      JUSER_ANC15 = SI(315,Sect_No,im_index)
      JUSER_ANC16 = SI(316,Sect_No,im_index)
      JUSER_ANC17 = SI(317,Sect_No,im_index)
      JUSER_ANC18 = SI(318,Sect_No,im_index)
      JUSER_ANC19 = SI(319,Sect_No,im_index)
      JUSER_ANC20 = SI(320,Sect_No,im_index)
      JUSER_MULT1(1)  = SI(321,Sect_No,im_index)
      JUSER_MULT2(1)  = SI(322,Sect_No,im_index)
      JUSER_MULT3(1)  = SI(323,Sect_No,im_index)
      JUSER_MULT4(1)  = SI(324,Sect_No,im_index)
      JUSER_MULT5(1)  = SI(325,Sect_No,im_index)
      JUSER_MULT6(1)  = SI(326,Sect_No,im_index)
      JUSER_MULT7(1)  = SI(327,Sect_No,im_index)
      JUSER_MULT8(1)  = SI(328,Sect_No,im_index)
      JUSER_MULT9(1)  = SI(329,Sect_No,im_index)
      JUSER_MULT10(1) = SI(330,Sect_No,im_index)
      JUSER_MULT11(1) = SI(331,Sect_No,im_index)
      JUSER_MULT12(1) = SI(332,Sect_No,im_index)
      JUSER_MULT13(1) = SI(333,Sect_No,im_index)
      JUSER_MULT14(1) = SI(334,Sect_No,im_index)
      JUSER_MULT15(1) = SI(335,Sect_No,im_index)
      JUSER_MULT16(1) = SI(336,Sect_No,im_index)
      JUSER_MULT17(1) = SI(337,Sect_No,im_index)
      JUSER_MULT18(1) = SI(338,Sect_No,im_index)
      JUSER_MULT19(1) = SI(339,Sect_No,im_index)
      JUSER_MULT20(1) = SI(340,Sect_No,im_index)

! Set for multi-level user ancillaries
      DO LEV=2,MODEL_LEVELS
        JUSER_MULT1(LEV)  = JUSER_MULT1(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT2(LEV)  = JUSER_MULT2(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT3(LEV)  = JUSER_MULT3(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT4(LEV)  = JUSER_MULT4(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT5(LEV)  = JUSER_MULT5(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT6(LEV)  = JUSER_MULT6(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT7(LEV)  = JUSER_MULT7(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT8(LEV)  = JUSER_MULT8(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT9(LEV)  = JUSER_MULT9(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT10(LEV) = JUSER_MULT10(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT11(LEV) = JUSER_MULT11(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT12(LEV) = JUSER_MULT12(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT13(LEV) = JUSER_MULT13(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT14(LEV) = JUSER_MULT14(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT15(LEV) = JUSER_MULT15(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT16(LEV) = JUSER_MULT16(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT17(LEV) = JUSER_MULT17(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT18(LEV) = JUSER_MULT18(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT19(LEV) = JUSER_MULT19(LEV-1)+THETA_FIELD_SIZE
        JUSER_MULT20(LEV) = JUSER_MULT20(LEV-1)+THETA_FIELD_SIZE
      END DO

! Tiled vegetation and triffid
      JFRAC_TYP     = SI(216,Sect_No,im_index) ! surface type fractions
      JFRAC_CON1    = SI(442,Sect_No,im_index) ! surface type fractions
      JFRAC_CON2    = SI(443,Sect_No,im_index) ! surface type fractions
      JFRAC_CON3    = SI(444,Sect_No,im_index) ! surface type fractions
      JFRAC_CON4    = SI(445,Sect_No,im_index) ! surface type fractions
      JFRAC_CON5    = SI(446,Sect_No,im_index) ! surface type fractions
      JFRAC_CON6    = SI(447,Sect_No,im_index) ! surface type fractions
      JFRAC_CON7    = SI(448,Sect_No,im_index) ! surface type fractions
      JFRAC_CON8    = SI(449,Sect_No,im_index) ! surface type fractions
      JFRAC_CON9    = SI(450,Sect_No,im_index) ! surface type fractions
      JLAI_PFT      = SI(217,Sect_No,im_index) ! leaf area index of PFTs
      JCANHT_PFT    = SI(218,Sect_No,im_index) ! canopy height of PFTs
      JDISTURB      = SI(219,Sect_No,im_index) ! Veg disturbed fraction
      JSOIL_ALB     = SI(220,Sect_No,im_index) ! Snow-free soil albedo
      JSOIL_CARB    = SI(223,Sect_No,im_index) ! Soil carbon content
      JSOIL_CARB1   = SI(466,Sect_No,im_index) ! Soil carbon content DPM
      JSOIL_CARB2   = SI(467,Sect_No,im_index) ! Soil carbon content RPM
      JSOIL_CARB3   = SI(468,Sect_No,im_index) ! Soil carbon content BIO
      JSOIL_CARB4   = SI(469,Sect_No,im_index) ! Soil carbon content HUM
      JNPP_PFT_ACC  = SI(224,Sect_No,im_index) ! Accumulated NPP on PFTs
      JG_LF_PFT_ACC = SI(225,Sect_No,im_index) ! Accumulated leaf
!                                              ! turnover rate on PFTs
      JG_PHLF_PFT_ACC=SI(226,Sect_No,im_index) ! Accumulat. phenological
!                                              ! leaf turnover rate PFTs
      JRSP_W_PFT_ACC= SI(227,Sect_No,im_index) ! Accum. wood resp PFTs
      JRSP_S_ACC    = SI(228,Sect_No,im_index) ! Accumulated soil resp
      JRSP_S_ACC1 = SI(470,Sect_No,im_index)  ! Soil respiration DPM
      JRSP_S_ACC2 = SI(471,Sect_No,im_index)  ! Soil respiration RPM
      JRSP_S_ACC3 = SI(472,Sect_No,im_index)  ! Soil respiration BIO
      JRSP_S_ACC4 = SI(473,Sect_No,im_index)  ! Soil respiration HUM
      JCAN_WATER_TILE=SI(229,Sect_No,im_index) ! Canopy water content
!                                              ! on tiles
      JCATCH_TILE   = SI(230,Sect_No,im_index) ! Canopy capacity on
!                                              ! tiles
      JRGRAIN_TILE  = SI(231,Sect_No,im_index) ! Snow grain size on
!                                              ! tiles
      JTSTAR_TILE   = SI(233,Sect_No,im_index) ! Tiled surface temp
      JZ0_TILE      = SI(234,Sect_No,im_index) ! Tiled surface roughness
! Stash number for snow tile not finalised yet
      JSNODEP_TILE  = SI(240,Sect_No,im_index) ! Tiled snow depth
      JINFIL_TILE   = SI(236,Sect_No,im_index) ! Max tile infilt rate
! Stash codes for DORL, LW_DOWN, SW_TILE not finalised yet
      JDOLR         = SI(239,Sect_No,im_index) ! TOA surface up LW
      JLW_DOWN      = SI(238,Sect_No,im_index) ! Surface down LW
      JSW_TILE      = SI(237,Sect_No,im_index) ! Surface net SW on tiles


! River routing fields
      JRIV_SEQUENCE  = SI(151,Sect_No,im_index) ! River sequence
      JRIV_DIRECTION = SI(152,Sect_No,im_index) ! River Direction
      JRIV_STORAGE   = SI(153,Sect_No,im_index) ! River Water Storage
      JTOT_SURFROFF  = SI(155,Sect_No,im_index) ! Acc. surface runoff
      JTOT_SUBROFF   = SI(156,Sect_No,im_index) ! Acc. sub-surf runoff
! Set pointer for inland basin outflow
      JRIV_INLANDATM    = SI(511,Sect_No,im_index)
       !Inland basin outflow



! required for energy correction
      JNET_FLUX=SI(222,Sect_No,im_index)    ! store for energy flux
      JNET_MFLUX=SI(235,Sect_No,im_index)   ! store for moisture flux

!Fields carried forward from previous version.
      JTSTAR_ANOM    = SI(39,Sect_No,im_index)

#if defined (SLAB)
!     Set up pointers for slab model fields
      im_ident = slab_im
      im_index = internal_model_index(im_ident)

      JTCLIM = SI(178,Sect_No,im_index)      !  Ref SST
      JHCLIM = SI(179,Sect_No,im_index)      !  Clim SeaIce Depth
      JTSLAB = SI(210,Sect_No,im_index)      !  Slab temperature
      JCHEAT = SI(142,Sect_No,im_index)      !  Caryheat I2O
      JOIFLX = SI(143,Sect_No,im_index)      !  OI heat flux
      JUICE = SI(292,Sect_No,im_index)       !  X-comp of ice velocity
      JVICE = SI(293,Sect_No,im_index)       !  Y-comp of ice velocity
      JSIG11NE = SI(280,Sect_No,im_index)    !  Internal stresses for
      JSIG11SE = SI(281,Sect_No,im_index)    !  EVP ice dynamics
      JSIG11SW = SI(282,Sect_No,im_index)
      JSIG11NW = SI(283,Sect_No,im_index)
      JSIG12NE = SI(284,Sect_No,im_index)
      JSIG12SE = SI(285,Sect_No,im_index)
      JSIG12SW = SI(286,Sect_No,im_index)
      JSIG12NW = SI(287,Sect_No,im_index)
      JSIG22NE = SI(288,Sect_No,im_index)
      JSIG22SE = SI(289,Sect_No,im_index)
      JSIG22SW = SI(290,Sect_No,im_index)
      JSIG22NW = SI(291,Sect_No,im_index)

      im_ident  = atmos_im
      im_index  = internal_model_index(im_ident)
#endif



! Set all pointers referencing Lateral Boundary Conditions

      Sect_No=31  ! LBC section

      JOROG_LBC     = SI(1,Sect_No,im_index)
      JU_LBC        = SI(2,Sect_No,im_index)
      JV_LBC        = SI(3,Sect_No,im_index)
      JW_LBC        = SI(4,Sect_No,im_index)
      JRHO_LBC      = SI(5,Sect_No,im_index)
      JTHETA_LBC    = SI(6,Sect_No,im_index)
      JQ_LBC        = SI(7,Sect_No,im_index)
      JQCL_LBC      = SI(8,Sect_No,im_index)
      JQCF_LBC      = SI(9,Sect_No,im_index)
      JEXNER_LBC    = SI(10,Sect_No,im_index)
      JU_ADV_LBC    = SI(11,Sect_No,im_index)
      JV_ADV_LBC    = SI(12,Sect_No,im_index)
      JW_ADV_LBC    = SI(13,Sect_No,im_index)
      JQCF2_LBC     = SI(14,Sect_No,im_index)
      JQRAIN_LBC    = SI(15,Sect_No,im_index)
      JQGRAUP_LBC   = SI(16,Sect_No,im_index)
      JCF_BULK_LBC  = SI(17,Sect_No,im_index)
      JCF_LIQUID_LBC= SI(18,Sect_No,im_index)
      JCF_FROZEN_LBC= SI(19,Sect_No,im_index)
      JMURK_LBC     = SI(20,Sect_No,im_index)

      JU_LBC_TEND     = SI(257,Sect_No,im_index)
      JV_LBC_TEND     = SI(258,Sect_No,im_index)
      JW_LBC_TEND     = SI(259,Sect_No,im_index)
      JRHO_LBC_TEND   = SI(260,Sect_No,im_index)
      JTHETA_LBC_TEND = SI(261,Sect_No,im_index)
      JQ_LBC_TEND     = SI(262,Sect_No,im_index)
      JQCL_LBC_TEND   = SI(263,Sect_No,im_index)
      JQCF_LBC_TEND   = SI(264,Sect_No,im_index)
      JEXNER_LBC_TEND = SI(265,Sect_No,im_index)
      JU_ADV_LBC_TEND = SI(266,Sect_No,im_index)
      JV_ADV_LBC_TEND = SI(267,Sect_No,im_index)
      JW_ADV_LBC_TEND = SI(268,Sect_No,im_index)
      JQCF2_LBC_TEND   = SI(269,Sect_No,im_index)
      JQRAIN_LBC_TEND  = SI(270,Sect_No,im_index)
      JQGRAUP_LBC_TEND = SI(271,Sect_No,im_index)
      JCF_BULK_LBC_TEND  = SI(272,Sect_No,im_index)
      JCF_LIQUID_LBC_TEND= SI(273,Sect_No,im_index)
      JCF_FROZEN_LBC_TEND= SI(274,Sect_No,im_index)
      JMURK_LBC_TEND   = SI(275,Sect_No,im_index)

      IF (TR_VARS  >   0) THEN
        ivar=1
        DO jvar= A_TRACER_FIRST,A_TRACER_LAST
          IF (SI(jvar,33,im_index)  /=  1) THEN
            JTRACER_LBC(ivar)=                                          &
     &        SI(100+jvar-A_TRACER_FIRST,Sect_no,im_index)
            JTRACER_LBC_TEND(ivar)=                                     &
     &        SI(356+jvar-A_TRACER_FIRST,Sect_no,im_index)
            ivar=ivar+1
          ENDIF
        ENDDO ! j
      ELSE
        JTRACER_LBC(1) = 1  ! ensure sensible address even if no tracers
        JTRACER_LBC_TEND(1) = 1
      ENDIF ! IF (TR_VARS  >   0)

! Tracer prognostics are now in section 33, not section 0
      sect_no = 33   ! tracers section
      JVAR=0         ! JVAR+1 is the current tracer to be found
      IF (TR_VARS >  0) THEN
        DO IVAR=A_TRACER_FIRST,A_TRACER_LAST
          IF(SI(IVAR,sect_no,im_index) /= 1) THEN ! tracer in use
            JVAR=JVAR+1
            JTRACER(1,JVAR) = SI(IVAR,sect_no,im_index)
            DO LEV=2,TR_LEVELS
              JTRACER(LEV,JVAR)=JTRACER(LEV-1,JVAR)+THETA_OFF_SIZE
            END DO
            A_TR_INDEX(IVAR-A_TRACER_FIRST+1)=JVAR
          END IF
        END DO
      ELSE
        JTRACER(1,1)=1   ! Ensure a sensible address even if no tracers
      ENDIF
      IF(JVAR /= TR_VARS) THEN
        WRITE(6,*) 'STATMPT: TR_VARS and SI are inconsistent'
        WRITE(6,*) 'TR_VARS=',TR_VARS,' .     But, SI implies :',JVAR
        CMESSAGE=  'STATMPT: TR_VARS and SI  inconsistent, see output'
        ICODE=100
        GOTO 9999 ! error return
      END IF
! Oxidant concentrations from UKCA for use in HadGEM sulphur
! cycle (these are in Section 33):
      IF(L_SULPC_ONLINE_OXIDANTS .AND. L_UKCA) THEN 
         JOH_UKCA(1)   =SI(251,sect_no,im_index) 
         JH2O2_UKCA(1) =SI(8,sect_no,im_index)   
         JHO2_UKCA(1)  =SI(252,sect_no,im_index) 
         JO3_UKCA(1)   =SI(1,sect_no,im_index)   
      END IF 

! UKCA tracer prognostics are in section 34.
      sect_no = 34   ! UKCA tracers section
      JVAR=0         ! JVAR+1 is the current tracer to be found
        DO IVAR=A_UKCA_FIRST,A_UKCA_LAST
          IF(SI(IVAR,sect_no,im_index) /= 1) THEN ! tracer in use
            JVAR=JVAR+1
            JTR_UKCA(1,JVAR) = SI(IVAR,sect_no,im_index)
            DO LEV=2,TR_LEVELS
              JTR_UKCA(LEV,JVAR)=JTR_UKCA(LEV-1,JVAR)+THETA_OFF_SIZE
            END DO
            A_UKCA_INDEX(IVAR-A_UKCA_FIRST+1)=JVAR
          END IF
        END DO
      IF (JVAR  ==  0) THEN
        JTR_UKCA(1,1)=1   ! Ensure a sensible address when no tracers
      ENDIF

! Set pointers to level dependent constants for atmosphere.
      JETATHETA      =1              ! eta_theta_levels(0:model_levels)
      JETARHO        =JETATHETA+model_levels+1 ! eta_rho_levels
      JRHCRIT        =JETARHO+model_levels+1   ! rhcrit
      JSOIL_THICKNESS=JRHCRIT+model_levels+1   ! soil level depths
! For height definition z(i,j,k) = zsea(k) + C(k)*zorog(i,j)
      Jzseak_theta   =JSOIL_THICKNESS+model_levels+1 ! zsea (theta levs)
      JCk_theta      =Jzseak_theta   +model_levels+1 ! C    (theta levs)
      Jzseak_rho     =JCk_theta      +model_levels+1 ! zsea (rho levs)
      JCk_rho        =Jzseak_rho     +model_levels+1 ! C    (rho levs)
! Set pointers to Row dependent constants for atmosphere.
      JPHI_INPUT_P       =1
      JPHI_INPUT_V       =JPHI_INPUT_P + global_rows
! Set pointers to Col dependent constants for atmosphere.      
      JLAMBDA_INPUT_P    =1
      JLAMBDA_INPUT_U    =JLAMBDA_INPUT_P + global_row_length

! Set pointers to Row and Col dependent constants for atmosphere.
      JPHI_INPUT_P       =1
      JPHI_INPUT_V       =JPHI_INPUT_P + global_rows
      JLAMBDA_INPUT_P    =1
      JLAMBDA_INPUT_U    =JLAMBDA_INPUT_P + global_row_length

! The following if block should only be required for test purposes
! during the transition of vn5.2. This ensures old values for
! setting blev/bulev/brlev on old style dumps (before 'smooth'
! algorithm introduced and new definitions introduced).
      if(A_LEN2_LEVDEPC  <=  4) then ! ie before 5.2 change
         JCk_rho=jetarho
         Jck_theta=jetatheta
         Jzseak_rho=jetarho
         Jzseak_theta=jetatheta
         write(6,*) 'SETCONA_CTL: WARNING. This dump has not been '//   &
     &   'reconfigured with new 5.2 set of level dependent constants:'
         write(6,*) 'PP headers will revert to 5.0/5.1 style '//        &
     &   'blev/bulev/brlev definitions'
      endif

 9999 CONTINUE ! ERROR GOTO point.
      RETURN
      END SUBROUTINE SET_ATM_POINTERS

#endif
