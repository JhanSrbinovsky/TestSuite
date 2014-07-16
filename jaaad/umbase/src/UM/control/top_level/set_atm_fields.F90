#if defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Points atmosphere fields to the appropriate sections of D1
!
! Subroutine Interface: 
SUBROUTINE Set_Atm_Fields ( &
! "argptra.h" contains jpointers
#include "argptra.h"
! "argsts.h" contains SI (STASH index array); used to check tracers
#include "argsts.h" 
      D1, LD1, ID1 )

       USE atm_fields_mod   ! atmosphere fields
       USE field_length_mod ! field_length function 

IMPLICIT NONE
!
! Description: 
!   Routine to point atmosphere fields to the appropriate sections of D1.
!   After calling this subroutine, the fields can be used directly without
!   referring to D1
!
! Method: 
!   Assumming SET_ATM_POINTERS has been called beforehand, this subroutine
!   points each field to an area of D1 starting at the corresponding
!   "jpointer" (at the first level) and ending at the "jpointer" (at the 
!   last level) plus the size of a single level of that field. 
!
!   Tracers are dealt with differently:   First the number of active tracers
!   is computed so that the correct sections of the corresponding tracer
!   jpointers can be used in pointing the tracer fields to D1.  If no tracers
!   are active then the fields are pointed to a dummy array
!
! Owner:  A. Treshansky
!
! Code description: 
!   Language:  Fortran 90.
!   This code is written to UM programming standards version 8.0.
!

! Subroutine arguments

#include "parvars.h"
#include "typsize.h"
#include "typptra.h"
#include "typsts.h"
! constants (N_INTERNAL_MODELS_MAX) in "csubmodl.h" are needed by "ccontrol.h"
#include "csubmodl.h"
! constants (NUNITS) in "chsunits.h" are needed by "ccontrol.h"
#include "chsunits.h"
! constants (L_3D_CCA) in "ccontrol.h" are needed to determine field sizes below
#include "ccontrol.h" 
! constants (A_TRACER_FIRST, etc.) in "ctracera.h" are needed for tracers
#include "ctracera.h"

      REAL,    TARGET, INTENT(IN) :: D1(LEN_TOT)
      LOGICAL, TARGET, INTENT(IN) :: LD1(LEN_TOT)
      INTEGER, TARGET, INTENT(IN) :: ID1(LEN_TOT)

! Local variables

      INTEGER :: nTracer ! loop counter over available tracers
      INTEGER :: nActiveTracers ! number of tracers actually being used

! End of header

! 1.0 Start of subroutine code; point fields to D1

!    atmospheric primary variables
     U         => D1(JU(1) : JU(model_levels)+u_off_size )
     V         => D1(JV(1) : JV(model_levels)+v_off_size )
     THETA     => D1(JTHETA(1) : JTHETA(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     Q         => D1(JQ(1) : JQ(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QCF       => D1(JQCF(1) : JQCF(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     TSTAR     => D1(JTSTAR : JTSTAR+field_length(theta_points,no_halo,1) )
     LAND      => D1(JLAND : JLAND+field_length(theta_points,no_halo,1) )
     OROGRAPHY => D1(JOROG : JOROG+field_length(theta_points,no_halo,1) )
     W         => D1(JW(0) : JW(model_levels)+theta_off_size)
     RHO       => D1(JRHO(1) : JRHO(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     QCL       => D1(JQCL(1) : JQCL(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QCF2      => D1(JQCF2(1) : JQCF2(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QRAIN     => D1(JQRAIN(1) : JQRAIN(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     QGRAUP    => D1(JQGRAUP(1) : JQGRAUP(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     EXNER_RHO_LEVELS => D1(JEXNER_RHO_LEVELS(1) : JEXNER_RHO_LEVELS(1)+ &
      field_length(theta_points,single_halo,model_levels+1) )

!    Coastal Tiling
     FRAC_LAND  => D1(JFRAC_LAND:JFRAC_LAND+field_length(land_points,no_halo,1))
     TSTAR_LAND => D1(JTSTAR_LAND:JTSTAR_LAND+ &
      field_length(theta_points,no_halo,1))
     TSTAR_SEA  => D1(JTSTAR_SEA:JTSTAR_SEA+ &
      field_length(theta_points,no_halo,1))
     TSTAR_SICE => D1(JTSTAR_SICE:JTSTAR_SICE+ &
      field_length(theta_points,no_halo,1))

!    SeaIce and Land albedos
     SICE_ALB => D1(JSICE_ALB : JSICE_ALB+field_length(theta_points,no_halo,1))
     LAND_ALB => D1(JLAND_ALB : JLAND_ALB+field_length(theta_points,no_halo,1))

!    Large-Scale hydrology
     TI_MEAN   => D1(JTI_MEAN:JTI_MEAN+field_length(land_points,no_halo,1))
     TI_SIG    => D1(JTI_SIG:JTI_SIG+field_length(land_points,no_halo,1))
     FEXP      => D1(JFEXP:JFEXP+field_length(land_points,no_halo,1))
     GAMMA_INT => D1(JGAMMA_INT:JGAMMA_INT+field_length(land_points,no_halo,1))
     FSFC_SAT  => D1(JFSFC_SAT:JFSFC_SAT+field_length(land_points,no_halo,1))
     F_WETLAND => D1(JF_WETLAND:JF_WETLAND+field_length(land_points,no_halo,1))
     WATER_TABLE => D1(JWATER_TABLE:JWATER_TABLE+ &
      field_length(land_points,no_halo,1))

     STHZW   => D1(JSTHZW  : JSTHZW  + field_length(land_points,no_halo,1) )
     A_FSAT  => D1(JA_FSAT : JA_FSAT + field_length(land_points,no_halo,1) )
     C_FSAT  => D1(JC_FSAT : JC_FSAT + field_length(land_points,no_halo,1) )
     A_FWET  => D1(JA_FWET : JA_FWET + field_length(land_points,no_halo,1) )
     C_FWET  => D1(JC_FWET : JC_FWET + field_length(land_points,no_halo,1) )

!    Optional atmospheric primary variables
     ZH        => D1(JZH : JZH+field_length(theta_points,no_halo,1) )
     U_ADV     => D1(JU_ADV(1) : JU_ADV(1)+ &
      field_length(u_points,extended_halo,model_levels) )
     V_ADV     => D1(JV_ADV(1) : JV_ADV(1)+ &
      field_length(v_points,extended_halo,model_levels) )
     W_ADV     => D1(JW_ADV(0) : JW_ADV(0)+ &
      field_length(theta_points,extended_halo,model_levels+1) )
     NTML      => D1(JNTML : JNTML+field_length(theta_points,no_halo,1) )
     NBDSC     => D1(JNBDSC : JNBDSC+field_length(theta_points,no_halo,1) )
     NTDSC     => D1(JNTDSC : JNTDSC+field_length(theta_points,no_halo,1) )
     CUMULUS   => D1(JCUMULUS : JCUMULUS+field_length(theta_points,no_halo,1) )
     T1_SD     => D1(JT1_SD : JT1_SD+field_length(theta_points,no_halo,1) )
     Q1_SD     => D1(JQ1_SD : JQ1_SD+field_length(theta_points,no_halo,1) )
     CF_AREA   => D1(JCF_AREA(1) : JCF_AREA(1)+ &
      field_length(theta_points,no_halo,wet_levels) )
     CF_BULK   => D1(JCF_BULK(1) : JCF_BULK(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     CF_LIQUID => D1(JCF_LIQUID(1) : JCF_LIQUID(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     CF_FROZEN => D1(JCF_FROZEN(1) : JCF_FROZEN(1)+ &
      field_length(theta_points,extended_halo,wet_levels) )
     ! size of cca varies according to L_3D_CCA
     CCA => D1(JCCA(1):JCCA(1)+field_length(theta_points,no_halo,n_cca_lev))
     CCB          => D1(JCCB : JCCB+field_length(theta_points,no_halo,1) )
     CCT          => D1(JCCT : JCCT+field_length(theta_points,no_halo,1) )
     CCLWP        => D1(JCCLWP : JCCLWP+field_length(theta_points,no_halo,1) )
     CANOPY_WATER => D1(JCANOPY_WATER : JCANOPY_WATER+ &
      field_length(land_points,no_halo,1) )
     LCBASE       => D1(JLCBASE : JLCBASE+field_length(theta_points,no_halo,1))
     CCW_RAD      => D1(JCCW_RAD(1) : JCCW_RAD(wet_levels)+theta_field_size)

!    Secondary Fields in D1
     EXNER_THETA_LEVELS => D1(JEXNER_THETA_LEVELS(1):JEXNER_THETA_LEVELS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     P => D1(JP(1):JP(1)+field_length(theta_points,single_halo,model_levels+1))
     P_THETA_LEVELS => D1(JP_THETA_LEVELS(1):JP_THETA_LEVELS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     PSTAR => D1(JPSTAR : JPSTAR+field_length(theta_points,no_halo,1) )
     SW_INCS => D1(JSW_INCS(0) : JSW_INCS(MODEL_LEVELS+1)+theta_field_size)
     LW_INCS => D1(JLW_INCS(0) : JLW_INCS(MODEL_LEVELS)+theta_field_size)

!    Direct PAR flux for STOCHEM
     DIRPAR => D1(JDIRPAR : JDIRPAR+field_length(theta_points,no_halo,1) )

!    Soil Fields
     SMCL => D1(JSMCL(1):JSMCL(1)+field_length(land_points,no_halo,sm_levels))
     DEEP_SOIL_TEMP => D1(J_DEEP_SOIL_TEMP(1):J_DEEP_SOIL_TEMP(1)+ &
      field_length(land_points,no_halo,sm_levels))
     VOL_SMC_WILT   => D1(JVOL_SMC_WILT:JVOL_SMC_WILT+ &
      field_length(land_points,no_halo,1))
     VOL_SMC_CRIT   => D1(JVOL_SMC_CRIT:JVOL_SMC_CRIT+ &
      field_length(land_points,no_halo,1))
     VOL_SMC_SAT    => D1(JVOL_SMC_SAT:JVOL_SMC_SAT+ &
      field_length(land_points,no_halo,1))
     SAT_SOIL_COND  => D1(JSAT_SOIL_COND:JSAT_SOIL_COND+ &
      field_length(land_points,no_halo,1))
     THERM_CAP  =>D1(JTHERM_CAP:JTHERM_CAP+field_length(land_points,no_halo,1))
     THERM_COND =>D1(JTHERM_COND:JTHERM_COND+ &
      field_length(land_points,no_halo,1))
     CLAPP_HORN =>D1(JCLAPP_HORN:JCLAPP_HORN+ &
      field_length(land_points,no_halo,1))
     SAT_SOILW_SUCTION => D1(JSAT_SOILW_SUCTION:JSAT_SOILW_SUCTION+ &
      field_length(land_points,no_halo,1))
     STHU => D1(JSTHU(1):JSTHU(1)+field_length(land_points,no_halo,sm_levels))
     STHF => D1(JSTHF(1):JSTHF(1)+field_length(land_points,no_halo,sm_levels))

!    Vegetation Fields
     Z0          => D1(JZ0:JZ0+field_length(theta_points,no_halo,1))
     VEG_FRAC    => D1(JVEG_FRAC:JVEG_FRAC+field_length(land_points,no_halo,1))
     ROOT_DEPTH  => D1(JROOT_DEPTH:JROOT_DEPTH+ &
      field_length(land_points,no_halo,1))
     SFA         => D1(JSFA:JSFA+field_length(land_points,no_halo,1))
     MDSA        => D1(JMDSA:JMDSA+field_length(land_points,no_halo,1))
     SURF_RESIST => D1(JSURF_RESIST:JSURF_RESIST+ &
      field_length(land_points,no_halo,1))
     SURF_CAP    => D1(JSURF_CAP:JSURF_CAP+field_length(land_points,no_halo,1))
     INFILT      => D1(JINFILT:JINFILT+field_length(land_points,no_halo,1))
     LAI         => D1(JLAI:JLAI+field_length(land_points,no_halo,1))
     CANHT       => D1(JCANHT:JCANHT+field_length(land_points,no_halo,1))
     GS          => D1(JGS:JGS+field_length(land_points,no_halo,1))

!    CABLE
     TSOIL_TILE  => D1(JTSOIL_TILE(1):JTSOIL_TILE(1)+  &
                         field_length(land_points,no_halo,sm_levels*ntiles))
     SMCL_TILE   => D1(JSMCL_TILE(1):JSMCL_TILE(1)+  &
                          field_length(land_points,no_halo,sm_levels*ntiles))
!  Not needed MRD
!!$     STHU_TILE   => D1(JSTHU_TILE:JSTHU_TILE+  &
!!$                          field_length(land_points,no_halo,sm_levels*ntiles))
     STHF_TILE   => D1(JSTHF_TILE(1):JSTHF_TILE(1)+  &
                   field_length(land_points,no_halo,sm_levels*ntiles))
! MRD - should be a parameter for number of snow levels here rather than 3
     SNOW_DEPTH3L=> D1(JSNOW_DEPTH3L(1):JSNOW_DEPTH3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_MASS3L => D1(JSNOW_MASS3L(1):JSNOW_MASS3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_TMP3L  => D1(JSNOW_TMP3L(1):JSNOW_TMP3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_RHO3L  => D1(JSNOW_RHO3L(1):JSNOW_RHO3L(1)+  &
                          field_length(land_points,no_halo,3*ntiles))
     SNOW_RHO1L  => D1(JSNOW_RHO1L:JSNOW_RHO1L+  &
                          field_length(land_points,no_halo,ntiles))
     SNOW_AGE    => D1(JSNOW_AGE:JSNOW_AGE+  &
                          field_length(land_points,no_halo,ntiles))
     SNOW_FLG3L  => D1(JSNOW_FLG3L:JSNOW_FLG3L+  &
                          field_length(land_points,no_halo,ntiles))
! Lestevens Sept 2012 - CASA-CNP
     CPOOL_TILE  => D1(JCPOOL_TILE(1):JCPOOL_TILE(1)+ &
                          field_length(land_points,no_halo,10*ntiles))
     NPOOL_TILE  => D1(JNPOOL_TILE(1):JNPOOL_TILE(1)+ &
                          field_length(land_points,no_halo,10*ntiles))
     PPOOL_TILE  => D1(JPPOOL_TILE(1):JPPOOL_TILE(1)+ &
                          field_length(land_points,no_halo,12*ntiles))
     SOIL_ORDER  => D1(JSOIL_ORDER:JSOIL_ORDER+ &
                          field_length(land_points,no_halo,1))
     NIDEP       => D1(JNIDEP:JNIDEP+ &
                          field_length(land_points,no_halo,1))
     NIFIX       => D1(JNIFIX:JNIFIX+ &
                          field_length(land_points,no_halo,1))
     PWEA        => D1(JPWEA:JPWEA+ &
                          field_length(land_points,no_halo,1))
     PDUST       => D1(JPDUST:JPDUST+ &
                          field_length(land_points,no_halo,1))
     GLAI        => D1(JGLAI:JGLAI+ &
                          field_length(land_points,no_halo,ntiles))
     PHENPHASE   => D1(JPHENPHASE:JPHENPHASE+ &
                          field_length(land_points,no_halo,ntiles))
     CABLE_LAI => D1(JCABLE_LAI(1):JCABLE_LAI(1)+  &
                         field_length(land_points,no_halo,ntiles)) 

!    Orography Fields
     OROG_SIL => D1(JOROG_SIL : JOROG_SIL+field_length(land_points,no_halo,1) )
     OROG_HO2 => D1(JOROG_HO2 : JOROG_HO2+field_length(land_points,no_halo,1) )
     OROG_SD => D1(JOROG_SD : JOROG_SD+field_length(land_points,no_halo,1) )
     OROG_GRAD_X  => D1(JOROG_GRAD_X : JOROG_GRAD_X+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_Y  => D1(JOROG_GRAD_Y : JOROG_GRAD_Y+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_XX => D1(JOROG_GRAD_XX : JOROG_GRAD_XX+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_XY => D1(JOROG_GRAD_XY : JOROG_GRAD_XY+ &
      field_length(land_points,no_halo,1) )
     OROG_GRAD_YY => D1(JOROG_GRAD_YY : JOROG_GRAD_YY+ &
      field_length(land_points,no_halo,1) )

!    Sea/Sea Ice Fields
     U_SEA => D1(JU_SEA : JU_SEA+field_length(u_points,no_halo,1) )
     V_SEA => D1(JV_SEA : JV_SEA+field_length(v_points,no_halo,1) )
     ICE_FRACTION  => D1(JICE_FRACTION  : JICE_FRACTION + &
      field_length(theta_points_sea_only,no_halo,1) )
     ICE_THICKNESS => D1(JICE_THICKNESS : JICE_THICKNESS+ &
      field_length(theta_points_sea_only,no_halo,1) )
     TI => D1(JTI : JTI+field_length(theta_points_sea_only,no_halo,1) )
     ICE_FRACT_CAT => D1(JICE_FRACT_CAT : JICE_FRACT_CAT+ &
      field_length(theta_points,no_halo,nice) )
     ICE_THICK_CAT => D1(JICE_THICK_CAT : JICE_THICK_CAT+ &
      field_length(theta_points,no_halo,nice) )
     TI_CAT => D1(JTI_CAT : JTI_CAT+field_length(theta_points,no_halo,nice) )
     U_0_P => D1(JU_0_P : JU_0_P+field_length(theta_points,no_halo,1) )
     V_0_P => D1(JV_0_P : JV_0_P+field_length(theta_points,no_halo,1) )

!    Snow Fields
     SNODEP => D1(JSNODEP : JSNODEP+field_length(theta_points,no_halo,1) )
     SNODEP_SEA => D1(JSNODEP_SEA : JSNODEP_SEA+ &
      field_length(theta_points_sea_only,no_halo,1) )
     SNODEP_SEA_CAT => D1(JSNODEP_SEA_CAT : JSNODEP_SEA_CAT+ &
      field_length(theta_points,no_halo,nice) )
! SNSOOT may not be used as of vn6.6
     SNSOOT => D1(JSNSOOT : JSNSOOT+field_length(theta_points,no_halo,1) )
     CATCH_SNOW => D1(JCATCH_SNOW : JCATCH_SNOW+ &
      field_length(land_points,no_halo,1) )
     SNOW_GRND => D1(JSNOW_GRND : JSNOW_GRND+ &
      field_length(land_points,no_halo,1) )

!    OZONE
     O3 => D1(JOZONE(1) : JOZONE(1)+ &
      field_length(ozone_points,no_halo,ozone_levels) )
!    Tropopause-based Ozone
     IF (TPPS_OZONE_LEVELS > 0) THEN
       TPPSOZONE => D1(JTPPSOZONE(1) : JTPPSOZONE(1)+ &
        field_length(ozone_points,no_halo,ozone_levels) )
     ELSE
       TPPSOZONE => dummy_field
     END IF

!    Ozone tracer field and cariolle parameters
     OZONE_TRACER => D1(JOZONE_TRACER(1)   : JOZONE_TRACER(1)+       &
                 field_length(theta_points,single_halo,model_levels))
     O3_PROD_LOSS => D1(JO3_PROD_LOSS(1)   : JO3_PROD_LOSS(1)+       &
                             (ROWS*MODEL_LEVELS))
     O3_P_L_VMR   => D1(JO3_P_L_VMR(1)   : JO3_P_L_VMR(1)+           &
                             (ROWS*MODEL_LEVELS))
     O3_VMR       => D1(JO3_VMR (1)   : JO3_VMR(1)+                  &
                             (ROWS*MODEL_LEVELS))
     O3_P_L_TEMP  => D1(JO3_P_L_TEMP(1)   : JO3_P_L_TEMP(1)+         &
                             (ROWS*MODEL_LEVELS))
     O3_TEMP      => D1(JO3_TEMP(1)   : JO3_TEMP(1)+                 &
                             (ROWS*MODEL_LEVELS))    
     O3_P_L_COLO3 => D1(JO3_P_L_COLO3(1)   : JO3_P_L_COLO3(1)+       &
                             (ROWS*MODEL_LEVELS))
     O3_COLO3     => D1(JO3_COLO3(1)   : JO3_COLO3(1)+               &
                             (ROWS*MODEL_LEVELS))
!    STOCHEM fields
     CH4_STOCH => D1(JCH4_STOCH(1) : JCH4_STOCH(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     O3_STOCH  => D1(JO3_STOCH(1)  : JO3_STOCH(1) + &
      field_length(theta_points,no_halo,model_levels) )

!    Sources and Aerosol Ancillaries
     MURK_SOURCE => D1(JMURK_SOURCE(1) : JMURK_SOURCE(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     SO2_EM      => D1(JSO2_EM : JSO2_EM+field_length(theta_points,no_halo,1) )
     DMS_EM      => D1(JDMS_EM : JDMS_EM+field_length(theta_points,no_halo,1) )
     MURK        => D1(JMURK(1) : JMURK(1)+ &
      field_length(theta_points,single_halo,model_levels) )

!    Sulphur cycle
     SO2         => D1(JSO2(1) : JSO2(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DMS         => D1(JDMS(1) : JDMS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO4_AITKEN  => D1(JSO4_AITKEN(1) : JSO4_AITKEN(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO4_ACCU    => D1(JSO4_ACCU(1) : JSO4_ACCU(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO4_DISS    => D1(JSO4_DISS(1) : JSO4_DISS(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     H2O2        => D1(JH2O2(1) : JH2O2(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     NH3         => D1(JNH3(1) : JNH3(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SOOT_NEW    => D1(JSOOT_NEW(1) : JSOOT_NEW(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SOOT_AGD    => D1(JSOOT_AGD(1) : JSOOT_AGD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SOOT_CLD    => D1(JSOOT_CLD(1) : JSOOT_CLD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     BMASS_NEW   => D1(JBMASS_NEW(1) : JBMASS_NEW(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     BMASS_AGD   => D1(JBMASS_AGD(1) : JBMASS_AGD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     BMASS_CLD   => D1(JBMASS_CLD(1) : JBMASS_CLD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     OCFF_NEW    => D1(JOCFF_NEW(1) : JOCFF_NEW(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     OCFF_AGD    => D1(JOCFF_AGD(1) : JOCFF_AGD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     OCFF_CLD    => D1(JOCFF_CLD(1) : JOCFF_CLD(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     SO2_NATEM   => D1(JSO2_NATEM(1) : JSO2_NATEM(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     OH          => D1(JOH(1) : JOH(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     HO2         => D1(JHO2(1) : JHO2(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     H2O2_LIMIT  => D1(JH2O2_LIMIT(1) : JH2O2_LIMIT(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     O3_CHEM     => D1(JO3_CHEM(1) : JO3_CHEM(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     SO2_HILEM   => D1(JSO2_HILEM : JSO2_HILEM+ &
      field_length(theta_points,no_halo,1) )
     NH3_EM      => D1(JNH3_EM : JNH3_EM+ &
      field_length(theta_points,no_halo,1) )
     SOOT_EM     => D1(JSOOT_EM : JSOOT_EM+ &
      field_length(theta_points,no_halo,1) )
     SOOT_HILEM  => D1(JSOOT_HILEM : JSOOT_HILEM+ &
      field_length(theta_points,no_halo,1) )
     BMASS_EM    => D1(JBMASS_EM : JBMASS_EM+ &
      field_length(theta_points,no_halo,1) )
     BMASS_HILEM => D1(JBMASS_HILEM : JBMASS_HILEM+ &
      field_length(theta_points,no_halo,1) )
     OCFF_EM     => D1(JOCFF_EM : JOCFF_EM+ &
      field_length(theta_points,no_halo,1) )
     OCFF_HILEM  => D1(JOCFF_HILEM : JOCFF_HILEM+ &
      field_length(theta_points,no_halo,1) )
     DMS_CONC    => D1(JDMS_CONC : JDMS_CONC+ &
      field_length(theta_points,no_halo,1) )
     DMS_OFLUX   => D1(JDMS_OFLUX : JDMS_OFLUX+ &
      field_length(theta_points,no_halo,1) )

! Aerosol climatologies
     ARCLBIOG_BG => D1(JARCLBIOG_BG(1) : JARCLBIOG_BG(model_levels)+theta_field_size)
     ARCLBIOM_FR => D1(JARCLBIOM_FR(1) : JARCLBIOM_FR(model_levels)+theta_field_size)
     ARCLBIOM_AG => D1(JARCLBIOM_AG(1) : JARCLBIOM_AG(model_levels)+theta_field_size)
     ARCLBIOM_IC => D1(JARCLBIOM_IC(1) : JARCLBIOM_IC(model_levels)+theta_field_size)
     ARCLBLCK_FR => D1(JARCLBLCK_FR(1) : JARCLBLCK_FR(model_levels)+theta_field_size)
     ARCLBLCK_AG => D1(JARCLBLCK_AG(1) : JARCLBLCK_AG(model_levels)+theta_field_size)
     ARCLSSLT_FI => D1(JARCLSSLT_FI(1) : JARCLSSLT_FI(model_levels)+theta_field_size)
     ARCLSSLT_JT => D1(JARCLSSLT_JT(1) : JARCLSSLT_JT(model_levels)+theta_field_size)
     ARCLSULP_AC => D1(JARCLSULP_AC(1) : JARCLSULP_AC(model_levels)+theta_field_size)
     ARCLSULP_AK => D1(JARCLSULP_AK(1) : JARCLSULP_AK(model_levels)+theta_field_size)
     ARCLSULP_DI => D1(JARCLSULP_DI(1) : JARCLSULP_DI(model_levels)+theta_field_size)
     ARCLDUST_B1 => D1(JARCLDUST_B1(1) : JARCLDUST_B1(model_levels)+theta_field_size)
     ARCLDUST_B2 => D1(JARCLDUST_B2(1) : JARCLDUST_B2(model_levels)+theta_field_size)
     ARCLDUST_B3 => D1(JARCLDUST_B3(1) : JARCLDUST_B3(model_levels)+theta_field_size)
     ARCLDUST_B4 => D1(JARCLDUST_B4(1) : JARCLDUST_B4(model_levels)+theta_field_size)
     ARCLDUST_B5 => D1(JARCLDUST_B5(1) : JARCLDUST_B5(model_levels)+theta_field_size)
     ARCLDUST_B6 => D1(JARCLDUST_B6(1) : JARCLDUST_B6(model_levels)+theta_field_size)
     ARCLOCFF_FR => D1(JARCLOCFF_FR(1) : JARCLOCFF_FR(model_levels)+theta_field_size)
     ARCLOCFF_AG => D1(JARCLOCFF_AG(1) : JARCLOCFF_AG(model_levels)+theta_field_size)
     ARCLOCFF_IC => D1(JARCLOCFF_IC(1) : JARCLOCFF_IC(model_levels)+theta_field_size)
     ARCLDLTA_DL => D1(JARCLDLTA_DL(1) : JARCLDLTA_DL(model_levels)+theta_field_size)
     
!    Mineral Dust Schema
     SOIL_CLAY  => D1(JSOIL_CLAY : JSOIL_CLAY+ &
      field_length(theta_points,no_halo,1) )
     SOIL_SILT  => D1(JSOIL_SILT : JSOIL_SILT+ &
      field_length(theta_points,no_halo,1) )
     SOIL_SAND  => D1(JSOIL_SAND : JSOIL_SAND+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL1 => D1(JDUST_MREL1 : JDUST_MREL1+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL2 => D1(JDUST_MREL2 : JDUST_MREL2+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL3 => D1(JDUST_MREL3 : JDUST_MREL3+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL4 => D1(JDUST_MREL4 : JDUST_MREL4+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL5 => D1(JDUST_MREL5 : JDUST_MREL5+ &
      field_length(theta_points,no_halo,1) )
     DUST_MREL6 => D1(JDUST_MREL6 : JDUST_MREL6+ &
      field_length(theta_points,no_halo,1) )
     DUST_DIV1  => D1(JDUST_DIV1(1) : JDUST_DIV1(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV2  => D1(JDUST_DIV2(1) : JDUST_DIV2(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV3  => D1(JDUST_DIV3(1) : JDUST_DIV3(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV4  => D1(JDUST_DIV4(1) : JDUST_DIV4(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV5  => D1(JDUST_DIV5(1) : JDUST_DIV5(1)+ &
      field_length(theta_points,single_halo,model_levels) )
     DUST_DIV6  => D1(JDUST_DIV6(1) : JDUST_DIV6(1)+ &
      field_length(theta_points,single_halo,model_levels) )

!    Carbon Cycle
     CO2FLUX   => D1(J_CO2FLUX : J_CO2FLUX+ &
      field_length(theta_points,no_halo,1) )
     CO2_EMITS => D1(J_CO2_EMITS : J_CO2_EMITS+ &
      field_length(theta_points,no_halo,1) )
     CO2       => D1(JCO2(1) : JCO2(model_levels)+theta_off_size)

!    level dependent constants
     zseak_theta => D1(JZSEAK_THETA : JZSEAK_THETA+(model_levels+1) )
     Ck_theta    => D1(JCK_THETA : JCK_THETA+(model_levels+1) )
     zseak_rho   => D1(JZSEAK_RHO : JZSEAK_RHO+(model_levels+1) )
     Ck_rho      => D1(JCK_RHO : JCK_RHO+(model_levels+1) )

!    Tracer Fluxes - kdcorbin, 05/10
     TRACER_FLUX1  => D1(JTRACER_FLUX1 : JTRACER_FLUX1 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX2  => D1(JTRACER_FLUX2 : JTRACER_FLUX2 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX3  => D1(JTRACER_FLUX3 : JTRACER_FLUX3 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX4  => D1(JTRACER_FLUX4 : JTRACER_FLUX4 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX5  => D1(JTRACER_FLUX5 : JTRACER_FLUX5 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX6  => D1(JTRACER_FLUX6 : JTRACER_FLUX6 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX7  => D1(JTRACER_FLUX7 : JTRACER_FLUX7 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX8  => D1(JTRACER_FLUX8 : JTRACER_FLUX8 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX9  => D1(JTRACER_FLUX9 : JTRACER_FLUX9 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX10 => D1(JTRACER_FLUX10: JTRACER_FLUX10 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX11 => D1(JTRACER_FLUX11: JTRACER_FLUX11 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX12 => D1(JTRACER_FLUX12: JTRACER_FLUX12 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX13 => D1(JTRACER_FLUX13: JTRACER_FLUX13 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX14 => D1(JTRACER_FLUX14: JTRACER_FLUX14 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX15 => D1(JTRACER_FLUX15: JTRACER_FLUX15 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX16 => D1(JTRACER_FLUX16: JTRACER_FLUX16 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX17 => D1(JTRACER_FLUX17: JTRACER_FLUX17 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX18 => D1(JTRACER_FLUX18: JTRACER_FLUX18 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX19 => D1(JTRACER_FLUX19: JTRACER_FLUX19 + &
      field_length(theta_points,no_halo,1) )
     TRACER_FLUX20 => D1(JTRACER_FLUX20: JTRACER_FLUX20 + &
      field_length(theta_points,no_halo,1) )

!    User ancillaries
     USER_ANC1   => D1(JUSER_ANC1  : JUSER_ANC1 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC2   => D1(JUSER_ANC2  : JUSER_ANC2 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC3   => D1(JUSER_ANC3  : JUSER_ANC3 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC4   => D1(JUSER_ANC4  : JUSER_ANC4 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC5   => D1(JUSER_ANC5  : JUSER_ANC5 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC6   => D1(JUSER_ANC6  : JUSER_ANC6 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC7   => D1(JUSER_ANC7  : JUSER_ANC7 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC8   => D1(JUSER_ANC8  : JUSER_ANC8 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC9   => D1(JUSER_ANC9  : JUSER_ANC9 + &
      field_length(theta_points,no_halo,1) )
     USER_ANC10  => D1(JUSER_ANC10 : JUSER_ANC10+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC11  => D1(JUSER_ANC11 : JUSER_ANC11+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC12  => D1(JUSER_ANC12 : JUSER_ANC12+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC13  => D1(JUSER_ANC13 : JUSER_ANC13+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC14  => D1(JUSER_ANC14 : JUSER_ANC14+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC15  => D1(JUSER_ANC15 : JUSER_ANC15+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC16  => D1(JUSER_ANC16 : JUSER_ANC16+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC17  => D1(JUSER_ANC17 : JUSER_ANC17+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC18  => D1(JUSER_ANC18 : JUSER_ANC18+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC19  => D1(JUSER_ANC19 : JUSER_ANC19+ &
      field_length(theta_points,no_halo,1) )
     USER_ANC20  => D1(JUSER_ANC20 : JUSER_ANC20+ &
      field_length(theta_points,no_halo,1) )
     USER_MULT1  => D1(JUSER_MULT1(1)  : JUSER_MULT1(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT2  => D1(JUSER_MULT2(1)  : JUSER_MULT2(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT3  => D1(JUSER_MULT3(1)  : JUSER_MULT3(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT4  => D1(JUSER_MULT4(1)  : JUSER_MULT4(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT5  => D1(JUSER_MULT5(1)  : JUSER_MULT5(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT6  => D1(JUSER_MULT6(1)  : JUSER_MULT6(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT7  => D1(JUSER_MULT7(1)  : JUSER_MULT7(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT8  => D1(JUSER_MULT8(1)  : JUSER_MULT8(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT9  => D1(JUSER_MULT9(1)  : JUSER_MULT9(1) + &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT10 => D1(JUSER_MULT10(1) : JUSER_MULT10(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT11 => D1(JUSER_MULT11(1) : JUSER_MULT11(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT12 => D1(JUSER_MULT12(1) : JUSER_MULT12(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT13 => D1(JUSER_MULT13(1) : JUSER_MULT13(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT14 => D1(JUSER_MULT14(1) : JUSER_MULT14(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT15 => D1(JUSER_MULT15(1) : JUSER_MULT15(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT16 => D1(JUSER_MULT16(1) : JUSER_MULT16(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT17 => D1(JUSER_MULT17(1) : JUSER_MULT17(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT18 => D1(JUSER_MULT18(1) : JUSER_MULT18(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT19 => D1(JUSER_MULT19(1) : JUSER_MULT19(1)+ &
      field_length(theta_points,no_halo,model_levels) )
     USER_MULT20 => D1(JUSER_MULT20(1) : JUSER_MULT20(1)+ &
      field_length(theta_points,no_halo,model_levels) )

!    Tiled vegetation and triffid
     FRAC_TYP  => D1(JFRAC_TYP:JFRAC_TYP+field_length(land_points,no_halo,1))
     FRAC_CON1 => D1(JFRAC_CON1:JFRAC_CON1+field_length(land_points,no_halo,1))
     FRAC_CON2 => D1(JFRAC_CON2:JFRAC_CON2+field_length(land_points,no_halo,1))
     FRAC_CON3 => D1(JFRAC_CON3:JFRAC_CON3+field_length(land_points,no_halo,1))
     FRAC_CON4 => D1(JFRAC_CON4:JFRAC_CON4+field_length(land_points,no_halo,1))
     FRAC_CON5 => D1(JFRAC_CON5:JFRAC_CON5+field_length(land_points,no_halo,1))
     FRAC_CON6 => D1(JFRAC_CON6:JFRAC_CON6+field_length(land_points,no_halo,1))
     FRAC_CON7 => D1(JFRAC_CON7:JFRAC_CON7+field_length(land_points,no_halo,1))
     FRAC_CON8 => D1(JFRAC_CON8:JFRAC_CON8+field_length(land_points,no_halo,1))
     FRAC_CON9 => D1(JFRAC_CON9:JFRAC_CON9+field_length(land_points,no_halo,1))
     LAI_PFT   => D1(JLAI_PFT:JLAI_PFT+field_length(land_points,no_halo,1))
     CANHT_PFT => D1(JCANHT_PFT:JCANHT_PFT+field_length(land_points,no_halo,1))
     DISTURB_VEG => D1(JDISTURB:JDISTURB+field_length(land_points,no_halo,1))
     SOIL_ALB  => D1(JSOIL_ALB:JSOIL_ALB+field_length(land_points,no_halo,1))
     SOIL_CARB => D1(JSOIL_CARB:JSOIL_CARB+field_length(land_points,no_halo,1))
     SOIL_CARB1 => D1(JSOIL_CARB1:JSOIL_CARB1+ &
      field_length(land_points,no_halo,1))
     SOIL_CARB2 => D1(JSOIL_CARB2:JSOIL_CARB2+ &
      field_length(land_points,no_halo,1))
     SOIL_CARB3 => D1(JSOIL_CARB3:JSOIL_CARB3+ &
      field_length(land_points,no_halo,1))
     SOIL_CARB4 => D1(JSOIL_CARB4:JSOIL_CARB4+ &
      field_length(land_points,no_halo,1))
     NPP_PFT_ACC    => D1(JNPP_PFT_ACC:JNPP_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     G_LF_PFT_ACC   => D1(JG_LF_PFT_ACC:JG_LF_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     G_PHLF_PFT_ACC => D1(JG_PHLF_PFT_ACC:JG_PHLF_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     RSP_W_PFT_ACC  => D1(JRSP_W_PFT_ACC:JRSP_W_PFT_ACC+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC  => D1(JRSP_S_ACC:JRSP_S_ACC+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC1 => D1(JRSP_S_ACC1:JRSP_S_ACC1+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC2 => D1(JRSP_S_ACC2:JRSP_S_ACC2+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC3 => D1(JRSP_S_ACC3:JRSP_S_ACC3+ &
      field_length(land_points,no_halo,1))
     RSP_S_ACC4 => D1(JRSP_S_ACC4:JRSP_S_ACC4+ &
      field_length(land_points,no_halo,1))
     CAN_WATER_TILE => D1(JCAN_WATER_TILE:JCAN_WATER_TILE+ &
      field_length(land_points,no_halo,1))
     CATCH_TILE  => D1(JCATCH_TILE:JCATCH_TILE+ &
      field_length(land_points,no_halo,1))
     RGRAIN_TILE => D1(JRGRAIN_TILE:JRGRAIN_TILE+ &
      field_length(land_points,no_halo,1))
     ! TSNOW => no longer used
     TSTAR_TILE  => D1(JTSTAR_TILE:JTSTAR_TILE+ &
      field_length(land_points,no_halo,1))
     Z0_TILE     => D1(JZ0_TILE:JZ0_TILE+field_length(land_points,no_halo,1))
     SNODEP_TILE => D1(JSNODEP_TILE:JSNODEP_TILE+ &
      field_length(land_points,no_halo,1))
     INFIL_TILE  => D1(JINFIL_TILE:JINFIL_TILE+ &
      field_length(land_points,no_halo,1))
     DOLR_FIELD  => D1(JDOLR:JDOLR+field_length(theta_points,no_halo,1))
     LW_DOWN     => D1(JLW_DOWN:JLW_DOWN+field_length(theta_points,no_halo,1))
     SW_TILE_RTS => D1(JSW_TILE:JSW_TILE+field_length(land_points,no_halo,1))

!    River routing fields
     RIV_SEQUENCE  => D1(JRIV_SEQUENCE : JRIV_SEQUENCE+ &
      field_length(river_points,no_halo,1) )
     RIV_DIRECTION => D1(JRIV_DIRECTION : JRIV_DIRECTION+ &
      field_length(river_points,no_halo,1) )
     RIV_STORAGE   => D1(JRIV_STORAGE : JRIV_STORAGE+ &
      field_length(river_points,no_halo,1) )
     TOT_SURFROFF  => D1(JTOT_SURFROFF : JTOT_SURFROFF+ &
      field_length(river_points,no_halo,1) )
     TOT_SUBROFF   => D1(JTOT_SUBROFF : JTOT_SUBROFF+ &
      field_length(river_points,no_halo,1) )
     RIV_INLANDATM => D1(JRIV_INLANDATM : JRIV_INLANDATM+ &
      field_length(land_points,no_halo,1) )
     ! these are unitialised upon entering ATM_STEP
     RIV_IAREA     => dummy_field !D1(1:1+row_length*rows)
     RIV_SLOPE     => dummy_field !D1(1:1+row_length*rows)
     RIV_FLOWOBS1  => dummy_field !D1(1:1+row_length*rows)
     RIV_INEXT     => dummy_field !D1(1:1+row_length*rows)
     RIV_JNEXT     => dummy_field !D1(1:1+row_length*rows)
     RIV_LAND      => dummy_field !D1(1:1+row_length*rows)
     RIV_SUBSTORE  => dummy_field !D1(1:1+row_length*rows)
     RIV_SURFSTORE => dummy_field !D1(1:1+row_length*rows)
     RIV_FLOWIN    => dummy_field !D1(1:1+row_length*rows)
     RIV_BFLOWIN   => dummy_field !D1(1:1+row_length*rows)

!    Fields to be tretained in dumps for coupled models using OASIS
#define ACCESS_SOLAR 1
#if ! defined(ACCESS_SOLAR)
!    gol124: because we restart coupling fields from oasis3
!    fiels we do not need to save them in um dump; this means
!    we don't need c_solar etc. to point into D1 - a regular
!    global array should suffice (advise from CH)
!    test it!: c_solar allocate in oasis_inita2o
     C_SOLAR => D1(JC_SOLAR : JC_SOLAR + &
                    field_length(theta_points,no_halo,1))

     C_BLUE =>  D1(JC_BLUE : JC_BLUE + &
                    field_length(theta_points,no_halo,1))
#endif

     C_DOWN =>  D1(JC_DOWN : JC_DOWN + &
                    field_length(theta_points,no_halo,1))

#if ! defined(ACCESS_SOLAR)
     C_LONGWAVE => D1(JC_LONGWAVE : JC_LONGWAVE + &
                    field_length(theta_points,no_halo,1))

     C_TAUX => D1(JC_TAUX : JC_TAUX + &
                    field_length(u_points,no_halo,1))

     C_TAUY => D1(JC_TAUY : JC_TAUY + &
                    field_length(v_points,no_halo,1))

     C_WINDMIX => D1(JC_WINDMIX : JC_WINDMIX + &
                    field_length(theta_points,no_halo,1))

     C_SENSIBLE => D1(JC_SENSIBLE : JC_SENSIBLE + &
                    field_length(theta_points,no_halo,1))

     C_SUBLIM =>  D1(JC_SUBLIM : JC_SUBLIM + &
                    field_length(theta_points,no_halo,1))

     C_EVAP =>  D1(JC_EVAP : JC_EVAP + &
                    field_length(theta_points,no_halo,1))

     C_BOTMELTN => D1(JC_BOTMELTN : JC_BOTMELTN + &
                    field_length(theta_points,no_halo,nice))

     C_TOPMELTN => D1(JC_TOPMELTN : JC_TOPMELTN + &
                    field_length(theta_points,no_halo,nice))

     C_LSRAIN =>  D1(JC_LSRAIN : JC_LSRAIN + &
                    field_length(theta_points,no_halo,1))

     C_LSSNOW =>  D1(JC_LSSNOW : JC_LSSNOW + &
                    field_length(theta_points,no_halo,1))

     C_CVRAIN =>  D1(JC_CVRAIN : JC_CVRAIN + &
                    field_length(theta_points,no_halo,1))

     C_CVSNOW =>  D1(JC_CVSNOW : JC_CVSNOW + &
                    field_length(theta_points,no_halo,1))

     C_RIVEROUT => D1(JC_RIVEROUT : JC_RIVEROUT + &
                    field_length(theta_points,no_halo,1))
#endif
 
#if defined(ACCESS)    
! gol124: auscom coupling
! probably not needed: coupling code seems to use def
! in include/argument/arg_atm_fields.h
#if ! defined(ACCESS_SOLAR)
     C_PRESS => D1(JC_PRESS : JC_PRESS + &
                    field_length(theta_points,no_halo,1))
#endif
#endif

!    Required for energy correction
     NET_FLUX  => D1(JNET_FLUX:JNET_FLUX+field_length(theta_points,no_halo,1))
     NET_MFLUX => D1(JNET_MFLUX:JNET_MFLUX+field_length(theta_points,no_halo,1))

!    Fields carried forward from previous version
     TSTAR_ANOM => D1(JTSTAR_ANOM : JTSTAR_ANOM+field_length(theta_points,no_halo,1) )

! SLAB model removed as of vn7.0
!#if defined (SLAB)
!
!     TCLIM => D1(JTCLIM : JTCLIM+field_length(theta_points,no_halo,1))
!     HCLIM => D1(JHCLIM : JHCLIM+field_length(theta_points,no_halo,1))
!     TSLAB => D1(JTSLAB : JTSLAB+field_length(theta_points_sea_only,no_halo,1))
!     CHEAT => D1(JCHEAT : JCHEAT+field_length(theta_points_sea_only,no_halo,1))
!     OIFLX => D1(JOIFLX : JOIFLX+field_length(theta_points_sea_only,no_halo,1))
!     UICE  => D1(JUICE  : JUICE+field_length(u_points,no_halo,1))
!     VICE  => D1(JVICE  : JVICE+field_length(v_points,no_halo,1))
!     SIG11NE => D1(JSIG11NE : JSIG11NE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG11SE => D1(JSIG11SE : JSIG11SE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG11SW => D1(JSIG11SW : JSIG11SW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG11NW => D1(JSIG11NW : JSIG11NW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12NE => D1(JSIG12NE : JSIG12NE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12SE => D1(JSIG12SE : JSIG12SE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12SW => D1(JSIG12SW : JSIG12SW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG12NW => D1(JSIG12NW : JSIG12NW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22NE => D1(JSIG22NE : JSIG22NE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22SE => D1(JSIG22SE : JSIG22SE+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22SW => D1(JSIG22SW : JSIG22SW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!     SIG22NW => D1(JSIG22NW : JSIG22NW+ &
!      field_length(theta_points_sea_only,no_halo,1) )
!
!#endif ! defined (SLAB)

!    lateral boundary conditions

     OROG_LBC  => D1(JOROG_LBC : JOROG_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*1) )
     U_LBC     => D1(JU_LBC : JU_LBC + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_LBC     => D1(JV_LBC : JU_LBC +  &
     (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_LBC     => D1(JW_LBC : JW_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)*(model_levels+1)) )
     RHO_LBC   => D1(JRHO_LBC : JRHO_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     THETA_LBC => D1(JTHETA_LBC : JTHETA_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     Q_LBC     => D1(JQ_LBC : JQ_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCL_LBC   => D1(JQCL_LBC : JQCL_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCF_LBC   => D1(JQCF_LBC : JQCF_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     EXNER_LBC => D1(JEXNER_LBC : JEXNER_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     U_ADV_LBC => D1(JU_ADV_LBC : JU_ADV_LBC + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_ADV_LBC => D1(JV_ADV_LBC : JV_ADV_LBC + &
      (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_ADV_LBC => D1(JW_ADV_LBC : JW_ADV_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     QCF2_LBC  => D1(JQCF2_LBC : JQCF2_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QRAIN_LBC => D1(JQRAIN_LBC : JQRAIN_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QGRAUP_LBC  => D1(JQGRAUP_LBC : JQGRAUP_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_BULK_LBC => D1(JCF_BULK_LBC : JCF_BULK_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_LIQUID_LBC => D1(JCF_LIQUID_LBC : JCF_LIQUID_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_FROZEN_LBC => D1(JCF_FROZEN_LBC : JCF_FROZEN_LBC + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     MURK_LBC  => D1(JMURK_LBC : JMURK_LBC + &
      (LENRIMA(fld_type_p,halo_type_single,1)*model_levels) )

     U_LBC_TEND     => D1(JU_LBC_TEND : JU_LBC_TEND + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_LBC_TEND     => D1(JV_LBC_TEND : JV_LBC_TEND + &
      (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_LBC_TEND     => D1(JW_LBC_TEND : JW_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,rima_type_norm)*(model_levels+1)) )
     RHO_LBC_TEND   => D1(JRHO_LBC_TEND : JRHO_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     THETA_LBC_TEND => D1(JTHETA_LBC_TEND : JTHETA_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*model_levels) )
     Q_LBC_TEND     => D1(JQ_LBC_TEND : JQ_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCL_LBC_TEND   => D1(JQCL_LBC_TEND : JQCL_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QCF_LBC_TEND   => D1(JQCF_LBC_TEND : JQCF_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     EXNER_LBC_TEND => D1(JEXNER_LBC_TEND : JEXNER_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     U_ADV_LBC_TEND => D1(JU_ADV_LBC_TEND : JU_ADV_LBC_TEND + &
      (LENRIMA(fld_type_u,halo_type_extended,1)*model_levels) )
     V_ADV_LBC_TEND => D1(JV_ADV_LBC_TEND : JV_ADV_LBC_TEND + &
      (LENRIMA(fld_type_v,halo_type_extended,1)*model_levels) )
     W_ADV_LBC_TEND => D1(JW_ADV_LBC_TEND : JW_ADV_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*(model_levels+1)) )
     QCF2_LBC_TEND => D1(JQCF2_LBC_TEND : JQCF2_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QRAIN_LBC_TEND => D1(JQRAIN_LBC_TEND : JQRAIN_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     QGRAUP_LBC_TEND => D1(JQGRAUP_LBC_TEND : JQGRAUP_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_BULK_LBC_TEND => D1(JCF_BULK_LBC_TEND : JCF_BULK_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_LIQUID_LBC_TEND => D1(JCF_LIQUID_LBC_TEND : JCF_LIQUID_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     CF_FROZEN_LBC_TEND => D1(JCF_FROZEN_LBC_TEND : JCF_FROZEN_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_extended,1)*wet_levels) )
     MURK_LBC_TEND => D1(JMURK_LBC_TEND : JMURK_LBC_TEND + &
      (LENRIMA(fld_type_p,halo_type_single,1)*model_levels) )

! Oxidant concentrations from UKCA for use in HadGEM sulphur
! cycle (these are in Section 33):
      IF(L_SULPC_ONLINE_OXIDANTS .AND. L_UKCA) THEN 
         OH_UKCA   => D1(JOH_UKCA(1)   : JOH_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
         H2O2_UKCA => D1(JH2O2_UKCA(1) : JH2O2_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
         HO2_UKCA  => D1(JHO2_UKCA(1)  : JHO2_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
         O3_UKCA   => D1(JO3_UKCA(1)   : JO3_UKCA(1)+ &
          field_length(theta_points,single_halo,model_levels))
      ELSE 
        OH_UKCA   => dummy_field
        H2O2_UKCA => dummy_field
        HO2_UKCA  => dummy_field
        O3_UKCA   => dummy_field
      END IF  
! 1.1 point tracer fields to D1

      ! find out how many tracers are active
      nActiveTracers=0
      DO nTracer=A_TRACER_FIRST,A_TRACER_LAST
        IF (SI(nTracer,33,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer
      IF (nActiveTracers /= 0) THEN
        ! set the pointer to the appropriate section of D1
        TRACER => D1(JTRACER(1,A_TRACER_FIRST) : &
         JTRACER(tr_levels,nActiveTracers)+theta_off_size)
        TRACER_LBC => D1(JTRACER_LBC(A_TRACER_FIRST) : &
         JTRACER_LBC(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels) )
        TRACER_LBC_TEND => D1(JTRACER_LBC_TEND(A_TRACER_FIRST) : &
         JTRACER_LBC_TEND(nActiveTracers) + &
         (LENRIMA(fld_type_p,halo_type_extended,1)*tr_levels) )
      ELSE
        ! or set it to something non-null if there are no active tracers
        TRACER => dummy_field
        TRACER_LBC => dummy_field
        TRACER_LBC_TEND => dummy_field
      END IF
     
      ! do the same for section 34 (UKCA) tracers     
      nActiveTracers=0
      DO nTracer=A_UKCA_FIRST,A_UKCA_LAST
        IF (SI(nTracer,34,atmos_im) /= 1) THEN
          nActiveTracers = nActiveTracers+1
        END IF
      END DO ! nTracer          
      IF (nActiveTracers /= 0) THEN
        TRACER_UKCA => D1(JTR_UKCA(1,A_UKCA_FIRST) : &
         JTR_UKCA(tr_levels,nActiveTracers)+theta_off_size)
      ELSE
        TRACER_UKCA => dummy_field
      END IF

END SUBROUTINE Set_Atm_Fields

#endif
