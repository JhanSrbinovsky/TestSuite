#if defined(ATMOS)
#if defined(CONTROL) || defined(SCMA)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Global data module for switches concerned with radiation.

MODULE rad_switches_mod

  IMPLICIT NONE
  SAVE

! Description:
!   Module containing logical switches used by the radiation code.
!
! Method:
!   Switches are initialised to false and then set in atm_step to
!   the corresponding switch used in cntlatm.h and read in from the 
!   UMUI. The module may then be used directly where the switches
!   are needed within the radiation code.
!
! Current Code Owner: James Manners
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

  LOGICAL :: lrad_ctile_fix = .false.
!             Needed in r2_set_surface_field_lw-lwrad3a.
!             Fix the selection of fractional sea points so that
!             machine tolerance is not used. This then agrees with
!             the criteria used later in diff_planck_source-dfpln3a.

  LOGICAL :: lrad_cldtop_t_fix = .false.
!             Corrects convective cloud top temperature. Existing code
!             (=<Vn6.5) incorrectly selects temperature at the top of
!             atmosphere to partition cloud liquid/ice.


  LOGICAL :: lrad_quad_src_fix = .false.
!             Needed in trans_source_coeff-trsfc3b.
!             Use asymptotic form of quadratic correction to the LW 
!             source term when tau is less than the quad root of 
!             EPSILON. (Fix's instability in LW source function near 
!             the top of the atmosphere.)


  LOGICAL :: lrad_ccrad = .false.
!             Allows access to ccrad code and the logicals
!             lrad_3d_ccw, lrad_ovrlap and lrad_ccw_scav 


  LOGICAL :: lrad_3d_ccw = .false.
!             (Only active if l_ccrad from convection is .TRUE.)
!             Used in set_cloud_field. If TRUE, the radiation code 
!             will use the convective cloud water (CCW) profile passed
!             to it by convection. Otherwise it will construct a CCW 
!             profile using cclwp, ccb and cct.


  LOGICAL :: lrad_ovrlap = .false.
!             (Only active if l_ccrad from convection is .TRUE.)
!             (THIS IS EXPERIMENTAL AND USED FOR DEVELOPMENT ONLY).
!             Current convective/large-scale cloud fractions in the
!             radiation scheme are mutally exclusive. This assumes CCA
!             and the large-scale to overlap with CCA taking dominance.
!             I.E. Large-scale cloud fraction must exceed the convective
!             cloud fraction before having any presence.


  LOGICAL :: lrad_ccw_scav = .false.
!             (Only active if l_ccrad from convection is .TRUE.
!                       .AND. lrad_ovrlap is  .TRUE.)
!             (THIS IS EXPERIMENTAL AND USED FOR DEVELOPMENT ONLY)
!             Allowing the CCA to negate large-scale fractions of lower
!             values means that the large-scale cloud water in the
!             overlapping fraction is lost. This switch will scavenge 
!             the large-scale cloud water from the overlapping fraction
!             and combine it with the convective cloud water to
!             conpensate.

  LOGICAL :: LRAD_TRIPLECLOUDS = .FALSE.
!             If this switch is set to TRUE, the Tripleclouds method of
!             introducing horizontal cloud inhomogeneity is invoked.

  LOGICAL :: LRAD_EXPRAND = .TRUE.
!             If this switch is set to TRUE, exponential random overlap
!             is invoked.
 
  LOGICAL  :: LRAD_EMIS_LAND_GEN = .false.
!             Needed to preserve bit reproducibility when introducing the 
!             functionality to allow non-unity values for aggregate land 
!             surface emissivity. 

  REAL     :: RAD_EMIS_LAND_GEN = 1.0
!             Namelist value for aggregate land surface emissivity. 

!- End of header

END MODULE rad_switches_mod
#endif
#endif

