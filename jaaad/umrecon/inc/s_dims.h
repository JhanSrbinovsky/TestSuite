!include file: s_dims.h
! Description:
!   This deck defines the dimensions used in the SCM.
!
! Current Code Owner: Ricky Wong
!
! Declarations:

      INTEGER                                                           &
     &  model_levels                                                    &
     &, wet_levels                                                      &
     &, bl_levels                                                       &
     &, cloud_levels                                                    &
     &, tr_levels                                                       &
     &, tr_vars, tr_ukca                                                &
     &, st_levels                                                       &
     &, sm_levels                                                       &
     &, ozone_levels                                                    &
     &, ntiles

      ! Max. no. of STASH sections  per internal model (44 in practice)
      ! Note: not used for STASH purposes in the SCM.
      ! Copied from version.h
      INTEGER,PARAMETER :: NSECTP=99
      ! String array to hold version choices for various physics
      ! schemes. Only radiation choices (indices 1 and 2) are used in
      ! SCM.  Copied from model.h
      CHARACTER(LEN=2) :: ATMOS_SR(0:NSECTP)

!- End of COMDECK declaration

      NAMELIST/NLSIZES/ CLOUD_LEVELS, MODEL_LEVELS, WET_LEVELS,         &
     &  TR_LEVELS, TR_VARS, TR_UKCA, ST_LEVELS, SM_LEVELS, BL_LEVELS,   &
     &  OZONE_LEVELS, NTILES, ATMOS_SR
      ! Note that it is a hack to include ATMOS_SR from the STSHCOMP
      ! UMUI namelist here, but it enables the radiation options to
      ! work for the SCM.

      COMMON/NLSIZES/ CLOUD_LEVELS, MODEL_LEVELS, WET_LEVELS,           &
     &  TR_LEVELS, TR_VARS, TR_UKCA, ST_LEVELS, SM_LEVELS, BL_LEVELS,   &
     &  OZONE_LEVELS, NTILES, ATMOS_SR

!---------------------------------------------------------------------
