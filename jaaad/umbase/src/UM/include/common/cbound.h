! CBOUND start

!  History:
!  Date      Vn     Modification
!  31/10/01  5.3    Remove RIMWEIGHTS_OROG. D. Robinson
!  17/09/02  5.4    Variables for controlling 2nd bndy file. A Clayton
!  17/02/03  5.5    Allow Wave model to use boundary code. D.Holmes-Bell
!  20/01/06  6.2    Add Current_LBC_Step. Dave Robinson
!  01/03/06  6.2    Remove RIMWEIGHTSW. Dave Robinson

      ! These 3 arrays are set by namelist read in IN_BOUND hence
      ! cannot be in argument list nor in COMMON if array lengths are
      ! passed variables. Only way seems to be to set MAX allowed
      ! lengths consistent with User Interface so that can be in COMMON
      INTEGER, PARAMETER :: MAX_BND_FLDS=4
      INTEGER, PARAMETER :: MAX_RIMWIDTH=10
      INTEGER :: BOUND_FIELDCODE(MAX_BND_FLDS)  ! Set by NAMELIST
      REAL :: RIMWEIGHTSA(MAX_RIMWIDTH)      ! Set by NAMELIST
      REAL :: RIMWEIGHTSO(MAX_RIMWIDTH)      ! Set by NAMELIST

      ! Variable for controlling LBC updating
      ! - initialised in INBOUNDA and updated in BOUNDVAL

      Integer  :: Current_LBC_Step     ! Timestep at which LBCs were
                                       ! last updated.

      ! Variables for controlling 2nd atmos boundary file. All are
      ! calculated within the code from other data.
      INTEGER :: ALBC_num              ! Number of atmos boundary file
                                       ! currently in use
      INTEGER :: ALBC2_StartTime_steps ! VT of first block of data in
                                       ! 2nd atmos boundary file, in
                                       ! steps from start of run
      INTEGER :: ALBC_SwapStep         ! Step on which to swap to 2nd
                                       ! atmos boundary file

      COMMON /BOUND_CT/ BOUND_FIELDCODE,                                &
     &                  RIMWEIGHTSA, RIMWEIGHTSO,                       &
     &                  ALBC_num, ALBC2_StartTime_steps, ALBC_SwapStep, &
     &                  Current_LBC_Step

! CBOUND end
