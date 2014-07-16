! Start AODTYPE

! Description:
!   Sets aerosol type numbers for aerosol optical depth diags.
!
! Current Code Owner: Nicolas Bellouin
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   6.2    02/03/06 Original code. N. Bellouin
!
! A given aerosol type may gather several aerosol components
! (see AERCMP3A)
!
      INTEGER, PARAMETER :: IP_TYPE_SULPHATE = 1
      INTEGER, PARAMETER :: IP_TYPE_DUST     = 2
      INTEGER, PARAMETER :: IP_TYPE_SEASALT  = 3
      INTEGER, PARAMETER :: IP_TYPE_SOOT     = 4
      INTEGER, PARAMETER :: IP_TYPE_BIOMASS  = 5
      INTEGER, PARAMETER :: IP_TYPE_BIOGENIC = 6
      INTEGER, PARAMETER :: IP_TYPE_OCFF     = 7
      INTEGER, PARAMETER :: IP_TYPE_DELTA    = 8
! End AODTYPE
