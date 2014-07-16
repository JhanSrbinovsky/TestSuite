! Start c_st_bdy
!
! Contains resistance factors for dry deposition of 3 modes of soot
!
! Current Code Owner: Paul Davison
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.4    10/07/02 Original code, based on parameters used in the soot
!                   scheme at v4.5.                 Paul Davison.
!   5.5    28/02/03 Revise deposition parameters to be consistent
!                   with new size distributions.   P. Davison
!
! Rb(soot modes)/Rb(H2O)
       REAL, PARAMETER :: RESB_FreshSoot = 1840.77
       REAL, PARAMETER :: RESB_AgedSoot = 1840.77
       REAL, PARAMETER :: RESB_SootInCloud = 0.0

!       Stomatal resistance set to zero for particles
       REAL, PARAMETER :: RESS_Soot = 0.0

! End c_st_bdy
