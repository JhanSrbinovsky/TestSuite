! Start c_bm_bdy
!
! Contains resistance factors for dry deposition of biomass smoke
!
! Current Code Owner: Paul Davison
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.5    05/02/03 Original code                   Paul Davison.
!   6.2    01/11/05 Add modified values for version 2B of the aerosol
!                   scheme.        Andy Jones & Nicolas Bellouin.
!
! Rb(smoke modes)/Rb(H2O)
#if !defined(A17_2B)
       REAL, PARAMETER :: RESB_FreshBmass = 2261.0
       REAL, PARAMETER :: RESB_AgedBmass = 5598.0
#endif
#if defined(A17_2B)
       REAL, PARAMETER :: RESB_FreshBmass = 2492.4
       REAL, PARAMETER :: RESB_AgedBmass = 2970.7
#endif
       REAL, PARAMETER :: RESB_BmassInCloud = 0.0

!       Stomatal resistance set to zero for particles
       REAL, PARAMETER :: RESS_Bmass = 0.0

! End c_bm_bdy
