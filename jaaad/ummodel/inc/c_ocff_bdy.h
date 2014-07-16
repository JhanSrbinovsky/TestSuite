! Start c_ocff_bdy
!
! Contains resistance factors for dry deposition of fossil fuel
! organic carbon aerosols
!
! Current Code Owner: M S Reddy
!
! Rb(smoke modes)/Rb(H2O)
       REAL, PARAMETER :: RESB_Freshocff = 2492.0
       REAL, PARAMETER :: RESB_Agedocff  = 2970.0
       REAL, PARAMETER :: RESB_ocffInCloud = 0.0

!       Stomatal resistance set to zero for particles
       REAL, PARAMETER :: RESS_ocff = 0.0

! End c_ocff_bdy
