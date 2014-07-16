! CANCFTNA start
!
! Purpose : Store Fortran Unit numbers for Atmosphere Ancillary Files
!
! Added units 188-207 - kdcorbin, 05/10
! -------------------------------------------------------------------
      ! Fortran Unit Numbers
      INTEGER :: FTN_ANCIL_A (NANCIL_DATASETSA) =                       &
     &  (/30,  31,  32,  33,  34,  35,  36,  38,  39,  96,              &
     &     1, 110,  48, 109, 111, 112, 115, 116, 117, 135,              &
     &   136, 137, 139, 118, 119, 120,  75,  76,  95,  77,              &
     &    78,  79,   1,   1,   1,   1,   1, 154, 155, 156,              &
     &   157, 158, 159, 160, 161, 148, 128, 188, 189, 190,              &
     &   191, 192, 193, 194, 195, 196, 197, 198, 199, 200,              &
     &   201, 202, 203, 204, 205, 206, 207/)

!     FTN_ANCIL_A (33-37) : Currently not used in UM
!     FTN_ANCIL_A (11)    : No longer used (slab model)

! CANCFTNA end
