#if defined(OCEAN)
! TYPOSTFL start
      INTEGER :: DT_SIZE  !}
      INTEGER :: DS_SIZE  !} Dimensions for STASH work space for DT, DS,
      INTEGER :: BIO_SIZE !} BIOLOGY and Air-Sea Flux diagnostics.
      INTEGER :: ASF_SIZE !}
      INTEGER :: SI_DT(*)  !}
      INTEGER :: SI_DS(*)  !} Index pointers to individual items in DT,
      INTEGER :: SI_BIO(*) !} DS, BIOLOGY and AS Flux STASH work arrays.
      INTEGER :: SI_ASF(*) !}

      REAL :: swrk_dt(*)   ! STASH work space - DT diagnostics
      REAL :: swrk_ds(*)   ! STASH work space - DS diagnostics
      REAL :: swrk_bio(*)  ! STASH work space - BIOLOGY diagnostics
      REAL :: swrk_asf(*)  ! STASH work space - Air-Sea Flux diagnostics

      LOGICAL :: SF_DT(*)  !}
      LOGICAL :: SF_DS(*)  !} STASH flags for DT, DS, BIOLOGY
      LOGICAL :: SF_BIO(*) !} and Air-Sea Flux diagnostics
      LOGICAL :: SF_ASF(*) !}
! TYPOSTFL end
#endif
