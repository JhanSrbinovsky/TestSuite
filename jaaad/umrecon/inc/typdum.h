! TYPDUM start. TYPSIZE must be included first

      ! Dump headers
      INTEGER :: FIXHD(LEN_FIXHD)   ! fixed length header
      INTEGER :: INTHD(LEN_INTHD)   ! integer header
      INTEGER :: CFI1(LEN_CFI1+1)   ! compress field index
      INTEGER :: CFI2(LEN_CFI2+1)   ! compress field index
      INTEGER :: CFI3(LEN_CFI3+1)   ! compress field index
      REAL :: REALHD(LEN_REALHD)                   ! real header
      REAL :: LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC+1) ! level  dep const
      REAL :: ROWDEPC(LEN1_ROWDEPC*LEN2_ROWDEPC+1) ! row    dep const
      REAL :: COLDEPC(LEN1_COLDEPC*LEN2_COLDEPC+1) ! column dep const
      REAL :: FLDDEPC(LEN1_FLDDEPC*LEN2_FLDDEPC+1) ! field  dep const
      REAL :: EXTCNST(LEN_EXTCNST+1)               ! extra constants
      REAL :: DUMPHIST(LEN_DUMPHIST+1)             ! temporary hist file

      ! PP headers
      INTEGER :: LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! lookup heads
      INTEGER :: MPP_LOOKUP(MPP_LEN1_LOOKUP,LEN2_LOOKUP)
      INTEGER :: ixsts(len_ixsts)                ! stash index array
      INTEGER :: spsts(len_spsts)                      ! stash array
! TYPDUM end
