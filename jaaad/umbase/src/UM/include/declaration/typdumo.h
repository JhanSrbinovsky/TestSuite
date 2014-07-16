#if defined(OCEAN)
! TYPDUMO needs TYPSIZE included first
      ! Dump headers (ocean)
!L --------------- Dump headers (ocean) -----------------
#if defined(OCEAN)
      INTEGER :: O_FIXHD(LEN_FIXHD)     ! fixed length header
      INTEGER :: O_INTHD(O_LEN_INTHD)   ! integer header
      INTEGER :: O_CFI1(O_LEN_CFI1+1)   ! compress field index
      INTEGER :: O_CFI2(O_LEN_CFI2+1)   ! compress field index
      INTEGER :: O_CFI3(O_LEN_CFI3+1)   ! compress field index

      REAL::O_REALHD(O_LEN_REALHD)                    ! real header
      REAL::O_LEVDEPC(O_LEN1_LEVDEPC*O_LEN2_LEVDEPC+1)! level  dep const
      REAL::O_ROWDEPC(O_LEN1_ROWDEPC*O_LEN2_ROWDEPC+1)! row    dep const
      REAL::O_COLDEPC(O_LEN1_COLDEPC*O_LEN2_COLDEPC+1)! column dep const
      REAL::O_FLDDEPC(O_LEN1_FLDDEPC*O_LEN2_FLDDEPC+1)! field  dep const
      REAL::O_EXTCNST(O_LEN_EXTCNST+1)                ! extra constants
      REAL::O_DUMPHIST(LEN_DUMPHIST+1)                ! temp hist file

      ! PP headers
      INTEGER :: O_LOOKUP(LEN1_LOOKUP,O_LEN2_LOOKUP)  ! lookup heads
      INTEGER :: O_MPP_LOOKUP(MPP_LEN1_LOOKUP,O_LEN2_LOOKUP)
      INTEGER :: o_ixsts(len_o_ixsts)              ! stash index array

      INTEGER :: o_spsts(len_o_spsts)              ! ocean stash array
! TYPDUMO end
#endif
#endif
