#if defined(OCEAN)
! TYPCONO start
      ! Array containing pre-calculated ocean fields - carried down to
      ! ocean routines as a single array and decomposed to constituent
      ! routines only within the lower ocean routines.This is a special
      ! case in which the array size is also passed down as an argument.
      INTEGER :: O_SPCON_LEN
      REAL    :: O_SPCON(O_SPCON_LEN)
! TYPCONO end
#endif
