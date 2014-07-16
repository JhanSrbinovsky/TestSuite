! TYPSPD1 super array of D1 (NOTE: only single component array)
#if defined(T3E)
      ! Pad to ensure on a cache line boundary
      real :: spd1(spd1_len+8)
#else
      REAL :: SPD1(SPD1_LEN)
#endif
! TYPSPD1 end
