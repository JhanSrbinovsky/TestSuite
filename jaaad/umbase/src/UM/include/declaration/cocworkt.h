#if defined(OCEAN)
! Reflects structure of COCTWRKA except master slab arrays TW and TBW
! are used rather than the individual TM,T,TP,TPP etc.
!     COMDECK COCWORKT
!     ----------------
      REAL                                                              &
     & TA(NSLAB),TF(NSLAB)                                              &
     &,TW(NSLAB,4)                                                      &
     &,TBW(NSLAB,4)
#include "octwkgen.h"
#include "typocloc.h"
#include "typocinc.h"
#include "coctmom.h"
#endif
