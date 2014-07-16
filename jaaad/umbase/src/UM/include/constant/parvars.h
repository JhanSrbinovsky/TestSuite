! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the MPP-UM
! Added new comdeck AMAXSIZE required for new arrays in PARCOMM
! Add non-MPP option
!                                                      P.Burton
#include "parparm.h"
#include "amaxsize.h"
#include "parcomm.h"
#if defined(UTILIO) || defined(FLDIO) || defined(UTILHIST) \
 || defined(FLUXPROC)
      INTEGER mype,nproc
      PARAMETER (mype=0,nproc=1)
#endif
