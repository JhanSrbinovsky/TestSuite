!LL  Comdeck: CMEANCTL -------------------------------------------------
!LL
!LL  Purpose: Contains control information local to the climate meaning
!LL           control subroutines INITMEAN, MEANCTL etc.
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL  Logical components covered: C500
!LL
!LL  Project task: C5
!LL
!LLEND -------------------------------------------------------------
      ! Length of largest IO buffer needed (A/O)
      INTEGER :: IBUFLEN(2)
!
      COMMON /CMEANCTL/ IBUFLEN
!
