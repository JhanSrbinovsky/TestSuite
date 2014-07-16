#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine IN_BOUND
!LL
!LL Purpose : Takes as input,the code defining whether updates of
!LL  boundary data are required.  The physical files required are
!LL  identified, and the headers lookup tables are read into common
!LL  blocks.  Reads the update intervals from the boundary datasets.
!LL  Where the update interval is in months or years, the check will be
!LL  made daily.
!LL
!LL Control routine for CRAY YMP
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version no. 1, dated 15/01/90
!LL
!LL Logical components covered : C720
!LL
!LL System task : C7
!LL
!LL Documentation : Unified Model Documentation Paper No C7
!LLEND
!
!*L  Arguments
      SUBROUTINE IN_BOUND (                                             &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "argbnd.h"
     &   A_LEN1_LEVDEPCDA,A_LEN2_LEVDEPCDA,                             &
     &   A_LEN1_ROWDEPCDA,A_LEN2_ROWDEPCDA,                             &
     &   A_LEN1_COLDEPCDA,A_LEN2_COLDEPCDA,                             &         
#include "argppx.h"
     &           ICODE,CMESSAGE)

      IMPLICIT NONE

#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "typduma.h"
#include "typsts.h"
#include "typptra.h"
#include "typbnd.h"

      INTEGER                                                           &
     &         A_LEN1_LEVDEPCDA,                                        &
                                   ! IN : copy of A_LEN1_LEVDEPC
     &         A_LEN2_LEVDEPCDA,                                        &
                                   ! IN : copy of A_LEN2_LEVDEPC
     &         A_LEN1_ROWDEPCDA,                                        &
                                   ! IN : copy of A_LEN1_ROWDEPC
     &         A_LEN2_ROWDEPCDA,                                        &
                                   ! IN : copy of A_LEN2_ROWDEPC
     &         A_LEN1_COLDEPCDA,                                        &
                                   ! IN : copy of A_LEN1_COLDEPC
     &         A_LEN2_COLDEPCDA,                                        &
                                   ! IN : copy of A_LEN2_COLDEPC                                   
     &        ICODE            ! Return code = 0 Normal Exit
!                              !    "     "  > 0 Error Exit

      CHARACTER*(80) CMESSAGE  ! Error message if ICODE > 0
!*

#include "ppxlook.h"

!L      Internal Structure

#if (defined(ATMOS) && !defined(GLOBAL))

! Call inbound for the atmosphere

! DEPENDS ON: inbounda
      CALL INBOUNDA(                                                    &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "argbnd.h"
#include "argppx.h"
     &   A_LEN1_LEVDEPCDA,A_LEN2_LEVDEPCDA,                             &
     &   A_LEN1_ROWDEPCDA,A_LEN2_ROWDEPCDA,                             &
     &   A_LEN1_COLDEPCDA,A_LEN2_COLDEPCDA)
#endif

!L  4   End of routine

      RETURN
      END SUBROUTINE IN_BOUND
!-----------------------------------------------------------------------

#endif
