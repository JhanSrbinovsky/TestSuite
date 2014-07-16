#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: UM_INDEX-----------------------------------------------
!LL
!LL  Purpose: Calculate addresses of component arrays within a
!LL           series of super arrays, made up of combinations of arrays
!LL           that require dynamic allocation. Lengths of the super
!LL           arrays are calculated and passed into U_MODEL for
!LL           dynamic allocation, reducing the no. of arguments needed
!LL           to be passed between top-level routines.
!LL
!LL  Tested under compiler:   cft77
!LL  Tested under OS version: UNICOS 6.1.5A
!LL
!LL  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!LL
!LL  Logical components covered: C0
!LL
!LL  Project task: C0
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
      SUBROUTINE UM_INDEX(                                              &
#include "argszsp.h"
#include "argszspa.h"
#include "argszspc.h"
     &              ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------
      IMPLICIT NONE
!
!  Local variables
!
      INTEGER ICODE             ! Work - Internal return code
      CHARACTER*80  CMESSAGE    ! Work - Internal error message
!
!  Configuration-dependent sizes for dynamic arrays
!
! Provides N_INTERNAL_MODEL parameter for STASH array
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
! For STASH sizes
#include "typstsz.h"
!
!  ppxref parameters needed for a STASH array
#include "cppxref.h"
!
!  Super array sizes for dynamic allocation in U_MODEL
!
#include "typszsp.h"
#include "typszspa.h"
#include "typszspc.h"
#include "chsunits.h"
#include "ccontrol.h"
!
! Holds D1_LIST_LEN needed to calculate array size
#include "d1_addr.h"
!

#include "decomptp.h"
#include "decompdb.h"
!
!  Addresses of arrays in super arrays.
!
#include "spindex.h"
! to use PrintStatus:
#include "cprintst.h"
!
      INTEGER LEN_D1_ADDR
!L----------------------------------------------------------------------
!L 0. Start Timer running
!L
      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('UM_SETUP',4)
! DEPENDS ON: timer
        CALL TIMER('UM_INDEX',3)
      END IF
      ICODE=0

! 0.1 Initialise variables in TYPSIZE
      LEN_A_IXSTS = A_IXSTS_LEN
      LEN_A_SPSTS = 0
!L----------------------------------------------------------------------
!L 1. Calculate     addresses in super array and each super array length
!L
!L 1.1   D1       super array
!L
!L          super array addresses
      LEN_D1_ADDR=D1_LIST_LEN*N_OBJ_D1_MAX*N_SUBMODEL_PARTITION
#if defined(T3E)
!--make sure the increment aligns on to an Scache boundary,
!  if the initial location of D1 is set properly
      len_d1_addr=((len_d1_addr+7)/8)*8
#endif
      IXD1(1)=1
      IXD1(2)=LEN_D1_ADDR+1
      IXD1(3)=LEN_D1_ADDR+1         ! array not used
      IXD1(4)=LEN_D1_ADDR+1         ! array not used
!L
!L          super array length
      SPD1_LEN=LEN_TOT+LEN_D1_ADDR
!L
!L
!  1.2   STASH super array

!        super array addresses
!     SF
      IXSTS( 1)=1
!     STINDEX
      IXSTS( 2)=IXSTS( 1)+ (NITEMS+1)*(NSECTS+1)
!     STLIST
      IXSTS( 3)=IXSTS( 2)+2*NITEMS   *(NSECTS+1)                        &
     &                     *N_INTERNAL_MODEL
!     SI
      IXSTS( 4)=IXSTS( 3)+  LEN_STLIST*TOTITEMS
!     STTABL
      IXSTS( 5)=IXSTS( 4)+  NITEMS   *(NSECTS+1)                        &
     &                     *N_INTERNAL_MODEL
!     STASH_MAXLEN
      IXSTS( 6)=IXSTS( 5)+  NSTTIMS  *NSTTABL
!     PPINDEX
      IXSTS( 7)=IXSTS( 6)+ (NSECTS+1)*N_INTERNAL_MODEL
!     STASH_LEVELS
      IXSTS( 8)=IXSTS( 7)+  NITEMS   *N_INTERNAL_MODEL
!     STASH_PSEUDO_LEVELS
      IXSTS( 9)=IXSTS( 8)+ (NUM_STASH_LEVELS+1)*NUM_LEVEL_LISTS
!     STASH_SERIES
      IXSTS(10)=IXSTS( 9)+ (NUM_STASH_PSEUDO+1)*NUM_PSEUDO_LISTS
!     STASH_SERIES_INDEX
      IXSTS(11)=IXSTS(10)+  TIME_SERIES_REC_LEN*NSTASH_SERIES_RECORDS
!     MOS_MASK
      IXSTS(12)=IXSTS(11)+2*NSTASH_SERIES_BLOCK

!L
!L          super array length
      SPSTS_LEN    =IXSTS(12)+ MOS_MASK_LEN
      SPSTS_LEN    =SPSTS_LEN -1
!L
!L
!L 1.3   Input boundary conditions   super array
!L
!L          super array addresses
      IXBND(1) =1
!L
!L          super array length
      SPBND_LEN    =IXBND(1)  ! References to this superarray are now
                              ! redundant (at 5.0) and can be removed.
!L
!L----------------------------------------------------------------------
!L 2.    atmosphere super arrays
#if defined(ATMOS)
!L
! DEPENDS ON: um_index_a
      CALL UM_INDEX_A(                                                  &
#include "argszspa.h"
     &              ICODE,CMESSAGE)
#endif
!L
!L----------------------------------------------------------------------
!L 5.    coupled    super arrays  for river routing only
!L
!L Get global sizes (all PEs) for river routing regridding routines
      AOCPL_ROW_LENGTH                                                  &
     &            =decomp_db_glsize(1,fld_type_p,decomp_standard_atmos)
      AOCPL_P_ROWS=decomp_db_glsize(2,fld_type_p,decomp_standard_atmos)
!L
!L          super array addresses
      AO_IXCPL(1) =1
      AO_IXCPL(2) =AO_IXCPL(1) + 0
      AO_IXCPL(3) =AO_IXCPL(2) + 0
      AO_IXCPL(4) =AO_IXCPL(3) + 0
      AO_IXCPL(5) =AO_IXCPL(4) + 0
      AO_IXCPL(6) =AO_IXCPL(5) + AOCPL_ROW_LENGTH+1
      AO_IXCPL(7) =AO_IXCPL(6) + AOCPL_ROW_LENGTH+1
      AO_IXCPL(8) =AO_IXCPL(7) + AOCPL_ROW_LENGTH+1
      AO_IXCPL(9) =AO_IXCPL(8) + AOCPL_P_ROWS
      AO_IXCPL(10)=AO_IXCPL(9) + AOCPL_P_ROWS
!L
!L          super array length
      AO_SPCPL_LEN  =AO_IXCPL(10)+ AOCPL_P_ROWS+1
      AO_SPCPL_LEN  =AO_SPCPL_LEN -1
!L----------------------------------------------------------------------
!L 6. Exit processing
!L

      IF(PrintStatus  ==  PrStatus_Diag) THEN
      write(6,*) 'UM_INDEX: superarray lengths:'

      write(6,*) 'Super array lengths:',                                &
#include "argszsp.h"
     &'   not sub-model'
      write(6,*) 'Super array lengths:',                                &
#include "argszspa.h"
     &'   atmosphere'
      write(6,*) 'Super array lengths:',                                &
#include "argszspc.h"
     &'   coupled for river routing'

      ENDIF ! PrintStatus

      IF(LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('UM_INDEX',4)
      END IF
!
      RETURN
      END SUBROUTINE UM_INDEX
#endif
