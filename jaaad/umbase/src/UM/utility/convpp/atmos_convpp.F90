#if defined(CONVPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  PROGRAM MAIN_CONVPP --------------------------------------------
!LL
!LL  Purpose: Converts a UM file into PP format.
!LL
!LL  Written by A. Dickinson 05/07/93
!LL
!LL  Model            Modification history:
!LL version  Date
!LL
!LL  3.3   31/10/93   Dimension of data array set to maximum value
!LL                   Author: A. Dickinson      Reviewer: P.Burton
!LL
!LL   3.3   15/12/93  Rename subroutine PRINTDUMP to CONVPP. D Robinson
!LL
!LL   3.3   08/12/93  Extra argument for READFLDS. D. Robinson
!LL
!LL   3.4   23/09/94  Extended to process ocean dumps. Alternative
!LL                   subroutine introduced
!LL                   Author D.M.Goddard
!LL   3.5  24/03/95    Changed OPEN to FILE_OPEN  P.Burton
!LL   4.4  24/10/96   Ocean data is written out without wrap points.
!LL                   Catherine Jones
!LL   4.4  23/04/97   Compressed fields are uncompressed using the
!LL                   subroutine UNPACK
!LL                   Catherine Jones


!LL   4.2  Oct. 96    DEF CRAY replaced by DEF T3E
!LL                             S.J.Swarbrick
!LL  4.4   Oct. 1997 Changed error handling from routine HDPPXRF
!LL                  so only fatal (+ve) errors are handled.
!LL                                             Shaun de Witt
!   4.4  23/04/97   Corrections to processing of Land-sea mask and
!                   Land compressed fields
!                   D.M. Goddard
!   4.4  24/10/97   Initialise ICODE as it is no longer
!                   initialised in HDPPXRF
!                   Author D.M. Goddard
!     4.5  14/10/97   Sets most significant number in packing indicator
!                     to zero (Native format) to enable PP package to be
!                     used on fieldsfile output
!                     Author D.M. Goddard
!     5.1  11/05/00   Added call to InitPrintStatus. D.P.Matthews
!LL  5.2  7/11/99  Enable run length encoding of ocean fieldsfiles to
!LL                compress the sequences of mdi values that represent
!LL                the land points. Ian Edmond
!    5.3  22/11/01  Enable MPP as the only option for
!                   small executables         E.Leung
!    5.4  03/09/02  Arguments are added to DECOMPOSE_SMEXE
!                                                  E.Leung
!    5.5  31/01/03  - Use INTHD 6,7 instead of LOOKUP 18,19 as GLSIZE
!                   - Enable convpp to be built on HP
!                                                             E.Leung
!    6.0  11/09/03  Add call to gc_exit                       P.Dando
!    6.2  10/04/05  Removed calls to ABORT.  J. Gill
!LL
!LL  Programming standards:
!LL
!LL  Logical components covered:
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  -----------------------------------------------------------------
!LL  SUBROUTINE ATMOS_CONVPP -----------------------------------------
!LL
!LL Purpose: Converts UM file to PP format.
!LL
!LL  Written by A. Dickinson
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL
!LL   3.4   23/09/94  New output lookup table array introduced because
!LL                   element 21 set to  0 for output file would need
!LL                   to be reset before attempting to read next record
!LL                   Routine renamed because separate routine for
!LL                   ocean dump introduced.
!LL                   Author D.M.Goddard
!     4.1  18/06/96  Changes to cope with changes in STASH addressing
!                    Author D.M. Goddard.
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
       SUBROUTINE ATMOS_CONVPP                                          &
     &  (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE)
!L
!L

      IMPLICIT NONE

      INTEGER                                                           &

     & LEN_FIXHD                                                        &
                    !IN Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                    !IN Length of integer header on input file
     &,LEN_REALHD                                                       &
                    !IN Length of real header on input file
     &,LEN1_LEVDEPC                                                     &
                    !IN 1st dim of lev dependent consts on input file
     &,LEN2_LEVDEPC                                                     &
                    !IN 2nd dim of lev dependent consts on input file
     &,LEN1_ROWDEPC                                                     &
                    !IN 1st dim of row dependent consts on input file
     &,LEN2_ROWDEPC                                                     &
                    !IN 2nd dim of row dependent consts on input file
     &,LEN1_COLDEPC                                                     &
                    !IN 1st dim of col dependent consts on input file
     &,LEN2_COLDEPC                                                     &
                    !IN 2nd dim of col dependent consts on input file
     &,LEN1_FLDDEPC                                                     &
                    !IN 1st dim of field dependent consts on input fi
     &,LEN2_FLDDEPC                                                     &
                    !IN 2nd dim of field dependent consts on input fi
     &,LEN_EXTCNST                                                      &
                    !IN Length of extra consts on input file
     &,LEN_DUMPHIST                                                     &
                    !IN Length of history header on input file
     &,LEN_CFI1                                                         &
                    !IN Length of index1 on input file
     &,LEN_CFI2                                                         &
                    !IN Length of index2 on input file
     &,LEN_CFI3                                                         &
                    !IN Length of index3 on input file
     &,LEN1_LOOKUP                                                      &
                    !IN 1st dim of LOOKUP on input file
     &,LEN2_LOOKUP                                                      &
                    !IN 2nd dim of LOOKUP on input file
     &,LEN_DATA                                                         &
                    !IN Length of data on input file
     &,P_FIELD                                                          &
                    !IN No of p-points per level on input file
     &,MAX_FIELD_SIZE !Maximum field size on file

      INTEGER                                                           &
     & NFTIN


! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & FIXHD(LEN_FIXHD),                                                &
                                                 !
     & INTHD(LEN_INTHD),                                                &
                                                 !\  integer
     & CFI1(LEN_CFI1+1),CFI2(LEN_CFI2+1),                               &
                                                 ! > file headers
     & CFI3(LEN_CFI3+1),                                                &
                                                 !/
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                  &
                                                 !
     &,LOOKUP_OUT(LEN1_LOOKUP) ! Output lookup table

! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows,                                                  &
                          ! IN  :number of P rows of entire mode
     &  tot_levels        ! IN  :total number of levels

      REAL                                                              &
     & REALHD(LEN_REALHD),                                              &
     & LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC),                            &
                                                 !
     & ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC),                            &
                                                 !
     & COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC),                            &
                                                 !\  real
     & FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC),                            &
                                                 ! > file headers
     & EXTCNST(LEN_EXTCNST+1),                                          &
                                                 !/
     & DUMPHIST(LEN_DUMPHIST+1),                                        &
                                                 !
     & D1(MAX_FIELD_SIZE)  ! Data array used to read in each field
      REAL     D1_TMP(MAX_FIELD_SIZE)

      LOGICAL  LAND_SEA_MASK(MAX_FIELD_SIZE)

! External subroutines called:------------------------------------------
      EXTERNAL ABORT,ABORT_IO,READHEAD,READFLDS,HDPPXRF,GETPPX,         &
     &         FROM_LAND_POINTS
!*----------------------------------------------------------------------
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "decomptp.h"
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------

      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,L      ! Loop indices
      INTEGER  ROWNUMBER     ! Row number
      INTEGER  IROWS         ! Number of points north-south
      INTEGER  ICOLS         ! Number of points east-west
      REAL     NROWS         ! Number of points north-south
      REAL     NCOLS         ! Number of points east-west
      REAL     RMDI          ! Real missing data indicator


      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,STRING*20    ! Format control for header printout
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)
!*----------------------------------------------------------------------

!L 0. Read in PPXREF

      ppxRecs=1
      RowNumber=0
      cmessage = ' '
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_A')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_O')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_S'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_S')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_W'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_W')

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_A',RowNumber,                  &
#include "argppx.h"
     &           ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_O',RowNumber,                  &
#include "argppx.h"
     &           ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_S',RowNumber,                  &
#include "argppx.h"
     &           ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,'             ',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,'             ',RowNumber,                     &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF

!L 1. Read in file header

! DEPENDS ON: readhead
      CALL READHEAD(NFTIN,FIXHD,LEN_FIXHD,                              &
     &                INTHD,LEN_INTHD,                                  &
     &                REALHD,LEN_REALHD,                                &
     &                LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                &
     &                ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                &
     &                COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                &
     &                FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                &
     &                EXTCNST,LEN_EXTCNST,                              &
     &                DUMPHIST,LEN_DUMPHIST,                            &
     &                CFI1,LEN_CFI1,                                    &
     &                CFI2,LEN_CFI2,                                    &
     &                CFI3,LEN_CFI3,                                    &
     &                LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                   &
     &                LEN_DATA,                                         &
#include "argppx.h"
     &                START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE
! DEPENDS ON: ereport
        CALL EREPORT('ATMOS_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF

      NROWS        = INTHD(7)
      NCOLS        = INTHD(6)
      IROWS        = INTHD(7)
      ICOLS        = INTHD(6)

      IF(LEN_REALHD >= 29)THEN
        RMDI         = REALHD(29)
      ENDIF

! Get decomposition information
! -----------------------------
      TOT_LEVELS=-99

! DEPENDS ON: decompose_smexe
      CALL DECOMPOSE_SMEXE(ICOLS,IROWS,0,0,TOT_LEVELS)
! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)

      DO I=1,LEN2_LOOKUP
        IF(LOOKUP(42,I) == 30)THEN

! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                   &
     &                  LAND_SEA_MASK,MAX_FIELD_SIZE,FIXHD,             &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVPP',CMESSAGE,ICODE,NFTIN)
        END IF
      END DO

!L  Print out individual fields
      DO I=1,LEN2_LOOKUP
        IF(LOOKUP(1,I) == -99)GOTO 100

! Fill output lookup table
        DO K=1,LEN1_LOOKUP
          LOOKUP_OUT(K)=LOOKUP(K,I)
        ENDDO

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                     &
     &                D1,MAX_FIELD_SIZE,FIXHD,                          &
#include "argppx.h"
     &                ICODE,CMESSAGE)

      LOOKUP_OUT(21)=MOD(LOOKUP_OUT(21),1000)
      IF((LOOKUP_OUT(21)/10)*10 == 120)THEN
!Data compressed on to land points.
!Copy data to temporary array
         DO K=1,LOOKUP_OUT(15)
           D1_TMP(K)=D1(K)
        END DO
!Set unpacked array to RMDI
         DO K=1,NROWS*NCOLS
           D1(K)=RMDI
         END DO

!Uncompress data
! DEPENDS ON: from_land_points
         CALL FROM_LAND_POINTS(D1,D1_TMP,LAND_SEA_MASK,                 &
     &                         MAX_FIELD_SIZE,LOOKUP_OUT(15))
         LOOKUP_OUT(15)=NROWS*NCOLS
         LOOKUP_OUT(18)=NROWS
         LOOKUP_OUT(19)=NCOLS
         LOOKUP_OUT(21)=0
       END IF
        WRITE(10)(LOOKUP_OUT(K),K=1,64)
        WRITE(10) (D1(K),K=1,LOOKUP_OUT(15))
      ENDDO

 100  CONTINUE
      WRITE(6,*)I-1,' pp fields written out'

      RETURN
      END SUBROUTINE ATMOS_CONVPP
!LL  SUBROUTINE OCEAN_CONVPP-----------------------------------------
!LL
!LL Purpose: Converts UM ocean file to PP format.
!LL
!LL  Written by D.M. Goddard
!LL
!LL  Model            Modification history from model version 3.4:
!LL version  Date
!LL
!LL   3.4   23/09/94  New routine at version 3.4
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------

#endif
