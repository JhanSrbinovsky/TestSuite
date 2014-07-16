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
      SUBROUTINE OCEAN_CONVPP                                           &
     &  (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE,JOC_NO_SEAPTS,LEN_OCFLD)
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
     &,MAX_FIELD_SIZE                                                   &
                      !IN Maximum field size on file
     &,JOC_NO_SEAPTS                                                    &
                      !IN Number of points in compressed array
     &,LEN_OCFLD      !IN Length of uncompressed ocean field

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
     & D1(MAX_FIELD_SIZE),                                              &
                           ! Array used to read in non-compressed fields
     & E1(MAX_FIELD_SIZE),                                              &
                           ! Array used to read in non-compressed fields
                           ! without wrap points
     & C1(JOC_NO_SEAPTS),                                               &
                           ! Array used to read in compressed fields
     & U1(LEN_OCFLD)       ! Array used to hold  uncompressed fields


! External subroutines called:------------------------------------------
       EXTERNAL ABORT,ABORT_IO,READHEAD,READFLDS,HDPPXRF,GETPPX,UNPACK
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
     &,I,J,K,L                                                          &
                     ! Loop indices
     &,IJ_IN,IJ_OUT                                                     &
                     ! More loop indices
     &,NROWS                                                            &
                     ! Number of points north-south
     &,IROWS                                                            &
     &,NROWS_FIELD                                                      &
                     ! Number of rows in a field
     &,NCOLS_IN                                                         &
                     ! Number of points east-west
     &,ICOLS_IN                                                         &
     &,NCOLS_OUT                                                        &
                     ! Number of points east-west for pp fields
     &,NLEVS                                                            &
                     ! Number of points in vertical
     &,NT                                                               &
                     ! Number of tracers
     &,NCOMP                                                            &
                     ! Number of compressed fields
     &,RECNUM                                                           &
                     ! Record number of field in lookup table
     &,POSIN                                                            &
                     ! Start position of field within C1
     &,POSU1                                                            &
                     ! Start position of field within U1
     &,FIELD_CODE                                                       &
                     ! field code for this field
     &,LBPACK        ! packing indicator from lookup table

      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,STRING*20    ! Format control for header printout

      REAL                                                              &
     & RMDI         ! Real missing data indicator

      LOGICAL                                                           &
     & LL_CYCLIC_IN    ! T => cyclic ; F => not cyclic

      INTEGER RowNumber

      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)

!*----------------------------------------------------------------------

!L 0. Read in PPXREF

      ppxRecs=1
      RowNumber=0
      CMESSAGE=''
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_A')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_O')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_S'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_S')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) 'Error reading STASHmaster_W'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_W')

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_A',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_O',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_S',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_W',RowNumber,                  &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,'             ',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,'             ',RowNumber,                     &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
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
        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)

      ENDIF


      NROWS        = INTHD(7)
      NCOLS_IN     = INTHD(6)
      IROWS        = INTHD(7)
      ICOLS_IN     = INTHD(6)
      LBPACK       = 21
      NLEVS        = INTHD(8)
      NT           = INTHD(14)

      IF(LEN_REALHD >= 29)THEN
        RMDI         = REALHD(29)
      ENDIF

! Get decomposition information
! -----------------------------
      TOT_LEVELS=-99

! DEPENDS ON: decompose_smexe
      CALL DECOMPOSE_SMEXE(ICOLS_IN,IROWS,0,0,TOT_LEVELS)
! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)

! Determine whether input data is cyclic and number of columns to output

      IF ( MOD ( FIXHD(4), 100 )  /=  3 ) THEN
        LL_CYCLIC_IN = .TRUE.
      ELSE
        LL_CYCLIC_IN = .FALSE.
      ENDIF

!L 2. Read in compressed data

      RECNUM=1

! Decide whether there are any compressed fields and on number of
! compressed fields. Use LBPACK to work out whether the first field
! contains sea points only.

      IF ( MOD(LOOKUP(LBPACK,1)/10,10)  ==  0) THEN

       NCOMP = 0

      ELSE

       NCOMP = NT + 2

      DO L=1,NCOMP

! Loop over levels storing all levels in one 1-D array
        POSIN=1
        DO K=1,NLEVS

! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,RECNUM,LOOKUP,LEN1_LOOKUP,C1(POSIN),    &
     &                  MAX_FIELD_SIZE,FIXHD,                           &
#include "argppx.h"
     &                  ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVPP',CMESSAGE,ICODE,NFTIN)
          POSIN=POSIN+LOOKUP(15,K+(L-1)*NLEVS)
          RECNUM=RECNUM+1

        ENDDO

!L 3. Uncompress 3-D field
! DEPENDS ON: unpack
        CALL UNPACK(1,NROWS,1,NLEVS,NROWS,NLEVS,NCOLS_IN,NROWS,NLEVS,   &
     &          CFI1,CFI2,LEN_CFI1,CFI3,JOC_NO_SEAPTS,                  &
     &          C1,U1,RMDI,LL_CYCLIC_IN)

!L 4. Output data level by level
        DO K=1,NLEVS

! Fill output lookup table
          DO I=1,LEN1_LOOKUP
            LOOKUP_OUT(I)=LOOKUP(I,K+(L-1)*NLEVS)
          ENDDO

          FIELD_CODE = LOOKUP_OUT(23)

          IF (  FIELD_CODE  >   600 .AND. FIELD_CODE  <   700) THEN
            NROWS_FIELD = NROWS
          ELSE IF ( FIELD_CODE  >   699 .AND. FIELD_CODE  <   800 ) THEN
            NROWS_FIELD = NROWS - 1
          ELSE
            write(6,*) ' unknown field code : exiting '
            go to 9999
          END IF

! Determine number of columns to output

          IF ( LL_CYCLIC_IN ) THEN
            NCOLS_OUT = NCOLS_IN - 2
          ELSE
            NCOLS_OUT = NCOLS_IN
          ENDIF

          LOOKUP_OUT(15)=NROWS_FIELD*NCOLS_OUT
          LOOKUP_OUT(18)=NROWS_FIELD
          LOOKUP_OUT(19)=NCOLS_OUT
          LOOKUP_OUT(21)=0
          WRITE(10)(LOOKUP_OUT(I),I=1,64)

          POSU1=(K-1)*NROWS*NCOLS_IN
          DO J=1,NROWS_FIELD
            DO I=1,NCOLS_OUT
              IJ_IN = I + (J-1) * NCOLS_IN
              IJ_OUT   = I + (J-1) * NCOLS_OUT
              E1(IJ_OUT) = U1(IJ_IN+POSU1)
        ENDDO
          ENDDO


          WRITE(10) (E1(I),I=1,NROWS_FIELD*NCOLS_OUT)

      ENDDO

      ENDDO

      END IF  ! LBPACK


!L 5.  Now processing non compressed fields
!  Print out individual fields
      DO L=RECNUM,LEN2_LOOKUP

        IF(LOOKUP(1,L) == -99)GOTO 100

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,L,LOOKUP,LEN1_LOOKUP,                     &
     &                D1,MAX_FIELD_SIZE,FIXHD,                          &
#include "argppx.h"
     &                ICODE,CMESSAGE)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('CONVPP',CMESSAGE,ICODE,NFTIN)
        IF(FIXHD(5) /= 3)THEN

! Take off the extra columns if the dump is cyclic using the E1 array

          IF ( LL_CYCLIC_IN ) THEN
            NCOLS_OUT = NCOLS_IN - 2
          ELSE
            NCOLS_OUT = NCOLS_IN
          ENDIF

          DO J=1,NROWS
            DO I=1,NCOLS_OUT
              IJ_IN = I + (J-1) * NCOLS_IN
              IJ_OUT   = I + (J-1) * NCOLS_OUT
              E1(IJ_OUT) = D1(IJ_IN)
            ENDDO
          ENDDO

        ELSE

! Fieldsfile. NO cyclic columns

          DO I=1,LOOKUP(15,L)
            E1(I) = D1(I)
          ENDDO

        ENDIF
        DO K=1,LEN1_LOOKUP
          LOOKUP_OUT(K)=LOOKUP(K,L)
        ENDDO

        IF(FIXHD(5) /= 3)THEN
          LOOKUP_OUT(15)=NROWS*NCOLS_OUT
          LOOKUP_OUT(19)=NCOLS_OUT
        ENDIF
        LOOKUP_OUT(21)=MOD(LOOKUP_OUT(21),1000)
        WRITE(10)(LOOKUP_OUT(K),K=1,64)
        WRITE(10) (E1(K),K=1,LOOKUP_OUT(15))

      ENDDO

 100  CONTINUE
      WRITE(6,*)L-1,' pp fields written out'

9999  continue

      RETURN
      END SUBROUTINE OCEAN_CONVPP

#endif
