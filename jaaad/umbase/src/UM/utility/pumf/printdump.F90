#if defined(PUMF)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PRINTDUMP---------------------------------------------
!  
!   Purpose: Prints a summary of contents of atmosphere, ocean or
!            ancillary file.
!            PRINTDUMP reads in headers and data fields from unit NFTIN
!            printing out a summary of their contents.
!            Printout of headers is written to unit 7.
!            Printout of data fileds is written to unit 6.
!  
!    Written by A. Dickinson
!  
!    Documentation: UM Doc Paper F5
!  
!    System Tasks: F3,F4,F6
!  
!    -----------------------------------------------------------------
!    Arguments:-------------------------------------------------------
      SUBROUTINE PRINTDUMP                                              &
     &  (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  P_ROWS,ROW_LENGTH,                                              &
     &  nftin,max_field_size,wgdos_expand)
! 
! 

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
     &,P_ROWS                                                           &
     &,ROW_LENGTH                                                       &
     &,MAX_FIELD_SIZE                                                   &
                      !IN Maximum field size on file
     &,wgdos_expand ! IN set to 1 to exapnd WGDOS Fields for comparison

      integer lblrec_1

      INTEGER                                                           &
     &  NFTIN,                                                          &
     &  NFTOUT

      INTEGER EXPPXI
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "cstash.h"
#include "clookadd.h"
#include "c_mdi.h"
#include "parvars.h"
#include "decomptp.h"
#include "sx_size.h"
#include "atm_lsm.h"

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
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)           !

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
     & RLOOKUP(LEN1_LOOKUP,LEN2_LOOKUP),                                &
     & D1(MAX_FIELD_SIZE)  ! Data array used to read in each field

      LOGICAL LAND_MASK(MAX_FIELD_SIZE)  ! land-sea mask
      INTEGER RowNumber
      CHARACTER*100 PAREXE_ENV  ! hold name of the // exec script
      INTEGER ME_GC,NPROC_GC

!*----------------------------------------------------------------------
!    Local variables:---------------------------------------------------

      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,MASK_ROW                                                         &
                    ! row length of land sea mask
     &,MASK_COL                                                         &
                    ! no. of col of land sea mask
     &,FIELD_ITEM                                                       &
     &,FIELD_SECT                                                       &
     &,FIELD_MODEL                                                      &
     &,GRID_TYPE                                                        &
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,L      ! Loop indices

      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,STRING*20                                                        &
                    ! Format control for header printout
     &,FILENAME*80  !Name of user preSTASH master file
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)

      character(len=80), parameter :: routinename='PrintDump'

      INTEGER POS_MAX,POS_MIN
      REAL    F_MAX,F_MIN
      REAL RCODE   ! return code from READ_LAND_SEA
! Variables for printing out observation fields
      INTEGER IIII             !loop counter
      INTEGER BLK, NBLK, PBEGIN, PEND, OBS_LEFT
      INTEGER NumMeta, NumLevs, Lev, NumObs, NumObs_print, N
      INTEGER NumXdata_print
      INTEGER NumItem,NumObVariables, Variable, Shift, IPGE

      character*5 CNPRINT,CXPRINT

      real PGE(8),PGE1(8),PGE2(8)

! Dummy variables used for call to F_TYPE
      INTEGER                                                           &
     & N_TYPES                                                          &
     &,PP_NUM(LEN2_LOOKUP)                                              &
                             ! No of successive fields with same co
     &,PP_LEN(LEN2_LOOKUP)                                              &
                             ! Length of field
     &,PP_ITEMC(LEN2_LOOKUP)                                            &
                             ! PP code of field
     &,PP_TYPE(LEN2_LOOKUP)                                             &
                             ! Integer/real/timeseries
     &,PP_POS(LEN2_LOOKUP)                                              &
                             ! Pointer to number of PP field
     &,PP_LS(LEN2_LOOKUP)    ! Data stored on land or sea pts

! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  HALOX,HALOY       ! EW and NS halo

      CHARACTER*80 F_TYPE_TITLE
      logical :: summary
      character(len=4) :: csummary

!*----------------------------------------------------------------------
! Get the environment variable SUMMARY to decide whether summary or
! full output is required
      summary = .false.
      call fort_get_env("SUMMARY",7,csummary,4,icode)
      if(icode  >   0)then
        cmessage='Problem picking up value of "summary"                 &
     & from wrapper script'
! DEPENDS ON: ereport
        call ereport(RoutineName,icode,cmessage)
      endif
      if(csummary == 'true')then
        summary = .true.
      endif


      NFTOUT=7                  ! Out file on unit 7.
      cmessage = ' '

!  0. Read in PPXREF

      ppxRecs=1
      RowNumber=0
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   'Error reading STASHmaster_A')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   'Error reading STASHmaster_O')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_S',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_S'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   'Error reading STASHmaster_S')

      END IF
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_W',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_W'
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('PRINTDUMP', ICODE,                                &
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
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   CMESSAGE)

      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,'             ',ppxRecs,ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   CMESSAGE)

      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,'             ',RowNumber,                     &
#include "argppx.h"
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE
! DEPENDS ON: ereport
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   CMESSAGE)

      ENDIF

!  1. Read in file header

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
        CALL EREPORT('PRINTDUMP', ICODE,                                &
     &   CMESSAGE)

      ENDIF

!     Write contents summary for model dumps or for fields files
!     and lbcs only when requested
      IF ((FIXHD(5) == 1) .or.                                          &
     &    (fixhd(5) == 3 .and. summary) .or.                            &
     &    (fixhd(5) == 5 .and. summary)) THEN
        WRITE (6,*) '=================================================='
        select case(fixhd(5))
          case(3)
            F_TYPE_TITLE='Contents of Fields File'
          case(5)
            F_TYPE_TITLE='Contents of LBC file'
          case default
            F_TYPE_TITLE='Contents of dump or other file'
        end select
! DEPENDS ON: f_type
        CALL F_TYPE(LOOKUP,LEN2_LOOKUP,PP_NUM,N_TYPES,                  &
     &              PP_LEN,PP_ITEMC,PP_TYPE,PP_POS,                     &
     &              PP_LS,FIXHD,                                        &
#include "argppx.h"
     &              F_TYPE_TITLE)
        WRITE (6,*)
        WRITE (6,*) '=================================================='
      ENDIF

!       Open up unit NFTOUT (N.Farnon)

      CALL GET_FILE(NFTOUT,FILENAME,80,ICODE)
      OPEN(NFTOUT,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
      ELSE
        WRITE(6,*) 'OPEN: ',NFTOUT,':',FILENAME,'has been created'
      ENDIF


!  2. Print out Fixed Length Header
      WRITE(NFTOUT,*)
      WRITE(NFTOUT,*)'                FIXED LENGTH HEADER'
      WRITE(NFTOUT,*)'                -------------------'
      WRITE(NFTOUT,*)
! DEPENDS ON: print_inte
      CALL PRINT_INTE(FIXHD,LEN_FIXHD,LEN_FIXHD,1,NFTOUT)

!  3. Print out Integer Header
      IF(LEN_INTHD >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               INTEGER HEADER'
        WRITE(NFTOUT,*)'               --------------'
        WRITE(NFTOUT,*)
! DEPENDS ON: print_inte
        CALL PRINT_INTE(INTHD,LEN_INTHD,LEN_INTHD,1,NFTOUT)
      ENDIF

!  4. Print out Real Header
      IF(LEN_REALHD >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'                REAL HEADER'
        WRITE(NFTOUT,*)'                -----------'
        WRITE(NFTOUT,*)
! DEPENDS ON: print_real
        CALL PRINT_REAL(REALHD,LEN_REALHD,LEN_REALHD,1,NFTOUT)
      ENDIF

!  5. Print out Level Dependent Constants
      IF(FIXHD(110) >  0 .AND. LEN2_LEVDEPC >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'                LEVEL DEPENDENT CONSTS'
        WRITE(NFTOUT,*)'                ----------------------'
        WRITE(NFTOUT,*)
        DO K=1,LEN2_LEVDEPC
        WRITE(NFTOUT,*)K,':'
! DEPENDS ON: print_real
        CALL PRINT_REAL(LEVDEPC,LEN1_LEVDEPC,LEN1_LEVDEPC,K,NFTOUT)
        ENDDO
      ENDIF

!  7. Print out Row Dependent Constants
      IF(FIXHD(115) >  0 .AND. LEN2_ROWDEPC >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'                  ROW DEPENDENT CONSTS'
        WRITE(NFTOUT,*)'                  --------------------'
        WRITE(NFTOUT,*)
        DO K=1,LEN2_ROWDEPC
        WRITE(NFTOUT,*)K,':'
! DEPENDS ON: print_real
        CALL PRINT_REAL(ROWDEPC,LEN1_ROWDEPC,LEN1_ROWDEPC,K,NFTOUT)
        ENDDO
      ENDIF

!  8. Print out Column Dependent Consts
      IF(FIXHD(120) >  0 .AND. LEN2_COLDEPC >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               COLUMN DEPENDENT CONSTS'
        WRITE(NFTOUT,*)'               -----------------------'
        WRITE(NFTOUT,*)
        DO K=1,LEN2_COLDEPC
        WRITE(NFTOUT,*)K,':'
! DEPENDS ON: print_real
        CALL PRINT_REAL(COLDEPC,LEN1_COLDEPC,LEN1_COLDEPC,K,NFTOUT)
        ENDDO
      ENDIF

!  9. Print out Field Dependent Consts
      IF(FIXHD(125) >  0 .AND. LEN2_FLDDEPC >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               FIELD DEPENDENT CONSTS'
        WRITE(NFTOUT,*)'               ----------------------'
        WRITE(NFTOUT,*)
        DO K=1,LEN2_FLDDEPC
        WRITE(NFTOUT,*)K,':'
! DEPENDS ON: print_real
        CALL PRINT_REAL(FLDDEPC,LEN1_FLDDEPC,LEN1_FLDDEPC,K,NFTOUT)
        ENDDO
      ENDIF

!  10. Print out Extra Constants
      IF(FIXHD(130) >  0 .AND. LEN_EXTCNST >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'                EXTRA CONSTS'
        WRITE(NFTOUT,*)'                ------------'
        WRITE(NFTOUT,*)
! DEPENDS ON: print_real
        CALL PRINT_REAL(EXTCNST,LEN_EXTCNST,LEN_EXTCNST,1,NFTOUT)
      ENDIF

!  11. Print out CFI1
      IF(FIXHD(140) >  0 .AND. LEN_CFI1 >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               COMPRESSED FIELD INDEX 1'
        WRITE(NFTOUT,*)'               ------------------------'
        WRITE(NFTOUT,*)
! DEPENDS ON: print_inte
        CALL PRINT_INTE(CFI1,LEN_CFI1,LEN_CFI1,1,NFTOUT)
      ENDIF

!  12. Print out CFI2
      IF(FIXHD(142) >  0 .AND. LEN_CFI2 >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               COMPRESSED FIELD INDEX 2'
        WRITE(NFTOUT,*)'               ------------------------'
        WRITE(NFTOUT,*)
! DEPENDS ON: print_inte
        CALL PRINT_INTE(CFI2,LEN_CFI2,LEN_CFI2,1,NFTOUT)
      ENDIF

!  12. Print out CFI3
      IF(FIXHD(144) >  0 .AND. LEN_CFI3 >  0)THEN
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               COMPRESSED FIELD INDEX 3'
        WRITE(NFTOUT,*)'               ------------------------'
        WRITE(NFTOUT,*)
! DEPENDS ON: print_inte
        CALL PRINT_INTE(CFI3,LEN_CFI3,LEN_CFI3,1,NFTOUT)
      ENDIF

!  13. Print out LOOKUP Headers
      IF(LEN2_LOOKUP >  0)THEN
      DO J=1,LEN2_LOOKUP
        DO I=1,LEN1_LOOKUP
          RLOOKUP(I,J)=LOOKUP(I,J)
        ENDDO
      ENDDO
        WRITE(NFTOUT,*)
        WRITE(NFTOUT,*)'               LOOKUP HEADERS'
        WRITE(NFTOUT,*)'               --------------'
        WRITE(NFTOUT,*)
        DO K=1,LEN2_LOOKUP
        IF (LOOKUP(1,K) /= -99) THEN
        WRITE(NFTOUT,*)K,':'
        WRITE(NFTOUT,*) 'Words 1-45'
! DEPENDS ON: print_inte
        CALL PRINT_INTE(LOOKUP(1,1),45,LEN1_LOOKUP,K,NFTOUT)
        WRITE(NFTOUT,*) 'Words 46-64'
! DEPENDS ON: print_real
        CALL PRINT_REAL(LOOKUP(46,1),19,LEN1_LOOKUP,K,NFTOUT)
        IF (LEN1_LOOKUP >  64) THEN
          WRITE(NFTOUT,*) 'Words 65-128'
! DEPENDS ON: print_inte
          CALL PRINT_INTE(LOOKUP(65,1),64,LEN1_LOOKUP,K,NFTOUT)
        ENDIF
        ENDIF
        ENDDO
      ENDIF

!  13. Print out individual fields
      if(.not. summary)then
      WRITE(6,*)
      WRITE(6,*)'               DATA FIELDS'
      WRITE(6,*)'               -----------'
      WRITE(6,*)

      IF (FIXHD(5) >= 8 .AND. FIXHD(5) <= 10) THEN !Cx/Cov/ObS

        WRITE (6,*)
        WRITE (6,*) 'Observation file : Observations not printed out'
        WRITE (6,*)

      ELSE


! Get decompostion information
! ----------------------------
! TOT_LEVELS not used in SX
        global_row_length=ROW_LENGTH
        global_rows=P_ROWS

! DEPENDS ON: decompose_smexe
        CALL DECOMPOSE_SMEXE(ROW_LENGTH, P_ROWS,                        &
     &                       0,0,TOT_LEVELS)
! DEPENDS ON: change_decomposition
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)

! READ_LAND_SEA call removed

      DO I=1,LEN2_LOOKUP

        IF (LOOKUP(1,I) /= -99) THEN
          FIELD_ITEM=MOD(LOOKUP(42,I),1000)
          FIELD_SECT=(LOOKUP(42,I)-FIELD_ITEM)/1000
          FIELD_MODEL=LOOKUP(45,I)

!         All diagnostics under model code of 10 are in section 20
!         of Atmos StashMaster file.
          IF (FIELD_MODEL == 10) THEN
            FIELD_MODEL = 1
          END IF

! DEPENDS ON: exppxi
          GRID_TYPE=EXPPXI(FIELD_MODEL,FIELD_SECT,FIELD_ITEM,           &
     &                            ppx_grid_type,                        &
#include "argppx.h"
     &                            ICODE,CMESSAGE)
          IF ((GRID_TYPE <  0).OR.(GRID_TYPE >  100)) THEN
            GRID_TYPE=1
!           WRITE(6,*)'PRINTDU: CANNOT GET GRID_TYPE INFO'
!           WRITE(6,*)'         FROM STASHMASTER FOR '
!           WRITE(6,*)'         SECTION ',FIELD_SECT,
!    &                ' ITEM ',FIELD_ITEM
!           WRITE(6,*)'         GRID_TYPE SET TO 1 (NORMAL GRID)'
          ENDIF
       lblrec_1=lookup(lblrec, i)

!        Since, for LBC data, header definition for pre Vn5.0
!        and Vn5.2- is different, a check is to be done.
         IF ((GRID_TYPE == ppx_atm_lbc_theta).OR.                       &
     &       (GRID_TYPE == ppx_atm_lbc_u).or.                           &
     &       (GRID_TYPE == ppx_atm_lbc_v).or.                           &
     &       (GRID_TYPE == ppx_atm_lbc_orog).or.                        &
     &       (GRID_TYPE == ppx_ocn_lbc_theta).or.                       &
     &       (GRID_TYPE == ppx_ocn_lbc_u)) then

           IF (FIXHD(12) >= 502) THEN
             IF (FIXHD(2) == 1) RIMWIDTHA=LOOKUP(41,i)/10000
             IF (FIXHD(2) == 2) RIMWIDTHO=LOOKUP(41,i)/10000
             haloY=MOD(LOOKUP(41,i),10000)
            haloY=haloY/100
             haloX=MOD(haloY,100)
             TOT_LEVELS=LOOKUP(17,i)-100
           ENDIF
         ENDIF

!       field is atmos/ocean LB field
        IF ((GRID_TYPE == ppx_atm_lbc_theta).OR.                        &
     &      (GRID_TYPE == ppx_atm_lbc_u).OR.                            &
     &      (GRID_TYPE == ppx_atm_lbc_v).OR.                            &
     &      (GRID_TYPE == ppx_atm_lbc_orog).OR.                         &
     &      (GRID_TYPE == ppx_ocn_lbc_theta).OR.                        &
     &      (GRID_TYPE == ppx_ocn_lbc_u)) THEN

          CURRENT_DECOMP_TYPE=-1

! DEPENDS ON: decompose_smexe
          CALL DECOMPOSE_SMEXE(ROW_LENGTH,P_ROWS,                       &
     &                       haloX, haloY,                              &
     &                       TOT_LEVELS)
! DEPENDS ON: change_decomposition
          CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)

        ENDIF

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                     &
     &                D1,MAX_FIELD_SIZE,FIXHD,                          &
#include "argppx.h"
     &               wgdos_expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('PRINTDUMP',CMESSAGE,ICODE,NFTIN)
        IF(FIXHD(5) == 5)THEN
! Boundary dataset. READFLDS does not write out max and min values
! for boundary datasets.
          F_MIN=D1(1)
          F_MAX=D1(1)
          POS_MAX=1
          POS_MIN=1
          DO J=1,LOOKUP(LBLREC,I)
            IF(D1(J) >  F_MAX)THEN
              F_MAX=D1(J)
              POS_MAX=J
            ENDIF
            IF(D1(J) <  F_MIN)THEN
              F_MIN=D1(J)
              POS_MIN=J
            END IF
          END DO

          WRITE(6,'('' MINIMUM='',E12.5,'' POSITION='',I8,              &
     &              '' MAXIMUM='',E12.5,'' POSITION='',I8)')            &
     &          F_MIN,POS_MIN,F_MAX,POS_MAX

          WRITE(6,'('' '')')
        END IF
        IF(LOOKUP(LBEXT,I) >  0) THEN ! got some extra data
          CALL FORT_GET_ENV("XPRINT",6,CXPRINT,5,ICODE)
          IF(ICODE /= 0) THEN
            NumXdata_print=5
          ELSE
            READ(CXPRINT,'(I5)') NumXdata_print
          ENDIF
          IF(NumXdata_print /= 0) THEN
! DEPENDS ON: print_extra
            CALL PRINT_EXTRA(LOOKUP,D1,NumXdata_print,I)
          ENDIF
        ENDIF

        IF (FIXHD(5)  ==  6) THEN ! ACOBS file
          NumMeta=5
          NumLevs=INT(LEVDEPC((I-1)*LEN1_LEVDEPC+2))
          NumObs=LOOKUP(66,I)
          CALL FORT_GET_ENV("NPRINT",6,CNPRINT,5,ICODE)
          IF (ICODE  /=  0) THEN
           WRITE(6,'(A33)') 'ERROR ENCOUNTERED IN FORT_GET_ENV'
           RETURN
          END IF
          READ(CNPRINT,'(I5)') NumObs_print
          NumObs_print=MIN(NumObs_print,NumObs)
          NBLK=INT(NumObs_print/8)
          WRITE (6,'(/,A10,I5,A6,I3,A13)') 'There are ',                &
     &           NumObs,' type ',                                       &
     &           LOOKUP(65,I),' observations'
          PEND=0
          DO BLK = 1, NBLK
           N=0
           PBEGIN=((BLK-1)*8)+1
           PEND=((BLK-1)*8)+8
           WRITE (6,*)
           WRITE (6,'(A12,8I12)')     'Observation:',                   &
     &           ((PBEGIN-1+IIII),IIII=1,8)
           WRITE (6,'(A12,8F12.2)') 'Latitude  : ',                     &
     &           (D1(0*NumObs+IIII),IIII=PBEGIN,PEND)
           WRITE (6,'(A12,8F12.2)') 'Longitude : ',                     &
     &           (D1(1*NumObs+IIII),IIII=PBEGIN,PEND)
           WRITE (6,'(A12,8F12.2)') 'Time      : ',                     &
     &           (D1(2*NumObs+IIII),IIII=PBEGIN,PEND)
           WRITE (6,'(A12,8F12.2)') 'MOT       : ',                     &
     &           (D1(3*NumObs+IIII),IIII=PBEGIN,PEND)
           IF (INT(LEVDEPC((I-1)*LEN1_LEVDEPC+1))  ==  1) THEN
            N=1
            WRITE (6,'(A12,8F12.2)') 'Pressure  : ',                    &
     &            (D1(5*NumObs+IIII)/100.,IIII=PBEGIN,PEND)
           END IF
           DO Lev = 1, NumLevs
            WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,                 &
     &            (D1((5+N)*NumObs+(NumObs*(Lev-1))+IIII),              &
     &                 IIII=PBEGIN,PEND)
           END DO
           DO Lev = 1, NumLevs
            WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,                 &
     &            (D1((5+N)*NumObs+(NumObs*(NumLevs+Lev-1))+IIII),      &
     &                 IIII=PBEGIN,PEND)
           END DO
           IF ( LOOKUP(65,I)  ==  301 .OR.                              &
     &          LOOKUP(65,I)  ==  302 .OR.                              &
     &          LOOKUP(65,I)  ==  303 .OR.                              &
     &          LOOKUP(65,I)  ==  304 .OR.                              &
     &          LOOKUP(65,I)  ==  305 .OR.                              &
     &          LOOKUP(65,I)  ==  306 ) THEN
             DO Lev = 1, NumLevs
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((5+N)*NumObs+(NumObs*(2*NumLevs+Lev-1))+IIII),  &
     &                   IIII=PBEGIN,PEND)
             END DO
           END IF
          END DO
          OBS_LEFT=NumObs_print - NBLK*8
          IF (OBS_LEFT  /=  0 .AND. NumObs  /=  0) THEN
           N=0
           WRITE (6,*)
           WRITE (6,'(A12,8I12)')     'Observation:',                   &
     &           (IIII,IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'Latitude  : ',                     &
     &           (D1(0*NumObs+IIII),IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'Longitude : ',                     &
     &           (D1(1*NumObs+IIII),IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'Time      : ',                     &
     &           (D1(2*NumObs+IIII),IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'MOT       : ',                     &
     &           (D1(3*NumObs+IIII),IIII=PEND+1,PEND+OBS_LEFT)
           IF (INT(LEVDEPC((I-1)*LEN1_LEVDEPC+1))  ==  1) THEN
            N=1
            WRITE (6,'(A12,8F12.2)') 'Pressure  : ',                    &
     &            (D1(5*NumObs+IIII)/100.,IIII=PEND+1,PEND+OBS_LEFT)
           END IF
           DO Lev = 1, NumLevs
             WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,                &
     &             (D1((5+N)*NumObs+(NumObs*(Lev-1))+IIII),             &
     &                  IIII=PEND+1,PEND+OBS_LEFT)
           END DO
           DO Lev = 1, NumLevs
            WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,                 &
     &            (D1((5+N)*NumObs+(NumObs*(NumLevs+Lev-1))+IIII),      &
     &                 IIII=PEND+1,PEND+OBS_LEFT)
           END DO
           IF ( LOOKUP(65,I)  ==  301 .OR.                              &
     &          LOOKUP(65,I)  ==  302 .OR.                              &
     &          LOOKUP(65,I)  ==  303 .OR.                              &
     &          LOOKUP(65,I)  ==  304 .OR.                              &
     &          LOOKUP(65,I)  ==  305 .OR.                              &
     &          LOOKUP(65,I)  ==  306 ) THEN
            DO Lev = 1, NumLevs
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((5+N)*NumObs+(NumObs*(2*NumLevs+Lev-1))+IIII),  &
     &                   IIII=PEND+1,PEND+OBS_LEFT)
            END DO
           END IF
          ENDIF
        ENDIF

        IF (FIXHD(5)  ==  7) THEN ! VAROBS file
          NumMeta=7
          NumItem=3
          NumObVariables=(LOOKUP(67,I)-NumMeta)/NumItem
          NumLevs=INT(LEVDEPC((I-1)*LEN1_LEVDEPC+2))
          Shift=NumMeta+(NumLevs*NumObVariables*NumItem)
          NumObs=LOOKUP(66,I)
          CALL FORT_GET_ENV("NPRINT",6,CNPRINT,5,ICODE)
          IF (ICODE  /=  0) THEN
           WRITE(6,'(A33)') 'ERROR ENCOUNTERED IN FORT_GET_ENV'
           RETURN
          END IF
          READ(CNPRINT,'(I5)') NumObs_print
          NumObs_print=MIN(NumObs_print,NumObs)
          NBLK=INT(NumObs_print/8)
          PEND=0
          DO BLK = 1, NBLK
           PBEGIN=((BLK-1)*8)+1
           PEND=((BLK-1)*8)+8
           WRITE (6,*)
           WRITE (6,'(A12,8I12)')     'Observation:',                   &
     &           ((PBEGIN-1+IIII),IIII=1,8)
           WRITE (6,'(A12,8F12.2)') 'Latitude  : ',                     &
     &           (D1((IIII-1)*Shift+1),IIII=PBEGIN,PEND)
           WRITE (6,'(A12,8F12.2)') 'Longitude : ',                     &
     &           (D1((IIII-1)*Shift+2),IIII=PBEGIN,PEND)
           WRITE (6,'(A12,8F12.2)') 'Time      : ',                     &
     &           (D1((IIII-1)*Shift+3),IIII=PBEGIN,PEND)
           WRITE (6,'(A12,8F12.2)') 'MOT       : ',                     &
     &           (D1((IIII-1)*Shift+4),IIII=PBEGIN,PEND)
           DO Variable = 1, NumObVariables
            WRITE (6,*)
            DO Lev = 1, NumLevs
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((IIII-1)*Shift+NumMeta+                         &
     &       (NumObVariables*NumItem*(Lev-1))+(Variable-1)*NumItem+1),  &
     &               IIII=PBEGIN,PEND)
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((IIII-1)*Shift+NumMeta+                         &
     &              (NumObVariables*NumItem*(Lev-1))+                   &
     &              (Variable-1)*NumItem+2),                            &
     &               IIII=PBEGIN,PEND)
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((IIII-1)*Shift+NumMeta+                         &
     &              (NumObVariables*NumItem*(Lev-1))+                   &
     &              (Variable-1)*NumItem+3),                            &
     &               IIII=PBEGIN,PEND)
              DO IPGE=1,8
                PGE(IPGE)=D1((PBEGIN+IPGE-2)*Shift+NumMeta+             &
     &            (NumObVariables*NumItem*(Lev-1))+                     &
     &            (Variable-1)*NumItem+3)
                IF (PGE(IPGE)  /=  RMDI) THEN
                  PGE2(IPGE)=INT(PGE(IPGE))
                  PGE1(IPGE)=PGE(IPGE)-PGE2(IPGE)
                  PGE2(IPGE)=PGE2(IPGE)/10000.0
                ELSE
                  PGE1(IPGE)=RMDI
                  PGE2(IPGE)=RMDI
                END IF
              END DO
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (PGE1(IIII),IIII=1,8)
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (PGE2(IIII),IIII=1,8)
            END DO
           END DO
          END DO
          OBS_LEFT=NumObs_print - NBLK*8
          IF (OBS_LEFT  /=  0 .AND. NumObs  /=  0) THEN
           WRITE (6,*)
           WRITE (6,'(A12,8I12)')     'Observation:',                   &
     &           (IIII,IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'Latitude  : ',                     &
     &           (D1((IIII-1)*Shift+1),IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'Longitude : ',                     &
     &           (D1((IIII-1)*Shift+2),IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'Time      : ',                     &
     &           (D1((IIII-1)*Shift+3),IIII=PEND+1,PEND+OBS_LEFT)
           WRITE (6,'(A12,8F12.2)') 'MOT       : ',                     &
     &           (D1((IIII-1)*Shift+4),IIII=PEND+1,PEND+OBS_LEFT)
           DO Variable = 1, NumObVariables
            WRITE (6,*)
            DO Lev = 1, NumLevs
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((IIII-1)*Shift+NumMeta+                         &
     &              (NumObVariables*NumItem*(Lev-1))+                   &
     &              (Variable-1)*NumItem+1),                            &
     &               IIII=PEND+1,PEND+OBS_LEFT)
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((IIII-1)*Shift+NumMeta+                         &
     &              (NumObVariables*NumItem*(Lev-1))+                   &
     &              (Variable-1)*NumItem+2),                            &
     &               IIII=PEND+1,PEND+OBS_LEFT)
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (D1((IIII-1)*Shift+NumMeta+                         &
     &              (NumObVariables*NumItem*(Lev-1))+                   &
     &              (Variable-1)*NumItem+3),                            &
     &               IIII=PEND+1,PEND+OBS_LEFT)
              DO IPGE=1,OBS_LEFT
                PGE(IPGE)=D1((IPGE+PEND-2)*Shift+NumMeta+               &
     &                   (NumObVariables*NumItem*(Lev-1))+              &
     &                   (Variable-1)*NumItem+3)
                IF (PGE(IPGE)  /=  RMDI) THEN
                  PGE2(IPGE)=INT(PGE(IPGE))
                  PGE1(IPGE)=PGE(IPGE)-PGE2(IPGE)
                  PGE2(IPGE)=PGE2(IPGE)/10000.0
                ELSE
                  PGE1(IPGE)=RMDI
                  PGE2(IPGE)=RMDI
                END IF
              END DO
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (PGE1(IIII),IIII=PEND+1,PEND+OBS_LEFT)
              WRITE (6,'(A6,I3,3X,8F12.3)') ' Level',Lev,               &
     &              (PGE2(IIII),IIII=PEND+1,PEND+OBS_LEFT)
            END DO
           END DO
          END IF
        END IF
       ENDIF
       lookup(lblrec, i)=lblrec_1

      ENDDO

      ENDIF
      endif

 9999 CONTINUE
      RETURN
      END SUBROUTINE PRINTDUMP
#endif
