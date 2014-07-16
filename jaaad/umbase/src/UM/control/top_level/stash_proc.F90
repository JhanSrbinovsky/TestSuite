#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Control routine for processing of basis library STASH file
!
! Subroutine Interface:

      SUBROUTINE STASH_PROC(NFTPPXREF,NFTSTMSTU  ,L_RECONF,             &
     &                                ppxRecs,                          &
     &                                ErrorStatus,CMESSAGE)
      IMPLICIT NONE

! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.  S.J.Swarbrick
!   4.1     Apr. 96    Generalise & incorporate
!                       wave model     S.J.Swarbrick
!   4.2     28/11/96   MPP code: Initialised new global_LPRIM and
!                      global_LDUMP variables.  P.Burton
!    vn4.4    9/4/97 Null string in argument list in call to GETPPX
!                    padded out to 13 characters (defined as
!                    CHARACTER*13 in GETPPX) to allow NAG f90 compiled
!                    code to run. IEdmond
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:

#include "lenfil.h"
#include "csubmodl.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "parparm.h"
#include "typsize.h"
#include "model.h"
#include "cstash.h"
#include "stextend.h"
#include "c_mdi.h"

! Subroutine arguments

!   Scalar arguments with intent(in):
      INTEGER NFTPPXREF   ! Unit no. for PPXREF file
      INTEGER NFTSTMSTU   ! Unit no. for user STASH master
      LOGICAL L_RECONF    ! Not used (but may be some day)

!   Array arguments with intent(out):
      CHARACTER*(80) CMESSAGE    ! Error return message

!   Error status:
      INTEGER        ErrorStatus ! +ve = fatal error

! Local scalars
      INTEGER I,J,L,II
      INTEGER RI           ! Row Index
      INTEGER Model
      INTEGER Section
      INTEGER Item
      INTEGER RowNumber    ! Row no. counter for PPXI, PPXC arrays
      INTEGER NLEVELS
      INTEGER NRECS
      INTEGER NTIMES
      INTEGER NPROFDP6
      PARAMETER(NPROFDP6=NPROFDP*6)

! Function and subroutine calls:
      EXTERNAL GETPPX,ADDRES,SETMODL                                    &
     &,PRELIM ,ORDER  ,SINDX ,DUPLIC,INACTR,TIMSER,                     &
     & POINTR ,OUTPTL ,INPUTL,WSTLST,RDBASIS

!- End of Header ---------------------------------------------------

!Initialisation
      NPSLISTS    =0   ! Counter for no. of pseudo levels lists
      NSERIES     =0   ! Time series block counter
      NSERREC_S   =0   ! Total no. of time series records
      NSERBLK_S   =0   ! Total no. of time series blocks
      ErrorStatus =0

      DO I=1,NPROFDP
        NRECS_TS(I)=0
         NPOS_TS(I)=0
      END DO

! Initialisation of data length arrays
      DO I=1,N_SUBMODEL_PARTITION_MAX
        LPRIM(I)=0   ! LENGTH OF PRIMARY DATA
        LDUMP(I)=0   ! LENGTH OF DUMP EXTENSION (DIAGNOSTIC)
#if defined(MPP)
        global_LPRIM(I)=0   ! LENGTH OF global PRIMARY DATA
        global_LDUMP(I)=0   ! LENGTH OF global DUMP EXTENSION
#endif
        LSECD(I)=0   ! LENGTH OF SECONDARY ATMOS
        LEXTRA(I)=0  ! LENGTH OF SPACE THAT IS ADDRESSED IN MODEL
        LWORK(I)=0   ! LENGTH OF WORK
        NHeadSub(I)=0! No. of pp headers for each submodel
      END DO
      DO I=1,N_INTERNAL_MODEL_MAX
        NHEAD  (I)=0 ! NUMBER OF PP HEADERS for each internal model
        LPrimIM(I)=0 ! Primary data length for each internal model
        LDumpIM(I)=0 ! Diagnostic do.
#if defined(MPP)
        global_LPrimIM(I)=0 ! GLOBAL Primary data length for each
!                           ! internal model
        global_LDumpIM(I)=0 ! GLOBAL Diagnostic do.
#endif
        LSecdIM(I)=0 ! Secondary do.
      END DO
      LPRIM_O2=0  ! LENGTH OF PRIMARY DATA OCEAN, SECOND TIME LEVEL

      DO I=OUTFILE_S,OUTFILE_E
        PPlen2LkUp(I)=4096
        FTOutUnit (I)=' '
      END DO


      DO I=1,NPROFDP6
        LLISTTY(I)=' '
      END DO


      DO J=1,NLEVP_S
      DO I=1,NLEVLSTSP
        RLEVLST_S(J,I)=RMDI
         LEVLST_S(J,I)=IMDI
      END DO
      END DO

      DO L=1,N_INTERNAL_MODEL_MAX
      DO J=0,NSECTP
      DO I=1,NITEMP
          IN_S(1,L,J,I)=0
          IN_S(2,L,J,I)=0
        INDX_S(1,L,J,I)=0
        INDX_S(2,L,J,I)=0
      END DO
      END DO
      END DO

      DO J=1,NELEMP+1
      DO I=1,NRECDP
        LIST_S(J,I)=0
      END DO
      END DO

      DO J=1,  NTIMEP
      DO I=1,2*NPROFTP+2
        ITIM_S(J,I)=-1
      END DO
      END DO

      DO J=1,N_INTERNAL_MODEL_MAX
      DO I=1,NITEMP
        PPIND_S(J,I)=0
      END DO
      END DO

      DO I=1,NPROFDP
        NRECS_TS(I)=0
         NPOS_TS(I)=0
      END DO

      DO I=1,NPSLISTP
        LENPLST(I)=0
      END DO

      DO I=OUTFILE_S,OUTFILE_E
        NHEAD_FILE(I)=0
      END DO

      DO I=1,N_INTERNAL_MODEL_MAX
      DO J=0,NSECTP
        H_VERS(I,J)=0
      END DO
      END DO

      DO I=1,MAX_AOBS
        AOBINC(I)=0
        AOBGRP(I)=0
      END DO

! Read stash basis file from job library
! DEPENDS ON: rdbasis
      CALL RDBASIS(4,CMESSAGE,ErrorStatus)
! Adjust stash time series records (if any)
      IF (NSERIES >  0) THEN
! DEPENDS ON: timser
        CALL TIMSER(CMESSAGE,ErrorStatus)
      END IF
      IF (ErrorStatus  /=  0) GO TO 9999

! Read STASHmaster files into look-up arrays PPXI, PPXC
      ErrorStatus = 0
      RowNumber   = 0
! Initialise arrays
      DO   I = 1,ppxRecs
        DO J = 1,PPXREF_CODELEN
          PPXI(I,J)   = 0
        END DO
        DO J = 1,PPXREF_CHARLEN
          PPXC(I,J)   = ' '
        END DO
      END DO
      DO I = 1,NDIAGP
        OriginFlag(I) =' '
        RowIndex  (I) = 0
      END DO
      DO II = 1,N_INTERNAL_MODEL
        DO   I = 0,PPXREF_SECTIONS
          DO J = 1,PPXREF_ITEMS
            PPXPTR(II,I,J) = 0
          END DO
        END DO
      END DO

      IF (INTERNAL_MODEL_INDEX(A_IM) >  0) THEN
! DEPENDS ON: getppx
        CALL GETPPX(NFTPPXREF,NFTSTMSTU,'STASHmaster_A',RowNumber,      &
#include "argppx.h"
     &                        ErrorStatus,CMESSAGE)
      END IF
      IF (INTERNAL_MODEL_INDEX(O_IM) >  0) THEN
! DEPENDS ON: getppx
        CALL GETPPX(NFTPPXREF,NFTSTMSTU,'STASHmaster_O',RowNumber,      &
#include "argppx.h"
     &                        ErrorStatus,CMESSAGE)
      END IF
      IF (INTERNAL_MODEL_INDEX(S_IM) >  0) THEN
! DEPENDS ON: getppx
        CALL GETPPX(NFTPPXREF,NFTSTMSTU,'STASHmaster_S',RowNumber,      &
#include "argppx.h"
     &                        ErrorStatus,CMESSAGE)
      END IF
      IF (INTERNAL_MODEL_INDEX(W_IM) >  0) THEN
! DEPENDS ON: getppx
        CALL GETPPX(NFTPPXREF,NFTSTMSTU,'STASHmaster_W',RowNumber,      &
#include "argppx.h"
     &                        ErrorStatus,CMESSAGE)
      END IF
!Read user STASHmaster files (which may be empty)
! DEPENDS ON: getppx
        CALL GETPPX(0,NFTSTMSTU,'             ',RowNumber,              &
#include "argppx.h"
     &                        ErrorStatus,CMESSAGE)

        IF (ErrorStatus  /=  0) GO TO 9999

! Define submodel and section/version configuration
! DEPENDS ON: setmodl
          CALL SETMODL(ErrorStatus,CMESSAGE)
          IF (ErrorStatus  /=  0) GO TO 9999

      NRECS=0
      NTIMES=0
      NLEVELS=0
      DO I=1,NPSLISTP
        LENPLST(I)=0
      END DO


! Construct preliminary STASH list
! DEPENDS ON: prelim
        CALL PRELIM(NRECS,                                              &
#include "argppx.h"
     &              NTIMES,NLEVELS,ErrorStatus,CMESSAGE)
        IF (ErrorStatus >  0) GO TO 9999

! REORDER STASH LIST & SET UP INDEX

! DEPENDS ON: order
        CALL ORDER(NRECS)
! DEPENDS ON: sindx
        CALL SINDX(NRECS)

! DELETE DUPLIC ENTRIES, CONCATONATE OVERLAP LEVELS ETC
! DELETE DUPLICATE STASH_TIMES, REORDER

! DEPENDS ON: duplic
        CALL DUPLIC(NRECS,NTIMES,NLEVELS)

! ADD INACTIVE AND IMPLIED RECORDS

! DEPENDS ON: inactr
        CALL INACTR(                                                    &
#include "argppx.h"
     &              NRECS,ErrorStatus,CMESSAGE)
          IF (ErrorStatus  /=  0) GO TO 9999

! REORDER STASH LIST & SET UP INDEX

! DEPENDS ON: order
        CALL ORDER(NRECS)
! DEPENDS ON: sindx
        CALL SINDX(NRECS)

! CHANGE POINTER SYSTEM, ADD ADDRESSES, LENGTHS AND INPUT LEVELS
! DEPENDS ON: pointr
        CALL POINTR(NRECS)

! OUTPUT LENGTH
! DEPENDS ON: outptl
        CALL OUTPTL(                                                    &
#include "argppx.h"
     &              NRECS,ErrorStatus,CMESSAGE)
          IF (ErrorStatus  /=  0) GO TO 9999

! INPUT LENGTH AND INPUT LEVELS, SET STLIST(NELEMP+1,I) TO MODEL_ST
! ALSO INPUT PSEUDO LEVELS
! DEPENDS ON: inputl
        CALL INPUTL(NRECS,                                              &
#include "argppx.h"
     &              NLEVELS,ErrorStatus,CMESSAGE)
          IF (ErrorStatus  /=  0) GO TO 9999


!     ADDRESSING
! DEPENDS ON: addres
      CALL ADDRES(                                                      &
#include "argppx.h"
     &            NRECS,ErrorStatus,CMESSAGE)
          IF (ErrorStatus  /=  0) GO TO 9999

! SET RETURN VALUES FOR OTHER FILES & write out STASH list

        NRECS_S=NRECS
        NTIMES_S=NTIMES
        NLEVL_S=NLEVELS
        ITEM_MAX_ALL=NITEMP
!       ITEM_MAX_REQ IS DONE IN WSTLIST
        NMAXLEV_S=1
        DO I =1,NLEVELS
          NMAXLEV_S=MAX(NMAXLEV_S,LEVLST_S(1,I))
        END DO

        NPSLISTS_S=NPSLISTS
        NMAXPSL_S=1
        DO I =1,NPSLISTS
          NMAXPSL_S=MAX(NMAXPSL_S,LENPLST(I))
        END DO
!       LSTUSER=NUSERD >= 1

! Assign values for STSIZES common block.
! Write output file (for checking purposes).

! DEPENDS ON: wstlst
        CALL WSTLST(NRECS,NTIMES,NLEVELS)

      DO I=1,NUM_DIAG_MAX
        IF (OriginFlag(I) == 'P'.OR.OriginFlag(I) == 'U') THEN
!   Determine model,section,item values for this row
          RI=RowIndex(I)
          Model  =     RI/100000
          Section=(RI-(RI/100000)*100000)/1000
          Item   =(RI-(RI/1000  )*1000  )
          IF (IN_S(1,Model,Section,Item)  ==  0) THEN
!   No entry in stash address list - overwrite this entry
            OriginFlag(I)=' '
          END IF
        END IF
      END DO
      DO J=1,NUM_DIAG_MAX
        IF (Originflag(J) /= ' ') THEN
          RI=RowIndex(J)
          Model  =     RI/100000
          Section=(RI-(RI/100000)*100000)/1000
          Item   =(RI-(RI/1000  )*1000  )
        END IF
      END DO
!   Close up gaps in OriginFlag array
      DO I=1,NUM_DIAG_MAX
        IF (OriginFlag(I) == ' ') THEN
          DO J=I+1,NUM_DIAG_MAX
            IF (OriginFlag(J) /= ' '.AND.                               &
     &          OriginFlag(I) == ' ') THEN
              OriginFlag(I) = OriginFlag(J)
              Originflag(J) = ' '
            END IF
          END DO
        END IF
      END DO

 9999 RETURN
      END SUBROUTINE STASH_PROC

!- End of Subroutine Code ----------------------------------------------
#endif
