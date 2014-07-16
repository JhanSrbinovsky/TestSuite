#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Subroutine interface:
      SUBROUTINE WSTLST(NRECS,NTIMES,NLEVELS)
      IMPLICIT NONE
! Description: Print STASH control arrays [but also initialises some
!              STASH control variables]
!
! Method:  Simple interception of control arrays for printing.
!          [ Note that modularity would be improved if the function of
!            simply printing out variables was separated from the
!            initialisation of a no. of derived STASH variable
!            interspersed through the routine.]
!
! Current code owner:  S.J.Swarbrick
!
! History:
! Version   Date       Comment
! =======   ====       =======
!   3.5     Mar. 95    Original code.     S.J.Swarbrick
!   4.0     Sept.95                       S.J.Swarbrick
!   4.1     Apr. 96    Mods relating to code
!                       generalisation.   S.J.Swarbrick
!   4.2     29/11/96   MPP code : Added global version of A_LEN_DATA
!                      P.Burton
!   4.2      11/10/96  Enable atmos-ocean coupling for MPP.
!                      (2): Swap D1 memory. Add copies of D1 for atmos
!                      and ocean.                         R.Rawlins
!                      Initialise O_LEN_DUALDATA.         S.Ineson
!   4.5      05/08/98  Remove redundant code. Ft_Output for
!                      boundary files now initialised in
!                      INTF_CTL. D. Robinson.
!   5.0      21/05/99  STASH sizes now held separately. R.Rawlins
!   5.1      13/07/00  Make write statements of control arrays to
!                      standard output conditional on PrintStatus.
!                      R Rawlins
!LL  5.1  22/02/00  Add PARPARM for TYPSIZE                 P.Burton
!   5.4      21/03/02  Remove comment on same line as #include
!                                                   S. Carroll
!   5.5      02/08/00  Modification for parallelisation of WAM.
!                      Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
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
#include "version.h"
#include "stparam.h"
#include "cstash.h"
#include "stextend.h"
#include "parparm.h"
#include "typsize.h"
! For STASH sizes
#include "typstsz.h"
#include "model.h"
! Print status information in CPRINTST:
#include "cprintst.h"

      INTEGER  NUMBER_TIMES
      INTEGER  I,J,K
      INTEGER  ft_unit
      INTEGER  LEND1A
      INTEGER  LEND1O
      INTEGER  LEND1S
      INTEGER  LEND1W
      INTEGER  NLEVELS
      INTEGER  NRECS
      INTEGER  NTIMES
      INTEGER  BlkId           !Time series block identifier
      INTEGER  BlkSt           !Start position of ts block data
      INTEGER  IDP             !Domain profile loop counter
      INTEGER  IPOS
      INTEGER  Nrecs_prev      !No of recs in previous time ser block

      NAMELIST/STSIZES/                                                 &
     &        A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                        &
     &        S_LEN2_LOOKUP,S_LEN_DATA,                                 &
     &        O_LEN2_LOOKUP,O_LEN_DATA,O_LEN_DUALDATA,O_LEN_D1,         &
     &        W_LEN2_LOOKUP,W_LEN_DATA,W_LEN_D1,                        &
     &        LEN_TOT,                                                  &
     &        NSECTS,N_REQ_ITEMS,NITEMS,TOTITEMS,                       &
     &        NSTTABL,NUM_STASH_LEVELS,NUM_LEVEL_LISTS,                 &
     &        NUM_STASH_PSEUDO,NUM_PSEUDO_LISTS,                        &
     &        NSTTIMS,NSTASH_SERIES_BLOCK,                              &
     &        NSTASH_SERIES_RECORDS,N_PPXRECS
!
!
      LEND1A=0
      LEND1O=0
      LEND1S=0
      LEND1W=0

      IF(PrintStatus >= PrStatus_Oper) THEN

      write(6,120) NRECS
  120 FORMAT(I5,' STASH LIST RECORDS')
      DO 100 J=1,NRECS
      WRITE(6,*) ' '
      write(6,110)(LIST_S(I,J),I=1,NELEMP)
  110 FORMAT(10I8)
  100 CONTINUE
!
      write(6,310) NTIMES
  310 FORMAT(I4,' STASH TIMES')
      DO I=1,NTIMES
        NUMBER_TIMES=0
        DO J=1,NTIMEP
          IF (ITIM_S(J,I) >  0) THEN
            NUMBER_TIMES=NUMBER_TIMES+1
          END IF
        END DO
      write(6,320) (ITIM_S(J,I),J=1,NUMBER_TIMES)
      END DO
  320 FORMAT(100I5)
!
      write(6,410) NLEVELS
  410 FORMAT(I4,' STASH LEVELS LIST(S)')
      DO 400 I=1,NLEVELS
      write(6,420) LEVLST_S(1,I),LLISTTY(I)
  420 FORMAT(I5,A1)
      IF(LLISTTY(I) == 'I') THEN
      write(6,430) (LEVLST_S(J,I),J=2,LEVLST_S(1,I)+1)
  430 FORMAT(16I4)
      ELSE
      write(6,440) (RLEVLST_S(J,I),J=2,LEVLST_S(1,I)+1)
  440 FORMAT(6F12.3)
      END IF
  400 CONTINUE
!
!
      write(6,710) NPSLISTS
  710 FORMAT(I4,' STASH PSEUDO DIMENSION LIST(S)')
      DO 700 I=1,NPSLISTS
      write(6,720) LENPLST(I)
  720 FORMAT(I5)
      write(6,730) (PSLIST_D(J,I),J=1,LENPLST(I))
  730 FORMAT(16I4)
  700 CONTINUE
!
      write(6,510) NSERBLK_S
  510 FORMAT(I4,' STASH TIME SERIES BLOCKS')
!
      BlkSt =1
      DO IDP=1,NDPROF
        IF(NPOS_TS(IDP) >  0) THEN
          BlkId = NPOS_TS (IDP)
          IF (BlkId >  1) THEN
            BlkSt=BlkSt+Nrecs_prev
          END IF
          WRITE(6,530) NPOS_TS(IDP),NRECS_TS(NPOS_TS(IDP))
  530     FORMAT('SERIES NUMBER',I4,' WITH ',I4,' RECORDS')
          WRITE(6,*) '   NORTH   SOUTH    EAST    WEST  BOTTOM     TOP'
          DO IPOS=BlkSt,BlkSt+NRECS_TS(NPOS_TS(IDP))-1
            WRITE(6,560) NLIM_TS(IPOS),SLIM_TS(IPOS),ELIM_TS(IPOS),     &
     &                   WLIM_TS(IPOS),BLIM_TS(IPOS),TLIM_TS(IPOS)
          END DO
  560     FORMAT(6I8)
          WRITE(6,*)
        Nrecs_prev=NRECS_TS(NPOS_TS(IDP)) ! For next TS block
        END IF
      END DO

      write(6,*)                                                        &
     &' MODL SECT ITEM   IN_S(1)   IN_S(2) INDX_S(1) INDX_S(2)',        &
     &'   PPIND_S'
      write(6,*)                                                        &
     &'                  St addr   St  len StListPos StListNum '

      ENDIF ! on PrintStatus


      N_PPXRECS=0
      ITEM_MAX_REQ=1

      DO K=1,N_INTERNAL_MODEL_MAX
      DO J=0,44
      DO I=1,NITEMP
        IF(IN_S(1,K,J,I) /= 0) THEN
          N_PPXRECS=N_PPXRECS+1
          ITEM_MAX_REQ=MAX(J,ITEM_MAX_REQ)

          IF(PrintStatus >= PrStatus_Oper) THEN
          IF(J == 0) THEN
            write(6,610)                                                &
     &        K,J,I,  IN_S(1,K,J,I),  IN_S(2,K,J,I),                    &
     &              INDX_S(1,K,J,I),INDX_S(2,K,J,I),PPIND_S(K,I)
  610       FORMAT(3I5,5I10)
          ELSE
            write(6,610)                                                &
     &        K,J,I,  IN_S(1,K,J,I),  IN_S(2,K,J,I),                    &
     &              INDX_S(1,K,J,I),INDX_S(2,K,J,I)
          END IF
          ENDIF ! on PrintStatus

        END IF
      END DO
      END DO
      END DO

      IF(PrintStatus >= PrStatus_Oper) THEN
         I=-1
         write(6,610) I,I,I,I,I,I
      ENDIF ! on PrintStatus

! Variables in COMMON STSIZES - for UMINDEX routine.
!          N_PPXRECS was obtained in the loop above.

      A_LEN2_LOOKUP=NHEAD(A_IM)
      A_LEN_DATA   =MAX(1,LPrimIM(A_IM)+LDumpIM(A_IM))
#if defined(MPP)
      global_A_LEN_DATA   =                                             &
     &  MAX(1,global_LPrimIM(A_IM)+global_LDumpIM(A_IM))
#endif
      LEND1A       =LPrimIM(A_IM)+LDumpIM(A_IM)+LSecdIM(A_IM)           &
     &                                         +LEXTRA (A_SM)
      S_LEN2_LOOKUP=NHEAD(S_IM)
      S_LEN_DATA   =MAX(1,LPrimIM(S_IM)+LDumpIM(S_IM))
      LEND1S       =LPrimIM(S_IM)+LDumpIM(S_IM)+LSecdIM(S_IM)
      O_LEN2_LOOKUP=NHEAD(O_IM)
      O_LEN_DATA   =MAX(1,LPrimIM(O_IM)+LDumpIM(O_IM))
      LEND1O       =LPrimIM(O_IM)+LDumpIM(O_IM)+LSecdIM(O_IM)           &
     &                           +LPRIM_O2
#if defined(MPP)
      global_O_LEN_DATA   =                                             &
     &  MAX(1,global_LPrimIM(O_IM)+global_LDumpIM(O_IM))
#endif
      W_LEN2_LOOKUP=NHEAD(W_IM)
      W_LEN_DATA   =MAX(1,LPrimIM(W_IM)+LDumpIM(W_IM))
#if defined(MPP)
      global_W_LEN_DATA   =                                             &
     &  MAX(1,global_LPrimIM(W_IM)+global_LDumpIM(W_IM))
#endif
      LEND1W       =LPrimIM(W_IM)+LDumpIM(W_IM)+LSecdIM(W_IM)           &
     &                                   +LEXTRA (W_SM)
      LEN_TOT      =MAX(LEND1A+LEND1S,LEND1O,LEND1W)

      O_LEN_DUALDATA=LPRIM_O2

      A_LEN_D1=LEND1A
      O_LEN_D1=LEND1O
      W_LEN_D1=LEND1W

#if defined(SLAB)
      a_len2_lookup = a_len2_lookup + s_len2_lookup
      a_len_data    = a_len_data    + s_len_data
!      len_tot       = len_tot       + s_len_data
      write (6,*) ' wtslst ; slab: a_len2_lookup updated to ',          &
     &              a_len2_lookup
      write (6,*) ' wstlst ; slab: a_len_data updated to',              &
     &              a_len_data
      write (6,*) ' wstlst ; slab: len_tot updated to',                 &
     &              len_tot
      write (6,*) ' slab: stsizes with updated values'
      write (6,stsizes)
#endif
      NSECTS       =NSECTP
      N_REQ_ITEMS  =ITEM_MAX_REQ
      NITEMS       =NITEMP
      TOTITEMS     =MAX(1,NRECS)
      NSTTABL      =MAX(1,NTIMES)

      NUM_STASH_LEVELS=MAX(1,NMAXLEV_S)
      NUM_LEVEL_LISTS =MAX(1,NLEVL_S)
      NUM_STASH_PSEUDO=MAX(1,NMAXPSL_S)
      NUM_PSEUDO_LISTS=MAX(1,NPSLISTS_S)
      NSTTIMS         =NTIMEP

      NSTASH_SERIES_BLOCK =MAX(1,NSERBLK_S)
      NSTASH_SERIES_RECORDS=MAX(1,NSERREC_S)

!Assign values to PPlen2LkUp, FTOutUnit
      DO I = OUTFILE_S,OUTFILE_E
        PPlen2LkUp(I) = MAX(4096,NHEAD_FILE(I))
        IF ((NHEAD_FILE(I) >  0).AND.(I /= 27)) THEN
          FTOutUnit(I)='Y'
        ELSE
          FTOutUnit(I)='N'
        END IF
      END DO
      PPlen2LkUp(27) = 4096
!Output unit numbers for pp files and macros
      DO I = 1,NRECS
        IF (LIST_S(st_output_code,I) <  0) THEN
!Note that for pp files
! -LIST_S(st_output_code,I)=LIST_S(st_output_addr,I)=FT unit no.
          ft_unit             =-LIST_S(st_output_code,I)
          IF (ft_unit  /=  27) THEN
          PPlen2LkUp(ft_unit) = 4096
          FTOutUnit (ft_unit) = 'Y'
          END IF
        END IF
      END DO

      IF(PrintStatus >= PrStatus_Oper) THEN
         write(6,STSIZES)
      ENDIF ! on PrintStatus

      RETURN
      END SUBROUTINE WSTLST
!- End of subroutine code --------------------------------------------
#endif
