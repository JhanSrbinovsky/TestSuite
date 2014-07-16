#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine INITDUMP -------------------------------------------
!LL
!LL Purpose:To read atmosphere or ocean dumps, and to calculate
!LL additional constants based on the dump header information.
!LL
!LL Extra constants needed for cloud types calulated within SETDCFLD
!LL
!LL Level 2 control routine for Cray YMP
!LL
!LL Programming Standard : UM documentation paper no. 3
!LL                        version no. 1, dated 15/01/90
!LL
!LL System components covered : R30,C26
!LL System task : P0
!LL
!LL Documentation : U.M. Documentation Paper no. P0.
!LL                 U.M. Documentation paper no F3,draft version
!LL                 number 3, dated 18/12/89
!LL
!LLEND--------------------------------------------------------------
!
!*L Arguments

      SUBROUTINE INITDUMP(                                              &
#include "argd1.h"
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
#include "argcona.h"
#include "arglndm.h"
#include "argppx.h"
     &             sm_ident,ICODE,CMESSAGE)

      IMPLICIT NONE

!*L Arguments
!L
#include "cmaxsize.h"
#include "csubmodl.h"
#include "parvars.h"
#include "typsize.h"
#include "nstypes.h"
#include "typd1.h"
#include "typduma.h"
! Contains *CALL CPPXREF
#include "typsts.h"
#include "typptra.h"
#include "typcona.h"
#include "typlndm.h"

      INTEGER   sm_ident  ! Sub-model indicator
      INTEGER   ICODE     ! Return code
      INTEGER   NFTIN     ! FTN number for read
      INTEGER   NFTSWAP   ! FTN number for swapping radiation incrs
      INTEGER   SEGSTART  ! Pointer to start of radiation incrs
      INTEGER   SEGEND    ! Pointer to end of radiation incrs
      INTEGER   LEN_IO    ! Length of data transferred
      INTEGER   ERROR     ! Error code returned by OPEN

      REAL      A_IO      ! IO completion code

      CHARACTER*80                                                      &
     &          CMESSAGE  ! Error message

#include "chsunits.h"
#include "ccontrol.h"
#include "cenvir.h"
#include "chistory.h"
#include "c_mdi.h"
#include "clookadd.h"
#include "c_global.h"
#include "c_writd.h"
#include "cruntimc.h"

#include "ppxlook.h"

      LOGICAL                                                           &
     &   L_A_DUMP                                                       &
     &  ,L_A_PROG_ONLY 
                            ! ) Switches set if only prognostic
                            ! ) fields to be read in.

      INTEGER   i,j,ij    ! Loop counters
      INTEGER                                                           &
     &   LEN2_LOOKUP                                                    &
     &  ,LEN_DATA                                                       &
     &  ,N_PROG_FLDS                                                    &
                            ! No of prognostic fields in dumps
     &  ,N_PROG_LOOKUP                                                  &
     &  ,LEN_PROG                                                       &
     &  ,TOT_LEN_DATA                                                   &
     &  ,D1_ADDR_SUBMODEL_ID  ! submodel id in D1_ADDR array
#if defined(ATMOS)
      INTEGER :: A_MPP_ADDR(A_LEN2_LOOKUP),                             &
     &           A_MPP_LEN(A_LEN2_LOOKUP)
#endif
#include "decomptp.h"
      INTEGER info  ! return code from GC operations

      REAL                                                              &
     & Dummyarg             ! Dummy argument needed for legal call


#include "comocflw.h"

      INTEGER                                                           &
     &        IFLD                                                      &
                            ! Loop variable
     &       ,cox_tracer                                                &
                            !
     &       ,IM_IDENT                                                  &
                            ! internal model identifier
     &       ,IM_INDEX      ! internal model index for STASH arrays

!*---------------------------------------------------------------------
!L   Internal Structure

#if defined(ATMOS)

!L 1.0 Read atmosphere dump and initialise atmosphere model.
      IF (sm_ident == atmos_sm) THEN

!L 1.1 Open unit for atmosphere dump, read fixed length header
!L     and set buffer length
        NFTIN = 21
! DEPENDS ON: file_open
      CALL FILE_OPEN(NFTIN,FT_ENVIRON(NFTIN),                           &
     &               LEN_FT_ENVIR(NFTIN),0,0,ERROR)

! DEPENDS ON: read_flh
        CALL READ_FLH (NFTIN,A_FIXHD,LEN_FIXHD,ICODE,CMESSAGE)
        IF (ICODE >  0) RETURN

! DEPENDS ON: setpos
        CALL SETPOS (NFTIN,0,ICODE)

!       Test if atmos dump.
        L_A_DUMP = A_FIXHD(5) == 1 .AND. A_FIXHD(2) == atmos_sm

!       Test if only prognostic fields to be read in
        L_A_PROG_ONLY = L_A_DUMP .AND. H_STEPim(a_im) == 0

!       Get no of prognostic fields in atmos dump
        N_PROG_FLDS = A_FIXHD(153)

!       Check N_PROG_FLDS has been set.
        IF (N_PROG_FLDS == IMDI) THEN
          WRITE (6,*) '  '
          WRITE (6,*) ' No of prognostic fields not set in FIXHD(153)'
          WRITE (6,*) ' Run RECONFIGURATION to set FIXHD(153)'
          CMESSAGE = 'INITDUMP: FIXHD(153) not set in atmos dump'
          ICODE = 101
          GO TO 9999  !  Return
        ENDIF

!       Check N_PROG_FLDS matches with A_PROG_LOOKUP set up by the UI
        N_PROG_LOOKUP = A_PROG_LOOKUP
        IF (N_PROG_FLDS /= N_PROG_LOOKUP) THEN
          WRITE (6,*) ' '
          WRITE (6,*) ' Mismatch in no of prognostic fields.'
          WRITE (6,*) ' No of prog fields in Atmos dump ',N_PROG_FLDS
          WRITE (6,*) ' No of prog fields expected      ',N_PROG_LOOKUP
          WRITE (6,*) ' '
          WRITE (6,*) ' Run RECONFIGURATION to get correct no of',      &
     &                ' prognostic fields in atmos dump'
          WRITE (6,*) ' or'
          WRITE (6,*) ' Check/Reset experiment in User Interface'
          WRITE (6,*) ' '
          CMESSAGE = 'INITDUMP: Wrong no of atmos prognostic fields'
          ICODE = 102
          GO TO 9999  !  Return
        ENDIF

! Initialise D1 to prevent uninitialised data in unused rows of U fields
! *DIR$ CACHE_BYPASS D1
          DO I = 1,LEN_TOT
            D1(I)=0.0
          ENDDO
!       Determine no of fields to be read in
        IF (L_A_PROG_ONLY) THEN

!         Prognostic fields only
          LEN2_LOOKUP = N_PROG_FLDS
          LEN_PROG = A_PROG_LEN
          TOT_LEN_DATA = A_LEN_DATA

          LEN_DATA = LEN_PROG

          WRITE (6,*) ' '
          WRITE (6,*) ' Read in ',N_PROG_FLDS,' prognostic fields.'


!      INITIALISE DIAGNOSTIC AREA OF D1 TO RMDI
       DO I = LEN_DATA+1, TOT_LEN_DATA
         D1(I)=RMDI
       END DO

        ELSE

!         All fields.
          LEN2_LOOKUP = A_LEN2_LOOKUP
          LEN_DATA    = A_LEN_DATA

          WRITE (6,*) ' '
          WRITE (6,*) ' Read in all ',LEN2_LOOKUP,' fields.'

        ENDIF

! Ensure that domain decomposition is consistent with submodel

! DEPENDS ON: change_decomposition
      CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)


!L 1.2 Call READDUMP to read atmosphere dump.
        IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('READDUMP',3)

        END IF


      D1_ADDR_SUBMODEL_ID = SUBMODEL_FOR_SM(atmos_sm)

! DEPENDS ON: um_readdump
      CALL UM_READDUMP(NFTIN, A_FIXHD, LEN_FIXHD,                       &
     &    A_INTHD, A_LEN_INTHD,                                         &
     &    A_REALHD, A_LEN_REALHD,                                       &
     &    A_LEVDEPC, A_LEN1_LEVDEPC, A_LEN2_LEVDEPC,                    &
     &    A_ROWDEPC, A_LEN1_ROWDEPC, A_LEN2_ROWDEPC,                    &
     &    A_COLDEPC, A_LEN1_COLDEPC, A_LEN2_COLDEPC,                    &
     &    A_FLDDEPC, A_LEN1_FLDDEPC, A_LEN2_FLDDEPC,                    &
     &    A_EXTCNST, A_LEN_EXTCNST,                                     &
     &    A_DUMPHIST, LEN_DUMPHIST,                                     &
     &    A_CFI1, A_LEN_CFI1,                                           &
     &    A_CFI2, A_LEN_CFI2,                                           &
     &    A_CFI3, A_LEN_CFI3,                                           &
     &    A_LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                             &
     &    A_MPP_LOOKUP,MPP_LEN1_LOOKUP,                                 &
     &         atmos_sm,                                                &
     &         NO_OBJ_D1(D1_ADDR_SUBMODEL_ID),                          &
     &         D1_ADDR(1,1,D1_ADDR_SUBMODEL_ID),                        &
     &         LEN_DATA,D1,                                             &
#include "argppx.h"
     &    .TRUE.)

        IF (LTIMER) THEN
! DEPENDS ON: timer
          CALL TIMER('READDUMP',4)

        END IF

! DEPENDS ON: file_close
      CALL FILE_CLOSE(NFTIN,FT_ENVIRON(NFTIN),                          &
     &                LEN_FT_ENVIR(NFTIN),0,0,ICODE)

        IF (ICODE  >   0) RETURN

! Check validity of integer header data and print out information.
! Pass through the global numbers so the validity check works
! glsize(1) is the global ROW_LENGTH
! glsize(2) is the global ROWS
! DEPENDS ON: pr_inhda
        CALL PR_INHDA (A_INTHD, A_LEN_INTHD, glsize(1,fld_type_p),      &
     &     glsize(2,fld_type_p),MODEL_LEVELS,WET_LEVELS,TR_LEVELS,      &
     &     ST_LEVELS, SM_LEVELS, BL_LEVELS,                             &
     & TR_VARS, ICODE, CMESSAGE)

        IF (ICODE >  0) RETURN

! Check validity of real header data and print out information.
! DEPENDS ON: pr_rehda
        CALL PR_REHDA (A_REALHD, A_LEN_REALHD)

! DEPENDS ON: check_idealise
        CALL CHECK_IDEALISE                                             &
     &  ( A_LEN_INTHD, A_LEN_REALHD, A_LEN1_LEVDEPC, A_LEN2_LEVDEPC,    &
     &    A_INTHD, A_REALHD, A_LEVDEPC )

        IF (L_A_PROG_ONLY) THEN

!         Need to pass field address and length info to ADDRESS_CHECK
          DO I=1,LEN2_LOOKUP
            A_MPP_ADDR(I) = A_MPP_LOOKUP(P_NADDR,I)
            A_MPP_LEN(I)  = A_MPP_LOOKUP(P_LBLREC,I)
          ENDDO

! DEPENDS ON: address_check
          CALL ADDRESS_CHECK (A_LOOKUP,A_MPP_ADDR,                      &
     &      A_MPP_LEN,LEN1_LOOKUP,LEN2_LOOKUP,                          &
     &                        SI,NITEMS,NSECTS,LEN_DATA,                &
#include "argppx.h"
     &                        ICODE,CMESSAGE)
          IF (ICODE >  0) RETURN

        ENDIF

! ------------------------------------------------------------
!       Check that packing codes in lookup table is consistent
!       with packing requested in DUMP_PACKim.
! --------------------------------------------------------------

! DEPENDS ON: check_dump_packing
        CALL CHECK_DUMP_PACKING (                                       &
     &       A_FIXHD, LEN_FIXHD,                                        &
     &       A_LOOKUP, LEN1_LOOKUP, LEN2_LOOKUP,                        &
#include "argppx.h"
     &       DUMP_PACKim(sm_ident), ATMOS_IM )

!       Reset A_FIXHD to correspond to Output Dump
        A_FIXHD(152) = A_LEN2_LOOKUP
        A_FIXHD(160) = A_FIXHD(150) + LEN1_LOOKUP*A_LEN2_LOOKUP
        A_FIXHD(161) = global_A_LEN_DATA

!L 1.3 Call SET_ATM_POINTERS to set integer pointers to data in
!L     atmosphere dump and secondary storage area in D1 array.
! DEPENDS ON: set_atm_pointers
      CALL SET_ATM_POINTERS(                                            &
#include "argduma.h"
#include "argsts.h"
#include "argptra.h"
     &                  ICODE,CMESSAGE)

!L Call READLSTA to read namelists to control atmosphere integration
!L and diagnostic point print.
        REWIND 5
! DEPENDS ON: readlsta
      CALL READLSTA(                                                    &
#include "arglndm.h"
     &Dummyarg)

      IF (ICODE >  0) RETURN
! SETCONA is now called from INITIAL
! Set ELF flag
        ELF=(A_FIXHD(4) == 3.OR.A_FIXHD(4) == 103)

      END IF

#endif

 9999 CONTINUE

      RETURN
      END SUBROUTINE INITDUMP


#endif
