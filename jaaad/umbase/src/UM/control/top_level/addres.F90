#if defined(CONTROL)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Set the STASH addresses for D1
! Subroutine Interface:
      SUBROUTINE ADDRES(                                                &
#include "argppx.h"
     &                  NRECS,                                          &
     &                  ErrorStatus,CMESSAGE)
      IMPLICIT NONE
! Description:
!
! Method:
!
! Current code owner:  S.J.Swarbrick
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
! Global variables:
#include "csubmodl.h"
#include "lenfil.h"
#include "cppxref.h"
#include "ppxlook.h"
#include "decomptp.h"
#include "parvars.h"
#include "typsize.h"
#include "model.h"
#include "cstash.h"
#include "stextend.h"
#include "stparam.h"



! Subroutine arguments:
!   Scalar arguments with intent(in):
      INTEGER NRECS

!   Scalar arguments with intent(out):
      CHARACTER*80 CMESSAGE

! ErrorStatus:
      INTEGER ErrorStatus

! Local scalars:
      INTEGER TOTIMP
      INTEGER Im_ident  !Internal model identifier (absolute - CSMID)
      INTEGER Im_index  !Internal model index (expt. dependent)
      INTEGER Sm_ident  !Submodel identifier (absolute)
      INTEGER ISEC
      INTEGER ISEC_loop
      INTEGER IITM
      INTEGER RLEVS
      INTEGER RADDRESS
      INTEGER PIrow
      INTEGER I,J
      INTEGER IFIRST
      INTEGER IFREQ
      INTEGER IHOURS
      INTEGER ILAST
      INTEGER IREC
      INTEGER IH,IL,IP,IT
      INTEGER LWORK_S(N_SUBMODEL_PARTITION_MAX)
#if defined(MPP)
      INTEGER ICODE ! return from CHANGE_DECOMPOSITION
#endif

! Local arrays:
!    Submodel definitions array: stores list of Im_index's
!     for each submodel partition
      INTEGER                                                           &
     & SM_def(N_SUBMODEL_PARTITION_MAX,N_INTERNAL_MODEL_MAX)

! Function and subroutine calls:
      INTEGER  EXPPXI
      EXTERNAL PRIMARY,EXPPXI

!- End of Header ----------------------------------------------------


! 1.  Set STASHIN addresses and input lengths for primary fields

!   The address loop for primary fields is performed for each
!   internal model in turn. Hence, each internal model's primary
!   data occupies a contiguous block in D1. The order of these blocks
!   is the same as the order of the internal models given in the
!   array INTERNAL_MODEL_LIST.
!   User-defined prognostics are included in this primary addressing
!   routine, since they are incorporated into the ppxref lookup
!   arrays PPXI, PPXC in routine GETPPX.

!   Initialisation
      N_OBJ_D1_MAX=0
      DO I = 1,N_SUBMODEL_PARTITION_MAX
        N_OBJ_D1(I)=0
      ENDDO

      DO I = 1,N_SUBMODEL_PARTITION_MAX
      DO J = 1,N_INTERNAL_MODEL_MAX
        SM_def(I,J) = 0
      END DO
      END DO

!   Obtain submodel definitions and store in SMdef array
      DO Im_index = 1,N_INTERNAL_MODEL
!   Submodel ident.
         Sm_ident =   SUBMODEL_FOR_IM(Im_index)
!   Internal model index
         SM_def(Sm_ident,Im_index) = Im_index
      END DO

!   Primary address loop

!     Loop over submodel partitions
      DO  Sm_ident = 1,N_SUBMODEL_PARTITION_MAX

!       Initialise LEXTRA
        LEXTRA(Sm_ident)=0

!     Initialise address for reconfiguration
      RADDRESS = 1

!      Loop over internal models for each SM partition
       DO Im_index = 1,N_INTERNAL_MODEL

!       Test whether current SM contains this IM
        IF (SM_def(Sm_ident,Im_index) >  0) THEN

!        Obtain internal model identifier
         Im_ident   = INTERNAL_MODEL_LIST(Im_index)

#if defined(MPP)
! Set the correct decomposition in PARVARS

      ICODE=0

      IF (Im_ident  ==  A_IM) THEN
        IF (current_decomp_type  /=  decomp_standard_atmos)             &
! DEPENDS ON: change_decomposition
     &  CALL CHANGE_DECOMPOSITION(decomp_standard_atmos,ICODE)

      ELSE  ! unsupported decomposition type
        WRITE(6,*) 'ADDRES1 : Error - Only atmosphere ',                &
     &             'submodel is currently supported for UM code.'
        ErrorStatus=-1
        CMESSAGE='Unsupported submodel for UM code'
        GOTO 9999
      ENDIF

      IF (ICODE  /=  0) THEN
        WRITE(6,*) 'ADDRES1 : Error - Could not set decomposition ',    &
     &             'for selected submodel.'
        ErrorStatus=-2
        CMESSAGE='Unsupported decomposition selected for UM code'
        GOTO 9999
      ENDIF

#endif

!        Initialise primary data lengths for TYPSIZE
         IF (Im_ident == A_IM) A_PROG_LEN=0
         IF (Im_ident == A_IM) A_PROG_LOOKUP=0
           PIrow  = 0
      !Added section 3 for tracer fluxes - kdcorbin, 05/10
      DO ISEC_loop = 1,6   ! Currently there are five sections that
                           ! contain "primary" type fields
        IF (ISEC_loop  ==  1) ISEC=0  ! section zero primary
        IF (ISEC_loop  ==  2) ISEC=3  ! Boundary Layer Variables
        IF (ISEC_loop  ==  3) ISEC=31 ! LBC Input (treated as primary)
        IF (ISEC_loop  ==  4) ISEC=32 ! LBC Output(treated as primary)
        IF (ISEC_loop  ==  5) ISEC=33 ! Tracers (treated as primary)
        IF (ISEC_loop  ==  6) ISEC=34 ! UKCA tracers (like primary)
!       Loop over section zero items
        DO IITM   = 1,PPXREF_ITEMS
!   Check whether there is a primary field corresponding
!         to this item number
        IF (PPXPTR(Im_index,ISEC,IITM) /= 0) THEN
! DEPENDS ON: exppxi
          VMSK    = EXPPXI(Im_ident,ISEC,IITM,ppx_version_mask ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          ISPACE  = EXPPXI(Im_ident,ISEC,IITM,ppx_space_code   ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          IGP     = EXPPXI(Im_ident,ISEC,IITM,ppx_grid_type    ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          ILEV    = EXPPXI(Im_ident,ISEC,IITM,ppx_lv_code      ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          IBOT    = EXPPXI(Im_ident,ISEC,IITM,ppx_lb_code      ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          ITOP    = EXPPXI(Im_ident,ISEC,IITM,ppx_lt_code      ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
          DO I=1,6
! DEPENDS ON: exppxi
          IOPN(I) = EXPPXI(Im_ident,ISEC,IITM,ppx_opt_code+I-1 ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
          END DO
! DEPENDS ON: exppxi
          IFLAG   = EXPPXI(Im_ident,ISEC,IITM,ppx_lev_flag      ,       &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          IPSEUDO = EXPPXI(Im_ident,ISEC,IITM,ppx_pt_code       ,       &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          IPFIRST = EXPPXI(Im_ident,ISEC,IITM,ppx_pf_code       ,       &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          IPLAST  = EXPPXI(Im_ident,ISEC,IITM,ppx_pl_code       ,       &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          halo_type = EXPPXI(Im_ident,ISEC,IITM,ppx_halo_type,          &
#include "argppx.h"
     &                       ErrorStatus,CMESSAGE)
          IF((ISPACE == 2).OR.(ISPACE == 3).OR.(ISPACE == 9)            &
     &                    .OR.(ISPACE == 4)                             &
     &                    .OR.(ISPACE == 5)                             &
     &                    .OR.(ISPACE == 10)                            &
     &                    .OR.(ISPACE == 8)) THEN ! Primary variable
! DEPENDS ON: primary
            CALL PRIMARY(ISEC,IITM,Im_index,Im_ident,Sm_ident,          &
     &                  RLEVS,RADDRESS,PIrow,ErrorStatus,CMESSAGE)
          END IF
        END IF    !  PPXPTR(m,s,i)  /=  0
        END DO    !  Loop over items
        END DO    !  ISEC_loop : Loop over sections
        END IF    !  test whether SM contains IM
       END DO     !  Loop over Im_index
      END DO      !  Loop over SM partitions

! LOOKUP array lengths for TYPSIZE
      A_PROG_LOOKUP = NHEAD(A_IM)
! Primary data lengths for TYPSIZE
      A_PROG_LEN = LPrimIM(A_IM)
      WRITE(6,*) ' ADDRES : A_PROG_LOOKUP = ',A_PROG_LOOKUP
      WRITE(6,*) ' ADDRES : A_PROG_LEN    = ',A_PROG_LEN
      WRITE(6,*) ' ADDRES : NHEAD,LPrimIM = ',NHEAD,LPrimIM

! 2. Loop through stash list to set output addresses and
!                 header positions for diagnostics
      DO IREC=1,NRECS

! Read internal model number from stash list. Stash list has already
! been ordered by internal model, section, item. Thus, all the atmos
! diagnostic addressing will be done first, followed by the slab
! addressing in the case of a slab model.
        Im_ident = LIST_S(st_model_code,IREC)
! Obtain submodel partition id.
        Sm_ident = SUBMODEL_PARTITION_INDEX(Im_ident)

! Set output address relative to D1
        IF(LIST_S(st_output_code,IREC) == 1) THEN

! Diagnostic output to dump rather than direct output pp file
!   Add the output length for this diag to LDUMP; total length of
!   dump so far = LPRIM + LDUMP; hence obtain the start address for
!   the output from the next diagnostic to be stored in dump.

          LIST_S(st_output_addr,IREC)                                   &
     &             = LPRIM(Sm_ident)+LDUMP(Sm_ident)+1
! Information for preliminary D1 addressing array
          N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
          IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN)THEN
            D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=diag
            D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
            D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IREC
          ENDIF
#if defined(MPP)
          LIST_S(st_dump_output_addr,IREC)=                             &
     &             global_LPRIM(Sm_ident)+global_LDUMP(Sm_ident)+1
#endif
          LDUMP(Sm_ident)                                               &
     &             = LDUMP(Sm_ident)+LIST_S(st_output_length,IREC)
          LDumpIM(Im_ident)                                             &
     &             = LDumpIM(Im_ident)+LIST_S(st_output_length,IREC)
#if defined(MPP)
          global_LDUMP(Sm_ident)=                                       &
     &      global_LDUMP(Sm_ident)+LIST_S(st_dump_output_length,IREC)
          global_LDUMPIM(Sm_ident)=                                     &
     &      global_LDUMPIM(Im_ident)+LIST_S(st_dump_output_length,IREC)
#endif

          IF(LIST_S(st_output_bottom,IREC) == 100) THEN
! Special levels
            RLEVS=1
          ELSE IF(LIST_S(st_series_ptr,IREC) /= 0) THEN
! Time series domain
            RLEVS=1
          ELSE IF(LIST_S(st_gridpoint_code,IREC) >= 10                  &
     &       .AND.LIST_S(st_gridpoint_code,IREC) <  20) THEN
! Vertical ave.
            RLEVS=1
          ELSE  IF(LIST_S(st_output_bottom,IREC) <  0) THEN
! Levels list
            RLEVS=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
          ELSE
! Range of model levels
            RLEVS=LIST_S(st_output_top   ,IREC)                         &
     &           -LIST_S(st_output_bottom,IREC)+1
          END IF

          IF (LIST_S(st_pseudo_out,IREC) >  0) THEN
! Pseudo levels
            RLEVS=RLEVS*LENPLST(LIST_S(st_pseudo_out,IREC))
          END IF

! Set position of pp lookup header in the dump
          LIST_S(st_lookup_ptr,IREC)=NHeadSub(Sm_ident)+1

! Increment NHEAD (there is one pp header for each level at
!  which a diagnostic is output
          NHEAD   (Im_ident)=NHEAD   (Im_ident)+RLEVS
          NHeadSub(Sm_ident)=NHeadSub(Sm_ident)+RLEVS

        ELSE IF(LIST_S(st_output_code,IREC) == 2) THEN

! Secondary data in D1.
! Compute and store secondary data lengths. Start address for
! secondary data is determined below, after total dump
! diagnostic length has been found.

          LIST_S(st_output_addr,IREC)=LSECD(Sm_ident)+1
          LSECD(Sm_ident)                                               &
     &   =LSECD(Sm_ident)+LIST_S(st_output_length,IREC)
          LSecdIM(Im_ident)                                             &
     &   =LSecdIM(Im_ident)+LIST_S(st_output_length,IREC)
! Set pointer for pp header
          LIST_S(st_lookup_ptr,IREC)=-1

        ELSE IF(LIST_S(st_output_code,IREC) <  0) THEN

! Diagnostic output to PP file

! Compute no. of pp headers for this diagnostic
!   = output levels * pseudo output levels * output times

!   No. of levels
          IF(LIST_S(st_output_bottom,IREC) == 100) THEN
! Special levels
            IL=1
          ELSE IF(LIST_S(st_series_ptr,IREC) /= 0) THEN
! Time series dom
            IL=1
          ELSE IF(LIST_S(st_gridpoint_code,IREC) >= 10                  &
     &      .AND.LIST_S(st_gridpoint_code,IREC) <  20) THEN
! Vertical average
            IL=1
          ELSE  IF(LIST_S(st_output_bottom,IREC) <  0) THEN
! Levels list
            IL=LEVLST_S(1,-LIST_S(st_output_bottom,IREC))
          ELSE
! Range of mod levs
            IL=LIST_S(st_output_top,IREC)                               &
     &       -LIST_S(st_output_bottom,IREC)+1
          END IF

!   No. of pseudo levels
          IF (LIST_S(st_pseudo_out,IREC) >  0) THEN
            IP=LENPLST(LIST_S(st_pseudo_out,IREC))
          ELSE
            IP=1
          END IF

!   No. of output times
          IF(LIST_S(st_freq_code,IREC) >  0) THEN
            IFIRST=LIST_S(st_start_time_code,IREC)
            IFREQ =LIST_S(st_freq_code      ,IREC)
            IF(LIST_S(st_end_time_code,IREC) == -1) THEN
! Output to continues to end of run
              IHOURS=1+8760*RUN_TARGET_END(1)                           &
     &                + 744*RUN_TARGET_END(2)                           &
     &                +  24*RUN_TARGET_END(3)                           &
     &                +     RUN_TARGET_END(4)
! DEPENDS ON: totimp
              ILAST=TOTIMP(IHOURS,'H ',Im_ident)
              if (ILAST  ==  -999) then
                 errorStatus = 1
                 cmessage = 'TOTIMP:UNEXPECTED TIME UNIT or             &
     &                IRREGULAR DUMPS FOR DUMP FREQUENCY'
                 GOTO 9999
              endif
            ELSE
! Last output time before end of run
              ILAST=LIST_S(st_end_time_code,IREC)
            END IF

            IT= 1 + (ILAST-IFIRST)/IFREQ
            IF (IT <  0) THEN
              IT=0
              WRITE(6,*)                                                &
     &      ' Output time error detected in routine ADDRESS:'
              WRITE(6,*)                                                &
     &      ' Output time starts after specified end of run'
              WRITE(6,*)                                                &
     &      ' STASH record no.,MODEL,SECTION,ITEM as follows: ',        &
     &                        IREC, LIST_S(st_model_code,IREC),         &
     &                              LIST_S(st_sect_code ,IREC),         &
     &                              LIST_S(st_item_code ,IREC)
              WRITE(6,*) 'OUTPUT CODE: ',                               &
     &                              LIST_S(st_output_code,IREC)
            END IF
          ELSE
! Times table in STASH_times array
            IT=1
            DO I=1,NTIMEP
              IF (ITIM_S(I,-LIST_S(st_freq_code,IREC)) == -1) THEN
                IT=I-1
                GOTO 260
              END IF
            END DO
 260        CONTINUE
          END IF
! No. of output "headers" - (levels)*(pseudo-levels)*(output times)
          IH=IL*IP*IT
! Assign output unit no. (nn) to (st_output_addr)
          LIST_S(st_output_addr,IREC)=-LIST_S(st_output_code,IREC)
! Assign no. of output headers to NHEAD_FILE(nn)
          NHEAD_FILE(LIST_S(st_output_addr,IREC))=                      &
     &    NHEAD_FILE(LIST_S(st_output_addr,IREC)) + IH
        ELSE IF (LIST_S(st_output_code,IREC) == 0) THEN
! Inactive record, not output
          LIST_S(st_output_addr,IREC)=-LIST_S(st_output_code,IREC)
        ELSE
          WRITE(6,*) 'ERROR detected in routine ADDRESS '
          WRITE(6,*) 'ILLEGAL OUTPUT CODE FOR STASH RECORD '
          WRITE(6,*)                                                    &
     &  ' STASH record no.,MODEL,SECTION,ITEM as follows: ',            &
     &                   IREC, LIST_S(st_model_code,IREC),              &
     &                         LIST_S(st_sect_code ,IREC),              &
     &                         LIST_S(st_item_code ,IREC)
        END IF

      END DO      ! End of loop over records for D1 addressing


!     Correct the addressing of SPACE=9 items from being relative
!     to start of LEXTRA space to being relative to start of dump

!     Loop over submodel partitions
      DO  Sm_ident = 1,N_SUBMODEL_PARTITION_MAX

!       Loop over internal models for each SM partition
        DO Im_index = 1,N_INTERNAL_MODEL

!         Test whether current SM contains this IM
          IF (SM_def(Sm_ident,Im_index) >  0) THEN

!           Obtain internal model identifier
            Im_ident   = INTERNAL_MODEL_LIST(Im_index)

            !Added tracer fluxes - kdcorbin, 05/10
            DO ISEC_loop=1,6

            IF (ISEC_loop  ==  1) ISEC=0  ! section zero primary
            IF (ISEC_loop  ==  2) ISEC=3  ! Boundary Layer Variables
            IF (ISEC_loop  ==  3) ISEC=31 ! LBC Input (primary)
            IF (ISEC_loop  ==  4) ISEC=32 ! LBC Output (primary)
            IF (ISEC_loop  ==  5) ISEC=33 ! Tracers (primary)
            IF (ISEC_loop  ==  6) ISEC=34 ! UKCA tracers (primary)
            DO IITM   = 1,PPXREF_ITEMS
!             Check whether there is a primary field corresponding
              IF (PPXPTR(Im_index,ISEC,IITM) /= 0) THEN
! DEPENDS ON: exppxi
                ISPACE  = EXPPXI(Im_ident,ISEC,IITM,ppx_space_code,     &
#include "argppx.h"
     &            ErrorStatus,CMESSAGE)
                IF (IN_S(1,Im_ident,ISEC,IITM) /= 0                     &
                                                      ! item is active
     &              .and. ISPACE == 9) THEN
                  IN_S(1,Im_ident,ISEC,IITM)=IN_S(1,Im_ident,ISEC,IITM)+&
     &              LPRIM(Sm_ident)+LDUMP(Sm_ident)
                  IF (Im_ident == O_IM) THEN
                    IN_S(1,Im_ident,ISEC,IITM)=                         &
     &                IN_S(1,Im_ident,ISEC,IITM)+LPRIM_O2
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
            ENDDO ! ISEC_loop
          ENDIF
        ENDDO
      ENDDO


! Set secondary data addresses relative to start of D1
      DO IREC=1,NRECS
        Im_ident = LIST_S(st_model_code,IREC)
        Sm_ident = SUBMODEL_PARTITION_INDEX(Im_ident)

        IF (LIST_S(st_output_code,IREC) == 2) THEN
          LIST_S(st_output_addr,IREC)  =LIST_S(st_output_addr,IREC)     &
     &  + LPRIM(Sm_ident)+LDUMP(Sm_ident)+LEXTRA(Sm_ident)
! Information for preliminary D1 addressing array
          N_OBJ_D1(Sm_ident)     =N_OBJ_D1(Sm_ident)+1
          IF (N_OBJ_D1(Sm_ident) <= MAX_D1_LEN)THEN
            D1_PADDR(d1_type,N_OBJ_D1(Sm_ident),Sm_ident)=seco
            D1_PADDR(d1_im,N_OBJ_D1(Sm_ident),Sm_ident)=Im_ident
            D1_PADDR(d1_extra_info,N_OBJ_D1(Sm_ident),Sm_ident)=IREC
          ENDIF
          IF (Im_ident == O_IM) THEN
            LIST_S(st_output_addr,IREC)=LIST_S(st_output_addr,IREC)     &
     &  +   LPRIM_O2
          END IF
        END IF
      END DO

! 3.  Set input addresses and work lengths for non-primary
!            fields (i.e., ISPACE=0,1,6 or 7)
      DO Im_ident=1,N_INTERNAL_MODEL_MAX
         Sm_ident=  SUBMODEL_PARTITION_INDEX(Im_ident)
      DO ISEC    =0,PPXREF_SECTIONS
! Re-initialise sectional work lengths
        DO I=1,N_SUBMODEL_PARTITION_MAX
          LWORK_S(I)=0
        END DO
        DO IITM  =1,PPXREF_ITEMS
          IF(INDX_S(2,Im_ident,ISEC,IITM) >  0) THEN
! Item in STASH list
! Obtain space code & section zero point-back code
!   from ppxref lookup array
! DEPENDS ON: exppxi
          ISPACE  = EXPPXI(Im_ident,ISEC,IITM,ppx_space_code   ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
! DEPENDS ON: exppxi
          PTR_PROG= EXPPXI(Im_ident,ISEC,IITM,ppx_ptr_code     ,        &
#include "argppx.h"
     &                                          ErrorStatus,CMESSAGE)
            IF ( (ISPACE == 0).OR.(ISPACE == 1).OR.                     &
     &           (ISPACE == 6).OR.(ISPACE == 7) ) THEN
! Compute length of work space required
              IF (ISPACE /= 7) THEN
! STASH_WORK address & length
                IN_S(1,Im_ident,ISEC,IITM)=LWORK_S(Sm_ident)+1
                LWORK_S(Sm_ident)=LWORK_S(Sm_ident)                     &
     &                           +IN_S(2,Im_ident,ISEC,IITM)
              ELSE
! Point-back to primary space in section 0
                IN_S(1,Im_ident,ISEC,IITM    )                          &
     &         =IN_S(1,Im_ident,0   ,PTR_PROG)
                IN_S(2,Im_ident,ISEC,IITM    )                          &
     &         =IN_S(2,Im_ident,0   ,PTR_PROG)
              END IF
            END IF
          END IF
        END DO   ! Items

! Find max sectional work length for each submodel partition
        DO I=1,N_SUBMODEL_PARTITION_MAX
          LWORK(I)=MAX(LWORK(I),LWORK_S(I))
        END DO

      END DO     ! Sections
      IF(Sm_ident /= 0)THEN
!       Save the maximum value for dimensioning full D1 address array
        N_OBJ_D1_MAX=MAX(N_OBJ_D1_MAX,N_OBJ_D1(Sm_ident))
        WRITE(6,*)N_OBJ_D1(Sm_ident),' D1 items in submodel ',Sm_ident
      ENDIF
      END DO     ! Models
      IF(N_OBJ_D1_MAX >  MAX_D1_LEN)THEN
        WRITE(6,*)'ADDRES1: No of items in D1 exceeds maximum allowed:'
        WRITE(6,*)'Number allowed ',MAX_D1_LEN,' Number requested '     &
     &    ,N_OBJ_D1_MAX
        WRITE(6,*)'Modify the COMDECK STEXTEND to increase'
        WRITE(6,*)'MAX_D1_LEN parameter as required'
        WRITE(6,*)'Such a change can be safely made'
        CMESSAGE='ADDRES1: No of D1 items exceeds max: See output'
        ErrorStatus=1
      ENDIF

 9999 CONTINUE
      RETURN
      END SUBROUTINE ADDRES
#endif
