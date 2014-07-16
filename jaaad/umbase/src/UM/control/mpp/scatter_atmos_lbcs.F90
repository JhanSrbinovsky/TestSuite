#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Passes out atmosphere LBCs to processors at boundaries

! Subroutine Interface
      SUBROUTINE SCATTER_ATMOS_LBCS(                                    &
     &  FULL_LBC,DECOMP_LBC,                                            &
     &  FULL_LBC_SIZE,DECOMP_LBC_SIZE,                                  &
     &  FULL_LBC_LEVELS,DECOMP_LBC_LEVELS,                              &
     &  FLD_TYPE,HALO_TYPE,RIM_TYPE,                                    &
     &  PE_FOR_LEVEL,                                                   &
     &  ICODE,CMESSAGE)

      IMPLICIT NONE

! Description:
! Scatters atmosphere LBCs to relevant processors at the grid boundaries

! Method:
! L--Loop over all processors : iproc
! L L--Loop over sides (North,East,South,West) : iside
! L L I--IF processor iproc has an extremity at iside THEN
! L L I    Calculate all the sizes and offsets of FULL_LBC data and
! L L I                                             DECOMP_LBC data
! L L I L--Loop over levels : k
! L L I L I--IF I'm processor iproc THEN
! L L I L I    Setup a recv_map entry
! L L I L I--ENDIF
! L L I L I--IF I'm processor containing level k of FULL_LBC THEN
! L L I L I    Setup a send_map entry
! L L I L I--ENDIF
! L L I L--Endlooop : k
! L L I--ENDIF
! L L--Endloop : iside
! L--Endloop : iproc
!
! Call GCG_RALLTOALLE with the send_map/recv_map set up above to send
! the data from FULL_LBC to DECOMP_LBC
!

! Subroutine Arguments:

      INTEGER                                                           &
     &  FULL_LBC_SIZE                                                   &
                           ! IN single level size of the FULL_LBC array
     &, DECOMP_LBC_SIZE                                                 &
                           ! IN single level size of the DECOMP_LBC
                           !    array
     &, FULL_LBC_LEVELS                                                 &
                           ! IN number of levels of FULL_LBC on this
                           !    processor
     &, DECOMP_LBC_LEVELS                                               &
                           ! IN number of levels of DECOMP_LBC
     &, FLD_TYPE                                                        &
                           ! IN Which fld_type is the LBC?
     &, HALO_TYPE                                                       &
                           ! IN Which halo_type is the LBC?
     &, RIM_TYPE                                                        &
                           ! IN Which rim_type is the LBC?
     &, PE_FOR_LEVEL(DECOMP_LBC_LEVELS)                                 &
                           ! IN which level of FULL_LBC is on
                           !    which processor
     &, ICODE              ! OUT Return code

      REAL                                                              &
     &  FULL_LBC(FULL_LBC_SIZE,FULL_LBC_LEVELS)                         &
                           ! IN Some levels of the full LBC
     &, DECOMP_LBC(DECOMP_LBC_SIZE,DECOMP_LBC_LEVELS)
                           ! OUT All levels of the decomposed LBC on
                           !     this processor

      CHARACTER*(80)                                                    &
     &  CMESSAGE           ! OUT Error message


! Parameters and COMMON

#include "parvars.h"
#include "typsize.h"
#include "gccom.h"

! Local variables

      INTEGER                                                           &
     &  iproc                                                           &
                           ! loop counter for loop over processors
     &, iside                                                           &
                           ! loop counter for loop over sides
     &, k                                                               &
                           ! loop counter for levels
     &, full_lbc_row_len                                                &
                           ! length of a row of data on the full LBC
     &, full_lbc_nrows                                                  &
                           ! number of rows of data on the full LBC
     &, decomp_lbc_row_len                                              &
                           ! length of a row of data on the
                           ! decomposed LBC
     &, decomp_lbc_nrows                                                &
                           ! number of rows of data on the
                           ! decomposed LBC
     &, first_lbc_pt                                                    &
                           ! first point in full LBC to start
                           ! copying from
     &, first_lbc_row                                                   &
                           ! first row in fill LBC to start
                           ! copying from
     &, level_index_pe(0:nproc-1)                                       &
                           ! How many levels on each PE
     &, level_index(DECOMP_LBC_LEVELS)                                  &
                           ! Which level full_LBC corresponds
                           ! to the real level
     &, full_lbc_start_pt                                               &
                           ! First point index on a level of the
                           ! full LBC to start sending
     &, decomp_lbc_start_pt                                             &
                           ! First point index on a level of the
                           ! decomposed LBC to start receiving
     &, n_recv                                                          &
                           ! Number of items of data to receive
     &, n_send                                                          &
                           ! Number of items of data to send
     &, flag                                                            &
                           ! GCOM input flag (ignored)
     &, info                                                            &
                           ! GCOM return code

     &, send_map(7,FULL_LBC_LEVELS*nproc*4)                             &
                           ! Send map for all the data this
                           ! processor will be sending
     &, recv_map(7,DECOMP_LBC_LEVELS*4)
                           ! Receive map for all the data this
                           ! processor will be receiving
      INTEGER g_DECOMP_LBC_SIZE(0:nproc-1)

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! 1.0 Set up indexing describing where each level of full LBC data is
!     held

      g_DECOMP_LBC_SIZE=0
      g_DECOMP_LBC_SIZE(mype)=DECOMP_LBC_SIZE
      CALL GC_IMAX(nproc,NPROC,info,g_DECOMP_LBC_SIZE(0))

      DO iproc=0,nproc-1
        level_index_pe=0
      ENDDO

! Calculate indexing of full LBC data on the processors.
! A given level "k" of full LBC data will be held on
! processor "PE_FOR_LEVEL(k)" and the index of the level
! on this processor will be "level_index(k)"

      DO k=1,DECOMP_LBC_LEVELS
        level_index_pe(PE_FOR_LEVEL(k))=                                &
     &    level_index_pe(PE_FOR_LEVEL(k))+1
        level_index(k)=level_index_pe(PE_FOR_LEVEL(k))
      ENDDO

! These variables will describe how many receive entries and how many
! send entries the send/receive maps will have on this processor.
      n_recv=0
      n_send=0

!--------------------------------------------------------------------
! 2.0 Main loop, calculating the distribution of data and setting up
!     the send/receive maps

      DO iproc=0,nproc-1  ! loop over all processors

        DO iside=1,4      ! loop over each side

          IF (g_at_extremity(iside,iproc)) THEN  ! If this processor is
                                                 ! at edge type iside,
                                                 ! and so needs LBC data

            IF ((iside  ==  PNorth) .OR. (iside  ==  PSouth)) THEN

              ! Calculate size of FULL_LBC

              IF (FLD_TYPE  ==  fld_type_u) THEN
                full_lbc_row_len=glsize(1,FLD_TYPE)-1
              ELSE
                full_lbc_row_len=glsize(1,FLD_TYPE)
              ENDIF

              full_lbc_row_len=full_lbc_row_len+2*halosize(1,HALO_TYPE)

              full_lbc_nrows=halosize(2,HALO_TYPE)+RIMWIDTHA(RIM_TYPE)

              ! Calculate size of DECOMP_LBC

              IF ((FLD_TYPE  ==  fld_type_u) .AND.                      &
     &            (g_at_extremity(PEast,iproc))) THEN
                decomp_lbc_row_len=                                     &
     &            g_lasize(1,FLD_TYPE,HALO_TYPE,iproc)-1
              ELSE
                decomp_lbc_row_len=                                     &
     &            g_lasize(1,FLD_TYPE,HALO_TYPE,iproc)
              ENDIF

              decomp_lbc_nrows=halosize(2,HALO_TYPE)+                   &
     &                         RIMWIDTHA(RIM_TYPE)

              ! Calculate first point of DECOMP_LBC in FULL_LBC

              first_lbc_pt=g_datastart(1,iproc)
              first_lbc_row=1

            ELSE ! East or West boundaries

              ! Calculate size of FULL_LBC

              full_lbc_row_len=halosize(1,HALO_TYPE)+                   &
     &                         RIMWIDTHA(RIM_TYPE)
              full_lbc_nrows=glsize(2,FLD_TYPE)-2*RIMWIDTHA(RIM_TYPE)

              ! Calculate size of DECOMP_LBC

              decomp_lbc_row_len=halosize(1,HALO_TYPE)+                 &
     &                             RIMWIDTHA(RIM_TYPE)
              decomp_lbc_nrows=g_lasize(2,FLD_TYPE,HALO_TYPE,iproc)
              IF (g_at_extremity(PNorth,iproc))                         &
     &          decomp_lbc_nrows=decomp_lbc_nrows-                      &
     &                             halosize(2,HALO_TYPE)-               &
     &                             RIMWIDTHA(RIM_TYPE)
              IF (g_at_extremity(PSouth,iproc))                         &
     &          decomp_lbc_nrows=decomp_lbc_nrows-                      &
     &                             halosize(2,HALO_TYPE)-               &
     &                             RIMWIDTHA(RIM_TYPE)

              ! Calculate first point of DECOMP_LBC in FULL_LBC

              first_lbc_pt=1
              first_lbc_row=g_datastart(2,iproc)
              IF (.NOT. g_at_extremity(PSouth,iproc))                   &
     &          first_lbc_row=first_lbc_row-                            &
     &                          RIMWIDTHA(RIM_TYPE)-                    &
     &                          halosize(2,HALO_TYPE)

            ENDIF ! North/South or East/West boundary

            full_lbc_start_pt=                                          &
     &        global_LBC_STARTA(iside,FLD_TYPE,HALO_TYPE,RIM_TYPE)+     &
     &        (first_lbc_row-1)*full_lbc_row_len +                      &
     &        first_lbc_pt-1

            decomp_lbc_start_pt=                                        &
     &        g_LBC_STARTA(iside,FLD_TYPE,HALO_TYPE,RIM_TYPE,iproc)

            ! Now loop over levels and build up the send/receive maps
            ! for the LBC of the side

            DO k=1,DECOMP_LBC_LEVELS

              IF (mype  ==  iproc) THEN
                ! Set up a receive map entry

                n_recv=n_recv+1

                recv_map(R_SOURCE_PE,n_recv)=PE_FOR_LEVEL(k)
                recv_map(R_BASE_ADDRESS_IN_RECV_ARRAY,n_recv)=          &
     &            (k-1)*DECOMP_LBC_SIZE + decomp_lbc_start_pt
                recv_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,n_recv)=          &
     &            decomp_lbc_nrows
                recv_map(R_STRIDE_IN_RECV_ARRAY,n_recv)=                &
     &            decomp_lbc_row_len
                recv_map(R_ELEMENT_LENGTH,n_recv)=                      &
     &            decomp_lbc_row_len
                recv_map(R_BASE_ADDRESS_IN_SEND_ARRAY,n_recv)=          &
     &            (level_index(k)-1)*FULL_LBC_SIZE +                    &
     &            full_lbc_start_pt
                recv_map(R_STRIDE_IN_SEND_ARRAY,n_recv)=                &
     &            full_lbc_row_len

              ENDIF ! IF ((mype  ==  iproc) .AND. at_edge)

              IF (mype  ==  PE_FOR_LEVEL(k)) THEN
                ! Set up a send map entry

                n_send=n_send+1

                send_map(S_DESTINATION_PE,n_send)=iproc
                send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,n_send)=          &
     &            (level_index(k)-1)*FULL_LBC_SIZE +                    &
     &            full_lbc_start_pt
                send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,n_send)=          &
     &            decomp_lbc_nrows
                send_map(S_STRIDE_IN_SEND_ARRAY,n_send)=                &
     &            full_lbc_row_len
                send_map(S_ELEMENT_LENGTH,n_send)=                      &
     &            decomp_lbc_row_len
                send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,n_send)=          &
     &            (k-1)*g_DECOMP_LBC_SIZE(iproc) + decomp_lbc_start_pt
                send_map(S_STRIDE_IN_RECV_ARRAY,n_send)=                &
     &            decomp_lbc_row_len

              ENDIF ! IF (mype  ==  PE_FOR_LEVEL(k))

            ENDDO ! k : loop over levels

          ENDIF ! (g_at_extremity(iside,iproc))

        ENDDO ! iside : loop over sides

      ENDDO ! iproc : loop over processors

!--------------------------------------------------------------------
! 3.0 Now use the send/receive maps to perform the communication

      flag=GC_NONE
      info=GC_NONE

      CALL GCG_RALLTOALLE(FULL_LBC,send_map,n_send,                     &
     &                    FULL_LBC_SIZE*FULL_LBC_LEVELS,                &
     &                    DECOMP_LBC,recv_map,n_recv,                   &
     &                    DECOMP_LBC_SIZE*DECOMP_LBC_LEVELS,            &
     &                    gc_all_proc_group,flag,info)

      IF (info  /=  GC_NONE) THEN
        ICODE=1
        CMESSAGE='SCATTER_ATMOS_LBCS : GCG_RALLTOALLE failed'
        GOTO 9999
      ENDIF

!--------------------------------------------------------------------
 9999 CONTINUE

      RETURN
      END SUBROUTINE SCATTER_ATMOS_LBCS
#endif
