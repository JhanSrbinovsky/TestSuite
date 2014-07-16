#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO) || defined(FLDIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! + Scatters zonal field from one processor to many processors
!
! Subroutine interface:

      SUBROUTINE GATHER_ZONAL_FIELD (                                   &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE  , GLOBAL_SIZE  ,                                    &
     &  LEVELS, GRID_CODE, GRID_TYPE ,HALO_TYPE,                        &
     &  GATHER_PE)

      IMPLICIT NONE

! Description:
! Takes a zonal field on a single processor, and decomposes it over
! many processors.
!
! Current code owner : P.Burton
!
! Subroutine arguments:

      INTEGER                                                           &
     &  LOCAL_SIZE                                                      &
                        ! IN: size of level of LOCAL_FIELD
     &, GLOBAL_SIZE                                                     &
                        ! IN: size of level of GLOBAL FIELD
     &, LEVELS                                                          &
                        ! IN: number of levels in field
     &, GRID_CODE                                                       &
                        ! IN: ppx grid code of field
     &, GRID_TYPE                                                       &
                        ! IN: P,U or V
     &, HALO_TYPE                                                       &
                        ! IN: halo type of field
     &, GATHER_PE       ! IN:  PE to gather GLOBAL_FIELD

      REAL                                                              &
     &  LOCAL_FIELD(LOCAL_SIZE,LEVELS)                                  &
                                         ! IN : field to gather
     &, GLOBAL_FIELD(GLOBAL_SIZE,LEVELS) ! OUT : gathered field


! Parameters and common blocks

#include "cppxref.h"
#include "gccom.h"

#include "parvars.h"

! Local variables

      INTEGER                                                           &
     &  info                                                            &
                  ! GCOM return code
     &, send_map(7,1)                                                   &
                       ! send map
     &, receive_map(7,MAXPROC)                                          &
                                   ! receive map
     &, flag                                                            &
                             ! GCOM dummy argument
     &, n_mess_to_receive                                               &
                               ! number of messages to receive
     &, n_mess_to_send                                                  &
                               ! number of messages to send
     &, k,iproc     ! loop counter



!====================================================================


!--------------------------------------------------------------------

      n_mess_to_receive=0

      IF (mype  ==  GATHER_PE) THEN
        DO iproc=0,nproc-1
          IF (g_gridpos(1,iproc)  ==  0) THEN
!           Only one processor per LPG row needs to send the data
!           as it will be the same for each processor along the
!           row.
            n_mess_to_receive=n_mess_to_receive+1
            receive_map(R_SOURCE_PE,n_mess_to_receive) = iproc
            receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,                   &
     &                  n_mess_to_receive) =                            &
     &                            g_datastart(2,iproc)
            receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,                   &
     &                  n_mess_to_receive) = 1
            receive_map(R_STRIDE_IN_RECV_ARRAY,                         &
     &                  n_mess_to_receive) = 0
            receive_map(R_ELEMENT_LENGTH,n_mess_to_receive) =           &
     &        g_blsize(2,GRID_TYPE,iproc)
            receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,                   &
     &          n_mess_to_receive) = halosize(2,HALO_TYPE)+1
            receive_map(R_STRIDE_IN_SEND_ARRAY,                         &
     &                  n_mess_to_receive) = 0
          ENDIF
        ENDDO
      ENDIF

      n_mess_to_send=0
        IF (at_extremity(PWest)) THEN
          ! only processors at the left of the LPG will  send anything

          send_map(S_DESTINATION_PE,1) = GATHER_PE
          send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =                    &
     &        halosize(2,HALO_TYPE)+1
          send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1) = 1
          send_map(S_STRIDE_IN_SEND_ARRAY,1) = 0
          send_map(S_ELEMENT_LENGTH,1) = blsize(2,GRID_TYPE)
          send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) = datastart(2)
          send_map(S_STRIDE_IN_RECV_ARRAY,1) = 0

          n_mess_to_send=1

        ENDIF


      DO k=1,LEVELS

        info=GC_NONE
        flag=GC_NONE

        CALL GCG_RALLTOALLE(                                            &
     &    LOCAL_FIELD(1,k),send_map,n_mess_to_send,                     &
     &    LOCAL_SIZE,                                                   &
     &    GLOBAL_FIELD(1,k),receive_map,n_mess_to_receive,              &
     &    GLOBAL_SIZE,GC_ALL_PROC_GROUP,flag,info)

      ENDDO

      RETURN

      END SUBROUTINE GATHER_ZONAL_FIELD
#endif
