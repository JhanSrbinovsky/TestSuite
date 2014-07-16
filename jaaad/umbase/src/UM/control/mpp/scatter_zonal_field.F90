#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
  || defined(UTILIO)
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

      SUBROUTINE SCATTER_ZONAL_FIELD (                                  &
     &  LOCAL_FIELD , GLOBAL_FIELD ,                                    &
     &  LOCAL_SIZE  , GLOBAL_SIZE  ,                                    &
     &  LEVELS, GRID_CODE, GRID_TYPE ,HALO_TYPE,SCATTER_PE)

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
     &, SCATTER_PE      ! IN:  PE on which GLOBAL_FIELD resides

      REAL                                                              &
     &  GLOBAL_FIELD(GLOBAL_SIZE,LEVELS)                                &
                                         ! IN : field to scatter
     &, LOCAL_FIELD(LOCAL_SIZE,LEVELS) ! OUT : local part of field

! Parameters and common blocks

#include "cppxref.h"
#include "gccom.h"

#include "parvars.h"

! Local variables

      INTEGER                                                           &
     &  fld_type                                                        &
                  ! type (P or U) of field
     &, info                                                            &
                  ! GCOM return code
     &, send_map(7,MAXPROC)                                             &
                             ! send map
     &, receive_map(7,1)                                                &
                             ! receive map
     &, flag                                                            &
                             ! GCOM dummy argument
     &, n_mess_to_send                                                  &
                             ! number of messages to send
     &, k,iproc     ! loop counter

      LOGICAL                                                           &
     &  mead_fld  ! indicates an ocean Mead diagnostic field

!====================================================================

      mead_fld=((GRID_CODE  ==  ppx_ocn_uzonal) .OR.                    &
     &          (GRID_CODE  ==  ppx_ocn_tzonal))

!--------------------------------------------------------------------

      n_mess_to_send=0

      IF (mype  ==  SCATTER_PE) THEN
        DO iproc=0,nproc-1
          send_map(S_DESTINATION_PE,iproc+1) = iproc
          send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =              &
     &      g_datastart(2,iproc)
          send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = 1
          send_map(S_STRIDE_IN_SEND_ARRAY,iproc+1) = 0
          send_map(S_ELEMENT_LENGTH,iproc+1) =                          &
     &      g_blsize(2,grid_type,iproc)
          IF (mead_fld) THEN
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) = 1
          ELSE
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =            &
     &        halosize(2,halo_type)+1
          ENDIF
          send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) = 0
        ENDDO
        n_mess_to_send=nproc
      ENDIF

      receive_map(R_SOURCE_PE,1) = SCATTER_PE
      IF (mead_fld) THEN
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) = 1
      ELSE
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) =                   &
     &    halosize(2,halo_type)+1
      ENDIF
      receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) = 1
      receive_map(R_STRIDE_IN_RECV_ARRAY,1) = 0
      receive_map(R_ELEMENT_LENGTH,1) = blsize(2,grid_type)
      receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) = datastart(2)
      receive_map(R_STRIDE_IN_SEND_ARRAY,1) = 0

      DO k=1,LEVELS

        info=GC_NONE
        flag=GC_NONE

          CALL GCG_RALLTOALLE(                                          &
     &      GLOBAL_FIELD(1,k),send_map,n_mess_to_send,                  &
     &      glsize(2,grid_type),                                        &
     &      LOCAL_FIELD(1,k), receive_map,1,                            &
     &      lasize(2,grid_type,halo_type),                              &
     &      GC_ALL_PROC_GROUP,flag,info)

      ENDDO

      RETURN

      END SUBROUTINE SCATTER_ZONAL_FIELD
#endif
