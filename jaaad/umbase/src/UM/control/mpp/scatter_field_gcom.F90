#if defined(C96_1A) || defined(C96_1B) || defined (C96_1C) \
 || defined(UTILIO)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Scatters a field from one processor to many processors
!
! Subroutine Interface:
      SUBROUTINE SCATTER_FIELD_GCOM(LOCAL_FIELD,GLOBAL_FIELD,           &
     &                              LOCAL_ROW_LEN,LOCAL_ROWS,           &
     &                              GLOBAL_ROW_LEN,GLOBAL_ROWS,         &
     &                              GRID_TYPE,HALO_TYPE,                &
     &                              SCATTER_PE,PROC_GROUP,              &
     &                              ICODE,CMESSAGE)

      IMPLICIT NONE

!
! Description:
!  Takes a model field which is stored entirely on one processor
!  and distributes it over a group of processors using the
!  standard UM decomposition.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a scatter to all processors in the
!  group from the SCATTER_PE
!
! Current Code Owner: Paul Selwood
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  LOCAL_ROW_LEN                                                   &
                         ! IN length of rows in local part of field
     &, LOCAL_ROWS                                                      &
                         ! IN number of rows in local part of field
     &, GLOBAL_ROW_LEN                                                  &
                         ! IN length of rows in global field
     &, GLOBAL_ROWS                                                     &
                         ! IN number of rows in global field
     &, GRID_TYPE                                                       &
                         ! IN type (P,U or V) of grid
     &, HALO_TYPE                                                       &
                         ! IN halo type (hence width) of grid
     &, SCATTER_PE                                                      &
                         ! IN processor to scatter global field from
     &, PROC_GROUP                                                      &
                         ! IN group ID of processors involved here
     &, ICODE            ! OUT return code

      REAL                                                              &
     &  LOCAL_FIELD(LOCAL_ROW_LEN*LOCAL_ROWS)                           &
!                        ! OUT local part of field
     &, GLOBAL_FIELD(GLOBAL_ROW_LEN*GLOBAL_ROWS)
!                        ! IN (on PE GATHER_PE) global field
      CHARACTER*80                                                      &
     &  CMESSAGE         ! OUT error message

! Parameters and Common blocks

#include "parvars.h"
#include "gccom.h"

! Local variables

      INTEGER                                                           &
     &   send_map(7,MAXPROC)                                            &
     &,  receive_map(7,1)                                               &
     &,  n_mess_to_send                                                 &
     &,  info                                                           &
     &,  i

      INTEGER                                                           &
     &  old_GLOBAL_ROW_LEN                                              &
                              ! value on last call
     &, old_GLOBAL_ROWS                                                 &
                              ! value on last call
     &, old_PROC_GROUP                                                  &
                              ! value on last call
     &, old_SCATTER_PE                                                  &
                              ! value on last call
     &, old_DECOMP                                                      &
                              ! value on last call
     &, old_GRID_TYPE                                                   &
                              ! value on last call
     &, old_HALO_TYPE         ! value on last call

      SAVE send_map,receive_map,n_mess_to_send,                         &
     &     old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
     &     old_SCATTER_PE,old_DECOMP,                                   &
     &     old_GRID_TYPE,old_HALO_TYPE
      DATA old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
     &     old_SCATTER_PE,old_DECOMP,                                   &
     &     old_GRID_TYPE,old_HALO_TYPE                                  &
     &   / -1234, -1234, -1234, -1234, -1234, -1234, -1234/

      INTEGER                                                           &
     &  iproc                                                           &
     &, flag

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

      IF ((GLOBAL_ROW_LEN  /=  old_GLOBAL_ROW_LEN) .OR.                 &
     &    (GLOBAL_ROWS     /=  old_GLOBAL_ROWS   ) .OR.                 &
     &    (PROC_GROUP      /=  old_PROC_GROUP    ) .OR.                 &
     &    (SCATTER_PE      /=  old_SCATTER_PE    ) .OR.                 &
     &    (GRID_TYPE       /=  old_GRID_TYPE     ) .OR.                 &
     &    (HALO_TYPE       /=  old_HALO_TYPE     ) .OR.                 &
     &    (current_decomp_type  /=  old_DECOMP  )) THEN
!       Different arguments from the last call so we need
!       to calculate a new send/receive map

        IF (GRID_TYPE  ==  fld_type_unknown) THEN
          WRITE(6,*) 'SCATTER_FIELD_GCOM : Bad field type'
          WRITE(6,*) 'Field will not be scattered.'
          CMESSAGE='SCATTER_FIELD_GCOM : Bad field type'
          ICODE=1
          GOTO 9999
        ENDIF
! 2.0 Set up the send map (for PE SCATTER_PE only)

! Assume here that this group consists of all processors
! We'll get some new GCG functionality soon to improve this

        n_mess_to_send=0

        IF (mype  ==  SCATTER_PE) THEN
          DO iproc=0,nproc-1
            send_map(S_DESTINATION_PE,iproc+1) = iproc
            send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =            &
     &        g_datastart_f(1,grid_type,iproc)+                         &
     &        (g_datastart_f(2,grid_type,iproc)-1)*GLOBAL_ROW_LEN
            send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) =            &
     &        g_blsize(2,grid_type,iproc)
            send_map(S_STRIDE_IN_SEND_ARRAY,iproc+1) = GLOBAL_ROW_LEN
            send_map(S_ELEMENT_LENGTH,iproc+1) =                        &
     &        g_blsize(1,grid_type,iproc)
            send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =            &
     &        halosize(2,halo_type)*                                    &
     &        g_lasize(1,grid_type,halo_type,iproc) +                   &
     &        halosize(1,halo_type) + 1
            send_map(S_STRIDE_IN_RECV_ARRAY,iproc+1) =                  &
     &        g_lasize(1,grid_type,halo_type,iproc)
          ENDDO
          n_mess_to_send=nproc
        ENDIF

! 3.0 Set up the receive map

        receive_map(R_SOURCE_PE,1) = SCATTER_PE
        receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,1) =                   &
     &   halosize(2,halo_type)*LOCAL_ROW_LEN+1+                         &
     &   halosize(1,halo_type)
        receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,1) =                   &
     &    blsize(2,GRID_TYPE)
        receive_map(R_STRIDE_IN_RECV_ARRAY,1) = LOCAL_ROW_LEN
        receive_map(R_ELEMENT_LENGTH,1) = LOCAL_ROW_LEN-                &
     &      2*halosize(1,halo_type)
        receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,1) =                   &
     &   datastart_f(1,grid_type) +                                     &
     &  (datastart_f(2,grid_type)-1)*GLOBAL_ROW_LEN
        receive_map(R_STRIDE_IN_SEND_ARRAY,1) = GLOBAL_ROW_LEN


        old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
        old_GLOBAL_ROWS=GLOBAL_ROWS
        old_PROC_GROUP=PROC_GROUP
        old_SCATTER_PE=SCATTER_PE
        old_DECOMP=current_decomp_type
        old_GRID_TYPE=GRID_TYPE
        old_HALO_TYPE=HALO_TYPE

      ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

      flag=GC_NONE  ! This is currently ignored at GCG v1.1
      info=GC_NONE

!For small exe, if a field contains odd number
!global_row_len (eg OZONE 1*73), simply copy global_field to
! local_field


      CALL GCG_RALLTOALLE(GLOBAL_FIELD,send_map,n_mess_to_send,         &
     &                    GLOBAL_ROW_LEN*GLOBAL_ROWS,                   &
     &                    LOCAL_FIELD,receive_map,1,                    &
     &                    LOCAL_ROW_LEN*LOCAL_ROWS,                     &
     &                    PROC_GROUP,flag,info)


 9999 CONTINUE

      RETURN
      END SUBROUTINE SCATTER_FIELD_GCOM
#endif
