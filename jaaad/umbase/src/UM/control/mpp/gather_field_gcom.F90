#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C)  \
 || defined(UTILIO) || defined(FLDIO) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers a field from many processors to one processor
!
! Subroutine Interface:
      SUBROUTINE GATHER_FIELD_GCOM(LOCAL_FIELD,GLOBAL_FIELD,       &
     &                             LOCAL_ROW_LEN,LOCAL_ROWS,       &
     &                             GLOBAL_ROW_LEN,GLOBAL_ROWS,     &
     &                             GRID_TYPE,HALO_TYPE,            &
     &                             GATHER_PE,PROC_GROUP,           &
     &                             ICODE,CMESSAGE)

      IMPLICIT NONE

!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field.
!
! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE
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
     &, GATHER_PE                                                       &
                         ! IN processor to gather global field to
     &, PROC_GROUP                                                      &
                         ! IN group ID of processors involved here
     &, ICODE            ! OUT return code

      REAL                                                              &
     &  LOCAL_FIELD(LOCAL_ROW_LEN*LOCAL_ROWS)                           &
!                        ! IN local part of field
     &, GLOBAL_FIELD(GLOBAL_ROW_LEN*GLOBAL_ROWS)
!                        ! OUT (on PE GATHER_PE) global field

      CHARACTER*80                                                      &
     &  CMESSAGE         ! OUT error message

! Parameters and Common blocks

#include "parvars.h"
#include "gccom.h"

! Local variables

      INTEGER                                                           &
     &   send_map(7,1)                                                  &
     &,  receive_map(7,MAXPROC)                                         &
     &,  n_mess_to_rec

      INTEGER                                                           &
     &  old_GLOBAL_ROW_LEN                                              &
                              ! value on last call
     &, old_GLOBAL_ROWS                                                 &
                              ! value on last call
     &, old_PROC_GROUP                                                  &
                              ! value on last call
     &, old_GATHER_PE                                                   &
                              ! value on last call
     &, old_DECOMP                                                      &
                              ! value on last call
     &, old_GRID_TYPE                                                   &
                              ! value on last call
     &, old_HALO_TYPE         ! value on last call

      SAVE send_map,receive_map,n_mess_to_rec,                          &
     &     old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
     &     old_GATHER_PE,old_DECOMP,                                    &
     &     old_GRID_TYPE,old_HALO_TYPE
      DATA old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
     &     old_GATHER_PE,old_DECOMP,                                    &
     &     old_GRID_TYPE,old_HALO_TYPE                                  &
     &   / -1234, -1234, -1234, -1234, -1234, -1234, -1234/

      INTEGER                                                           &
     &  iproc                                                           &
     &, info                                                            &
     &, flag

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

      IF ((GLOBAL_ROW_LEN  /=  old_GLOBAL_ROW_LEN) .OR.                 &
     &    (GLOBAL_ROWS     /=  old_GLOBAL_ROWS   ) .OR.                 &
     &    (PROC_GROUP      /=  old_PROC_GROUP    ) .OR.                 &
     &    (GATHER_PE       /=  old_GATHER_PE     ) .OR.                 &
     &    (GRID_TYPE       /=  old_GRID_TYPE     ) .OR.                 &
     &    (HALO_TYPE       /=  old_HALO_TYPE     ) .OR.                 &
     &    (current_decomp_type  /=  old_DECOMP  )) THEN
!       Different arguments from the last call so we need
!       to calculate a new send/receive map

        IF (GRID_TYPE  ==  fld_type_unknown) THEN
          WRITE(6,*) 'GATHER_FIELD_GCOM : Bad field type'
          WRITE(6,*) 'Field will not be gathered.'
          CMESSAGE='GATHER_FIELD_GCOM : Bad field type'
          ICODE=1
          GOTO 9999
        ENDIF


! 2.0 Set up send map

        send_map(S_DESTINATION_PE,1) = GATHER_PE

        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =                      &
     &    halosize(2,halo_type)*LOCAL_ROW_LEN+                          &
     &    1+halosize(1,halo_type)

        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=blsize(2,grid_type)

        send_map(S_STRIDE_IN_SEND_ARRAY,1) = LOCAL_ROW_LEN

        send_map(S_ELEMENT_LENGTH,1) =                                  &
     &    LOCAL_ROW_LEN-2*halosize(1,halo_type)

        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) =                      &
     &    datastart_f(1,grid_type) +                                    &
     &   (datastart_f(2,grid_type)-1)*GLOBAL_ROW_LEN

        send_map(S_STRIDE_IN_RECV_ARRAY,1) = GLOBAL_ROW_LEN

! 3.0 Set up the receive map (for PE GATHER_PE only)

        n_mess_to_rec=0

        IF (mype  ==  GATHER_PE) THEN
          DO iproc=0,nproc-1
            receive_map(R_SOURCE_PE,iproc+1) = iproc

            receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =         &
     &        g_datastart_f(1,grid_type,iproc)+                         &
     &        (g_datastart_f(2,grid_type,iproc)-1)*GLOBAL_ROW_LEN

            receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) =         &
     &        g_blsize(2,grid_type,iproc)

            receive_map(R_STRIDE_IN_RECV_ARRAY,iproc+1) =               &
     &        GLOBAL_ROW_LEN

            receive_map(R_ELEMENT_LENGTH,iproc+1) =                     &
     &        g_blsize(1,grid_type,iproc)

            receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =         &
     &        halosize(2,halo_type)*                                    &
     &        g_lasize(1,grid_type,halo_type,iproc)+                    &
     &        halosize(1,halo_type)+1

            receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) =               &
     &        g_lasize(1,grid_type,halo_type,iproc)

          ENDDO
          n_mess_to_rec=nproc
        ENDIF

        old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
        old_GLOBAL_ROWS=GLOBAL_ROWS
        old_PROC_GROUP=PROC_GROUP
        old_GATHER_PE=GATHER_PE
        old_DECOMP=current_decomp_type
        old_GRID_TYPE=GRID_TYPE
        old_HALO_TYPE=HALO_TYPE

      ENDIF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

      flag=0  ! This is currently ignored at GCG v1.1
      info=GC_NONE

      CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,1,                       &
     &                    LOCAL_ROW_LEN*LOCAL_ROWS,                     &
     &                    GLOBAL_FIELD,receive_map,n_mess_to_rec,       &
     &                    GLOBAL_ROW_LEN*GLOBAL_ROWS,                   &
     &                    PROC_GROUP,flag,info)

 9999 CONTINUE

      RETURN
      END SUBROUTINE GATHER_FIELD_GCOM
#endif
