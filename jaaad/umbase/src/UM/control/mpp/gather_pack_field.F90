#if defined(C96_1A) || defined(C96_1B) || defined(C96_1C) \
 || defined(UTILIO) || defined(FLDIO) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Gathers a field from many processors to one processor and packs
!
! Subroutine Interface:
        SUBROUTINE GATHER_PACK_FIELD(                                   &
                                LOCAL_FIELD,GLOBAL_FIELD,               &
                                LOCAL_ROW_LEN,LOCAL_ROWS,               &
                                GLOBAL_ROW_LEN,GLOBAL_ROWS,             &
                                GRID_TYPE,HALO_TYPE,                    &
                                GLOBAL_GATHER_PE,PROC_GROUP,            &
                                PACKING, IM_IDENT, LRLE, PACKING_TYPE,  &
                                NUM_OUT,                                &
                                COMP_ACCRCY, RMDI)

      IMPLICIT NONE

!
! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field. Optionally WGDOS packs the
!  gathered data as it is sent to get parallelism of the packing
!  process.
!
! Method:
!  In the most simple situation,
!  a send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE.
!
!  In the more ususal version, there is a 2 stage gather. The
!  first gathers a "block_factor" number of rows together. These
!  are packed (optionally) and then a final gather gets the full
!  field into place. This will require some adjustment to ensure
!  the packed data is correct.
!
!
! Current Code Owner: Paul Selwood
!
! Subroutine Arguments:

        INTEGER, INTENT(IN) ::                                          &
        LOCAL_ROW_LEN                                                   &
                         ! IN length of rows in local part of field
      , LOCAL_ROWS                                                      &
                           ! IN number of rows in local part of field
      , GLOBAL_ROW_LEN                                                  &
                         ! IN length of rows in global field
      , GLOBAL_ROWS                                                     &
                           ! IN number of rows in global field
      , GRID_TYPE                                                       &
                           ! IN type (P,U or V) of grid
      , HALO_TYPE                                                       &
                           ! IN halo type (hence width) of grid
      , GLOBAL_GATHER_PE                                                &
                         ! IN processor to gather global field to
      , PROC_GROUP         ! IN group ID of processors involved here

!
! Optional Arguments to handle the COEX packing if necessary
!
      LOGICAL, INTENT(IN) ::                                            &
        PACKING                                                         &
                           ! IN: Set .true. if packing of the input
                           !     field is to be packed
      , LRLE               ! IN: True if Run Length Encoding is required

      INTEGER, INTENT(IN) ::                                            &
        IM_IDENT           ! IN: Internal model identifier

      INTEGER, INTENT(INOUT) ::                                         &
        PACKING_TYPE       ! IN/OUT: This flag is zero on input,
                           !         then stash packing is selected,
                           !         and the routine returns the
                           !         packing flag.
                           !
                           !         If the variable is set to 2 on
                           !         input then 32-bit packing for
                           !         dumpfiles is selected

      INTEGER, INTENT(OUT) ::                                           &
        NUM_OUT            ! OUT: Number of 32-bit IBM words in the
                           !      Packed field for WDGOS packing

      INTEGER, INTENT(IN) ::                                            &
        COMP_ACCRCY        ! IN: Packing Accuracy in Power of 2

      REAL, INTENT(IN) ::                                               &
        RMDI               ! IN: Missing data indicator
!
! Remaining Non-Optional Arguments
!
      REAL, INTENT(IN) ::                                               &
        LOCAL_FIELD(LOCAL_ROW_LEN*LOCAL_ROWS)
                           ! IN local part of field

      REAL, INTENT(OUT) ::                                              &
        GLOBAL_FIELD(GLOBAL_ROW_LEN*GLOBAL_ROWS)
                           ! OUT (on PE GATHER_PE) global field


! Parameters and Common blocks

#include "csubmodl.h"
#include "parvars.h"
#include "gccom.h"

! Local variables

      INTEGER                                                           &
         send_map(7,1)                                                  &
      ,  receive_map(7,MAXPROC)                                         &
      ,  n_mess_to_recv                                                 &
      ,  n_mess_to_send                                                 &
      ,  send_map_2(7,1)                                                &
      ,  receive_map_2(7,MAXPROC)                                       &
      ,  n_mess_to_recv_2                                               &
      ,  n_mess_to_send_2

      INTEGER                                                           &
        old_GLOBAL_ROW_LEN                                              &
                              ! value on last call
      , old_GLOBAL_ROWS                                                 &
                              ! value on last call
      , old_PROC_GROUP                                                  &
                              ! value on last call
      , old_GATHER_PE                                                   &
                              ! value on last call
      , old_DECOMP                                                      &
                              ! value on last call
      , old_GRID_TYPE                                                   &
                              ! value on last call
      , old_HALO_TYPE         ! value on last call

        SAVE send_map,n_mess_to_send,receive_map,n_mess_to_recv,        &
           old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
           old_GATHER_PE,old_DECOMP,                                    &
           old_GRID_TYPE,old_HALO_TYPE
      DATA old_GLOBAL_ROW_LEN,old_GLOBAL_ROWS,old_PROC_GROUP,           &
           old_GATHER_PE,old_DECOMP,                                    &
           old_GRID_TYPE,old_HALO_TYPE                                  &
         / -1234, -1234, -1234, -1234, -1234, -1234, -1234/

      INTEGER                                                           &
        iproc                                                           &
      , jproc                                                           &
      , kproc                                                           &
      , lproc                                                           &
      , info                                                            &
      , flag                                                            &
      , gather_pe                                                       &
                             ! Local gather PE for a group/block of
                             ! rows
      , row_start_pe                                                    &
                             ! First PE in a block - may or may not
                             ! the gather PE
      , data_address                                                    &
      , data_size                                                       &
      , block_pe                                                        &
                             ! The PE holding the current block of
                             ! rows on the GLOBAL_GATHER_PE
      , pes_per_block                                                   & 
                             ! Number of PE's in a block (nproc_x*
                             ! block_factor)
      , n_local_rows                                                    &
                             ! Number of rows per block
      , packed_buffer_size                                              &
                             ! Size of the buffer to hold packed data
      , length_fullwrd                                                  &
                             ! Size of full word on this machine 64-bit
                             ! for Cray
      , local_packing_type   ! Local copy of the packing_type, which
                             ! maybe not present

      INTEGER :: icode       ! error code
 
      INTEGER, ALLOCATABLE :: icomp(:)
                             ! Local array to hold compressed data

! The block_factor is the number of rows that are gathered together
! for packing purposes before the final gather. This is best set to
! 1 where maximum parallelism can be exploited (eg on the NEC SX-6).
! Massively parallel systems may benefit from a higher value.
      INTEGER, PARAMETER   :: block_factor = 1
 
      REAL  :: my_global_buffer(global_row_len*global_rows+2)
 
      CHARACTER (LEN=*), PARAMETER :: ROUTINENAME='GATHER_PACK_FIELD'
      CHARACTER (LEN=80)           :: CMESSAGE

!-------------------------------------------------------

! Compute the local gather PE - head of row PE
      pes_per_block=block_factor*nproc_x
      gather_pe=(mype/pes_per_block)*pes_per_block
      row_start_pe=gather_pe
 
! Check if the global_gather_pe is in the block - if so
! make sure it is the gather_pe
      IF(global_gather_pe >= gather_pe .AND.                            &
         global_gather_pe <= MIN(nproc, gather_pe+pes_per_block)-1) THEN
        gather_pe=global_gather_pe
      END IF
 
! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

      IF ((GLOBAL_ROW_LEN  /=  old_GLOBAL_ROW_LEN) .OR.                 &
          (GLOBAL_ROWS     /=  old_GLOBAL_ROWS   ) .OR.                 &
          (PROC_GROUP      /=  old_PROC_GROUP    ) .OR.                 &
          (GATHER_PE       /=  old_GATHER_PE     ) .OR.                 &
          (GRID_TYPE       /=  old_GRID_TYPE     ) .OR.                 &
          (HALO_TYPE       /=  old_HALO_TYPE     ) .OR.                 &
          (current_decomp_type  /=  old_DECOMP  )) THEN

!       Different arguments from the last call so we need
!       to calculate a new send/receive map
        IF (GRID_TYPE  ==  fld_type_unknown) THEN
          WRITE(6,*) 'GATHER_PACK_FIELD : Bad field type'
          WRITE(6,*) 'Field will not be scattered.'
          CMESSAGE='GATHER_PACK_FIELD : Bad field type'
          ICODE=1

! DEPENDS ON: ereport
          CALL EREPORT(ROUTINENAME,ICODE,CMESSAGE)
        END IF


! 2.0 Set up send map

        send_map(S_DESTINATION_PE,1) = GATHER_PE

        send_map(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =                      &
          halosize(2,halo_type)*LOCAL_ROW_LEN+                          &
          1+halosize(1,halo_type)

        send_map(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=blsize(2,grid_type)

        send_map(S_STRIDE_IN_SEND_ARRAY,1) = LOCAL_ROW_LEN

        send_map(S_ELEMENT_LENGTH,1) =                                  &
          LOCAL_ROW_LEN-2*halosize(1,halo_type)

        send_map(S_BASE_ADDRESS_IN_RECV_ARRAY,1) =                      &
          datastart_f(1,grid_type) +                                    &
         (datastart_f(2,grid_type)-1)*GLOBAL_ROW_LEN

        send_map(S_STRIDE_IN_RECV_ARRAY,1) = GLOBAL_ROW_LEN

        n_mess_to_send=1

! 3.0 Set up the receive map (for PE GATHER_PE only)

        n_mess_to_recv=0

        IF (mype  ==  GATHER_PE) THEN
 
! Loop over PE's in this block
          DO jproc=row_start_pe,                                        &
           min(nproc, row_start_pe+pes_per_block)-1
            iproc=jproc-row_start_pe

            receive_map(R_SOURCE_PE,iproc+1) = jproc

            receive_map(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =         &
                g_datastart_f(1,grid_type,jproc)+                       &
                (g_datastart_f(2,grid_type,jproc)-1)                    &
                 *GLOBAL_ROW_LEN

            receive_map(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) =         &
                g_blsize(2,grid_type,jproc)

            receive_map(R_STRIDE_IN_RECV_ARRAY,iproc+1) =               &
              GLOBAL_ROW_LEN

            receive_map(R_ELEMENT_LENGTH,iproc+1) =                     &
                g_blsize(1,grid_type,jproc)

            receive_map(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =         &
              halosize(2,halo_type)*                                    &
                g_lasize(1,grid_type,halo_type,jproc)+                  &
              halosize(1,halo_type)+1

            receive_map(R_STRIDE_IN_SEND_ARRAY,iproc+1) =               &
                g_lasize(1,grid_type,halo_type,jproc)

            n_mess_to_recv=n_mess_to_recv+1

          END DO
        END IF

        old_GLOBAL_ROW_LEN=GLOBAL_ROW_LEN
        old_GLOBAL_ROWS=GLOBAL_ROWS
        old_PROC_GROUP=PROC_GROUP
        old_GATHER_PE=GATHER_PE
        old_DECOMP=current_decomp_type
        old_GRID_TYPE=GRID_TYPE
        old_HALO_TYPE=HALO_TYPE

      END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

      flag=0  ! This is currently ignored at GCG v1.1
      info=GC_NONE

      local_packing_type = 0
 
! Only the gather PE's need to do anything now
!
!      Run length encoding applies only to unpacked ocean fields, but
!      most ocean fields remain unpacked even when packing profiles
!      5 or 6 are set. Hence selecting, for example, both packing
!      profile 5 and run length encoding makes sense.
 
      IF(PACKING .AND. GLOBAL_ROWS  >=  2) THEN
        ! 1. Climate wgdos packing has been selected for current
        !    file stream via UMUI
        IF(COMP_ACCRCY  <=  -99 .AND. LRLE .AND.                        &
          im_ident  ==  ocean_im )THEN
           ! 2. STASH packing profile for the field is -99
           ! 3. Run Length Encoding has been selected
           ! 4. Submodel is Ocean
           PACKING_TYPE = 4
        ELSE IF(COMP_ACCRCY  >   -99) THEN
           ! 2. STASH packing profile for the field is set.
           PACKING_TYPE = 1
        END IF
      ELSE
         ! 1. Packing may or may not have been selected for current
         !    file stream. This section of code is not executed when
         !    GRIB packing selected.
         IF (LRLE .AND. im_ident  ==  ocean_im )THEN
           ! 2. Run Length Encoding has been selected
           ! 3. Submodel is Ocean
           PACKING_TYPE = 4
         END IF
      END IF

      num_out=0
      local_packing_type = packing_type

! If there is no packing, then use the output buffer
! if I am the gather PE, or if it has been reduced to a single
! stage gather.  Otherwise, use a temporary buffer.
      IF(local_packing_type /= 1 .AND.                                  &
       ((gather_pe == global_gather_pe) .OR.                            &
        (block_factor >= nproc_y))) THEN

        CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,n_mess_to_send,        &
                            LOCAL_ROW_LEN*LOCAL_ROWS,                   &
                            GLOBAL_FIELD,receive_map,n_mess_to_recv,    &
                            GLOBAL_ROW_LEN*GLOBAL_ROWS,                 &
                            PROC_GROUP,flag,icode)
      ELSE

        CALL GCG_RALLTOALLE(LOCAL_FIELD,send_map,                       &
                            n_mess_to_send,                             &
                            LOCAL_ROW_LEN*LOCAL_ROWS,                   &
                            my_global_buffer,receive_map,               &
                            n_mess_to_recv,                             &
                            GLOBAL_ROW_LEN*GLOBAL_ROWS,                 &
                            PROC_GROUP,flag,icode)
      END IF
 
      IF (ICODE /= 0) THEN
        CMESSAGE = ' Error in GCG_RALLTOALLE'

! DEPENDS ON: ereport
        CALL EREPORT(ROUTINENAME,ICODE,CMESSAGE)
      END IF

 
! Now we must check if the packing_type is 1
      IF(local_packing_type == 1) THEN
 
! Only gather_pe's need to do anything
        IF(mype == gather_pe) THEN
 
! Work out how much local data we have
          n_local_rows=0
          DO kproc=1, block_factor
            IF(row_start_pe+(kproc-1)*nproc_x <  nproc) THEN
              n_local_rows=n_local_rows+                                &
               g_blsize(2,grid_type,row_start_pe+(kproc-1)*nproc_x)
            END IF
          END DO
 
! Setup a buffer for the packed data
          packed_buffer_size=n_local_rows*global_row_len+2
          ALLOCATE (icomp(packed_buffer_size))
 
! Pack the data
          icode=0
          length_fullwrd=64
 
! Check if the2 stage gather has been eliminated - if so
! put the data straight into the output array
          IF(block_factor >= nproc_y) THEN
! DEPENDS ON: coex
            CALL coex(                                                  &
               my_global_buffer(                                        &
                         g_datastart_f(1,grid_type,row_start_pe)+       &
                         (g_datastart_f(2,grid_type,row_start_pe)-1)*   &
                                global_row_len),                        &
                                global_row_len*global_rows,             &
                                global_field, packed_buffer_size,       &
                                global_row_len, n_local_rows,           &
                                num_out,                                &
                                comp_accrcy, .true., rmdi,              &
                                length_fullwrd,                         &
                                icode, cmessage)
 
! Still doing 2 stage gather
          ELSE
 
! DEPENDS ON: coex
            CALL coex(                                                  &
               my_global_buffer(g_datastart_f(1,grid_type,row_start_pe)+&
                         (g_datastart_f(2,grid_type,row_start_pe)-1)*   &
                                global_row_len),                        &
                                global_row_len*global_rows,             &
                                icomp, packed_buffer_size,              &
                                global_row_len, n_local_rows,           &
                                num_out,                                &
                                comp_accrcy, .true., rmdi,              &
                                length_fullwrd,                         &
                                icode, cmessage)
          END IF

          IF (ICODE /= 0) THEN
            cmessage='Error in COEX'

! DEPENDS ON: ereport
            CALL EREPORT(ROUTINENAME, ICODE, CMESSAGE)
          END IF


        END IF ! mype == gather_pe
 
      END IF !packing_type is 1
 
      IF(local_packing_type == 1) THEN

        IF(mype == gather_pe) THEN

          IF(block_factor <  nproc_y) THEN

            CALL gc_rsend(1001+mype, (num_out+1)/2, global_gather_pe,   &
             info,                                                      &
             my_global_buffer(                                          &
             g_datastart_f(1,grid_type,row_start_pe)+                   &
            (g_datastart_f(2,grid_type,row_start_pe)-1)*global_row_len),&
             icomp)

          END IF ! are we doing 2 stage gather?

        END IF ! mype == gather_pe

      END IF !  packing_type is 1
 
      IF(local_packing_type == 1) THEN

! Now pack the into the output buffer on the global_gather_pe
        IF(mype == global_gather_pe) THEN

! Check if we are doing 2 stage gather
          IF(block_factor <  nproc_y) THEN
 
! Loop over each processors Contribution
            DO jproc=0, nproc_y-1, block_factor
              block_pe=jproc*nproc_x
 
! Preserve the row_start_pe for this block of data
              iproc=block_pe

! If the global_gather_pe is in the current block, then the block_pe
! needs to be set to global_gather_pe, not the row/block leader
!
              IF(row_start_pe == block_pe) THEN
                block_pe=global_gather_pe
              END IF
 
              n_local_rows=0
              DO kproc=1, block_factor
                IF(iproc+(kproc-1)*nproc_x <  nproc) THEN
                  n_local_rows=n_local_rows+                            &
                   g_blsize(2,grid_type,iproc+(kproc-1)*nproc_x)
                END IF
              END DO
 
              CALL gc_rrecv(1001+block_pe,                              &
               n_local_rows*global_row_len+2,                           &
               block_pe, info,                                          &
               my_global_buffer(                                        &
               g_datastart_f(1,grid_type,iproc)+                        &
               (g_datastart_f(2,grid_type,iproc)-1)*global_row_len),    &
               icomp)

! DEPENDS ON: unite_coex_files
              CALL unite_coex_files(                                    &
                my_global_buffer(g_datastart_f(1,grid_type,iproc)+      &
               (g_datastart_f(2,grid_type,iproc)-1)*global_row_len),    &
                global_field, num_out, iproc)
            END DO

          END IF ! are we doing 2 stage gather?

        END IF ! mype == global_gather_pe

      ELSE

! Normal Gather Operation, without Packing
! 5.0 Set up the second send map

        IF(block_factor <  nproc_y) THEN

          IF(mype == gather_pe .AND. mype /= global_gather_pe) THEN

            send_map_2(S_DESTINATION_PE,1) = global_gather_pe

            send_map_2(S_BASE_ADDRESS_IN_SEND_ARRAY,1) =                &
                  g_datastart_f(1,grid_type,mype)+                      &
                  (g_datastart_f(2,grid_type,mype)-1)*GLOBAL_ROW_LEN

            send_map_2(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=0
            DO kproc=1, block_factor
              IF(row_start_pe+(kproc-1)*nproc_x <  nproc) THEN
                send_map_2(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)=             &
                  send_map_2(S_NUMBER_OF_ELEMENTS_IN_ITEM,1)+           &
                  g_blsize(2,grid_type,mype+(kproc-1)*nproc_x)
              END IF
            END DO

            send_map_2(S_STRIDE_IN_SEND_ARRAY,1) =                      &
                  GLOBAL_ROW_LEN

            send_map_2(S_ELEMENT_LENGTH,1) =                            &
                  GLOBAL_ROW_LEN

            send_map_2(S_BASE_ADDRESS_IN_RECV_ARRAY,1) =                &
                  g_datastart_f(1,grid_type,mype)+                      &
                  (g_datastart_f(2,grid_type,mype)-1)*GLOBAL_ROW_LEN

            send_map_2(S_STRIDE_IN_RECV_ARRAY,1) =                      &
                  GLOBAL_ROW_LEN

            n_mess_to_send_2=1

          ELSE

            n_mess_to_send_2=0

          END IF

! 6.0 Set up the second receive map (for PE GLOBAL_GATHER_PE only)

          n_mess_to_recv_2=0

          IF(mype == global_gather_pe) THEN

            iproc=0
            DO jproc=0, nproc_y-1, block_factor
 
! Compute the block PE for this group of rows, and check it is not
! not the global_gather_pe
              block_pe=jproc*nproc_x
              lproc=block_pe
 
! Check if the global_gather_pe is in the block - if so
! make sure it is the block_pe
              IF(global_gather_pe >= block_pe .AND.                     &
                global_gather_pe <=                                     &
                MIN(nproc, block_pe+pes_per_block)-1) THEN
                 block_pe=global_gather_pe
              END IF
 
              IF (block_pe  /=  global_gather_pe) THEN

                receive_map_2(R_SOURCE_PE,iproc+1) = block_pe

                receive_map_2(R_BASE_ADDRESS_IN_RECV_ARRAY,iproc+1) =   &
                  g_datastart_f(1,grid_type,block_pe)+                  &
                  (g_datastart_f(2,grid_type,block_pe)-1)*GLOBAL_ROW_LEN

                receive_map_2(R_NUMBER_OF_ELEMENTS_IN_ITEM,iproc+1) = 0
                DO kproc=1, block_factor
                  IF(lproc+(kproc-1)*nproc_x <  nproc) THEN
                    receive_map_2(R_NUMBER_OF_ELEMENTS_IN_ITEM,         &
                                  iproc+1) =                            &
                     receive_map_2(R_NUMBER_OF_ELEMENTS_IN_ITEM,        &
                                  iproc+1) +                            &
                     g_blsize(2,grid_type,lproc+(kproc-1)*              &
                              nproc_x)
                  END IF
                END DO

                receive_map_2(R_STRIDE_IN_RECV_ARRAY,iproc+1) =         &
                  GLOBAL_ROW_LEN

                receive_map_2(R_ELEMENT_LENGTH,iproc+1) =               &
                  GLOBAL_ROW_LEN

                receive_map_2(R_BASE_ADDRESS_IN_SEND_ARRAY,iproc+1) =   &
                  g_datastart_f(1,grid_type,block_pe)+                  &
                  (g_datastart_f(2,grid_type,block_pe)-1)*GLOBAL_ROW_LEN

                receive_map_2(R_STRIDE_IN_SEND_ARRAY,iproc+1) =         &
                  GLOBAL_ROW_LEN

                iproc=iproc+1
                n_mess_to_recv_2=n_mess_to_recv_2+1

              END IF ! block_pe  /=  global_gather_pe

            END DO ! gather_pe's

          END IF ! mype == global_gather_pe

! 7.0 Do the exchange of data

          flag=0  ! This is currently ignored at GCG v1.1
          info=GC_NONE

          CALL GCG_RALLTOALLE(my_global_buffer,send_map_2,              &
                              n_mess_to_send_2,                         &
                              GLOBAL_ROW_LEN*GLOBAL_ROWS,               &
                              GLOBAL_FIELD,receive_map_2,               &
                              n_mess_to_recv_2,                         &
                              GLOBAL_ROW_LEN*GLOBAL_ROWS,               &
                              PROC_GROUP,flag,icode)

        END IF ! packing_type is 1

      END IF ! are we doing 2 stage gather
 
! Deallocate the temporary buffer for coex processed data
      IF(ALLOCATED(icomp)) DEALLOCATE(icomp)
 

      RETURN
      END SUBROUTINE GATHER_PACK_FIELD
#endif
