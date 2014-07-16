#if defined(C96_1B) 
#if defined(T3E) && defined(MPP)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Takes 1 or more levels of a model field that have been decomposed
!  over a group of processors, and gathers the data together so that
!  one complete global level is contained on one processor.  Processors
!  can hold one or more such levels, as determined by the 'LOCAL_LEVEL'
!  array, which gives the index for each level on the processor to
!  which it is sent.  For one level per PE, the setting of the values
!  LOCAL_LEVEL(1...GLOBAL_LEVELS) should all be one.  Successive ones
!  obviously range from 1 upwards.
!
! Method:
!  This routine uses SHMEM_PUT directly for each row, and staggers the
!  target PE's, based on the identity of the sending PE.
!
! Current Code Owner: Paul Burton
!
! History:
!  Model    Date      Modification history:
!  version
!    4.5    18/09/97  New code optimised for the T3E
!                       Author: P.Burton
!    5.0    14/09/99  - Added HALO_TYPE argument
!                     - Removed unused PROC_GROUP argument
!                     - Removed INFO argument
!                     - Changed to new format PARVARS variables
!    5.3    20/09/01  Adds FLD_TYPE argument and ensures number of
!                     rows sent is correct                P.Burton
!    6.0    05/09/03  Add new defs for use with makebc. R.Sempers
!    6.3    09/10/06  Remove MAKEBC DEF. D Robinson.
!
! Subroutine Interface:
      SUBROUTINE GATHER_FIELD_ML(                                       &
     &    LOCAL_FIELD,GLOBAL_FIELD,                                     &
     &    LOCAL_ROW_LEN,LOCAL_ROWS,LOCAL_LEVS,                          &
     &    GLOBAL_ROW_LEN,GLOBAL_ROWS,GLOBAL_LEVS,                       &
     &    PE_FOR_LEVEL,LOCAL_LEVEL,FLD_TYPE,HALO_TYPE)

      IMPLICIT NONE
!
! Subroutine Arguments:

      INTEGER                                                           &
     &  LOCAL_ROW_LEN                                                   &
                         ! IN length of rows in local part of field
     &, LOCAL_ROWS                                                      &
                         ! IN number of rows in local part of field
     &, LOCAL_LEVS                                                      &
                         ! IN number of levels in local field
     &, GLOBAL_ROW_LEN                                                  &
                         ! IN length of rows in global field
     &, GLOBAL_ROWS                                                     &
                         ! IN number of rows in global field
     &, GLOBAL_LEVS                                                     &
                         ! IN number of levels in global field
     &, PE_FOR_LEVEL(LOCAL_LEVS)                                        &
                         ! IN PE to gather each level to
     &, LOCAL_LEVEL(LOCAL_LEVS)                                         &
                         ! IN level index of level on gather pe
     &, FLD_TYPE                                                        &
                         ! IN field type of grid
     &, HALO_TYPE        ! IN halo type (hence width) of grid

      REAL                                                              &
     &  LOCAL_FIELD(LOCAL_ROW_LEN,LOCAL_ROWS,LOCAL_LEVS)                &
!                        ! IN local part of field
     &, GLOBAL_FIELD(GLOBAL_ROW_LEN,GLOBAL_ROWS,GLOBAL_LEVS)
!                        ! OUT (on PE GATHER_PE) global field

! Parameters and Common blocks

#include "parvars.h"

! Local variables

      INTEGER                                                           &
     &  level            ! loop index for levels

      INTEGER                                                           &
     &  n_rows_to_put                                                   &
                         ! Number of rows to send
     &, j                                                               &
                         ! loop index for rows
     &, start_row                                                       &
                         ! first row to send
     &, end_row                                                         &
                         ! last row to send
     &, real_level       ! The actual level to send - PE's
                         ! traverse the levels using their
                         ! PE numbers as an offset to reduce
                         ! network contention

! array to hold the address of the global fields on each PE
      INTEGER                                                           &
     &  address_global_field(0:MAXPROC)

      COMMON /shmem_align_address/                                      &
     &  address_global_field

! remote global field on other PE's, whose address is set up
! by a CRAY type pointer, after exchanging remote addresses
      REAL                                                              &
     & remote_GLOBAL_FIELD(GLOBAL_ROW_LEN, GLOBAL_ROWS, GLOBAL_LEVS)

      POINTER (PTR_remote_GLOBAL_FIELD, remote_GLOBAL_FIELD)


!-------------------------------------------------------

      n_rows_to_put=blsize(2,FLD_TYPE)
      start_row=halosize(2,halo_type)+1
      end_row=start_row+n_rows_to_put-1

! store the address of the global field I am collecting
! (sending PE's must get this before they send data)
      address_global_field(mype)=LOC(GLOBAL_FIELD)

! ensure everyone has set the address of their global fields
      CALL barrier()

! loop over the number of rows to put
      DO j=start_row,end_row

! loop over the levels to send, using our PE number as an
! offset to reduce network contention and spread the work
! out of different PE's
        DO level=1+mype,LOCAL_LEVS+mype

! compute the real level to send
          real_level=MOD(level-1,LOCAL_LEVS)+1

! first row for this level?  If so, we must find the remote
! address of the global field into which the data is to be sent
          IF (j  ==  start_row) THEN

            CALL shmem_get(                                             &
     &        address_GLOBAL_FIELD(PE_FOR_LEVEL(real_level)),           &
     &        address_GLOBAL_FIELD(PE_FOR_LEVEL(real_level)),           &
     &        1,PE_FOR_LEVEL(real_level))

          ENDIF

! set up the remote address of the global field
          PTR_remote_GLOBAL_FIELD=                                      &
     &      address_GLOBAL_FIELD(PE_FOR_LEVEL(real_level))

! send the data off to the collecting PE for this level
          CALL shmem_put(                                               &
     &      remote_GLOBAL_FIELD(datastart(1),                           &
     &                          datastart(2)+j-start_row,               &
     &                          LOCAL_LEVEL(real_level)),               &
     &      LOCAL_FIELD(halosize(1,halo_type)+1,j,real_level),          &
     &      LOCAL_ROW_LEN-2*halosize(1,halo_type),                      &
     &      PE_FOR_LEVEL(real_level))

        ENDDO ! level
      ENDDO ! j

! wait for everyone to finish sending their data
      CALL barrier()

9999  continue


      RETURN
      END SUBROUTINE GATHER_FIELD_ML
#endif
#endif
