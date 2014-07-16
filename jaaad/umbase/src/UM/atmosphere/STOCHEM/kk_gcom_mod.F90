#if defined(A25_1A)
#if defined(A25_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE KK_GCOM_MOD

! Module to substitue GCG_RALLTOALLE in subroutines stochem.swap
! and stochem.swap2.
! All data is collected in sendbuf, only one SEND/RECEIVE is done
! to other processes.

! KK_RALLTOALLE uses the same parameter list as GCG_RALLTOALLE,
! although not all paramter are used.

! NOTE: This is NOT a general substitution of GCG_RALLTOALLE. It works
!       only in context of stochem.swap and stochem.swap2
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.2  02/03/06   Tidied up code. Placed names of interface blocks on
!                  end statements behind comment for portability.
!                  T. Edwards / M.G. Sanderson

      IMPLICIT  NONE
      PRIVATE
#include "gccom.h"
#include "parvars.h"

! data area

      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nr_send,  nr_recv
      REAL, DIMENSION(:, :), ALLOCATABLE, SAVE  :: sendbuf,  recvbuf

! Interfaces

      INTERFACE KK_GCOM_INI
        MODULE PROCEDURE KK_GCOM_INI
      END INTERFACE  ! KK_GCOM_INI

      INTERFACE KK_RALLTOALLE
        MODULE PROCEDURE KK_RALLTOALLE_1
        MODULE PROCEDURE KK_RALLTOALLE_2
      END INTERFACE  ! KK_RALLTOALLE

      PUBLIC KK_GCOM_INI, KK_RALLTOALLE

      CONTAINS
        SUBROUTINE KK_GCOM_INI(max_num)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' at end of lines for F90 compatability.
!                  M.G. Sanderson
!  6.2  02/03/06   Changed variable maxval to max_num (maxval is a
!                  Fortran function name). M.G. Sanderson
!
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: max_num  ! Max Number of Elements to be
                                         ! sent to another PE

          ALLOCATE(nr_send(0:nproc-1))
          ALLOCATE(nr_recv(0:nproc-1))
          ALLOCATE(sendbuf(max_num,0:nproc-1))
          ALLOCATE(recvbuf(max_num,0:nproc-1))

          RETURN
        END SUBROUTINE KK_GCOM_INI

        SUBROUTINE KK_RALLTOALLE_1                                      &
     &                   (DATA_IN,SEND_MAP,N_ITEMS_SEND,SARR_LEN,       &
     &                    DATA_OUT,RECV_MAP,N_ITEMS_RECV,RARR_LEN,      &
     &                    GID,FLAG,INFO)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' at end of lines for F90 compatability.
!                  M.G. Sanderson
!  6.2  02/03/06   Tidied code. M.G. sanderson

!         One dimensional input/output array

          IMPLICIT NONE

          REAL,DIMENSION(:),INTENT(IN)      :: data_in
          REAL,DIMENSION(:),INTENT(OUT)     :: data_out
          INTEGER,DIMENSION(:,:),INTENT(IN) :: send_map,recv_map
          INTEGER,INTENT(IN)                :: n_items_send
          INTEGER,INTENT(IN)                :: n_items_recv
          INTEGER,INTENT(IN)                :: sarr_len,rarr_len
          INTEGER,INTENT(IN)                :: gid,flag
          INTEGER,INTENT(OUT)               :: info

          INTEGER                           :: i,j,k
          INTEGER                           :: ih,ns,nr
          INTEGER                           :: tag,istat
          CHARACTER(LEN=72)                 :: message

! fill send buffer

          nr_send = 1
          DO k=0,nproc-1
            j = nr_send(k)
            DO i=1,n_items_send
              IF (send_map(s_destination_pe,i) == k) THEN
                ih = send_map(s_base_address_in_send_array,i)
                sendbuf(j,k)   = data_in(ih)
                j = j+1
              END IF
            END DO
            nr_send(k) = j
          END DO
          nr_send = nr_send-1

! get number of receive Values

          nr_recv = 0
          DO k=0,nproc-1
            j = nr_recv(k)
            DO i=1,n_items_recv
              IF (recv_map(r_source_pe,i) == k) THEN
                j = j+1
              END IF
            END DO
            nr_recv(k) = j
          END DO

! exchange data between PEs

          j = 0
          tag = 1234
          DO k=0,nproc-1
             IF (nr_send(k) > 0) THEN
                ns = nr_send(k)
                CALL GC_RSEND(tag,ns,k,istat,recvbuf(1,k),sendbuf(1,k))
                IF (istat /= 0) THEN
                  WRITE(message,'(''GC_RSEND: istat = '',i8)') istat
! DEPENDS ON: ereport
                  CALL EREPORT('KK_RALLTOALLE_1',1,message)
                END IF
             END IF
          END DO
          CALL GC_SSYNC(nproc,istat)
          DO k=0,nproc-1
            IF (nr_recv(k) > 0) THEN
              nr = nr_recv(k)
              CALL GC_RRECV(tag,nr,k,istat,recvbuf(1,k),sendbuf(1,k))
              IF (istat /= 0) THEN
                WRITE(message,'(''GC_RRECV: istat = '',i8)') istat
! DEPENDS ON: ereport
                CALL EREPORT('KK_RALLTOALLE_1',1,message)
              END IF
            END IF
          END DO
          CALL GC_SSYNC(nproc,istat)

! Move data to data_out

          nr_recv = 1
          DO k=0,nproc-1
            j = nr_recv(k)
            DO i=1,n_items_recv
              IF (recv_map(r_source_pe,i) == k) THEN
                ih = recv_map(r_base_address_in_recv_array,i)
                data_out(ih) = recvbuf(j,k)
                j = j + 1
              END IF
            END DO
            nr_recv(k) = j
          END DO

          RETURN
        END SUBROUTINE KK_RALLTOALLE_1

        SUBROUTINE KK_RALLTOALLE_2                                      &
     &                   (data_in,SEND_MAP,N_ITEMS_SEND,SARR_LEN,       &
     &                    data_out,RECV_MAP,N_ITEMS_RECV,RARR_LEN,      &
     &                    GID,FLAG,INFO)
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  16/02/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' at end of lines for F90 compatability.
!                  M.G. Sanderson
!  6.2  02/03/06   Tidied code. M.G. sanderson

! Two dimensional input/output array

          IMPLICIT NONE

          REAL,DIMENSION(:,:),INTENT(IN)      :: data_in
          REAL,DIMENSION(:,:),INTENT(OUT)     :: data_out
          INTEGER,DIMENSION(:,:),INTENT(IN)   :: send_map,recv_map
          INTEGER,INTENT(IN)                  :: n_items_send
          INTEGER,INTENT(IN)                  :: n_items_recv
          INTEGER,INTENT(IN)                  :: sarr_len,rarr_len
          INTEGER,INTENT(IN)                  :: gid,flag
          INTEGER,INTENT(OUT)                 :: info

          INTEGER                             :: i,j,k
          INTEGER                             :: ih,ns,nr
          INTEGER                             :: nbl,tag,istat
          CHARACTER(LEN=72)                   :: message

! fill send buffer

          nbl = SIZE(data_in,1)

          nr_send = 1
          DO k=0,nproc-1
            j = nr_send(k)
            IF (nbl == 4) THEN                ! Special case POS array
              DO i=1,n_items_send
                IF (send_map(s_destination_pe,i) == k) THEN
                  ih = send_map(s_base_address_in_send_array,i)/4+1
                  sendbuf(j,k)   = data_in(1,ih)
                  sendbuf(j+1,k) = data_in(2,ih)
                  sendbuf(j+2,k) = data_in(3,ih)
                  sendbuf(j+3,k) = data_in(4,ih)
                  j = j+4
                END IF
              END DO
            ELSE
              DO i=1,n_items_send
                IF (send_map(s_destination_pe,i) == k) THEN
                  ih = send_map(s_base_address_in_send_array,i)/nbl+1
                  sendbuf(j:j+nbl-1,k) = data_in(1:nbl,ih)
                  j = j + nbl
                END IF
              END DO
            END IF
            nr_send(k) = j
          END DO
          nr_send = nr_send-1

! get number of receive Values

          nr_recv = 0
          DO k=0,nproc-1
            j = nr_recv(k)
            DO i=1,n_items_recv
              IF (recv_map(r_source_pe,i) == k) THEN
                j = j + nbl
              END IF
            END DO
            nr_recv(k) = j
          END DO

! exchange data between PEs

          j = 0
          tag = 2345
          DO k=0,nproc-1
            IF (nr_send(k) > 0) THEN
              ns = nr_send(k)
              CALL GC_RSEND(tag,ns,k,istat,recvbuf(1,k),sendbuf(1,k))
              IF (istat /= 0) THEN
                WRITE(message,'(''GC_RSEND: istat = '',i8)') istat
! DEPENDS ON: ereport
                CALL EREPORT('KK_RALLTOALLE_2',1,message)
              END IF
            END IF
          END DO
          CALL GC_SSYNC(nproc,istat)
          DO k=0,nproc-1
            IF (nr_recv(k) > 0) THEN
              nr = nr_recv(k)
              CALL GC_RRECV(tag,nr,k,istat,recvbuf(1,k),sendbuf(1,k))
              IF (istat /= 0) THEN
                WRITE(message,'(''GC_RRECV: istat = '',i8)') istat
! DEPENDS ON: ereport
                CALL EREPORT('KK_RALLTOALLE_2',1,message)
              END IF
            END IF
          END DO
          CALL GC_SSYNC(nproc,istat)

! Move data to data_out

          nr_recv = 1
          DO k=0,nproc-1
            j = nr_recv(k)
            IF (nbl == 4) THEN
              DO i=1,n_items_recv
                IF (recv_map(r_source_pe,i) == k) THEN
                  ih = recv_map(r_base_address_in_recv_array,i)/4+1
                  data_out(1,ih) = recvbuf(j,k)
                  data_out(2,ih) = recvbuf(j+1,k)
                  data_out(3,ih) = recvbuf(j+2,k)
                  data_out(4,ih) = recvbuf(j+3,k)
                  j = j + 4
                END IF
              END DO
            ELSE
              DO i=1,n_items_recv
                IF (recv_map(r_source_pe,i) == k) THEN
                   ih = recv_map(r_base_address_in_recv_array,i)/nbl+1
                   data_out(1:nbl,ih) = recvbuf(j:j+nbl-1,k)
                   j = j + nbl
                END IF
              END DO
            END IF
            nr_recv(k) = j
          END DO

          RETURN
        END SUBROUTINE KK_RALLTOALLE_2
      END MODULE KK_GCOM_MOD
#endif
#endif
