#if defined(C96_1C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Direct MPI version of swap_bounds to deal with multiple variables
! at once

SUBROUTINE SWAP_BOUNDS_MV(  &
   input_fields,            &        ! Fields to be swapped
   n_multi,                 &        ! The number of Fields to swap
   row_length,              &        ! field size
   halo_x, halo_y)                   ! halos

USE MPL, ONLY :             &
         MPL_REAL,          &
         MPL_STATUS_SIZE

USE SWAPABLE_FIELD_MOD, ONLY : &
    SWAPABLE_FIELD_POINTER_TYPE

IMPLICIT NONE

!  This code performs the same process as swap_bounds but
!  can swap multiple variables at once and hence utilise bandwidth
!  better. This is the MPL version.
!
!  Note that it can only deal with multiple variables with the same
!  row_length and same halo size.
!
! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. Data is swapped
!   across the poles for any non-zero halo size in the y direction
!   if it is a global model.
!
! Implementation
!   The logic flow is non-trivial!
!   The across pole differencing in particular must be handled carefully!
!   The basic idea is to copy the data to be transferred into a
!   send_buffer array, which is sent to the receiving processor
!   where it arrives in receive_buffer. The receiving processor
!   then copies this data into the appropriate halo region.
!   The East/West halos are done first, then the North/South
!   halos.
!
!   Note that due to the fact that pointers to the data are used,
!   addressing is (1:row_length+2*halo_x, 1:rows+2*halo_y) rather
!   than (1-halo_x:row_length+halo_x, 1-halo_y:rows+halo_y) as used
!   elsewhere. This is unavoidable as pointers change the addressing
!   mode.
!
! Author: Paul Selwood
! Current code owner: Paul Selwood
!
! Comdecks
#include "domtyp.h"
#include "parvars.h"

! Arguments:

INTEGER, INTENT(IN)  :: n_multi      ! number of fields to be swapped
INTEGER, INTENT(IN)  :: row_length   ! number of points on row (no halos)
INTEGER, INTENT(IN)  :: halo_x       ! sixze of "i" halo
INTEGER, INTENT(IN)  :: halo_y       ! sixze of "j" halo

! Fields to swap
TYPE(swapable_field_pointer_type), TARGET :: input_fields(n_multi)

! Local scalar variables
INTEGER :: levels             ! number of levels in field
INTEGER :: rows               ! number of rows in field (no halos)
INTEGER :: field_type    ! The grid type of the field (u,v,p)
INTEGER :: i,j,k         ! Spatial loop counters
INTEGER :: info          ! GCOM return code
INTEGER :: length
INTEGER :: i_field
INTEGER :: ierror        ! MPI return code
INTEGER :: nreq_r_ew     ! number of MPI recv requests ew
INTEGER :: nreq_s_ew     ! number of MPI send requests ew
INTEGER :: nreq_r_ns     ! number of MPI recv requests ns
INTEGER :: nreq_s_ns     ! number of MPI send requests ns
INTEGER :: full_row_length       ! length of row including halos
INTEGER :: max_full_rows         ! max no of rows including halo rows
INTEGER :: EW_Halo_Size          ! size of EW halo region
INTEGER :: NS_Halo_Size          ! size of NS halo region
INTEGER :: west_halo_source      ! first column of data for west halo
INTEGER :: east_halo_source      ! first column of data for east halo
INTEGER :: south_halo_source
INTEGER :: north_halo_source
INTEGER :: half_full_row_length  ! half of the total EW dimension
INTEGER :: half_row_length       ! half of the data (no halo) EW dim
INTEGER :: index_2_start         ! start address of index2 in buffers
INTEGER :: max_levels            ! maximum levels used over all fields
INTEGER :: max_rows              ! maximum rows used over all fields
INTEGER :: buffer_size           ! size of send/recv buffers
INTEGER :: MY_COMM               ! Communicator

LOGICAL :: l_vector      ! TRUE if a vector field


! Local arrays
INTEGER :: istat(MPL_STATUS_SIZE,4)  ! MPI status
INTEGER :: ireq_r_ew(4)              ! MPI requests
INTEGER :: ireq_s_ew(4)              ! MPI requests
INTEGER :: ireq_r_ns(4)              ! MPI requests
INTEGER :: ireq_s_ns(4)              ! MPI requests
INTEGER :: north_off(n_multi)        ! Offsets to use when copying data
INTEGER :: south_off(n_multi)        ! to send around poles

REAL, POINTER :: field(:, :, :)

! Send and receive buffers to be allocated dynamically
REAL, ALLOCATABLE :: send_buffer(:)
REAL, ALLOCATABLE :: receive_buffer(:)

LOGICAL :: change_sign(n_multi)     ! .TRUE. if sign change across pole


! Statement functions for addressing the buffer arrays
INTEGER :: EW_address
INTEGER :: NS_address
 
! Variables used in statement functions

INTEGER :: row
INTEGER :: point
INTEGER :: halo
INTEGER :: level
INTEGER :: index
INTEGER :: j_field

EW_address(row,halo,level,index,j_field)=                                 &
   n_multi*(index-1)*(max_levels*EW_Halo_Size) +                          &
   (j_field-1)*(max_levels*EW_Halo_Size) +                                &
   (level-1)*EW_Halo_Size +                                               &
   (halo-1)*max_full_rows +                                               &
   row

NS_address(point,halo,level,index,j_field)=                               &
   n_multi*(index-1)*(max_levels*NS_Halo_Size) +                          &
   (j_field-1)*(max_levels*NS_Halo_Size) +                                &
   (level-1)*NS_Halo_Size +                                               &
   (halo-1)*full_row_length +                                             &
   point


!------------------------------------------------------------------
! 0.0 Check if there is anything to do
IF (((halo_x == 0) .AND. (halo_y == 0)) .OR. n_multi == 0) GOTO 9999

!------------------------------------------------------------------
! 1.0 Initialise variables

! Maximum rows and levels
max_rows   = 0
max_levels = 0
DO i_field=1, n_multi
  max_levels= max(input_fields(i_field) % levels, max_levels)
  max_rows  = max(input_fields(i_field) % rows,   max_rows)
END DO

full_row_length = row_length + 2*halo_x
max_full_rows   = max_rows + 2*halo_y

half_full_row_length = full_row_length/2
half_row_length      = row_length/2

EW_Halo_Size         = max_full_rows   * halo_x
NS_Halo_Size         = full_row_length * halo_y
 
 
! Allocate buffers for communication
buffer_size =  n_multi * 2*max_levels *                                   &
               (max(max_rows*halo_x, row_length*halo_y) +                 &
               (2 * halo_x * halo_y))

ALLOCATE( send_buffer(buffer_size) )
ALLOCATE( receive_buffer(buffer_size) )

! Get communicator we will be using from GCOM
CALL GC_GET_COMMUNICATOR(MY_COMM,IERROR)

!------------------------------------------------------------------
! 2.0 East-West communications
!---------------------------------------
! 2.1 Simple case of only one processor
!     in the East-West direction

IF (halo_x > 0) THEN

nreq_s_ew = 0
nreq_s_ns = 0
IF (nproc_x == 1) THEN ! only 1 processor East-West

  IF (bound(1) == BC_CYCLIC) THEN        ! cyclic boundary conditions
    west_halo_source=row_length+1        ! copy from opposite end
    east_halo_source=halo_x+1            ! of each row

    DO i_field=1, n_multi
      field       => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows
!CDIR NOVECTOR
      DO I=1,halo_x
        DO K=1,levels
          DO J=1,rows + (2 * halo_y)

            ! Fill Western halo
            field(I,J,K) = field(west_halo_source+I-1,J,K)

            ! Fill Eastern halo
            field(row_length+halo_x+I,J,K) =                              &
                            field(east_halo_source+I-1,J,K)

          END DO ! J
        END DO ! K
      END DO ! I

    END DO ! loop over fields

  END IF   !  bound(1) == BC_CYCLIC

!---------------------------------------
! 2.1 Now the more common case of having
!     a number of processors in the
!     East-West direction

ELSE ! If there is more than 1 processor East-West

!---------------------------------------
! 2.1.1 Copy the data into send_buffer

  DO i_field=1, n_multi
    field       => input_fields(i_field) % field
    levels      = input_fields(i_field) % levels
    rows        = input_fields(i_field) % rows
  
    DO K=1,levels
      DO J=1,rows + (2*halo_y)
        DO I=1,halo_x

          ! Copy stuff from the Western side of the grid
          send_buffer(EW_address(J,I,K,1,i_field)) =                      &
                      field(halo_x+I,J,K)

          ! Copy stuff from the Eastern side of the grid
          send_buffer(EW_address(J,I,K,2,i_field)) =                      &
                      field(row_length+I,J,K)

        END DO ! I
      END DO ! J
    END DO ! K

  END DO ! loop over fields

!---------------------------------------
! 2.1.2 Send and receive the data

!---------------------------------------
! 2.1.2.1 Special case of 2 processors
!         East-West - both sides are
!         sent to the same processor
!         (if cyclic BC)

  IF ((nproc_x == 2) .AND. (bound(1) == BC_CYCLIC)) THEN

    length = 2 * n_multi * EW_halo_size * max_levels

    CALL GC_RSEND(1,length,neighbour(PEast),info,                         &
                   receive_buffer, send_buffer)

    CALL GC_RRECV(1,length,neighbour(PWest),info,                         &
                   receive_buffer, send_buffer)

!---------------------------------------
! 2.1.2.2 More general case when there
!         are more than 2 processors
!         in the EW direction. Each
!         halo can be sent seperately
!         as there is no danger of them
!         both being sent to the same
!         processor.

  ELSE ! more than 2 processors East-West

    ! index_2_start points to the start of the second index
    ! within the buffer arrays. The first index contains
    ! data sent from the Western side, the second index
    ! contains data sent from the Eastern side

    index_2_start=EW_address(1,1,1,2,1)

    length = n_multi * EW_halo_size * max_levels

    nreq_r_ew=0
    IF (neighbour(PEast) /= NoDomain) THEN

      ! Receive from East
      nreq_r_ew=nreq_r_ew+1
      CALL MPL_IRECV(receive_buffer,                                      &
                     length, MPL_REAL, neighbour(PEast), 2,               &
                     MY_COMM, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    IF (neighbour(PWest) /= NoDomain) THEN

      ! Receive from West
      nreq_r_ew=nreq_r_ew+1
      CALL MPL_IRECV(receive_buffer(index_2_start),                       &
                     length, MPL_REAL, neighbour(PWest), 3,               &
                     MY_COMM, ireq_r_ew(nreq_r_ew), ierror)

    END IF

    nreq_s_ew=0
    IF (neighbour(PEast) /= NoDomain) THEN
      
      ! Send East
      nreq_s_ew=nreq_s_ew+1
      CALL MPL_ISEND(send_buffer(index_2_start),                          &
                     length, MPL_REAL, neighbour(PEast), 3,               &
                     MY_COMM, ireq_s_ew(nreq_s_ew), ierror)
    END IF

    IF (neighbour(PWest) /= NoDomain) THEN
    
      ! Send West
      nreq_s_ew=nreq_s_ew+1
      CALL MPL_ISEND(send_buffer,                                         &
                     length, MPL_REAL, neighbour(PWest), 2,               &
                     MY_COMM, ireq_s_ew(nreq_s_ew), ierror)
    
    END IF
      
    CALL MPL_WAITALL( nreq_r_ew, ireq_r_ew, istat, ierror )

  END IF ! test on numbers of processors East-West

  CALL MPL_WAITALL( nreq_s_ew, ireq_s_ew, istat, ierror )

!---------------------------------------
! 2.1.2 Fill the halos with data

  IF (neighbour(PEast) /= NoDomain) THEN

    ! unpack data from receive_buffer into field

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows 

      DO K=1,levels
        DO J=1,rows + (2*halo_y)
          DO I=1,halo_x
            field(row_length+halo_x+I,J,K)=                               &
               receive_buffer(EW_address(J,I,K,1,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields

  ELSEIF (sb_Model_domain == mt_global) THEN
       ! No neighbour to my East (ie. at edge of the domain
       ! and it's not a cyclic boundary condition)
       ! Just copy data from last column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows 

      DO K=1,levels
        DO J=1,rows+ (2*halo_y)
          DO I=1,halo_x
            field(row_length+halo_x+I,J,K)=field(row_length+halo_x,J,K)
          END DO
        END DO
      END DO
    END DO ! loop over fields

  END IF ! IF (neighbour(PEast) /= NoDomain)

  IF (neighbour(PWest) /= NoDomain) THEN

    ! unpack data from receive_buffer into field

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO K=1,levels
        DO J=1,rows+ (2*halo_y)
          DO I=1,halo_x
            field(I,J,K)=                                                 &
               receive_buffer(EW_address(J,I,K,2,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields

  ELSEIF (sb_Model_domain == mt_global) THEN
       ! No neighbour to my West (ie. at edge of the domain
       ! and it's not a cyclic boundary condition)
       ! Just copy data from first column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO K=1,levels
        DO J=1,rows+ (2*halo_y)
          DO I=1,halo_x
            field(I,J,K)=field(1+halo_x,J,K)
          END DO
        END DO
      END DO
    END DO ! loop over fields
 
  END IF ! IF (neighbour(PWest) /= NoDomain)

END IF ! IF (nproc_x == 1)

END IF ! halo_x > 0

!------------------------------------------------------------------
! 3.0 North-South communications
!  section of code for cyclic N-S boundary conditions
IF (halo_y > 0) THEN

IF(bound(2) == BC_CYCLIC) THEN
  
  IF (nproc_y == 1) THEN ! only 1 processor north-south
  ! 2.1 Simple case of only one processor
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      south_halo_source=rows+1              ! copy from opposite end
      north_halo_source=1                   ! of each column

      DO K=1,levels
        DO J=1,halo_y
          DO I=1 ,row_length+ (2*halo_x)
  
            ! Fill southern halo
            field(I,J,K)=field(I,south_halo_source+J,K)
  
            ! Fill northern halo
            field(I,rows+halo_y+J,K)=field(I,north_halo_source+J,K)
  
          END DO ! I
        END DO ! J
      END DO ! K
    END DO ! loop over fields
  
  !---------------------------------------
  ! 2.1 Now the more common case of having
  !     a number of processors in the
  !     North-South direction
  
  ELSE ! If there is more than 1 processor north_south
  
  !---------------------------------------
  ! 2.1.1 Copy the data into buf_send
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows

      DO K=1,levels
        DO J=1,halo_y
          DO I=1, row_length + (2*halo_x)
  
            ! Copy stuff from the southern side of the grid
            send_buffer(NS_address(I,J,K,1,i_field)) =                    &
               field(I,halo_y+J,K)
  
            ! Copy stuff from the northern side of the grid
            send_buffer(NS_address(I,J,K,2,i_field)) =                    &
               field(I,rows+J,K)
  
          END DO ! I
        END DO ! J
      END DO ! K
    END DO ! loop over fields
  
  !---------------------------------------
  ! 2.1.2 Send and receive the data
  
  !---------------------------------------
  ! 2.1.2.1 Special case of 2 processors
  !         north-south - both sides are
  !         sent to the same processor
  !         as cyclic BC
  
   IF ( nproc_y == 2 ) THEN
  
     length=2*n_multi*NS_halo_size*max_levels
  
      CALL GC_RSEND(1,length,neighbour(PNorth),info,                      &
                     receive_buffer, send_buffer)
  
      CALL GC_RRECV(1,length,neighbour(PSouth),info,                      &
                     receive_buffer, send_buffer)
  
  !---------------------------------------
  ! 2.1.2.2 More general case when there
  !         are more than 2 processors
  !         in the NS direction. Each
  !         halo can be sent seperately
  !         as there is no danger of them
  !         both being sent to the same
  !         processor.
  
    ELSE ! more than 2 processors North-South
  
      ! index_2_start points to the start of the second index
      ! within the buffer arrays. The first index contains
      ! data sent from the southern side, the second index
      ! contains data sent from the northern side
  
      index_2_start=NS_address(1,1,1,2,1)
  
      length = n_multi * ns_halo_size * max_levels
  
      IF (neighbour(PSouth) /= NoDomain) THEN
  
        ! Send south
        CALL GC_RSEND(2,length,neighbour(Psouth),info,                    &
                       receive_buffer, send_buffer)
      END IF
  
      IF (neighbour(Pnorth) /= NoDomain) THEN
  
        ! Send north
        CALL GC_RSEND(3,length,neighbour(Pnorth),info,                    &
                       receive_buffer(index_2_start),                     &
                       send_buffer(index_2_start))
      END IF
  
      IF (neighbour(Pnorth) /= NoDomain) THEN
  
        ! Receive from north
        CALL GC_RRECV(2,length,neighbour(Pnorth),info,                    &
                       receive_buffer, send_buffer)
      END IF
  
      IF (neighbour(Psouth) /= NoDomain) THEN
  
        ! Receive from south
        CALL GC_RRECV(3,length,neighbour(Psouth),info,                    &
                       receive_buffer(index_2_start),                     &
                       send_buffer(index_2_start))
  
      END IF
  
    END IF ! test on numbers of processors north-south
  
  !---------------------------------------
  ! 2.1.2 Fill the halos with data
  
    IF (neighbour(Pnorth) /= NoDomain) THEN
  
      ! unpack data from receive_buffer into field
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows

        DO K=1,levels
          DO J=1,halo_y
            DO I=1,row_length+ (2*halo_x)
              field(I,J+halo_y+rows,K)=                                   &
              receive_buffer(NS_address(I,J,K,1,i_field))
                             
            END DO
          END DO
        END DO
      END DO ! loop over fields
  
    END IF ! IF (neighbour(Pnorth) /= NoDomain)
  
    IF (neighbour(Psouth) /= NoDomain) THEN
  
      ! unpack data from receive_buffer into field
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        DO K=1,levels
          DO J=1,halo_y
            DO I=1,row_length+ (2*halo_x)
              field(I,J,K)=                                               &
              receive_buffer(NS_address(I,J,K,2,i_field))
            END DO
          END DO
        END DO
      END DO ! loop over fields
  
  
    END IF ! IF (neighbour(Psouth) /= NoDomain)
  
  END IF ! IF (nproc_y == 1)
  
ELSE                 !!! bc_cyclic in NS

  ! Set up some variables
  
  ! Set up the offsets. When copying data that is to be passed over
  ! the pole, on wind (u or v) grid, then copy data one row away
  ! from the pole
  
  DO i_field=1, n_multi
    north_off(i_field)=0
    south_off(i_field)=0
    field_type = input_fields(i_field) % field_type
    IF ((sb_Model_domain == mt_global) .AND.                              &
            ((field_type == fld_type_p) .OR.                              &
             (field_type == fld_type_u))) THEN
  
      IF (at_extremity(PNorth)) north_off(i_field)=1
      IF (at_extremity(PSouth)) south_off(i_field)=1
   
    END IF
  
  ! Set up the sign factor. If l_vector is true and data has been passed
  ! over the poles in a global model, then the variables must change
  ! sign
  
    l_vector = input_fields(i_field) % vector
    IF (.NOT. ((l_vector) .AND. (sb_Model_domain == mt_global))) THEN
      change_sign(i_field)=.FALSE.
    ELSE
      change_sign(i_field)=.TRUE.
    END IF
  END DO
  
  !---------------------------------------
  ! 3.1 Copy data into the send_buffer
  !     But not if:
  !       - Not a global model and the data is at the North/South edge
  !       - A global model at the North/South edge but only 1 processor
  !         in the East-West direction.
  
  IF (.NOT. (at_extremity(PSouth) .AND.                                   &
              ((nproc_x == 1)  .OR.                                       &
              (sb_Model_domain /= mt_global))))                           &
     THEN
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+ (2*halo_x)
  
            send_buffer(NS_address(I,J,K,1,i_field)) =                    &
               field(I,J+halo_y+south_off(i_field),K)
  
          END DO
        END DO
      END DO
    END DO ! loop over fields
  END IF
  
  IF (.NOT. (at_extremity(PNorth) .AND.                                   &
              ((nproc_x == 1) .OR.                                        &
              (sb_Model_domain /= mt_global))))                           &
     THEN
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+ (2*halo_x)
  
            send_buffer(NS_address(I,J,K,2,i_field)) =                    &
               field(I,rows-north_off(i_field)+J,K)
  
          END DO
        END DO
      END DO
    END DO ! loop over fields
  END IF
  
  !---------------------------------------
  ! 3.2 Send and receive the data
  
  !---------------------------------------
  ! 3.2.1 The special case where nproc_y=1
  !       Both buffers are sent to the
  !       same processor as one message
  
  IF ((nproc_y == 1) .AND. (sb_Model_domain == mt_global) .AND.           &
       (nproc_x > 1)) THEN
  
    length = 2 * n_multi * ns_halo_size * max_levels
  
    CALL GC_RSEND(10,length,neighbour(PSouth),info,                       &
                     receive_buffer, send_buffer)
  
    CALL GC_RRECV(10,length,neighbour(PSouth),info,                       &
                   receive_buffer, send_buffer)
  END IF
  
  !---------------------------------------
  ! 3.2.2 The more general case, where
  !       each buffer is sent to a
  !       different processor
  
  index_2_start=NS_address(1,1,1,2,1)
  
  nreq_r_ns=0
  
  length = n_multi * ns_halo_size * max_levels
  
  IF (at_extremity(PSouth)) THEN
  
    IF ((neighbour(PSouth) /= NoDomain) .AND.                             &
         (neighbour(PSouth) /= mype)) THEN
  
       nreq_r_ns=nreq_r_ns+1
       CALL MPL_IRECV(receive_buffer(index_2_start),                      &
         length, MPL_REAL, neighbour(PSouth), 11, MY_COMM,                &
         ireq_r_ns(nreq_r_ns), ierror)
  
    END IF
  
  ELSE ! not at the South
  
     nreq_r_ns=nreq_r_ns+1
     CALL MPL_IRECV(receive_buffer(index_2_start),                        &
       length, MPL_REAL, neighbour(PSouth), 14, MY_COMM,                  &
       ireq_r_ns(nreq_r_ns), ierror)
  
  END IF
  
  IF (at_extremity(PNorth)) THEN
  
    IF ((neighbour(PNorth) /= NoDomain) .AND.                             &
         (neighbour(PNorth) /= mype)) THEN
  
       nreq_r_ns=nreq_r_ns+1
       CALL MPL_IRECV(receive_buffer,                                     &
         length, MPL_REAL, neighbour(PNorth), 13, MY_COMM,                &
         ireq_r_ns(nreq_r_ns), ierror)
  
    END IF
  
  ELSE
  
       nreq_r_ns=nreq_r_ns+1
       CALL MPL_IRECV(receive_buffer,                                     &
         length, MPL_REAL, neighbour(PNorth), 12, MY_COMM,                &
         ireq_r_ns(nreq_r_ns), ierror)
  
  END IF
  
  nreq_s_ns=0
  IF (at_extremity(PSouth)) THEN
  
    IF ((neighbour(PSouth) /= NoDomain) .AND.                             &
         (neighbour(PSouth) /= mype)) THEN
  
       nreq_s_ns=nreq_s_ns+1
       CALL MPL_ISEND(send_buffer,                                        &
         length, MPL_REAL, neighbour(PSouth), 11, MY_COMM,                &
         ireq_s_ns(nreq_s_ns),ierror)
  
    END IF
  
  ELSE ! not at the South
  
       nreq_s_ns=nreq_s_ns+1
       CALL MPL_ISEND(send_buffer,                                        &
         length, MPL_REAL, neighbour(PSouth), 12, MY_COMM,                &
         ireq_s_ns(nreq_s_ns),ierror)
  
  END IF
  
  IF (at_extremity(PNorth)) THEN
  
    IF ((neighbour(PNorth) /= NoDomain) .AND.                             &
         (neighbour(PNorth) /= mype)) THEN
  
     nreq_s_ns=nreq_s_ns+1
     CALL MPL_ISEND(send_buffer(index_2_start),                           &
       length, MPL_REAL, neighbour(PNorth), 13, MY_COMM,                  &
       ireq_s_ns(nreq_s_ns),ierror)
  
    END IF
  
  ELSE ! not at the North
  
     nreq_s_ns=nreq_s_ns+1
     CALL MPL_ISEND(send_buffer(index_2_start),                           &
       length, MPL_REAL, neighbour(PNorth), 14, MY_COMM,                  &
       ireq_s_ns(nreq_s_ns),ierror)
  
  END IF
  
  CALL MPL_WAITALL( nreq_r_ns, ireq_r_ns, istat, ierror )
  !---------------------------------------
  ! 3.3 Fill the halos with data
  
  !---------------------------------------
  ! 3.3.1 Southern halo
  
  IF (at_extremity(PSouth)) THEN
  
    IF (neighbour(PSouth) == NoDomain) THEN
      IF (sb_Model_domain /= mt_lam) THEN
  
      ! Just copy adjacent rows into halo area
  ! NOTE: This block of code should never be executed. It's being left
  !       in so that it can be used in the future by changing the logic
  !       in the IF statement above.
  
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,J,K)=field(I,halo_y+1,K)
              END DO
            END DO
          END DO
        END DO ! loop over fields
      END IF ! IF (sb_Model_domain /= mt_lam)
  
    ELSEIF (neighbour(PSouth) == mype) THEN
    ! Local across pole difference
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
                field(half_row_length+halo_x+I,1-J+halo_y,K)=             &
                   -field(I+halo_x,J+halo_y+south_off(i_field),K)
                field(I,1-J+halo_y,K)=                                    &
                   -field(half_row_length+I,J+halo_y+                     &
                   south_off(i_field),K)
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
                field(half_row_length+halo_x+I,1-J+halo_y,K)=             &
                   field(I+halo_x,J+halo_y+south_off(i_field),K)
                field(I,1-J+halo_y,K)=                                    &
                   field(half_row_length+I,J+halo_y+                      &
                   south_off(i_field),K)
              END DO
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields
  
    ELSE ! data is receive_buffer(index 2)
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,1-J+halo_y,K)=                                    &
                -receive_buffer(NS_address(I,J,K,2,i_field))
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,1-J+halo_y,K)=                                    &
                 receive_buffer(NS_address(I,J,K,2,i_field))
              END DO
            END DO
          END DO
        END IF !  IF (change_sign(i_field))
      END DO ! loop over fields
  
    END IF ! What type of South extremity
  
  ELSE ! IF (at_extremity(PSouth)
  
    ! not at a South extremity
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+(2*halo_x)
            field(I,J,K)=                                                 &
               receive_buffer(NS_address(I,J,K,2,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields
  
  END IF ! IF (at_extremity(PSouth)
  
  
  !---------------------------------------
  ! 3.3.2 Northern halo
  
  IF (at_extremity(PNorth)) THEN
  
    IF (neighbour(PNorth) == NoDomain) THEN
      IF (sb_Model_domain /= mt_lam) THEN
      ! Just copy adjacent rows into halo area
  ! NOTE: This block of code should never be executed. It's being left
  !       in so that it can be used in the future by changing the logic
  !       in the IF statement above.
  
        DO i_field=1, n_multi
          field      => input_fields(i_field) % field
          levels      = input_fields(i_field) % levels
          rows        = input_fields(i_field) % rows
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,rows+halo_y+J,K)=field(I,rows+halo_y,K)
              END DO
            END DO
          END DO
        END DO ! loop over fields
      END IF
  
    ELSEIF (neighbour(PNorth) == mype) THEN
    ! Local across pole difference
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
  
                field(half_row_length+halo_x+I,rows+halo_y+J,K)=          &
                   -field(I+halo_x,rows+halo_y-J+1-north_off(i_field),K)
                field(I,rows+J+halo_y,K)=                                 &
                   -field(half_row_length+I,                              &
                        rows-J+halo_y+1-north_off(i_field),K)
  
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,half_full_row_length
   
                field(half_row_length+halo_x+I,rows+halo_y+J,K)=          &
                   field(I+halo_x,rows+halo_y-J+1-north_off(i_field),K)
                field(I,rows+J+halo_y,K)=                                 &
                   field(half_row_length+I,                               &
                       rows-J+1+halo_y-north_off(i_field),K)
  
              END DO
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields
  
    ELSE ! data is receive_buffer(index 1)
  
      DO i_field=1, n_multi
        field      => input_fields(i_field) % field
        levels      = input_fields(i_field) % levels
        rows        = input_fields(i_field) % rows
        IF (change_sign(i_field)) THEN
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,rows+halo_y+J,K)=                                 &
                   -receive_buffer(                                       &
                   NS_address(I,halo_y-J+1,K,1,i_field))
              END DO
            END DO
          END DO
        ELSE ! don't change sign
          DO K=1,levels
            DO J=1,halo_y
              DO I=1,row_length+(2*halo_x)
                field(I,rows+halo_y+J,K)=                                 &
                   receive_buffer(                                        &
                   NS_address(I,halo_y-J+1,K,1,i_field))
              END DO
            END DO
          END DO
        END IF ! IF (change_sign(i_field))
      END DO ! loop over fields
  
    END IF ! What type of North extremity
  
  ELSE ! IF (at_extremity(PNorth)
  
    ! not at a North extremity
  
    DO i_field=1, n_multi
      field      => input_fields(i_field) % field
      levels      = input_fields(i_field) % levels
      rows        = input_fields(i_field) % rows
      DO K=1,levels
        DO J=1,halo_y
          DO I=1,row_length+ (2*halo_x)
            field(I,rows+halo_y+J,K)=                                     &
               receive_buffer(NS_address(I,J,K,1,i_field))
          END DO
        END DO
      END DO
    END DO ! loop over fields
  END IF ! IF (at_extremity(PNorth)
END IF ! BC_CYCLIC  for NS
  
CALL MPL_WAITALL( nreq_s_ns, ireq_s_ns, istat, ierror )

END IF ! halo_y > 0
  
9999 CONTINUE

IF (ALLOCATED(send_buffer))    DEALLOCATE(send_buffer)
IF (ALLOCATED(receive_buffer)) DEALLOCATE(receive_buffer)

RETURN
END SUBROUTINE SWAP_BOUNDS_MV
#endif
