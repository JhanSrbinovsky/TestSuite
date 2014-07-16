#if defined(C96_1A) || defined(MAKEBC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Portable GCOM based version SWAP_BOUNDS

      SUBROUTINE SWAP_BOUNDS(                                           &

     &  FIELD, ROW_LENGTH, ROWS, LEVELS,                                &
                                          ! field
     &  HALO_X, HALO_Y,                                                 &
                                          ! halos
     &  FIELD_TYPE, L_VECTOR                                            &
                                          ! supporting information
     &  )

      IMPLICIT NONE

!  This code replaces the new dynamics swap_bounds. It performs an
!  identical function, but:
!   - a number of arguments removed, information brought in via
!     COMMON blocks instead
!   - optimised to minimize the communication and synchronisation
!     overhead
!
! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. Data is swapped
!   across the poles for any non-zero halo size in the y direction
!   if it is a global model.
!
! Implementation
!   The logic flow is non-trivial! Special care has to be taken as
!   GCOM (under some message passing libraries) does not allow more
!   than one outstanding message from a given processor to another
!   given processor. Also the across pole differencing must be
!   handled carefully!
!   The basic idea is to copy the data to be transferred into a
!   send_buffer array, which is sent to the receiving processor
!   where it arrives in receive_buffer. The receiving processor
!   then copies this data into the appropriate halo region.
!   The East/West halos are done first, then the North/South
!   halos.
!
! Author: Paul Burton
! Current code owner: Paul Burton
!
! History
! Date       Version    Comment
! ----       -------    -------
! 6/7/99     5.0        New deck created.                P.Burton
! 17/3/00    5.1        Correction to some GCOM calls    P.Burton
! 24/5/00    5.1        Don't fill in "empty" extremity halos for
!                       non Global models which will have LBC data
!                       in these halos.  P.Burton
! 05/01/01   5.2     -  Corrected error in special case of 2 PEs in
!                       East-West - check now done against nproc_x
!                    -  Corrected indexing of buffer in NS exchange
!                    -  Corrected across N-pole swap for 1PE E-W case
!                    -  Corrected comms to do multiple levels
!                    -  Corrected indexing into send_buffer of N data
!                    -  Corrected offset for p/v fields at poles
!                    -  Corrected ordering of poles for across pole swap
!                                                          P.Burton
! 22/08/01   5.3        Immediately exit if no halos  P.Burton
! 22/8/01    5.3        Removed automatic multiplication of all
!                       north/south halo values, and only change
!                       sign if explicitly requested by l_vector
!                                                       P.Burton
! 14/09/01   5.3        Use sb_model_domain variable    P.Burton
! 14/02/02   5.3        add bi_cyclic LAM code          A.Malcolm
! 20/05/02   5.4        Add in UTILIO def for use to build sx on HP
!                                                           E.Leung
! 23/08/00   5.5        Modification for parallelisation of WAM.
!                       Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
! 09/12/03   6.0        Correct logic bug for LAM code with
!                       1 pe East-West.             A. Malcolm
! 14/10/04   6.1        FIX for EW-cyclic LAM           A. Malcolm
! 23/11/05   6.2        Removed all references to the wavemodel.
!                       T.Edwards
! 29/11/05   6.2        Replaced fixed length common buffers by
!                       appropriate length automatics. T.Edwards
! 21/10/05   6.2        Replace GSYNC with SSYNC.       P.Selwood
! 22/08/05   6.2        Fix defs to match exec_xref. P.Selwood
!
!
! Arguments:

      INTEGER                                                           &
     &  ROW_LENGTH                                                      &
                         ! IN: number of points on a row
                         !     (not including halos)
     &, ROWS                                                            &
                         ! IN: number of rows in a theta field
                         !     (not including halos)
     &, LEVELS                                                          &
                         ! IN: number of model levels
     &, HALO_X                                                          &
                         ! IN: size of halo in "i" direction
     &, HALO_Y                                                          &
                         ! IN: size of halo in "j" direction

     &, FIELD_TYPE       ! IN: Defines the grid interpolation type
                         !     of the input FIELD (u,v or w)

      LOGICAL                                                           &
     &  L_VECTOR         ! IN: TRUE:  Data is a horizontal vector
                         !            component
                         !     FALSE: Data is a scalar

      REAL                                                              &
     &  FIELD(1-HALO_X:ROW_LENGTH+HALO_X,                               &
     &        1-HALO_Y:ROWS+HALO_Y,                                     &
     &        LEVELS)    ! IN/OUT : Field to have its halos updated

! Comdecks
#include "domtyp.h"
#include "parvars.h"


! Local variables

! The send and recieve buffers, replacing common blocks
      REAL :: send_buffer(2*LEVELS*(max(ROWS*HALO_X,ROW_LENGTH*HALO_Y)  &
     &  +(2*HALO_X*HALO_Y))),                                           &
     &     receive_buffer(2*LEVELS*(max(ROWS*HALO_X,ROW_LENGTH*HALO_Y)  &
     &  +(2*HALO_X*HALO_Y)))

      INTEGER                                                           &
     &  full_ROW_LENGTH                                                 &
                                 ! length of row including halos
     &, full_ROWS                                                       &
                                 ! number of rows including halo rows
     &, EW_Halo_Size                                                    &
                                 ! size of EW halo region
     &, NS_Halo_Size                                                    &
                                 ! size of NS halo region
     &, west_halo_source                                                &
                                 ! first column of data for west halo
     &, east_halo_source                                                &
                                 ! first column of data for east halo
     &, south_halo_source                                               &
     &, north_halo_source                                               &
     &, half_full_ROW_LENGTH                                            &
                                 ! half of the total EW dimension
     &, half_ROW_LENGTH                                                 &
                                 ! half of the data (no halo) EW dim
     &, index_2_start                                                   &
                                 ! start address of index2 in buffers
     &, north_off                                                       &
                                 ! Offsets to use when copying data
     &, south_off                ! to send around poles

      INTEGER                                                           &
     &  I,J,K                                                           &
                                 ! Spatial loop counters
     &, info                     ! GCOM return code

      LOGICAL                                                           &
     &  change_sign              ! .TRUE. if sign change across pole

! Statement functions for addressing the buffer arrays

      INTEGER                                                           &
     &  EW_address                                                      &
     &, NS_address

! Variables used in statment functions

      INTEGER                                                           &
     &  row,point,halo,level,index

      EW_address(row,halo,level,index)=                                 &
     &  (index-1)*(levels*EW_Halo_Size) +                               &
     &  (level-1)*EW_Halo_Size +                                        &
     &  (halo-1)*full_ROWS +                                            &
     &  row+HALO_Y

      NS_address(point,halo,level,index)=                               &
     &  (index-1)*(levels*NS_Halo_Size) +                               &
     &  (level-1)*NS_Halo_Size +                                        &
     &  (halo-1)*full_ROW_LENGTH +                                      &
     &  point+HALO_X

!------------------------------------------------------------------
! 0.0 Check if there is anything to do

      IF ((HALO_X  ==  0) .AND. (HALO_Y  ==  0)) GOTO 9999

!------------------------------------------------------------------
! 1.0 Initialise variables

      full_ROW_LENGTH = ROW_LENGTH + 2*HALO_X
      full_ROWS       = ROWS + 2*HALO_Y

      half_full_ROW_LENGTH=full_ROW_LENGTH/2
      half_ROW_LENGTH=ROW_LENGTH/2

      EW_Halo_Size=full_ROWS*HALO_X
      NS_Halo_Size=full_ROW_LENGTH*HALO_Y

!------------------------------------------------------------------
! 2.0 East-West communications

!---------------------------------------
! 2.1 Simple case of only one processor

      IF(HALO_X >  0) then

!     in the East-West direction

      IF (nproc_x  ==  1) THEN ! only 1 processor East-West

        IF (bound(1)  ==  BC_CYCLIC) THEN ! cyclic boundary conditions
          west_halo_source=ROW_LENGTH-HALO_X+1 ! copy from opposite end
          east_halo_source=1                   ! of each row

          DO K=1,LEVELS
            DO J=1-HALO_Y,ROWS+HALO_Y
              DO I=1,HALO_X

                ! Fill Western halo
                FIELD(1-HALO_X+I-1,J,K)=FIELD(west_halo_source+I-1,J,K)

                ! Fill Eastern halo
                FIELD(ROW_LENGTH+I,J,K)=FIELD(east_halo_source+I-1,J,K)

              ENDDO ! I
            ENDDO ! J
          ENDDO ! K

        ENDIF

!---------------------------------------
! 2.1 Now the more common case of having
!     a number of processors in the
!     East-West direction


      ELSE ! If there is more than 1 processor East-West

        CALL GC_SSYNC(nproc,info) ! synchronise processors

!---------------------------------------
! 2.1.1 Copy the data into buf_send

        DO K=1,LEVELS
          DO J=1-HALO_Y,ROWS+HALO_Y
            DO I=1,HALO_X

              ! Copy stuff from the Western side of the grid
              send_buffer(EW_address(J,I,K,1))=FIELD(1+I-1,J,K)

              ! Copy stuff from the Eastern side of the grid
              send_buffer(EW_address(J,I,K,2))=                         &
     &          FIELD(ROW_LENGTH-HALO_X+I,J,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

!---------------------------------------
! 2.1.2 Send and receive the data

!---------------------------------------
! 2.1.2.1 Special case of 2 processors
!         East-West - both sides are
!         sent to the same processor
!         (if cyclic BC)

        IF ((nproc_x  ==  2) .AND. (bound(1)  ==  BC_CYCLIC)) THEN

          CALL GC_RSEND(1,2*EW_Halo_Size*LEVELS,neighbour(PEast),info,  &
     &                  receive_buffer,send_buffer)

          CALL GC_SSYNC(nproc,info)

          CALL GC_RRECV(1,2*EW_Halo_Size*LEVELS,neighbour(PWest),info,  &
     &                  receive_buffer,send_buffer)

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

          index_2_start=EW_address(1-HALO_Y,1,1,2)

          IF (neighbour(PWest)  /=  NoDomain) THEN

            ! Send West
            CALL GC_RSEND(2,EW_Halo_Size*LEVELS,neighbour(PWest),info,  &
     &                    receive_buffer,send_buffer)

          ENDIF

          IF (neighbour(PEast)  /=  NoDomain) THEN

            ! Send East
            CALL GC_RSEND(3,EW_Halo_Size*LEVELS,neighbour(PEast),info,  &
     &                    receive_buffer(index_2_start),                &
     &                    send_buffer(index_2_start))
          ENDIF

          CALL GC_SSYNC(nproc,info)

          IF (neighbour(PEast)  /=  NoDomain) THEN

            ! Receive from East
            CALL GC_RRECV(2,EW_Halo_Size*LEVELS,neighbour(PEast),info,  &
     &                    receive_buffer,send_buffer)

          ENDIF

          IF (neighbour(PWest)  /=  NoDomain) THEN

            ! Receive from West
            CALL GC_RRECV(3,EW_Halo_Size*LEVELS,neighbour(PWest),info,  &
     &                    receive_buffer(index_2_start),                &
     &                    send_buffer(index_2_start))

          ENDIF

        ENDIF ! test on numbers of processors East-West

!---------------------------------------
! 2.1.2 Fill the halos with data

        IF (neighbour(PEast)  /=  NoDomain) THEN

          ! unpack data from receive_buffer into FIELD

          DO K=1,LEVELS
            DO J=1-HALO_Y,ROWS+HALO_Y
              DO I=1,HALO_X
                FIELD(ROW_LENGTH+I,J,K)=                                &
     &            receive_buffer(EW_address(J,I,K,1))
              ENDDO
            ENDDO
          ENDDO

        ELSEIF (sb_Model_domain  ==  mt_global) THEN
             ! No neighbour to my East (ie. at edge of the domain
             ! and it's not a cyclic boundary condition)
             ! Just copy data from last column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

          DO K=1,LEVELS
            DO J=1-HALO_Y,ROWS+HALO_Y
              DO I=1,HALO_X
                FIELD(ROW_LENGTH+I,J,K)=FIELD(ROW_LENGTH,J,K)
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (neighbour(PEast)  /=  NoDomain)

        IF (neighbour(PWest)  /=  NoDomain) THEN

          ! unpack data from receive_buffer into FIELD

          DO K=1,LEVELS
            DO J=1-HALO_Y,ROWS+HALO_Y
              DO I=1,HALO_X
                FIELD(I-HALO_X,J,K)=                                    &
     &            receive_buffer(EW_address(J,I,K,2))
              ENDDO
            ENDDO
          ENDDO

        ELSEIF (sb_Model_domain  ==  mt_global) THEN
             ! No neighbour to my West (ie. at edge of the domain
             ! and it's not a cyclic boundary condition)
             ! Just copy data from first column
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the ELSE statement above.

          DO K=1,LEVELS
            DO J=1-HALO_Y,ROWS+HALO_Y
              DO I=1,HALO_X
                FIELD(I-HALO_X,J,K)=FIELD(1,J,K)
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (neighbour(PWest)  /=  NoDomain)

        CALL GC_SSYNC(nproc,info)

      ENDIF ! IF (nproc_x  ==  1)

      ENDIF    ! HALO_X >  0

!------------------------------------------------------------------
! 3.0 North-South communications

      IF(HALO_Y >  0) then

!  section of code for cyclic N-S boundary conditions
!  copied from above
      IF(bound(2)  ==  BC_CYCLIC) then

      IF (nproc_y  ==  1) THEN ! only 1 processor north-south
          south_halo_source=ROWs-HALO_y+1       ! copy from opposite end
          north_halo_source=1                   ! of each column
! 2.1 Simple case of only one processor

        DO K=1,LEVELS
          DO J=1,HALO_Y
            DO I=1-HALO_X ,row_length+halo_x

              ! Fill southern halo
              FIELD(I,1-halo_Y+J-1,K)=FIELD(I,south_halo_source+J-1,K)

              ! Fill northern halo
              FIELD(I,rows+J,K)=FIELD(I,north_halo_source+J-1,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

!---------------------------------------
! 2.1 Now the more common case of having
!     a number of processors in the
!     North-South direction


      ELSE ! If there is more than 1 processor north_south

        CALL GC_SSYNC(nproc,info) ! synchronise processors

!---------------------------------------
! 2.1.1 Copy the data into buf_send

        DO K=1,LEVELS
          DO J=1,HALO_Y
            DO I=1-HALO_X  , row_length +halo_x

              ! Copy stuff from the southern side of the grid
              send_buffer(NS_address(I,J,K,1))=FIELD(I,1+J-1,K)

              ! Copy stuff from the northern side of the grid
              send_buffer(NS_address(I,J,K,2))=                         &
     &          FIELD(I,rows-halo_y+J,K)

            ENDDO ! I
          ENDDO ! J
        ENDDO ! K

!---------------------------------------
! 2.1.2 Send and receive the data

!---------------------------------------
! 2.1.2.1 Special case of 2 processors
!         north-south - both sides are
!         sent to the same processor
!         (if cyclic BC)

       IF ((nproc_y  ==  2) .AND. (bound(2)  ==  BC_CYCLIC)) THEN

          CALL GC_RSEND(1,2*NS_Halo_Size*LEVELS,neighbour(PNorth),info, &
     &                  receive_buffer,send_buffer)

          CALL GC_SSYNC(nproc,info)

          CALL GC_RRECV(1,2*NS_Halo_Size*LEVELS,neighbour(PSouth),info, &
     &                  receive_buffer,send_buffer)

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

          index_2_start=NS_address(1-halo_x,1,1,2)


          IF (neighbour(PSouth)  /=  NoDomain) THEN

            ! Send south
            CALL GC_RSEND(2,NS_Halo_Size*LEVELS,neighbour(Psouth),info, &
     &                    receive_buffer,send_buffer)

          ENDIF

          IF (neighbour(Pnorth)  /=  NoDomain) THEN

            ! Send north
            CALL GC_RSEND(3,NS_Halo_Size*LEVELS,neighbour(Pnorth),info, &
     &                    receive_buffer(index_2_start),                &
     &                    send_buffer(index_2_start))
          ENDIF

          CALL GC_SSYNC(nproc,info)

          IF (neighbour(Pnorth)  /=  NoDomain) THEN

            ! Receive from north
            CALL GC_RRECV(2,NS_Halo_Size*LEVELS,neighbour(Pnorth),info, &
     &                    receive_buffer,send_buffer)

          ENDIF

          IF (neighbour(Psouth)  /=  NoDomain) THEN

            ! Receive from south
            CALL GC_RRECV(3,NS_Halo_Size*LEVELS,neighbour(Psouth),info, &
     &                    receive_buffer(index_2_start),                &
     &                    send_buffer(index_2_start))

          ENDIF

        ENDIF ! test on numbers of processors north-south

!---------------------------------------
! 2.1.2 Fill the halos with data

        IF (neighbour(Pnorth)  /=  NoDomain) THEN

          ! unpack data from receive_buffer into FIELD

          DO K=1,LEVELS
            DO J=1,HALO_Y
              DO I=1-halo_x,row_length+HALO_X
                FIELD(I,J+rows,K)=                                      &
     &            receive_buffer(NS_address(I,J,K,1))
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (neighbour(Pnorth)  /=  NoDomain)

        IF (neighbour(Psouth)  /=  NoDomain) THEN

          ! unpack data from receive_buffer into FIELD

          DO K=1,LEVELS
            DO J=1,HALO_Y
              DO I=1-halo_x,row_length+HALO_X
                FIELD(I,J-halo_y,K)=                                    &
     &            receive_buffer(ns_address(I,J,K,2))
              ENDDO
            ENDDO
          ENDDO

        ENDIF ! IF (neighbour(Psouth)  /=  NoDomain)

        CALL GC_SSYNC(nproc,info)

      ENDIF ! IF (nproc_y  ==  1)

      ELSE                 !!! bc_cyclic in NS
! Set up some variables

! Set up the offsets. When copying data that is to be passed over
! the pole, on wind (u or v) grid, then copy data one row away
! from the pole

      north_off=0
      south_off=0
      IF ((sb_Model_domain  ==  mt_global) .AND.                        &
     &    ((FIELD_TYPE  ==  fld_type_p) .OR.                            &
     &     (FIELD_TYPE  ==  fld_type_u))) THEN

        IF (at_extremity(PNorth)) north_off=1
        IF (at_extremity(PSouth)) south_off=1

      ENDIF

! Set up the sign factor. If L_VECTOR is true and data has been passed
! over the poles in a global model, then the variables must change
! sign

      IF (.NOT. ((L_VECTOR) .AND.                                       &
     &           (sb_Model_domain  ==  mt_global))) THEN
        change_sign=.FALSE.
      ELSE
        change_sign=.TRUE.
      ENDIF

!---------------------------------------
! 3.1 Copy data into the send_buffer
!     But not if:
!       - Not a global model and the data is at the North/South edge
!       - A global model at the North/South edge but only 1 processor
!         in the East-West direction.

      IF (.NOT. (at_extremity(PSouth) .AND.                             &
     &           ((nproc_x  ==  1)  .OR.                                &
     &           (sb_Model_domain  /=  mt_global))))                    &
     &  THEN

        DO K=1,LEVELS
          DO J=1,HALO_Y
            DO I=1-HALO_X,ROW_LENGTH+HALO_X

              send_buffer(NS_address(I,J,K,1))=                         &
     &          FIELD(I,J+south_off,K)

            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF (.NOT. (at_extremity(PNorth) .AND.                             &
     &           ((nproc_x  ==  1) .OR.                                 &
     &           (sb_Model_domain  /=  mt_global))))                    &
     &  THEN

        DO K=1,LEVELS
          DO J=1,HALO_Y
            DO I=1-HALO_X,ROW_LENGTH+HALO_X

              send_buffer(NS_address(I,J,K,2))=                         &
     &          FIELD(I,ROWS-HALO_Y-north_off+J,K)

            ENDDO
          ENDDO
        ENDDO
      ENDIF

!---------------------------------------
! 3.2 Send and receive the data

!---------------------------------------
! 3.2.1 The special case where nproc_y=1
!       Both buffers are sent to the
!       same processor as one message

      IF ((nproc_y  ==  1) .AND. (sb_Model_domain  ==  mt_global) .AND. &
     &    (nproc_x  >   1)) THEN

        CALL GC_RSEND(10,2*NS_Halo_Size*LEVELS,neighbour(PSouth),info,  &
     &                  receive_buffer,send_buffer)

        CALL GC_SSYNC(nproc,info)

        CALL GC_RRECV(10,2*NS_Halo_Size*LEVELS,neighbour(PSouth),info,  &
     &                receive_buffer,send_buffer)

      ENDIF

!---------------------------------------
! 3.2.2 The more general case, where
!       each buffer is sent to a
!       different processor

      index_2_start=NS_address(1-HALO_X,1,1,2)

      IF (at_extremity(PSouth)) THEN

        IF ((neighbour(PSouth)  /=  NoDomain) .AND.                     &
     &      (neighbour(PSouth)  /=  mype)) THEN

          CALL GC_RSEND(11,NS_Halo_Size*LEVELS,neighbour(PSouth),info,  &
     &                  receive_buffer(index_2_start),                  &
     &                  send_buffer)

        ENDIF

      ELSE ! not at the South

        CALL GC_RSEND(12,NS_Halo_Size*LEVELS,neighbour(PSouth),info,    &
     &                receive_buffer,send_buffer)

      ENDIF

      IF (at_extremity(PNorth)) THEN

        IF ((neighbour(PNorth)  /=  NoDomain) .AND.                     &
     &      (neighbour(PNorth)  /=  mype)) THEN

          CALL GC_RSEND(13,NS_Halo_Size*LEVELS,neighbour(PNorth),info,  &
     &                  receive_buffer,                                 &
     &                  send_buffer(index_2_start))

        ENDIF

      ELSE ! not at the North

        CALL GC_RSEND(14,NS_Halo_Size*LEVELS,neighbour(PNorth),info,    &
     &               receive_buffer(index_2_start),                     &
     &               send_buffer(index_2_start))

      ENDIF

      CALL GC_SSYNC(nproc,info)

      IF (at_extremity(PSouth)) THEN

        IF ((neighbour(PSouth)  /=  NoDomain) .AND.                     &
     &      (neighbour(PSouth)  /=  mype)) THEN

          CALL GC_RRECV(11,NS_Halo_Size*LEVELS,neighbour(PSouth),info,  &
     &                  receive_buffer(index_2_start),                  &
     &                  send_buffer)

        ENDIF

      ELSE ! not at the South

        CALL GC_RRECV(14,NS_Halo_Size*LEVELS,neighbour(PSouth),info,    &
     &                  receive_buffer(index_2_start),                  &
     &                  send_buffer(index_2_start))

      ENDIF

      IF (at_extremity(PNorth)) THEN

        IF ((neighbour(PNorth)  /=  NoDomain) .AND.                     &
     &      (neighbour(PNorth)  /=  mype)) THEN

          CALL GC_RRECV(13,NS_Halo_Size*LEVELS,neighbour(PNorth),info,  &
     &                  receive_buffer,                                 &
     &                  send_buffer(index_2_start))

        ENDIF

      ELSE

        CALL GC_RRECV(12,NS_Halo_Size*LEVELS,neighbour(PNorth),info,    &
     &                receive_buffer,                                   &
     &                send_buffer)

      ENDIF

!---------------------------------------
! 3.3 Fill the halos with data

!---------------------------------------
! 3.3.1 Southern halo

      IF (at_extremity(PSouth)) THEN

        IF (neighbour(PSouth)  ==  NoDomain) THEN
          IF (sb_Model_domain  /=  mt_lam) THEN

          ! Just copy adjacent rows into halo area
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the IF statement above.

          DO K=1,LEVELS
            DO J=1,HALO_Y
              DO I=1-HALO_X,ROW_LENGTH+HALO_X
                FIELD(I,J-HALO_Y,K)=FIELD(I,1,K)
              ENDDO
            ENDDO
          ENDDO
          ENDIF ! IF (sb_Model_domain  /=  mt_lam)

        ELSEIF (neighbour(PSouth)  ==  mype) THEN
        ! Local across pole difference

          IF (change_sign) THEN
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1,half_full_ROW_LENGTH
                    FIELD(half_ROW_LENGTH+I,1-J,K)=                     &
     &                -FIELD(I,J+south_off,K)
                    FIELD(I-HALO_X,1-J,K)=                              &
     &                -FIELD(half_ROW_LENGTH+I-HALO_X,J+south_off,K)
                ENDDO
              ENDDO
            ENDDO
          ELSE ! don't change sign
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1,half_full_ROW_LENGTH
                    FIELD(half_ROW_LENGTH+I,1-J,K)=                     &
     &                FIELD(I,J+south_off,K)
                    FIELD(I-HALO_X,1-J,K)=                              &
     &                FIELD(half_ROW_LENGTH+I-HALO_X,J+south_off,K)
                ENDDO
              ENDDO
            ENDDO
          ENDIF ! IF (change_sign)

        ELSE ! data is receive_buffer(index 2)

          IF (change_sign) THEN
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1-HALO_X,ROW_LENGTH+HALO_X
                  FIELD(I,1-J,K)=                                       &
     &              -receive_buffer(NS_address(I,J,K,2))
                ENDDO
              ENDDO
            ENDDO
          ELSE ! don't change sign
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1-HALO_X,ROW_LENGTH+HALO_X
                  FIELD(I,1-J,K)=                                       &
     &              receive_buffer(NS_address(I,J,K,2))
                ENDDO
              ENDDO
            ENDDO
          ENDIF !  IF (change_sign)

        ENDIF ! What type of South extremity

      ELSE ! IF (at_extremity(PSouth)

        ! not at a South extremity

        DO K=1,LEVELS
          DO J=1,HALO_Y
            DO I=1-HALO_X,ROW_LENGTH+HALO_X
              FIELD(I,J-HALO_Y,K)=                                      &
     &          receive_buffer(NS_address(I,J,K,2))
            ENDDO
          ENDDO
        ENDDO

      ENDIF ! IF (at_extremity(PSouth)


!---------------------------------------
! 3.3.2 Northern halo

      IF (at_extremity(PNorth)) THEN

        IF (neighbour(PNorth)  ==  NoDomain) THEN
          IF (sb_Model_domain  /=  mt_lam) THEN
          ! Just copy adjacent rows into halo area
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the IF statement above.

          DO K=1,LEVELS
            DO J=1,HALO_Y
              DO I=1-HALO_X,ROW_LENGTH+HALO_X
                FIELD(I,ROWS+J,K)=FIELD(I,ROWS,K)
              ENDDO
            ENDDO
          ENDDO
          ENDIF ! IF (sb_Model_domain  /=  mt_lam)

        ELSEIF (neighbour(PNorth)  ==  mype) THEN
        ! Local across pole difference

          IF (change_sign) THEN
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1,half_full_ROW_LENGTH

                  FIELD(half_ROW_LENGTH+I,ROWS+J,K)=                    &
     &              -FIELD(I,ROWS-J+1-north_off,K)
                  FIELD(I-HALO_X,ROWS+J,K)=                             &
     &              -FIELD(half_ROW_LENGTH+I-HALO_X,                    &
     &                     ROWS-J+1-north_off,K)

                ENDDO
              ENDDO
            ENDDO
          ELSE ! don't change sign
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1,half_full_ROW_LENGTH

                  FIELD(half_ROW_LENGTH+I,ROWS+J,K)=                    &
     &              FIELD(I,ROWS-J+1-north_off,K)
                  FIELD(I-HALO_X,ROWS+J,K)=                             &
     &              FIELD(half_ROW_LENGTH+I-HALO_X,                     &
     &                    ROWS-J+1-north_off,K)

                ENDDO
              ENDDO
            ENDDO
          ENDIF ! IF (change_sign)

        ELSE ! data is receive_buffer(index 1)

          IF (change_sign) THEN
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1-HALO_X,ROW_LENGTH+HALO_X
                  FIELD(I,ROWS+J,K)=                                    &
     &              -receive_buffer(NS_address(I,HALO_Y-J+1,K,1))
                ENDDO
              ENDDO
            ENDDO
          ELSE ! don't change sign
            DO K=1,LEVELS
              DO J=1,HALO_Y
                DO I=1-HALO_X,ROW_LENGTH+HALO_X
                  FIELD(I,ROWS+J,K)=                                    &
     &              receive_buffer(NS_address(I,HALO_Y-J+1,K,1))
                ENDDO
              ENDDO
            ENDDO
          ENDIF ! IF (change_sign)

        ENDIF ! What type of North extremity

      ELSE ! IF (at_extremity(PNorth)

        ! not at a North extremity

        DO K=1,LEVELS
          DO J=1,HALO_Y
            DO I=1-HALO_X,ROW_LENGTH+HALO_X
              FIELD(I,ROWS+J,K)=                                        &
     &          receive_buffer(NS_address(I,J,K,1))
            ENDDO
          ENDDO
        ENDDO

      ENDIF ! IF (at_extremity(PNorth)
      ENDIF

      ENDIF    ! HALO_Y >  0

      CALL GC_SSYNC(nproc,info)

 9999 CONTINUE
      RETURN
      END SUBROUTINE SWAP_BOUNDS
#endif
