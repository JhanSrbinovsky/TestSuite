#if defined(C96_1B) && defined(T3E)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ T3E specific optimised version of SWAP_BOUNDS

      SUBROUTINE SWAP_BOUNDS(                                           &

     &  FIELD, ROW_LENGTH, ROWS, LEVELS,                                &
                                          ! field
     &  HALO_X, HALO_Y,                                                 &
                                          ! halos
     &  FIELD_TYPE, L_VECTOR                                            &
                                          ! supporting information
     &  )

      IMPLICIT NONE

!  This code is based on the T3E version of SWAP_BOUNDS developed
!  for the ND model. The argument list has been reduced by bringing
!  some variables in via a common block.
!
! Purpose:
!   This subroutine takes care of all boundary swapping and
!   extending of arrays at the global boundaries. Data is swapped
!   across the poles for any non-zero halo size in the y direction
!   if it is a global model.
!
! Implementation/Portability
!   This version of SWAP_BOUNDS calls CRAY shmem routines
!   directly, rather than using the portable GCOM routines, so
!   is only suitable for running on CRAY/SGI machines supporting
!   shmem.
!
! Author: Paul Burton (based on UM(4.5) SWPBND1B)
! Current code owner: Paul Burton
!
! History
! Date       Version    Comment
! ----       -------    -------
! 12/10/98              CRAY/SGI specific optimised version.
!                       Paul Burton
! 17/3/99               1xn bug removed : Paul Burton
! 6/7/99     5.0        Argument list revised for UM5.0
! 25/5/00    5.2        Don't fill in "empty" extremity halos for
!                       non Global models which will have LBC data
!                       in these halos.  P.Burton
! 19/09/00   5.2        South halo was missed out in previous GPB1F501
!                       So now we add similar code to stop the south
!                       extremity halo being field in non global models
!                                                              P.Burton
! 12/1/00    5.2        Further fix to stop Mes. model exterior halos
!                       being overwritten.                     P.Burton
! 17/9/01    5.3        Fix to prevent north-south external halos from
!                       being overwritten in the mes model.  Z. Gardner
! 22/08/01   5.3        Immediately exit if no halos  P.Burton
! 14/09/01   5.3        Changed to use sb_model_domain variable
!                                                      P.Burton
! 22/11/01    5.3        Enable MPP as the only option for
!                        small executables         E.Leung
! 15/01/02   5.3        Add bi_cyclic LAM code        A. Malcolm
! 23/08/00   5.5        Modification for parallelisation of WAM.
!                       Bob Carruthers, Cray UK Inc(D.Holmes-Bell)
! 23/11/05   6.2        Removed all references to the wavemodel.
!                       T.Edwards
! 22/08/05   6.2        Fix defs to match exec_xref. P.Selwood
!
! Code Description:
!  Language: Fortran77 + CRAY extensions
!            (NB: Contains CRAY specific code, such as CRAY POINTERS,
!                 and calls to the shmem libraries)
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

! Local Variables

      INTEGER                                                           &
                        ! loop indicies:
     &  i,j,k                                                           &
                        ! spatial indicies

                        ! other quantities:
     &, half_size_x                                                     &
                        ! size_x/2
     &, half_row_length                                                 &
                        ! row_length/2
     &, j_off           ! offset for the j index required
                        ! for U+W fields across pole swaps

! Local variables required for the non-aligned SHMEM communications

      INTEGER                                                           &
     &  array_address                                                   &
                        ! address of my FIELD array
     &, size_x                                                          &
                        ! my total (including halos) X size
     &, size_y                                                          &
                        ! my total (including halos) Y size
     &, align_vars(3)                                                   &
                        ! equivalenced copies of above 3 variables

     &, remote_address                                                  &
                        ! address of FIELD array on remote processor
     &, remote_size_x                                                   & 
                        ! remote array's size_x
     &, remote_size_y                                                   & 
                        ! remote array's size_y
     &, remote_vars(3)  ! equivalenced copies of above 3 variables

      EQUIVALENCE                                                       &
! These EQUIVALENCEs allow us to communicate the data easily
! (as it's all contained within one contiguous array), but
! access the individual elements with meaningful variable names
     &  (align_vars(1)  , array_address)                                &
     &, (align_vars(2)  , size_x)                                       &
     &, (align_vars(3)  , size_y)                                       &

     &, (remote_vars(1) , remote_address)                               &
     &, (remote_vars(2) , remote_size_x)                                &
     &, (remote_vars(3) , remote_size_y)

! align_vars is put on a COMMON block to ensure that it is memory
! aligned over processors.
      COMMON /T3E_SWAPBOUNDS/ align_vars

      INTEGER                                                           &
                           ! variables describing the communications
     &  element_length                                                  &
                           ! length of element to send
     &, number_of_elements                                              &
                           ! number of elements to put/get
     &, local_stride                                                    &
                           ! stride through local array (1D addressing)
     &, remote_start                                                    &
                           ! start in remote array (1D addressing)
     &, remote_stride      ! stride in remote array (1D addressing)

      REAL                                                              &
     &  remote_array(                                                   &
     &  (ROW_LENGTH+1+2*HALO_X)*                                        &
     &  (ROWS+1+2*HALO_Y)*                                              &
     &  LEVELS)
      POINTER                                                           &
     & (ptr,remote_array)
! remote_array is used to generate the address of the FIELD array on
! the remote processor - so the size is made such that no out of bounds
! references will be made - nothing is physically written or read from
! this array. This size is defined assuming that the decomposition is
! such that the maximum difference between X_SIZE on any two processors
! is no more than 1 (and similarly for Y_SIZE).
! The Cray POINTER ptr allows the address of remote_array to be set by
! changing the value of ptr. ptr is set to the address of FIELD on the
! remote processor - so effectively allowing the address of any element
! of the remote FIELD array to be generated


!------------------------------------------------------------------
! 0.0 Check if there is anything to do

      IF ((HALO_X  ==  0) .AND. (HALO_Y  ==  0)) GOTO 9999

!------------------------------------------------------------------
! 1.0 Set up the align_vars COMMON block with my local information
!     so that other processors can access it.

      array_address=LOC(FIELD) ! set to the start address of my
                               ! FIELD array
      size_x=ROW_LENGTH+2*HALO_X
      size_y=ROWS+2*HALO_Y

      CALL barrier()

! Once this barrier is passed we know that we can safely access
! remote align_vars


       IF (nproc_x  ==  1) THEN
         half_size_x=size_x/2
         half_row_length=ROW_LENGTH/2
       ENDIF

!------------------------------------------------------------------
! 2.0 East-West Communications
!     We'll use a shmem_get as this means we don't have to do
!     another barrier before doing the North-South communications,
!     as the shmem_get operation ensures we've updated our EW halos
!     before we starting shmem_putting data in the North-South
!     direction. This is important for the corners, as we'll
!     effectively be shifting data from our EW neighbours to our
!     NS neighbours.

      IF (HALO_X  >   0) THEN  ! EW Halos exist

        local_stride=size_x
        number_of_elements=ROWS

!----------------------------
! 2.1 Get halo data from our Western neighbour

        IF (neighbour(PWest)  /=  NoDomain) THEN

! Get address and size information from my Western neighbour

          CALL shmem_get(remote_vars,align_vars,3,                      &
     &                   neighbour(PWest))

          ptr=remote_address ! address of FIELD on PE to West
          remote_stride=remote_size_x

          DO k=1,LEVELS
            DO i=1,HALO_X

              remote_start=(k-1)*remote_size_x*remote_size_y +          &
     &                     (HALO_Y+1)*remote_size_x -                   &
     &                      2*HALO_X + i

! The shmem_iget call is a strided shmem_get - the data for a
! whole column of one level of Western halo is got from the
! neighbouring processor by the single call

              CALL shmem_iget(                                          &
     &          FIELD(i-HALO_X,1,k),remote_array(remote_start),         &
     &          local_stride,remote_stride,                             &
     &          number_of_elements,                                     &
     &          neighbour(PWest))

            ENDDO ! i
          ENDDO ! k

        ELSEIF (sb_Model_domain  ==  mt_global) THEN
             ! No processor to my West, so just fill the
             ! halos with valid numbers

        ENDIF  ! IF (neighbour(PWest)  /=  NoDomain)

!----------------------------
! 2.2 Get halo data from our Eastern neighbour

        IF (neighbour(PEast)  /=  NoDomain) THEN

! Get address and size information from my Western neighbour

          CALL shmem_get(remote_vars,align_vars,3,                      &
     &                   neighbour(PEast))

          ptr=remote_address ! address of FIELD on PE to East
          remote_stride=remote_size_x

          DO k=1,LEVELS
            DO i=1,HALO_X

              remote_start=(k-1)*remote_size_x*remote_size_y +          &
     &                     HALO_Y*remote_size_x + HALO_X + i

! The shmem_iget call is a strided shmem_get - the data for a
! whole column of one level of Eastern halo is got from the
! neighbouring processor by the single call

              CALL shmem_iget(                                          &
     &          FIELD(ROW_LENGTH+i,1,k),                                &
     &          remote_array(remote_start),                             &
     &          local_stride,remote_stride,                             &
     &          number_of_elements,                                     &
     &          neighbour(PEast))

            ENDDO ! i
          ENDDO ! k

        ELSEIF (sb_Model_domain  ==  mt_global) THEN
             ! No processor to my East, so just fill the
             ! halos with valid numbers

        ENDIF  ! IF (neighbour(East)  /=  NoDomain)

      ENDIF ! IF (HALO_X  >   0)

!------------------------------------------------------------------
! 3.0 North-South Communications

      IF (HALO_Y  >   0) THEN  ! NS Halos exist

!----------------------------
! 3.1 Put data to our Northern neighbour

        IF (neighbour(PNorth)  /=  NoDomain) THEN

! Get address and size information from my Northern neighbour

          CALL shmem_get(remote_vars,align_vars,3,                      &
     &                   neighbour(PNorth))

          ptr=remote_address ! address of FIELD on PE to North
          remote_stride=remote_size_x*remote_size_y


          IF (.NOT. at_extremity(PNorth) .OR.                           &
     &                  Bound(2) == bc_cyclic) THEN

            remote_start=1
            element_length=size_x*HALO_Y

            DO k=1,LEVELS

              CALL shmem_put(                                           &
     &          remote_array(remote_start+(k-1)*remote_stride),         &
     &          FIELD(1-HALO_X,ROWS-HALO_Y+1,k),                        &
     &          element_length,neighbour(PNorth))

            ENDDO ! k

          ELSE ! this processor is at the Northern edge

            IF (sb_Model_domain  ==  mt_global) THEN
! Do across pole swap

              IF (FIELD_TYPE  ==  fld_type_v) THEN
                j_off=0
              ELSE ! Type_U or Type_W
                j_off=1
              ENDIF

              IF (nproc_x  >   1) THEN

                element_length=size_x

                DO k=1,LEVELS
                  DO j=1,HALO_Y

                    remote_start=(ROWS + HALO_Y + j - 1)*size_x + 1

                    CALL shmem_put(                                     &
     &                remote_array(remote_start+(k-1)*remote_stride),   &
     &                FIELD(1-HALO_X,ROWS-j+1-j_off,k),                 &
     &                element_length,neighbour(PNorth))

                  ENDDO ! j
                ENDDO ! k

              ELSE ! only one processor East-West

                DO k=1,LEVELS
                  DO j=1,HALO_Y
                    DO i=1,half_size_x

                      FIELD(half_ROW_LENGTH+i,ROWS+j,k)=                &
     &                  FIELD(i,ROWS-j+1-j_off,k)

                      FIELD(i-HALO_X,ROWS+j,k)=                         &
     &                  FIELD(half_ROW_LENGTH+i-HALO_X,                 &
     &                  ROWS-j+1-j_off,k)
                    ENDDO ! i
                  ENDDO ! j
                ENDDO ! k

              ENDIF ! IF (nproc_x  >   1)

            ELSE ! If this isn't a Global model domain
           IF (sb_Model_domain  ==  mt_global) THEN
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the IF statement above.
! Put some sensible numbers into the halo

              element_length=size_x*HALO_Y

              DO k=1,LEVELS

! We use a shmem_put to do a local copy, as it is faster
! than a straight Fortran copy loop which would go through
! cache.

                CALL shmem_put(                                         &
     &            FIELD(1-HALO_X,ROWS+1,k),                             &
     &            FIELD(1-HALO_X,ROWS-HALO_Y+1,k),                      &
     &            element_length,mype)

              ENDDO ! k
           ENDIF ! IF (Model_domain  ==  mt_global)

            ENDIF ! IF (Model_domain  ==  mt_Global)

          ENDIF ! IF (.NOT. at_extremity(North) or bound=cyclic )

        ELSEIF (bound(2)  ==  BC_STATIC .AND.                           &
     &          sb_Model_domain  /=  mt_lam) THEN
          ! If we have no neighbour to our North and it's static
          ! boundary conditions.
! Put some sensible numbers into the halo

          element_length=size_x*HALO_Y

          DO k=1,LEVELS

! We use a shmem_put to do a local copy, as it is faster
! than a straight Fortran copy loop which would go through
! cache.

            CALL shmem_put(                                             &
     &        FIELD(1-HALO_X,ROWS+1,k),                                 &
     &        FIELD(1-HALO_X,ROWS-HALO_Y+1,k),                          &
     &        element_length,mype)

          ENDDO ! k

        ENDIF ! IF (neighbour(North)  /=  NoDomain)


!----------------------------
! 3.3 Put data to our Southern neighbour

        IF (neighbour(PSouth)  /=  NoDomain) THEN

! Get address and size information from my Southern neighbour

          CALL shmem_get(remote_vars,align_vars,3,                      &
     &                   neighbour(PSouth))

          ptr=remote_address ! address of FIELD on PE to South
          remote_stride=remote_size_x*remote_size_y


          IF (.NOT. at_extremity(PSouth) .OR.                           &
     &                  Bound(2) == bc_cyclic) THEN

            remote_start=(remote_size_y-HALO_Y)*size_x + 1
            element_length=size_x*HALO_Y

            DO k=1,LEVELS

              CALL shmem_put(                                           &
     &          remote_array(remote_start+(k-1)*remote_stride),         &
     &          FIELD(1-HALO_X,1,k),                                    &
     &          element_length,neighbour(PSouth))

            ENDDO ! k

          ELSE ! this processor is at the Southern edge

            IF (sb_Model_domain  ==  mt_Global) THEN
! Do across pole swap

              IF (FIELD_TYPE  ==  fld_type_v) THEN
                j_off=0
              ELSE ! U or theta
                j_off=1
              ENDIF

              IF (nproc_x  >   1) THEN

                element_length=size_x

                DO k=1,LEVELS
                  DO j=1,HALO_Y

                    remote_start=(HALO_Y-j)*size_x + 1

                    CALL shmem_put(                                     &
     &                remote_array(remote_start+(k-1)*remote_stride),   &
     &                FIELD(1-HALO_X,j+j_off,k),                        &
     &                element_length,neighbour(PSouth))

                  ENDDO ! j
                ENDDO ! k

              ELSE ! only one processor East-West

                DO k=1,LEVELS
                  DO j=1,HALO_Y
                    DO i=1,half_size_x

                      FIELD(half_ROW_LENGTH+i,1-j,k)=                   &
     &                  FIELD(i,j+j_off,k)

                      FIELD(i-HALO_X,1-j,k)=                            &
     &                  FIELD(half_ROW_LENGTH+i-HALO_X,j+j_off,k)

                    ENDDO ! i
                  ENDDO ! j
                ENDDO ! k

              ENDIF ! IF (nproc_x  >   1)

            ELSE ! If this isn't a Global model domain
           IF (sb_Model_domain  ==  mt_global) THEN
! NOTE: This block of code should never be executed. It's being left
!       in so that it can be used in the future by changing the logic
!       in the IF statement above.
! Put some sensible numbers into the halo

              element_length=size_x*HALO_Y

              DO k=1,LEVELS

! We use a shmem_put to do a local copy, as it is faster
! than a straight Fortran copy loop which would go through
! cache.

                CALL shmem_put(                                         &
     &            FIELD(1-HALO_X,1-HALO_Y,k),                           &
     &            FIELD(1-HALO_X,1,k),                                  &
     &            element_length,mype)

              ENDDO ! k

           ENDIF ! IF (Model_domain  ==  mt_global)
            ENDIF ! IF (Model_domain  ==  Global)

          ENDIF ! IF (.NOT. at_extremity(North) or bound=cyclic )

        ELSEIF (bound(2)  ==  BC_STATIC .AND.                           &
     &          sb_Model_domain  /=  mt_lam) THEN
          ! If we have no neighbour to our South and it's static
          ! boundary conditions.
! Put some sensible numbers into the halo

          element_length=size_x*HALO_Y

          DO k=1,LEVELS

! We use a shmem_put to do a local copy, as it is faster
! than a straight Fortran copy loop which would go through
! cache.

            CALL shmem_put(                                             &
     &        FIELD(1-HALO_X,1-HALO_Y,k),                               &
     &        FIELD(1-HALO_X,1,k),                                      &
     &        element_length,mype)

          ENDDO ! k

        ENDIF ! IF (neighbour(PSouth)  /=  NoDomain)

      ENDIF ! IF (HALO_Y  >   0)

!------------------------------------------------------------------
! 4.0 Finished Communications
!     Do a barrier so we know that all the data has arrived

      CALL barrier()

!------------------------------------------------------------------
! 5.0 Tidy up Poles
!     Change sign of across pole swaps if L_VECTOR is set

      IF ((L_VECTOR) .AND. (sb_Model_domain  ==  mt_Global)) THEN

        IF ( at_extremity(PNorth) ) THEN

          DO k=1,LEVELS
            DO j=ROWS+1,ROWS+HALO_Y
              DO i=1-HALO_X,ROW_LENGTH+HALO_X

                FIELD(i,j,k)=-FIELD(i,j,k)

              ENDDO ! i
            ENDDO ! j
          ENDDO ! k

        ENDIF ! IF ( at_extremity(PNorth) )

        IF ( at_extremity(PSouth) ) THEN

          DO k=1,LEVELS
            DO j=1-HALO_Y,0
              DO i=1-HALO_X,ROW_LENGTH+HALO_X

                FIELD(i,j,k)=-FIELD(i,j,k)

              ENDDO ! i
            ENDDO ! j
          ENDDO ! k

        ENDIF ! IF ( at_extremity(PSouth) )

      ENDIF ! IF ((L_VECTOR) .AND. (Model_domain  ==  Global))

!------------------------------------------------------------------
! End of Routine

 9999 CONTINUE
      RETURN
      END SUBROUTINE SWAP_BOUNDS
#endif
