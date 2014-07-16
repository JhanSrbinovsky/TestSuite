#if defined(FLUXPROC) || defined(FLXPLPR)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
! History:
! version  date         change
! 4.5      03/09/98     New code
!
! Author:     M. J. Bell
!----------------------------------------------------------------------
! deck: FIELDADD
!
! contains routines: FieldAdd, etc. ScalarAdd etc.
!
! Purpose: Flux processing routines.
!          Simple arithmetic operations on pp fields
!----------------------------------------------------------------------






      SUBROUTINE FieldDiv                                               &
     &           (nx, ny, rmdi,                                         &
     &            in_field1, in_field2,                                 &
     &            out_field,                                            &
     &            icode, cmessage)


!     FieldDiv: subroutine to divide elements of one array by
!     --------- corresponding elements of another.
!               Takes account of missing data.

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  loop over all elements, dividing unless either is
!              missing or divisor is zero, in which case result
!              is missing.
!              in_field1 / in_field2.

      IMPLICIT NONE


!     Input:
!     ------

      INTEGER                                                           &
     &       nx                                                         &
                        ! IN number of columns of array
     &      ,ny         ! IN number of rows in array

       REAL                                                             &
     &       rmdi                                                       &
                        ! IN value of REAL missing data indicator

     &      ,in_field1(nx,ny)                                           &
                               ! IN array of input values (numerator)
     &      ,in_field2(nx,ny)  ! IN array of input values (denominator)

!     Output
!     ------

      REAL                                                              &
     &       out_field(nx,ny)  ! OUT results of division

      INTEGER                                                           &
     &       icode      ! OUT Completion code

      CHARACTER *(*)                                                    &
     &       cmessage   ! OUT Error message


!    Local variables
!    ---------------

      INTEGER                                                           &
     &        ix                                                        &
                        ! Loop counter over columns
     &       ,iy        ! Loop counter over rows


! ------------------------------------------------------------------

      icode =0
      cmessage='FieldDiv: division successful.'

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns

          IF (        in_field1(ix,iy)  /=  rmdi                        &
     &          .AND. in_field2(ix,iy)  /=  rmdi )                      &
     &      THEN

              IF ( in_field2(ix,iy)  /=  0.0 )                          &

     &          THEN  ! divide
                  out_field(ix,iy) = in_field1(ix,iy)/in_field2(ix,iy)

                ELSE  ! avoid division by zero
                  icode = 4
                  out_field(ix,iy) = rmdi
                END IF

            ELSE      ! set missing data if either input is missing

              out_field(ix,iy) = rmdi

            END IF

         END DO  ! ix: loop over columns
      END DO     ! iy: loop over columns

!     Set appropriate error message if division by zero has occurred

      IF (icode  /=  0) THEN
         cmessage = 'FieldDiv: division by zero'
      END IF

      RETURN
      END SUBROUTINE FieldDiv






























!----------------------------------------------------------------------
#endif
