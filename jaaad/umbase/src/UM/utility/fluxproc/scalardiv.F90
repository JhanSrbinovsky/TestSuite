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




























      SUBROUTINE ScalarDiv                                              &
     &           (nx, ny, rmdi,                                         &
     &            scalar, in_field,                                     &
     &            out_field,                                            &
     &            icode, cmessage)


!     ScalarDiv: subroutine to divide a scalar by
!     ---------- an array.
!                Takes account of missing data.

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  loop over all elements, dividing unless either is
!              missing or divisor is zero, in which case result
!              is missing.

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

     &      ,scalar                                                     &
                        ! IN scalar (numerator)
     &      ,in_field(nx,ny)   ! IN array of input values (denominator)

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
      cmessage='ScalarDiv: division successful.'

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns

          IF (    scalar  /=  rmdi                                      &
     &        .AND. in_field(ix,iy)  /=  rmdi )                         &
     &      THEN

              IF ( in_field(ix,iy)  /=  0.0 )                           &

     &          THEN  ! divide
                  out_field(ix,iy) = scalar/in_field(ix,iy)

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
         cmessage = 'ScalarDiv: division by zero'
      END IF

      RETURN
      END SUBROUTINE ScalarDiv








!----------------------------------------------------------------------
#endif
