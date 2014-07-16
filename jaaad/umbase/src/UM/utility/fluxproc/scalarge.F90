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


































      SUBROUTINE ScalarGE                                               &
     &           (nx, ny, rmdi,                                         &
     &            scalar, in_field,                                     &
     &            out_field,                                            &
     &            icode, cmessage)


!     ScalarGE: Set MDI if scalar >= field
!     -----------
!                 Takes account of missing data.

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  loop over all elements
!

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
                        ! IN scalar
     &      ,in_field(nx,ny)   ! IN array of input values

!     Output
!     ------

      REAL                                                              &
     &       out_field(nx,ny)  ! OUT results of test

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
      cmessage='ScalarGE: test successful.'

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns

          IF (    scalar  <   in_field(ix,iy) ) THEN

             out_field(ix,iy) = in_field(ix,iy)

           ELSE      ! set missing data if either input is missing

             out_field(ix,iy) = rmdi

           END IF

         END DO  ! ix: loop over columns
      END DO     ! iy: loop over columns

      RETURN
      END SUBROUTINE ScalarGE


!----------------------------------------------------------------------
#endif
