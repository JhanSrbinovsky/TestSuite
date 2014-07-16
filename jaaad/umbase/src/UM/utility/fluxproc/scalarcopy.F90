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





























      SUBROUTINE ScalarCopy                                             &
     &           (nx, ny, rmdi,                                         &
     &            scalar,                                               &
     &            out_field,                                            &
     &            icode, cmessage)


!     ScalarCopy: subroutine to copy a scalar to
!     ---------- an array.

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  loop over all elements.

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

     &      ,scalar     ! IN scalar

!     Output
!     ------

      REAL                                                              &
     &       out_field(nx,ny)  ! OUT results of addition

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
      cmessage='ScalarCopy: copy successful.'

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns


           out_field(ix,iy) = scalar

         END DO  ! ix: loop over columns
      END DO     ! iy: loop over columns

      RETURN
      END SUBROUTINE ScalarCopy







!----------------------------------------------------------------------
#endif
