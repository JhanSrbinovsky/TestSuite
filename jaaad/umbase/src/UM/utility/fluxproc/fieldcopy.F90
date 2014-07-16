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








      SUBROUTINE FieldCopy                                              &
     &           (nx, ny, rmdi,                                         &
     &            in_field,                                             &
     &            out_field,                                            &
     &            icode, cmessage)


!     FieldCopy: subroutine to copy one array to another.
!     ----------

!     Programmer: S J Foreman
!     -----------

!     Method:
!     -------  copy one array to another.

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

     &      ,in_field(nx,ny)  ! IN array of input values

!     Output
!     ------

      REAL                                                              &
     &       out_field(nx,ny)  ! OUT copy of values in input

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
      cmessage='FieldCopy: copy successful.'

      DO iy = 1, ny     ! iy: loop over rows
        DO ix = 1, nx  ! ix: loop over columns

           out_field(ix,iy) = in_field(ix,iy)

         END DO  ! ix: loop over columns
      END DO     ! iy: loop over columns

      RETURN
      END SUBROUTINE FieldCopy




























!----------------------------------------------------------------------
#endif
