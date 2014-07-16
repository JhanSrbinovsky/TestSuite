! Description: Define valid extra data vector types for
!              use in pp fields.
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.4   11/04/02   Original code.  R. Hill
!===============================================================
      INTEGER,PARAMETER :: NO_EXTRA_VECTORS=14 ! Number of valid extra
                                               ! data vector types
      INTEGER :: EXTRA_VECTOR(NO_EXTRA_VECTORS)

      ! Vector type 10 is not currently supported!
      DATA EXTRA_VECTOR/1,2,3,4,5,6,7,8,9,11,12,13,14,15/
