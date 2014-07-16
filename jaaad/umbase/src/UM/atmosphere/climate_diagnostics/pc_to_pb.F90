#if defined(A30_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Reformates a field from rho points to uv position on 'C' grid
!
! Subroutine Interface:
      SUBROUTINE pc_to_pb(pc,                                           &
     & row_length,rows,n_rows,levels,offx,offy,                         &
     & pb)

      IMPLICIT NONE
!
! Description:
!   pC_to_pB performs a simple horizontal interpolation of a field pC
!   at the p position for Arakawa 'C' grid staggering onto the uv
!   position of the 'B' grid.
! Method:
!   1. Copy field (without halos) to an array with halos.
!   2. Populate halos using swapbounds.
!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).
!
! Current Code Owner: S Wilson
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.1    25/01/00  Original code. S Wilson.
!  6.1   21/05/04    correct bit reproducibility problem     A. Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
#include "fldtype.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & row_length,rows                                                  &
                        ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows for last (N) row of pes
     &,levels                                                           &
                        ! vertical levels
     &,offx,offy        ! fixed halo sizes

!   Array  arguments with intent(in):
      REAL                                                              &
     & pc(row_length,rows,levels) ! Field on u points on C grid

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     & pb(row_length,n_rows,levels) ! Field on uv points on B grid

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='pc_to_pb')

! Local scalars:
      INTEGER                                                           &
     & i,j,k     !  loop indices

! Local dynamic arrays:
      REAL                                                              &
     & pc_halo(1-offx:row_length+offx,1-offy:rows+offy,levels)

! Function & Subroutine calls:
      External Swap_bounds

!- End of header

!   1. Copy field (without halos) to an array with halos

      DO k=1,levels
        DO j=1,rows
          DO i=1,row_length
            pc_halo(i,j,k)=pC(i,j,k)
          ENDDO                 ! i
        ENDDO                   ! j
      ENDDO                     ! k

!   2. Populate halos using swapbounds.

! DEPENDS ON: swap_bounds
      CALL Swap_bounds(pc_halo,row_length,rows,levels,                  &
     &                 offx,offy,fld_type_p,.true.)

!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).

      DO k=1,levels
!CDIR NOUNROLL
        DO j=1,n_rows
          DO i=1,row_length
            pb(i,j,k)=(pc_halo(i,j,k)+pc_halo(i,j+1,k)+                 &
     &        pc_halo(i+1,j,k)+pc_halo(i+1,j+1,k)) * 0.25
          ENDDO                 ! i
        ENDDO                   ! j
      ENDDO                     ! k

      RETURN
      END SUBROUTINE pc_to_pb
#endif
