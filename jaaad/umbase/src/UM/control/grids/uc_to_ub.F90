#if defined(C92_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate diagnostic quantities from the initial atmosphere dump
!
! Subroutine Interface:
      SUBROUTINE uC_to_uB(uC,                                           &
     & row_length,rows,n_rows,levels,offx,offy,                         &
     & uB)

      IMPLICIT NONE
!
! Description:
!   uC_to_uB performs a simple horizontal interpolation of a field uC
!   (without halos) at the u position for Arakawa 'C' grid staggering
!   onto the uv position of the 'B' grid.
! Method:
!   1. Copy field (without halos) to an array with halos.
!   2. Populate halos using swapbounds.
!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).
!
! Current Code Owner: R Rawlins
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.1    25/01/00  Original code. R Rawlins.
! 5.2    28/07/00  Fill bottom halo for LAM case   P.Burton
! 6.2    03/02/06  Moved to C92_2A. P.Selwood
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
     & uC(row_length,rows,levels) ! Field on u points on C grid

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     & uB(row_length,n_rows,levels) ! Field on uv points on B grid

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='uC_to_uB')

! Local scalars:
      INTEGER                                                           &
     & i,j,k     !  loop indices

! Local dynamic arrays:
      REAL                                                              &
     & uC_halo(1-offx:row_length+offx,1-offy:rows+offy,levels)

! Function & Subroutine calls:
      External Swap_bounds

!- End of header

!   1. Copy field (without halos) to an array with halos

      DO k=1,levels
        DO j=1,rows
        DO i=1,row_length
          uC_halo(i,j,k)=uC(i,j,k)
        ENDDO ! i
        ENDDO ! j
! The following lines are required for LAM where the bottom halo
! at the extremity of the model is not filled in by SWAP_BOUNDS
        DO i=1,row_length
          uc_halo(i,rows+1,k)=uc(i,rows,k)
        ENDDO
      ENDDO ! k

!   2. Populate halos using swapbounds.

! DEPENDS ON: swap_bounds
      CALL Swap_bounds(uC_halo,row_length,rows,levels,                  &
     &                 offx,offy,fld_type_u,.true.)

!   3. Make a simple average of adjacent rows and copy to output field
!     (without halos).

      DO k=1,levels
        DO j=1,n_rows
        DO i=1,row_length
          uB(i,j,k)=(uC_halo(i,j,k)+uC_halo(i,j+1,k)) * 0.5
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      RETURN
      END SUBROUTINE uC_to_uB
#endif
