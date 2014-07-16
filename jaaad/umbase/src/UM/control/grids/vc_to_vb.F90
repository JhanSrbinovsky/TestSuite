#if defined(C92_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate diagnostic quantities from the initial atmosphere dump
!
! Subroutine Interface:
      SUBROUTINE vC_to_vB(vC,                                           &
     & row_length,n_rows,levels,offx,offy,                              &
     & vB)

      IMPLICIT NONE
!
! Description:
!   vC_to_vB performs a simple horizontal interpolation of a field vC
!   (without halos) at the v position for Arakawa 'C' grid staggering
!   onto the uv position of the 'B' grid.
! Method:
!   1. Copy field (without halos) to an array with halos.
!   2. Populate halos using swapbounds.
!   3. Make a simple average of adjacent columns and copy to output
!      field (without halos).
!
! Current Code Owner: R Rawlins
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 5.1    25/01/00  Original code. R Rawlins.
! 5.2    28/07/00  Fill RHS halo for LAM case   P.Burton
! 5.4    11/03/02  Remove comment on same line as #include
!                                               S. Carroll
! 6.2    03/02/06  Moved to C92_2A. P.Selwood
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
! field_type meaningful names
#include "fldtype.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & row_length                                                       &
                        ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows for last (N) row of pes
     &,levels                                                           &
                        ! vertical levels
     &,offx,offy        ! fixed halo sizes

!   Array  arguments with intent(in):
      REAL                                                              &
     & vC(row_length,n_rows,levels) ! Field on v points on C grid

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     & vB(row_length,n_rows,levels) ! Field on uv points on B grid

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='vC_to_vB')

! Local scalars:
      INTEGER                                                           &
     & i,j,k     !  loop indices

! Local dynamic arrays:
      REAL                                                              &
     & vC_halo(1-offx:row_length+offx,1-offy:n_rows+offy,levels)

! Function & Subroutine calls:
      External Swap_bounds

!- End of header

!   1. Copy field (without halos) to an array with halos

      DO k=1,levels
        DO j=1,n_rows
        DO i=1,row_length
          vC_halo(i,j,k)=vC(i,j,k)
        ENDDO ! i
! The following line is required for LAM where the RHS halo
! at the extremity of the model is not filled in by SWAP_BOUNDS
        vc_halo(row_length+1,j,k)=vc(row_length,j,k)
        ENDDO ! j
      ENDDO ! k

!   2. Populate halos using swapbounds.

! DEPENDS ON: swap_bounds
      CALL Swap_bounds(vC_halo,row_length,n_rows,levels,                &
     &                 offx,offy,fld_type_v,.true.)

!   3. Make a simple average of adjacent columns and copy to output
!      field (without halos).

      DO k=1,levels
        DO j=1,n_rows
        DO i=1,row_length
          vB(i,j,k)=(vC_halo(i,j,k)+vC_halo(i+1,j,k)) * 0.5
        ENDDO ! i
        ENDDO ! j
      ENDDO ! k

      RETURN
      END SUBROUTINE vC_to_vB
#endif
