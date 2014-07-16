#if defined(C96_0A)

SUBROUTINE swap_bounds_mv(   &
        input_fields,        &
        n_multi,             &
        row_length,          &
        halo_x,              &
        halo_y )

! Description: 
!  Dummy stub version of routine SWAP_BOUNDS_MV to be used when 
!  MPP service routines are not selected.
!
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!
! Author: Paul Selwood
! Current code owner: Paul Selwood

  Use swapable_field_mod, Only: &
      swapable_field_pointer_type

  IMPLICIT NONE

! Comdecks

! Arguments:

INTEGER, INTENT(IN)  :: n_multi      ! number of fields to be swapped
INTEGER, INTENT(IN)  :: row_length   ! number of points on row (no halos)
INTEGER, INTENT(IN)  :: halo_x       ! sixze of "i" halo
INTEGER, INTENT(IN)  :: halo_y       ! sixze of "j" halo

! Fields to swap
TYPE(swapable_field_pointer_type), TARGET :: input_fields(n_multi)


! Local Variables:
INTEGER             :: ErrorStat  
CHARACTER(Len=14)   :: RoutineName
CHARACTER(Len=80)   :: Message    


!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------

ErrorStat   = 1
RoutineName = 'SWAP_BOUNDS_MV'
Message     = 'MPP service routines unavailable - see output.'

WRITE (6,*) '**ERROR**: SWAP_BOUNDS_MV called but is unavailable.'
WRITE (6,*) '  Sections C96_1A, C96_1B or C96_1C are required '
WRITE (6,*) '  (MPP Service routines)                         '

! DEPENDS ON: ereport
CALL ereport(RoutineName, ErrorStat, Message)


!-------------------------------------------------------------------------------

RETURN
END SUBROUTINE SWAP_BOUNDS
#endif
