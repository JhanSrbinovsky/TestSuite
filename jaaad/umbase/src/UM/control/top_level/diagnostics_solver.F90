#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Interface:
      SUBROUTINE Diagnostics_solver(                                    &
     & row_length,rows,n_rows,model_levels,                             &
! wind field increments after difusion:
     & R_u,R_v,R_w,                                                     &
! wind field increments before diffusion (on stashflag):
     & u_incr_diagnostic,v_incr_diagnostic,                             &
     & w_incr_diagnostic,                                               &
#include "argsts.h"
     & STASHwork)

      IMPLICIT NONE
!
! Description:
!   Diagnostics_solver calculates diagnostics for the Helmholtz solver
!   for output by STASH routines for UM section 10.
!
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays passed via argsts.h.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (185,10) u increment          = delta(R_u) across solver
!   (186,10) v increment          = delta(R_v) across solver
!   (187,10) w increment          = delta(R_v) across solver
!
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Original code author : R A Stratton
! Current  Code Owner  : A Malcolm
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5   03/01/03  Original code. R A Stratton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable  ! Description of variable
!
! Global variables (#include statements etc):
#include "parvars.h"
#include "csubmodl.h"
#include "typsts.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN) :: row_length,rows                            &
                                              ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows size for v-arrays
     &,model_levels     ! vertical levels

!   Array  arguments with intent(in):
      REAL, INTENT(IN) ::                                               &
     & R_u(1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,R_v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)    &
     &,R_w(row_length, rows, model_levels)

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
      REAL, INTENT(INOUT) ::                                            &
     & u_incr_diagnostic(row_length,        rows      , model_levels)   &
     &,v_incr_diagnostic(row_length,      n_rows      , model_levels)   &
     &,w_incr_diagnostic(row_length,        rows      , model_levels)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
      CHARACTER (LEN=*), PARAMETER :: RoutineName='Diagnostics_dif'

      INTEGER, PARAMETER :: sect =10  ! STASH section for diagnostics

! Local scalars:
      INTEGER  ::                                                       &
     & i,j,k                                                            &
                   !  loop indices
     &,im_index                                                         &
                   !  internal model index for STASH arrays
     &,item        !  STASH item of diagnostic
      INTEGER :: Errorstatus = 0  ! initial value for error code

      CHARACTER (LEN=80) :: CMessage !  Error message

! Local dynamic arrays: NONE at present
!      REAL  ::

! Function & Subroutine calls:
      External Copydiag_3D,Ereport

!- End of header

!
! 1. Initialisation
!
      im_index    = internal_model_index(atmos_im)
      Cmessage    = ''
!
! 2. Extract diagnostic fields dependent on STASHflags sf
!
! u wind increment
      item = 185           ! u increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                   &
     &                                        u_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

! v wind increment
      item = 186           ! v increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,n_rows
            DO i=1,row_length
              v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                   &
     &                                        v_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

! w  increment
      item = 187           ! w increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              w_incr_diagnostic(i,j,k) = R_w(i,j,k) -                   &
     &                             w_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        w_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

!
! 3. Error handling
!
      IF(Errorstatus /= 0) THEN
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,Errorstatus,Cmessage)
      ENDIF

      RETURN
      END SUBROUTINE Diagnostics_solver
#endif
