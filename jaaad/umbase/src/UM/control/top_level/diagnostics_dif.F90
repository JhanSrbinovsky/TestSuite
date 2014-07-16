#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate quantities from diffusion and divergence damping (sect 13)
!
! Subroutine Interface:
      SUBROUTINE Diagnostics_dif(                                       &
     & row_length,rows,n_rows,model_levels,wet_levels,                  &
! primary fields:
     & theta,q,                                                         &
! wind field increments after difusion:
     & R_u,R_v,                                                         &
! wind field increments before diffusion (on stashflag):
     & u_incr_diagnostic,v_incr_diagnostic,                             &
     & T_incr_diagnostic,q_incr_diagnostic,                             &
     & L_subfilter_horiz, L_subfilter_vert, visc_m, visc_h, shear,      &
     & RNEUTML, FM_3D, FH_3D, BL_LEVELS, BL_COEF_KM,BL_COEF_KH,         &
     & w_local_mask,                                                    &
! Current theta+dtheta and q+dq values
     & theta_star,q_star,                                               &
     & exner_theta_levels,                                              &
#include "argsts.h"
     & STASHwork)

      IMPLICIT NONE
!
! Description:
!   Diagnostics_fildif calculates diagnostics for divergence damping
!   and diffusion for output by STASH routines for UM section 13.
!
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays passed via argsts.h.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (181,13) T increment          = delta(theta)/exner
!   (182,13) q   increment        = delta(q)   across advection
!   (185,13) u increment          = delta(R_u) across advection
!   (186,13) v increment          = delta(R_v) across advection
!
!   _star fields are estimates of N+1 time-level quantities; R_u/v are
!   wind increments from all routines.
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Current Code Owner: A Malcolm
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.4   05/08/02  Original code. R A Stratton
!  5.5  23/02/03  Local q diffusion points    Terry Davies
!  6.2  10/02/06  Add 3D subgrid turbulence scheme.  Carol Halliwell
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
                        ! rows for last (N) row of pes
     &,model_levels                                                     &
                        ! vertical levels
     &,wet_levels       ! vertical levels with moisture

!   Array  arguments with intent(in):
      REAL, INTENT(IN) ::                                               &
     & theta(1-offx:row_length+offx,1-offy: rows+offy, model_levels)    &
     &,q(1-halo_i:row_length+halo_i,1-halo_j: rows+halo_j, wet_levels)  &
     &,R_u(1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,R_v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)    &
     &,theta_star                                                       &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,q_star                                                           &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,exner_theta_levels                                               &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)

       Logical                                                          &
     & L_subfilter_horiz                                                &
     &,L_subfilter_vert

       Integer                                                          &
     & BL_LEVELS

!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
      REAL, INTENT(INOUT) ::                                            &
     & u_incr_diagnostic(row_length,        rows      , model_levels)   &
     &,v_incr_diagnostic(row_length,      n_rows      , model_levels)   &
     &,T_incr_diagnostic(row_length,        rows      , model_levels)   &
     &,q_incr_diagnostic(row_length,        rows      , wet_levels)     &
     &,visc_m(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j          &
     &                                              ,model_levels)      &
!            ! diffusion coefficient for momentum from subgrid
!            ! turbulence scheme
     &,visc_h(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j          &
     &                                              ,model_levels)      &
!            ! diffusion coefficient for heat and moisture from
!            ! subgrid turbulence scheme
!     &,shear(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &,shear(row_length, rows, model_levels)                            &
!               ! S from subgrid turbulence scheme
     &,RNEUTML(1:row_length, 1:rows,model_levels)                       &
!               ! mixing length scale (lambda)
     &,FM_3D(row_length,rows,BL_LEVELS)                                 &
!               ! stability function for momentum transport.
!               ! level 1 value is dummy for use in diagnostics
     &,FH_3D(row_length,rows,BL_LEVELS)                                 &
!               ! stability function for heat and moisture.
!               ! level 1 value is dummy for use in diagnostics
     &,BL_COEF_KM(1:row_length, 1:rows, bl_levels-1)                    &
!               ! RHOKM from BL scheme
     &,BL_COEF_KH(1:row_length, 1:rows, bl_levels-1)                    &
!               ! RHOKH from BL scheme
     &,w_local_mask(row_length,        rows)

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
      CHARACTER (LEN=*), PARAMETER :: RoutineName='Diagnostics_dif'

      INTEGER, PARAMETER :: sect =13  ! STASH section for diagnostics

! Local scalars:
      INTEGER  ::                                                       &
     & i,j,k                                                            &
                   !  loop indices
     &,im_index                                                         &
                   !  internal model index for STASH arrays
     &,item        !  STASH item of diagnostic
      INTEGER :: Errorstatus = 0  ! initial value for error code

      CHARACTER (LEN=80) :: CMessage !  Error message

      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     & work_visc


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

! T increment
      item = 181           ! T increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              T_incr_diagnostic(i,j,k) = (theta_star(i,j,k) -           &
     &                        T_incr_diagnostic(i,j,k))                 &
     &                                  *exner_theta_levels(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        T_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

! q  increment
      item = 182           ! q increment
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

        DO k=1,wet_levels
          DO j=1,rows
            DO i=1,row_length
              q_incr_diagnostic(i,j,k) = q_star(i,j,k) -                &
     &                             q_incr_diagnostic(i,j,k)
            ENDDO  ! i
          ENDDO  ! j
        ENDDO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        q_incr_diagnostic,                                        &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      ENDIF ! sf(item,sect)

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
      ALLOCATE (work_visc(row_length, rows, model_levels))

      If (.not. L_subfilter_vert .and. .not. L_subfilter_horiz) then
        sf(190,sect)=.false.
        sf(191,sect)=.false.
        sf(192,sect)=.false.
        sf(193,sect)=.false.
        sf(194,sect)=.false.
        sf(195,sect)=.false.
        sf(196,sect)=.false.
        sf(197,sect)=.false.
      Else If (.not. L_subfilter_vert) then
        sf(196,sect)=.false.
        sf(197,sect)=.false.  
      End if

      item = 190           ! momentum viscosity coeff
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_visc(i,j,k) = visc_m(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_visc,                                                &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 191           ! scalar viscosity coeff
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_visc(i,j,k) = visc_h(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_visc,                                                &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If
      item = 192 ! shear
      If(sf(item,sect).AND.Errorstatus == 0) THEN

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              work_visc(i,j,k) = shear(i,j,k)
            End Do
          End Do
        End Do

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        work_visc,                                                &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 193 ! mixing length
      If(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        RNEUTML,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,                     &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 194 ! Ri dependent function FM_3D
      If(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        FM_3D,                                                    &
     &        row_length,rows,BL_levels,0,0,0,0,                        &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)
      End If

      item = 195 ! Ri dependent function FH_3D
      If(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        FH_3D,                                                    &
     &        row_length,rows,BL_levels,0,0,0,0,                        &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 196 ! Diffusion coeff from BL scheme
      If(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_COEF_KM,                                               &
     &        row_length,rows,BL_levels-1,0,0,0,0,                      &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      item = 197 ! Diffusion coeff from BL scheme
      If(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        Call copydiag_3d (stashwork(si(item,sect,im_index)),            &
     &        BL_COEF_KH,                                               &
     &        row_length,rows,BL_levels-1,0,0,0,0,                      &
     &        at_extremity,                                             &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      End If

      DEALLOCATE (work_visc)

! Counter for occurances of local q diffusion
      item = 201           ! local q diffusion at a point
      IF(sf(item,sect).AND.Errorstatus == 0) THEN

! DEPENDS ON: copydiag
        CALL copydiag(STASHwork(si(item,sect,im_index)),                &
     &        w_local_mask,                                             &
     &        row_length,rows,0,0,0,0,at_extremity,                     &
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
      END SUBROUTINE Diagnostics_dif
#endif
