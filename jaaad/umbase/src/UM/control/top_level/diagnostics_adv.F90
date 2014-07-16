#if defined(CONTROL) && defined(ATMOS)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate diagnostic quantities from the initial atmosphere dump
!
! Subroutine Interface:
      SUBROUTINE Diagnostics_adv(                                       &
     & row_length,rows,n_rows,model_levels,wet_levels,                  &
! primary wind fields:
     & u,v,theta,q,qcl,qcf,cf,cfl,cff,                                  &
! wind field increments after advection:
     & R_u,R_v,R_W,                                                     &
! wind field increments before advection (on stashflag):
     & u_incr_diagnostic,v_incr_diagnostic,                             &
     & T_incr_diagnostic,q_incr_diagnostic,                             &
     & qcl_incr_diagnostic,qcf_incr_diagnostic,                         &
     & cf_incr_diagnostic,cfl_incr_diagnostic,cff_incr_diagnostic,      &
     & theta_star,q_star,qcl_star,qcf_star,                             &
     & cf_star,cfl_star,cff_star,                                       &
     & exner_theta_levels,                                              &
! w departure point information
     & depart_lambda,depart_phi,depart_r,r_theta_levels,                &
#include "argsts.h"
     & STASHwork)

Use ac_diagnostics_mod, Only :                                          &
    qcl_adv

      IMPLICIT NONE
!
! Description:
!   Diagnostics_adv extracts diagnostics of (N+1) time-level estimates
!   of primary fields after advection has been called, to be processed
!   by STASH routines for UM section 12 (advection).
! Method:
!   Simple copying of relevant fields into STASHwork array for output.
!   Vertical compression of levels is handled by copydiag_3D routine
!   using STASH control arrays passed via argsts.h.
!   Sequential processing of each diagnostic is performed, dependent
!   upon STASH flags being set.
!   List of diagnostics: (item,section)
!   (  2,12) u wind               = u + R_u
!   (  3,12) v wind               = v + R_v
!   (  4,12) temperature          = theta_star/exner_theta_levels
!   ( 10,12) specific humidity    = q_star
!   (254,12) qcl                  = qcl_star
!   ( 12,12) qcf                  = qcf_star
!   (185,12) u increment          = delta(R_u) across advection
!   (186,12) v increment          = delta(R_v) across advection
!   (187,12) w increment          = delta(w  ) across advection
!   (181,12) T increment          = delta(theta)/exner
!   (182,12) q   increment        = delta(q)   across advection
!   (183,12) qcl increment        = delta(qcl) across advection
!   (184,12) qcf increment        = delta(qcf) across advection
!   (188,12) cf  increment        = delta(cf)   across advection
!   (199,12) cfl increment        = delta(cfl) across advection
!   (190,12) cff increment        = delta(cff) across advection
!   (204,12) departure point (w)  = depart_lambda
!   (203,12) departure point (w)  = depart_phi
!   (205,12) model height diff    = depart_r - r_theta_levels
!
!   _star fields are estimates of N+1 time-level quantities; R_u/v are
!   wind increments from physics1.
!
!   Note: no error trapping - cosmetic Errorstatus/cmessage needed for
!   compilation - no checks performed at lower levels.
!
! Current Code Owner: R Rawlins
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
#include "parvars.h"
#include "csubmodl.h"
#include "typsts.h"

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & row_length,rows                                                  &
                        ! horizontal dimensions
     &,n_rows                                                           &
                        ! rows for last (N) row of pes
     &,model_levels                                                     &
                        ! vertical levels
     &,wet_levels       ! vertical levels with moisture

!   Array  arguments with intent(in):
      REAL                                                              &
     &   u(1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,  v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)    &
     &,theta(1-offx:row_length+offx, 1-offy: rows+offy, model_levels)   &
     &,   q(1-halo_i:row_length+halo_i,                                 &
     &      1-halo_j:rows+halo_j, wet_levels)                           &
     &, qcl(1-halo_i:row_length+halo_i,                                 &
     &      1-halo_j:rows+halo_j, wet_levels)                           &
     &, qcf(1-halo_i:row_length+halo_i,                                 &
     &      1-halo_j:rows+halo_j, wet_levels)                           &
     &, cf (1-halo_i:row_length+halo_i,                                 &
     &      1-halo_j:rows+halo_j, wet_levels)                           &
     &, cfl(1-halo_i:row_length+halo_i,                                 &
     &      1-halo_j:rows+halo_j, wet_levels)                           &
     &, cff(1-halo_i:row_length+halo_i,                                 &
     &      1-halo_j:rows+halo_j, wet_levels)                           &
     &,R_u(1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,R_v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)    &
     &,R_w(row_length, rows, model_levels)                              &
     &,u_incr_diagnostic(row_length,        rows      , model_levels)   &
     &,v_incr_diagnostic(row_length,      n_rows      , model_levels)   &
     &,T_incr_diagnostic(row_length,        rows      , model_levels)   &
     &,q_incr_diagnostic(row_length,        rows      , wet_levels)     &
     &,qcl_incr_diagnostic(row_length,      rows      , wet_levels)     &
     &,qcf_incr_diagnostic(row_length,      rows      , wet_levels)     &
     &,cf_incr_diagnostic (row_length,      rows      , wet_levels)     &
     &,cfl_incr_diagnostic(row_length,      rows      , wet_levels)     &
     &,cff_incr_diagnostic(row_length,      rows      , wet_levels)     &
     &,theta_star                                                       &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,q_star                                                           &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,qcl_star                                                         &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,qcf_star                                                         &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,cf_star                                                          &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,cfl_star                                                         &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,cff_star                                                         &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy,   wet_levels)    &
     &,exner_theta_levels                                               &
     &    (1-offx:row_length+offx, 1-offy:  rows+offy, model_levels)    &
     &,depart_lambda(row_length, rows, model_levels)                    &
     &,depart_phi(row_length, rows, model_levels)                       &
     &,depart_r(row_length, rows, model_levels)                         &
     &,  R_THETA_LEVELS(1-halo_i:row_length+halo_i,                     &
     &            1-halo_j:rows+halo_j,0:model_levels)


!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      REAL                                                              &
     & STASHwork(*)   ! Output array holding diagnostic fields

! Local parameters:
      CHARACTER*(*) RoutineName
      PARAMETER (   RoutineName='Diagnostics_adv')

      INTEGER sect      ! STASH section for diagnostics
      PARAMETER ( sect = 12 )

! Local scalars:
      INTEGER                                                           &
     & i,j,k,ji                                                         &
                   !  loop indices
     &,im_index                                                         &
                   !  internal model index for STASH arrays
     &,item                                                             &
                   !  STASH item of diagnostic
     &,Errorstatus !  Error status

      CHARACTER*80 CMessage !  Error message

! Local dynamic arrays:
      REAL                                                              &
     & work_1(row_length,  rows,model_levels)                           &
     &,work_2(row_length,n_rows,model_levels)                           &
     &,w_incr_diagnostic(row_length, rows, model_levels)

! Function & Subroutine calls:
      External Copydiag_3D,Ereport

!- End of header

!
! 1. Initialisation
!
      im_index    = internal_model_index(atmos_im)
      Errorstatus = 0
      Cmessage    = ''
!
! 2. Extract diagnostic fields dependent on STASHflags sf
!

! u wind estimate = u + physics1 increment
      item = 2             ! u
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              work_1(i,j,k) = u(i,j,k) + R_u(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_1,      &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! v wind estimate = v + physics1 increment
      item = 3             ! v
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,n_rows
            DO i=1,row_length
              work_2(i,j,k) = v(i,j,k) + R_v(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_2,      &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)


! temperature estimate = theta / exner pressure
      item = 4             ! v
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              work_1(i,j,k)=theta_star(i,j,k)*exner_theta_levels(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),work_1,      &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! specific humidity
      item = 10            ! q
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),q_star,      &
     &        row_length,rows,wet_levels,0,0,1,1,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcl
      item = 254           ! qcl
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

      If (.not.allocated(qcl_adv)) Then
        allocate ( qcl_adv(row_length*rows,wet_levels) )
      End If

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),qcl_star,    &
     &        row_length,rows,wet_levels,0,0,1,1,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

! Copy non-halo area of qcl_star into qcl_adv in module ac_diagnostics_mod
      do k = 1,wet_levels
        do j = 1,rows
          do i = 1,row_length
            ji = (j-1)*row_length+i
            qcl_adv(ji,k) = qcl_star(i,j,k)
          end do
        end do
      end do

      END IF ! sf(item,sect)

! qcf
      item = 12            ! qcf
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),qcf_star,    &
     &        row_length,rows,wet_levels,0,0,1,1,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! u wind increment
      item = 185           ! u increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              u_incr_diagnostic(i,j,k) = R_u(i,j,k) -                   &
     &                                        u_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        u_incr_diagnostic,                                        &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! v wind increment
      item = 186           ! v increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,model_levels
          DO j=1,n_rows
            DO i=1,row_length
              v_incr_diagnostic(i,j,k) = R_v(i,j,k) -                   &
     &                                        v_incr_diagnostic(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        v_incr_diagnostic,                                        &
     &        row_length,n_rows,model_levels,0,0,0,0,at_extremity,      &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! w wind increment
      item = 187           ! w increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              w_incr_diagnostic(i,j,k) = R_w(i,j,k)
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

      END IF ! sf(item,sect)

! T wind increment
! theta_star now holds theta+dtheta whereas t_incr holds dtheta before
! advection
      item = 181           ! T increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,model_levels
          DO j=1,rows
!CDIR NOUNROLL
            DO i=1,row_length
              T_incr_diagnostic(i,j,k) = (theta_star(i,j,k) -           &
     &                 (theta(i,j,k) + T_incr_diagnostic(i,j,k)))       &
     &                                 *exner_theta_levels(i,j,k)
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

      END IF ! sf(item,sect)

! q wind increment
      item = 182           ! q increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,wet_levels
!CDIR NOUNROLL
          DO j=1,rows
            DO i=1,row_length
              q_incr_diagnostic(i,j,k) = q_star(i,j,k) -                &
     &                     (q(i,j,k) + q_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        q_incr_diagnostic,                                        &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcl wind increment
      item = 183           ! qcl increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,wet_levels
!CDIR NOUNROLL
          DO j=1,rows
            DO i=1,row_length
              qcl_incr_diagnostic(i,j,k) = qcl_star(i,j,k) -            &
     &                     (qcl(i,j,k) + qcl_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        qcl_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! qcf increment
      item = 184           ! qcf increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,wet_levels
!CDIR NOUNROLL
          DO j=1,rows
            DO i=1,row_length
              qcf_incr_diagnostic(i,j,k) = qcf_star(i,j,k) -            &
     &                     (qcf(i,j,k) + qcf_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        qcf_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cf increment
      item = 192           ! cf increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,wet_levels
!CDIR NOUNROLL
          DO j=1,rows
            DO i=1,row_length
              cf_incr_diagnostic(i,j,k) = cf_star(i,j,k) -              &
     &                     (cf(i,j,k) + cf_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cf_incr_diagnostic,                                       &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cfl increment
      item = 193           ! cfl increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,wet_levels
!CDIR NOUNROLL
          DO j=1,rows
            DO i=1,row_length
              cfl_incr_diagnostic(i,j,k) = cfl_star(i,j,k) -            &
     &                     (cfl(i,j,k) + cfl_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cfl_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! cf increment
      item = 194           ! cff increment
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,wet_levels
!CDIR NOUNROLL
          DO j=1,rows
            DO i=1,row_length
              cff_incr_diagnostic(i,j,k) = cff_star(i,j,k) -            &
     &                     (cff(i,j,k) + cff_incr_diagnostic(i,j,k))
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &        cff_incr_diagnostic,                                      &
     &        row_length,rows,wet_levels,0,0,0,0,at_extremity,          &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! Departure point diagnostics for w
! (a) lambda
      item = 204
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &         depart_lambda,                                           &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)
! (b) phi
      item = 205
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &         depart_phi,                                              &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)
! (c) dr  difference from model height of departure point
      item = 203
      IF (sf(item,sect) .AND. Errorstatus == 0) THEN
        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              work_1(i,j,k) = depart_r(i,j,k) -                         &
     &                                        r_theta_levels(i,j,k)
            END DO  ! i
          END DO  ! j
        END DO  ! k

! DEPENDS ON: copydiag_3d
        CALL copydiag_3d(STASHwork(si(item,sect,im_index)),             &
     &         work_1,                                                  &
     &        row_length,rows,model_levels,0,0,0,0,at_extremity,        &
     &        stlist(1,stindex(1,item,sect,im_index)),len_stlist,       &
     &        stash_levels,num_stash_levels+1,                          &
     &        atmos_im,sect,item,                                       &
     &        Errorstatus,cmessage)

      END IF ! sf(item,sect)

! 3. Error handling
!
      IF (Errorstatus /= 0) THEN
! DEPENDS ON: ereport
        CALL Ereport(RoutineName,Errorstatus,Cmessage)
      END IF

      RETURN
      END SUBROUTINE Diagnostics_adv
#endif
