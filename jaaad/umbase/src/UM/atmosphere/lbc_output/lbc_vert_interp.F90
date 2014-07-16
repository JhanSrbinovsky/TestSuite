#if defined(A32_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Vertical Interpolation of LBC data
!
! Subroutine Interface:

      Subroutine LBC_Vert_Interp (                                      &
     &           lbc_vi_data_in                                         &
     &,          lbc_vi_orog                                            &
     &,          lbc_vi_data_out                                        &
     &,          lbc_seg_size                                           &
     &,          level_type                                             &
! src
     &,          src_model_levels                                       &
     &,          src_levels                                             &
     &,          src_first_level                                        &
     &,          src_last_level                                         &
     &,          src_ht_gen_method                                      &
     &,          src_first_r_rho                                        &
     &,          src_z_top_model                                        &
     &,          src_eta_theta                                          &
     &,          src_eta_rho                                            &
! lbc
     &,          lbc_model_levels                                       &
     &,          lbc_levels                                             &
     &,          lbc_first_level                                        &
     &,          lbc_last_level                                         &
     &,          lbc_ht_gen_method                                      &
     &,          lbc_first_r_rho                                        &
     &,          lbc_z_top_model                                        &
     &,          lbc_eta_theta                                          &
     &,          lbc_eta_rho                                            &
     & )

      IMPLICIT NONE

!
! Description:
!   Performs Vertical interpolation of Lateral Boundary Conditions (LBC)
!   from model to LBC levels.
!
! Method:
!   1. Heights of rho & theta levels are calculated for both model and
!      LBC levels at the LBC points.
!   2. VERT_INTERP is called to interpolate the LBCs. Linear
!      interpolation used for all variables.
!
! Current Code Owner: Dave Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2    13/11/00   Original code. Dave Robinson
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!

      Integer  ::  lbc_seg_size
      Integer  ::  level_type     !  Theta or Rho levels

! source vertical grid
      Integer  :: src_model_levels
      Integer  :: src_levels
      Integer  :: src_first_level
      Integer  :: src_last_level
      Integer  :: src_ht_gen_method
      Integer  :: src_first_r_rho
      Real     :: src_z_top_model
      Real     :: src_eta_theta (0:src_model_levels)
      Real     :: src_eta_rho   (  src_model_levels)
! lbc vertical grid
      Integer  :: lbc_model_levels
      Integer  :: lbc_levels
      Integer  :: lbc_first_level
      Integer  :: lbc_last_level
      Integer  :: lbc_ht_gen_method
      Integer  :: lbc_first_r_rho
      Real     :: lbc_z_top_model
      Real     :: lbc_eta_theta (0:lbc_model_levels)
      Real     :: lbc_eta_rho   (  lbc_model_levels)

! data
      Real  ::  lbc_vi_data_in  (lbc_seg_size,                          &
     &                           src_first_level:src_last_level)
      Real  ::  lbc_vi_data_out (lbc_seg_size,                          &
     &                           lbc_first_level:lbc_last_level)
      Real  ::  lbc_vi_orog     (lbc_seg_size)

! Local variables
      Integer, Parameter :: Rho_Levels   = 1
      Integer, Parameter :: Theta_Levels = 2
      Integer, Parameter :: Interp_Order = 1 ! Linear Interpolation

      Character (Len=*), Parameter ::  RoutineName = 'LBC_Vert_Interp'

      Integer :: Level
      Integer :: ErrorStatus
      Character (Len=80) :: CMessage

! Height fields to be calculated
      Real, dimension (:,:), allocatable :: src_r_theta_levels
      Real, dimension (:,:), allocatable :: src_r_rho_levels
      Real, dimension (:,:), allocatable :: lbc_r_theta_levels
      Real, dimension (:,:), allocatable :: lbc_r_rho_levels

! ----------------------------------
! Check validity of height algorithm
! ----------------------------------

!     Height algorithm for model and lbcs must be the same.
      If (src_ht_gen_method /= lbc_ht_gen_method) Then
          Write (CMessage,*) 'Mismatch in height algorithm for ',       &
     &    'model ',src_ht_gen_method,' and LBCs ',lbc_ht_gen_method
          ErrorStatus = 10
! DEPENDS ON: ereport
          Call Ereport (RoutineName, ErrorStatus, CMessage)
      End If

! ------------------------------------
! Allocate space for the height fields
! ------------------------------------

      allocate (src_r_theta_levels(lbc_seg_size,0:src_model_levels)  )
      allocate (src_r_rho_levels  (lbc_seg_size,  src_model_levels+1))

      allocate (lbc_r_theta_levels(lbc_seg_size,0:lbc_model_levels)  )
      allocate (lbc_r_rho_levels  (lbc_seg_size,  lbc_model_levels+1))

! ---------------------------------------
! Calculate heights for the source levels
! ---------------------------------------

! DEPENDS ON: lbc_calc_heights
      Call lbc_calc_heights ( lbc_seg_size,                             &
     &                        src_model_levels, src_ht_gen_method,      &
     &                        src_first_r_rho, src_z_top_model,         &
     &                        src_eta_theta, src_eta_rho,               &
     &                        lbc_vi_orog,                              &
     &                        src_r_theta_levels, src_r_rho_levels )

! ------------------------------------
! Calculate heights for the lbc levels
! ------------------------------------

! DEPENDS ON: lbc_calc_heights
      Call lbc_calc_heights ( lbc_seg_size,                             &
     &                        lbc_model_levels, lbc_ht_gen_method,      &
     &                        lbc_first_r_rho, lbc_z_top_model,         &
     &                        lbc_eta_theta, lbc_eta_rho,               &
     &                        lbc_vi_orog,                              &
     &                        lbc_r_theta_levels, lbc_r_rho_levels )

! -----------------------------
! Do the Vertical Interpolation
! -----------------------------

      Select Case (Level_Type)

        Case (Rho_Levels)

          Do level = lbc_first_level, lbc_last_level

!           write (6,*) ' VI on Rho levels for level ',level

! DEPENDS ON: vert_interp
            Call Vert_Interp (lbc_vi_data_in, lbc_seg_size,             &
     &                        src_levels, lbc_r_rho_levels(1,level),    &
     &                        src_r_rho_levels(1,src_first_level),      &
     &                        interp_order,                             &
     &                        lbc_vi_data_out(1,level) )

          End Do

        Case (Theta_Levels)

          Do level = lbc_first_level, lbc_last_level

!           write (6,*) ' VI on Theta levels for level ',level

! DEPENDS ON: vert_interp
            Call Vert_Interp (lbc_vi_data_in, lbc_seg_size,             &
     &                        src_levels, lbc_r_theta_levels(1,level),  &
     &                        src_r_theta_levels(1,src_first_level),    &
     &                        interp_order,                             &
     &                        lbc_vi_data_out(1,level) )

          End Do

        Case Default

          Write (CMessage,*) 'LBC Vertical Interpolation not ',         &
     &                       'catered for Level Type ',Level_Type
          ErrorStatus = 20
! DEPENDS ON: ereport
          Call Ereport (RoutineName, ErrorStatus, CMessage)

      End Select

! --------------------
! Deallocate workspace
! --------------------

      deallocate ( src_r_theta_levels )
      deallocate ( src_r_rho_levels   )
      deallocate ( lbc_r_theta_levels )
      deallocate ( lbc_r_rho_levels   )

      Return
      END SUBROUTINE LBC_Vert_Interp
#endif
