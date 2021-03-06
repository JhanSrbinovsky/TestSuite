#if defined(A05_4A) || defined(A05_5A)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!  Subroutine CALC_3D_CCA: Calculates a conv. cld amt on theta model levels.
!
!  Subroutine Interface:

subroutine calc_3d_cca(np_field, npnts, nlev, nbl, cld_base, cld_top,         &
           p_lyr_bnds, frz_lev, cca_2d, cca_3d, z_theta, z_rho, l_q_interact, &
           l_use_sh_mask, l_pc2_diag_sh_pts)

  Use cv_cntl_mod, Only:                                                      &
      lcv_ccrad

  Use cv_run_mod, Only:                                                       &
      l_cloud_deep, tower_factor, anvil_factor, anv_opt

  IMPLICIT NONE
!
! Description:
! ------------
! Calculates a 3D convective cloud amount (i.e. on theta model levels) from
! the 2D convective cloud amount array according to parameters specified in
! the umui and the position of cloud base, cloud top and freezing level.
!
! Method:
! -------
! The 2D convective cloud amount is expanded into the vertical by applying it 
! between the cloud base and top with the additional constraints that:
!
!         (ia)  If the cloud base is in the boundary layer (Original)
!         (ib)  If the cloud base is below the freezing level (CCRad)
!         (ii)  Cloud top is above the freezing level and
!         (iii) The cloud is more than 500mb deep
!
! Then the cloud below the freezing level will be multiplied by TOWER_FACTOR,
! and the cloud above the freezing level will be linearly
! (model level/height/pressure(default)) increased to cloud top where it will
! be equal to the 2D fraction * ANVIL_FACTOR.
!
! NOTE: ***** The above method needs to be rewritten if these mods are********
!       ***** Implemented*****************************************************
!
! Current Code Owner: Julie M. Gregory
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Global variables (*CALLed COMDECKs etc...):

!-----------------------------------------------------------------------------
!   Scalar arguments with intent(in):
!-----------------------------------------------------------------------------

  INTEGER, intent(in)  :: npnts        ! Number of points
  INTEGER, intent(in)  :: np_field     ! Full data length
  INTEGER, intent(in)  :: nlev         ! Number of levels
  INTEGER, intent(in)  :: nbl          ! Number of Boundary layer levels

  LOGICAL, intent(in) :: l_q_interact  ! .TRUE. : PC2 cloud scheme is in use

  LOGICAL, intent(in) :: l_pc2_diag_sh_pts(npnts) 
                                       ! .TRUE. : If replacing pc2 cloud with
                                       !          diagnostic shallow Cloud

  LOGICAL, intent(in) :: l_use_sh_mask ! Logical to indicate whether or not
                                       ! to use l_pc2_diag_sh_pts as a test
                                       ! criteria


!-----------------------------------------------------------------------------
!   Array  arguments with intent(in):
!-----------------------------------------------------------------------------

  REAL,    intent(in)  :: z_theta    (np_field,   nlev) ! z (th layer centres)
  REAL,    intent(in)  :: z_rho      (np_field,   nlev) ! z (rh level  bounds)
  REAL,    intent(in)  :: p_lyr_bnds (np_field, 0:nlev)
                                                  ! Pressure on layer
                                                  ! boundaries (rho levels-1)

  REAL,    intent(in)  :: cca_2d   (npnts) ! 2D convective cloud amount
  INTEGER, intent(in)  :: cld_top  (npnts) ! Conv. cloud top  (theta level)
  INTEGER, intent(in)  :: cld_base (npnts) ! Conv. cloud base (theta level)
  INTEGER, intent(in)  :: frz_lev  (npnts) ! Freezing level   (theta level)




!-----------------------------------------------------------------------------
!   Array  arguments with intent(out):
!-----------------------------------------------------------------------------

  REAL,    intent(out) :: cca_3d(np_field, nlev)  ! Convective cloud amount on
                                                  ! model levels
                                                  ! (theta levels)


! Local variables:
! -----------------

  REAL,    parameter  :: deep_dp=50000.0  ! Min. depth for anvil criteria(Pa)

  ! For all the options below except 3), the anvil base is at the freezing
  ! Level

  INTEGER, parameter  :: anv_pressure                = 0 
  INTEGER, parameter  :: anv_height                  = 1
  INTEGER, parameter  :: anv_model_levels            = 2
  INTEGER, parameter  :: anv_limited_pressure_depth  = 3
  INTEGER             :: i, k             ! Loop counters

                         
  INTEGER :: anv_lev     ! Base level of 'anvil'
  REAL    :: anv_dep     ! Anvil depth in model levels
  REAL    :: anv_p_dep   ! Anvil depth in pressure
  REAL    :: anv_z_dep   ! Anvil depth in metres
  REAL    :: anv_p_base  ! Anvil base pressure (rho-level)

  REAL    :: p_cld_base  ! Pressure at lowest  cloud layer BOUNDARY
  REAL    :: p_cld_top   ! Pressure at highest cloud layer BOUNDARY



  LOGICAL :: cbct_crit (npnts) ! .TRUE. if cloud base/top are sensible
  LOGICAL :: dep_crit  (npnts) ! .TRUE. if depth criteria met
  LOGICAL :: base_crit (npnts) ! .TRUE. if cloud base criteria met
  LOGICAL :: anv_on    (npnts) ! .TRUE. if all anvil criteria met
  INTEGER :: tp_of_lp  (npnts) ! Index of cloud top, required so that CCRad
                               ! correction can be reverted

!-----------------------------------------------------------------------------
! Code Statements
!-----------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  ! 1.0 Initialise local arrays
  !---------------------------------------------------------------------------
  Do i=1, npnts
    cbct_crit(i)  = .FALSE.
    dep_crit(i)   = .FALSE.
    base_crit(i)  = .FALSE.
    anv_on(i)     = .FALSE.
    tp_of_lp(i)   = 0
  End Do



  !---------------------------------------------------------------------------
  ! 2.0 Option changes/bugfixes specified by CCRad
  !---------------------------------------------------------------------------
  If (lcv_ccrad) Then

    Do i=1, npnts
      ! Check for sensible Cloud base/top
      cbct_crit(i)  = ((cld_top(i)  >= cld_base(i)) .AND.                     &
                       (cld_top(i)  /= 0)           .AND.                     &
                       (cld_base(i) /= 0))

      ! Check for cloud base/top above/below freezing level
      base_crit(i)  = ((cld_base(i) < frz_lev(i))   .AND.                     &
                       (cld_top(i)  > frz_lev(i)))

      ! Change index of cloud top (Bug fix)
      tp_of_lp(i)   = cld_top(i)
    End Do

  Else 
    ! Original test criteria 
    Do i=1, npnts
      cbct_crit(i) = .TRUE.
      base_crit(i) = ((cld_base(i) < nbl)           .AND.                     &
                      (cld_top(i)  > frz_lev(i)))
      tp_of_lp(i)  = cld_top(i)-1
    End Do

  End If


  ! Locate grid points where depth criteria is satisfied
  If (l_cloud_deep) Then
    Do i=1, npnts
      p_cld_base  = p_lyr_bnds(i,MAX(CLD_BASE(i)-1,0))
      p_cld_top   = p_lyr_bnds(i,CLD_TOP(i))
      dep_crit(i) = (p_cld_base - p_cld_top) >= deep_dp
    End Do
  Else
    Do i=1,npnts
      dep_crit(i) = .TRUE.
    End Do
  End If



  !---------------------------------------------------------------------------
  ! 3.0 Locate grid points where all anvil criteria are satisfied
  !---------------------------------------------------------------------------
  Do i=1, npnts
    If (base_crit(i) .AND. dep_crit(i)) Then
      anv_on(i) = .TRUE.
    End If
  End Do


  !---------------------------------------------------------------------------
  ! 4.0 Apply CCA Profiles 
  !---------------------------------------------------------------------------
  Do i=1, npnts
    If ( cbct_crit(i)       .AND.                                             &
        (cca_2d(i) > 0.0) ) Then

      If (anv_on(i)) Then 

        !---------------------------------------------------------------------
        ! 4.1a Cloud satisfies anvil criteria: Apply Anvil
        !---------------------------------------------------------------------
        Select Case(anv_opt)
          Case(anv_height)
            !-----------------------------------------------------------------
            ! CCA increases with height from freezing level to cloud-top
            !-----------------------------------------------------------------
            anv_lev   = MAX(cld_base(i), frz_lev(i))
            anv_z_dep = z_rho(i,cld_top(i)) - z_rho(i,anv_lev-1)

            Do k=anv_lev, cld_top(i)

              cca_3d(i,k) =                                                   &
                      (anvil_factor - tower_factor)*cca_2d(i)                 &
                    * (z_theta(i,k) - z_rho(i,anv_lev-1)) / anv_z_dep         &
                    + (cca_2d(i) * tower_factor)

              If (cca_3d(i,k) >= 1.0) Then
                cca_3d(i,k) = 0.99
              End If
  
            End Do


          Case(anv_model_levels)
            !-----------------------------------------------------------------
            ! CCA increases with model level from freezing level to cloud-top:
            ! (original code)
            !-----------------------------------------------------------------
            anv_lev = MAX(cld_base(i), frz_lev(i))
            anv_dep = cld_top(i) - anv_lev

            Do k=anv_lev, tp_of_lp(i)

              cca_3d(i,k) =                                                   &
                      (anvil_factor - tower_factor) * cca_2d(i)               &
                    * (k - anv_lev + 1) / anv_dep                             &
                    + (cca_2d(i) * tower_factor)

              If (cca_3d(i,k)  >=  1.0) Then
                cca_3d(i,k) = 0.99
              End If

            End Do


          Case(anv_pressure)
            !-----------------------------------------------------------------
            ! CCA increases with pressure from freezing level to cloud-top
            !-----------------------------------------------------------------
            anv_lev   = MAX(cld_base(i), frz_lev(i))
            anv_p_dep = p_lyr_bnds(i,anv_lev-1) - p_lyr_bnds(i,cld_top(i))

            Do k=anv_lev, tp_of_lp(i)
    
              cca_3d(i,k) =                                                   &
                    (anvil_factor - tower_factor)*cca_2d(i)                   &
                  * (p_lyr_bnds(i,anv_lev-1) -  p_lyr_bnds(i,k)) / anv_p_dep  &
                  + (cca_2d(i) * tower_factor)

              If (cca_3d(i,k) >= 1.0) Then
                cca_3d(i,k) = 0.99
              End If

            End Do

          Case(anv_limited_pressure_depth)
            !-----------------------------------------------------------------
            ! Pressure based, but limit anvil depth to 5000.0 pa
            !-----------------------------------------------------------------
            anv_p_base = p_lyr_bnds(i,cld_top(i)) + 5000.0

            ! Ensure that Anvil is at least 2 levels deep, so only loop to 
            ! layer below cloud top so that it will be 2 levels deep even
            ! if k is set at top of loop
            Do k=1, tp_of_lp(i)-1
              If (p_lyr_bnds(i,k) > anv_p_base) Then
                anv_lev = k
              End If  
            End Do

            anv_p_dep = p_lyr_bnds(i,anv_lev-1) - p_lyr_bnds(i,cld_top(i))

            Do k=anv_lev, tp_of_lp(i)
    
              cca_3d(i,k) =                                                   &
                    (anvil_factor - tower_factor)*cca_2d(i)                   &
                  * (p_lyr_bnds(i,anv_lev-1) -  p_lyr_bnds(i,k)) / anv_p_dep  &
                  + (cca_2d(i) * tower_factor)

              If (cca_3d(i,k) >= 1.0) Then
                cca_3d(i,k) = 0.99
              End If

            End Do

        End Select

        !---------------------------------------------------------------------
        ! 4.1b Cloud satisfies anvil criteria: Apply Tower below anvil base
        !---------------------------------------------------------------------
        Do k=cld_base(i), anv_lev-1
          cca_3d(i,k) = tower_factor * cca_2d(i)
        End Do

      Else

        !---------------------------------------------------------------------
        ! 4.2 No anvil: Apply cca_2d to all model levels from cloud base to
        !     cloud top
        !---------------------------------------------------------------------
        If (l_q_interact .and. .not. lcv_ccrad) Then

          If (l_use_sh_mask) Then

            If (l_pc2_diag_sh_pts(i)) Then
              ! PC2 does not reset cca of shallow convection
              Do k=cld_base(i), tp_of_lp(i)
                cca_3d(i,k) = cca_2d(i)
              End Do
            Else
              Do k=cld_base(i), tp_of_lp(i)
                cca_3d(i,k) = 0.0
              End Do
            End If      ! l_pc2_diag_sh_pts

          Else

            ! PC2 will reset non-anvil cloud to zero
            Do k=cld_base(i), tp_of_lp(i)
              cca_3d(i,k) = 0.0
            End Do

          End If      ! l_use_sh_mask

        Else
          ! Non-PC2 code 
          Do k=cld_base(i), tp_of_lp(i)
            cca_3d(i,k) = cca_2d(i)
          End Do
        End If      ! l_q_interact 

      End If      ! End test on anvil criteria

      !-----------------------------------------------------------------------
      ! Finally check there is no cloud below cloud base or above cloud top:
      ! (original code)
      !-----------------------------------------------------------------------
      If (.not. lcv_ccrad) Then
        Do k=1, (cld_base(i)-1)
          cca_3d(i,k) = 0.0
        End Do

        Do k=cld_top(i), nlev
          cca_3d(i,k) = 0.0
        End Do
      End If      ! lcv_ccrad

    End If      ! cca_2d > 0 and sensible ccb/cct
  End Do      ! loop over npnts


!
!=============================================================================
!  End of anvil calculation
!=============================================================================
!
  RETURN
END SUBROUTINE CALC_3D_CCA
#endif
