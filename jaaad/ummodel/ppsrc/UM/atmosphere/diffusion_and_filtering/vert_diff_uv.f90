
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine vert_diff_uv

      Subroutine vert_diff_uv(                                          &
     &                        u, v, w,                                  &
     &                        r_theta_levels, r_rho_levels,             &
     &                        r_at_u, r_at_v,                           &
     &                        off_x, off_y, halo_i, halo_j,             &
     &                        me, n_proc, model_domain, at_extremity,   &
     &                        timestep, rows, n_rows, row_length,       &
     &                        model_levels, levels,                     &
     &                        start_level, end_level,                   &
     &                        vdiffuv_factor, vdiffuv_test,             &
     &                        R_u, R_v)

! Purpose:
!          Vertical diffusion of u and v based on
!          vertical shear limit
!
! Method:
!          Is described in ;
!
!
! Original Programmer: Terry Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date      Comment
! -------  ------     -------
!  6.2  24/01/06    Code introduced        Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows v.
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, me                                                              &
                      ! Processor id
     &, n_proc                                                          &
                      ! Total number of processors
     &, model_domain                                                    &
                         ! holds integer code for model domain
     &, start_level                                                     &
     &, end_level                                                       &
     &, levels

      Real, Intent(In) ::                                               &
     &  timestep                                                        &
     &, vdiffuv_test                                                    &
                            ! Vertical diffusion test value
     &, vdiffuv_factor     ! Vertical diffusion coefficient scaling

      Real, Intent(In) ::                                               &
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &                 1-halo_j:rows+halo_j, model_levels)              &
     &, r_at_u (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                     &
     &, r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)

      Real, Intent(In) ::                                               &
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     model_levels)                                                &
     &, v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     &     model_levels)                                                &
     &, w (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     0:model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real, Intent(InOut) ::                                            &
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &        model_levels)                                             &
     &, R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &        model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka                                                     &
                                  ! Loop indices
     &, j_start_u, j_stop_u                                             &
                                  ! loop bounds
     &, k_up, k_down, k_switch                                          &
                                  ! level pointers
     &, countu, countv                                                  &
     &, levelsp1                                                        &
     &, info

      Real                                                              &
     &  dvdz_at_u                                                       &
     &, dudz_at_v                                                       &
     &, shear_squ

! Local arrays

      Integer                                                           &
     &  level_sumu(levels)                                              &
     &, level_sumv(levels)                                              &
     &, icu(200), jcu(200), kcu(200)                                    &
     &, icv(200), jcv(200), kcv(200)

      Real                                                              &
     &  r_theta_at_u(1-off_x:row_length+off_x,                          &
     &        1-off_y:rows+off_y, levels + 1)                           &
     &, r_theta_at_v(1-off_x:row_length+off_x,                          &
     &        1-off_y:n_rows+off_y, levels + 1)                         &
     &, fluxu(row_length, rows, 2)                                      &
     &, fluxv(row_length, n_rows, 2)                                    &
     &, dudz(row_length, rows)                                          &
     &, dvdz(row_length, n_rows)

!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the mpp-UM
!
!   Two sets of parameters are set up -
!     i)  for the mpp-UM itself.
!     ii) for the interface to the Message Passing Software.
!
      !=================================================================
      ! Parameters needed for the mpp-UM
      !=================================================================
      ! maximum number of spatial dimensions
      INTEGER,PARAMETER:: Ndim_max = 3 ! 3d data

      ! number of different halo types
      INTEGER,PARAMETER:: NHalo_max = 3 ! for N.D. atmos. model

      INTEGER,PARAMETER:: halo_type_single   = 1
      INTEGER,PARAMETER:: halo_type_extended = 2
      INTEGER,PARAMETER:: halo_type_no_halo  = 3

! FLDTYPE definitions for the different field types recognised on the
! decomposition
      INTEGER,PARAMETER:: Nfld_max=7 ! maximum number of field types
      INTEGER,PARAMETER:: fld_type_p=1       ! grid on P points
      INTEGER,PARAMETER:: fld_type_u=2       ! grid on U points
      INTEGER,PARAMETER:: fld_type_v=3       ! grid on V points
      INTEGER,PARAMETER:: fld_type_comp_wave  = 4
                              ! Compressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_full_wave  = 5
                              ! Uncompressed WAM Wave Field
      INTEGER,PARAMETER:: fld_type_rim_wave   = 6
                              ! Boundary data for WAM Wave Field
      INTEGER,PARAMETER:: fld_type_r=7       ! grid on river points
      INTEGER,PARAMETER:: fld_type_unknown=-1! non-standard grid
! FLDTYPE end

      ! Used in addressing to indicate if calculation is for a local or
      ! global (ie. disk dump) size

      INTEGER,PARAMETER:: local_data=1
      INTEGER,PARAMETER:: global_dump_data=2

      ! maximum permitted size of a halo
      INTEGER,PARAMETER:: Max_Halo_Size=10

      !=================================================================
      ! Parameters needed for the Message Passing Software
      !=================================================================
      INTEGER,PARAMETER:: Maxproc = 512 ! Max number of processors

      ! Processor addresses in the neighbour array
      INTEGER,PARAMETER:: PNorth   = 1
      INTEGER,PARAMETER:: PEast    = 2
      INTEGER,PARAMETER:: PSouth   = 3
      INTEGER,PARAMETER:: PWest    = 4

      ! Value in neighbour array if the domain has  no neighbour in this
      ! direction. Otherwise the value will be the tid of the neighbor
      INTEGER,PARAMETER:: NoDomain = -1

      INTEGER,PARAMETER:: BC_STATIC   = 1 ! Static boundary conditions
      INTEGER,PARAMETER:: BC_CYCLIC   = 2 ! Cyclic boundary conditions
      INTEGER,PARAMETER:: BC_OVERPOLE = 3 ! Transfer over pole
! PARPARM end
! DOMTYP contains different model domain types
!
! Author : P.Burton
! History:
! Version  Date      Comment.
! 5.0      15/04/99  New comdeck
! 5.2      15/11/00  add bi_cyclic_lam domain   A. Malcolm

      INTEGER,PARAMETER:: mt_global        = 1
      INTEGER,PARAMETER:: mt_lam           = 2
      INTEGER,PARAMETER:: mt_cyclic_lam    = 3
      INTEGER,PARAMETER:: mt_bi_cyclic_lam = 4
      INTEGER,PARAMETER:: mt_single_column = 5
! DOMTYP end

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Vertical diffusion of horizontal wind
! ----------------------------------------------------------------------

      levelsp1 = levels + 1
!  Initialise sweep boundaries for filtering and diffusion
      j_start_u = 1
      j_stop_u = rows
      If (at_extremity(PSouth)) then
        j_start_u = 2
      endIf ! at_extremity(PSouth))
      If (at_extremity(PNorth)) then
        j_stop_u = rows - 1
      endIf ! at_extremity(PNorth)

      Do k = start_level - 1, end_level
        ka = k - start_level + 2
        Do j = 1, rows
          Do i = 1, row_length
            r_theta_at_u(i,j,ka) = .5 * (r_theta_levels(i+1,j,k) +      &
     &                                    r_theta_levels(i  ,j,k) )
          End Do
        End Do
        Do j = 1, n_rows
          Do i = 1, row_length
            r_theta_at_v(i,j,ka) = .5 * (r_theta_levels(i,j,k) +        &
     &                                    r_theta_levels(i,j+1,k) )
          End Do
        End Do
      End Do  ! k = start_level - 1, end_level

! call swap_bounds to set halo points
! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
     &                (r_theta_at_u,                                    &
     &                 row_length, rows, levelsp1,                      &
     &                 off_x, off_y, fld_type_u, .false.)
! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
     &                (r_theta_at_v,                                    &
     &                 row_length, n_rows, levelsp1,                    &
     &                 off_x, off_y, fld_type_v, .false.)

! Set downward flux = 0 for start_level
      k_down = 1
      Do j = 1, rows
        Do i = 1, row_length
          fluxu(i,j,k_down) = 0.0
        End Do
      End Do
      Do j = 1, n_rows
        Do i = 1, row_length
          fluxv(i,j,k_down) = 0.0
        End Do
      End Do

      countu = 0
      countv = 0
      k_up = 2
      Do k = start_level, end_level
        ka = k - start_level + 2
        Do j = 1, rows
          Do i = 1-off_x, row_length+off_x
            dudz(i,j) = ( u(i,j,k) - u(i,j,k-1) ) /                     &
     &                ( r_at_u(i,j,k) - r_at_u(i,j,k-1) )
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            dvdz(i,j) = ( v(i,j,k) - v(i,j,k-1) ) /                     &
     &                ( r_at_v(i,j,k) - r_at_v(i,j,k-1) )
          End Do
        End Do


        level_sumu(ka) = 0
        Do j = j_start_u, j_stop_u
          Do i = 1, row_length
            dvdz_at_u = 0.25 * ( dvdz(i,j) + dvdz(i+1,j) +              &
     &                           dvdz(i,j-1) + dvdz(i+1,j-1) )
            shear_squ = dudz(i,j) * dudz(i,j) +                         &
     &                   dvdz_at_u * dvdz_at_u
            if ( shear_squ > vdiffuv_test ) then
              if( countu < 200 ) then
                countu = countu + 1
                level_sumu(ka) = level_sumu(ka) + 1
                icu(countu) = i
                jcu(countu) = j
                kcu(countu) = k
              endif ! countu < 200
              fluxu(i,j,k_up) = vdiffuv_factor * dudz(i,j) *            &
     &                   r_theta_at_u(i,j,ka) * r_theta_at_u(i,j,ka)
            else
              fluxu(i,j,k_up) = 0.0
            endif ! shear_squ > vdiffuv_test
            R_u(i,j,k) = R_u(i,j,k) +                                   &
     &                       ( fluxu(i,j,k_up) - fluxu(i,j,k_down) ) /  &
     &                               ( r_at_u(i,j,k) * r_at_u(i,j,k) *  &
     &                (r_theta_at_u(i,j,ka) - r_theta_at_u(i,j,ka-1)) )
          End Do
        End Do

        level_sumv(ka) = 0
        Do j = 1, n_rows
          Do i = 1, row_length
            dudz_at_v = 0.25 * ( dudz(i,j) + dudz(i+1,j) +              &
     &                            dudz(i,j+1) + dudz(i+1,j+1) )
            shear_squ = dvdz(i,j) * dvdz(i,j) +                         &
     &                   dudz_at_v * dudz_at_v
            if ( shear_squ > vdiffuv_test ) then
              if( countv < 200 ) then
                countv = countv + 1
                level_sumv(ka) = level_sumv(ka) + 1
                icv(countv) = i
                jcv(countv) = j
                kcv(countv) = k
              endif ! countv < 200
              fluxv(i,j,k_up) = vdiffuv_factor * dvdz(i,j) *            &
     &                   r_theta_at_v(i,j,ka) * r_theta_at_v(i,j,ka)
            else
              fluxv(i,j,k_up) = 0.0
            endif ! shear_squ > vdiffuv_test
            R_v(i,j,k) = R_v(i,j,k) +                                   &
     &                       ( fluxv(i,j,k_up) - fluxv(i,j,k_down) ) /  &
     &                               ( r_at_v(i,j,k) * r_at_v(i,j,k) *  &
     &                (r_theta_at_v(i,j,ka) - r_theta_at_v(i,j,ka-1)) )
          End Do
        End Do

      k_switch = k_up
      k_up = k_down
      k_down = k_switch

      EndDo  ! k = start_level, end_level

      If( countu > 0 ) then
        write(6,*)' Large vertical u shear found at ', countu           &
     &        ,' points on processor ', me
      else   ! countu = 0
      write(6,*)'No points found having large vertical u shear on proc '&
     &        , me
      endIf  !  countu > 0
      If( countv > 0 ) then
        write(6,*)' Large vertical v shear found at ', countv           &
     &        ,' points on processor ', me
      else   ! countv = 0
      write(6,*)'No points found having large vertical v shear on proc '&
     &        , me
      endIf  !  countv > 0

! check count > 0 somewhere
      call gc_isum( levels, n_proc, info, level_sumu )
      call gc_isum( levels, n_proc, info, level_sumv )

      countu = 0
      countv = 0
      Do k = start_level, end_level
        ka = k - start_level + 2
        if( level_sumu(ka) > 0) then
          countu = countu + 1
          if( me == 0 ) then
            write(6,*)'  Large vertical u shear found on level ', k
          endif ! me == 0
        endif  !  level_sumu(ka) > 0
        if( level_sumv(ka) > 0) then
          countv = countv + 1
          if( me == 0 ) then
            write(6,*)'  Large vertical u shear found on level ', k
          endif ! me == 0
        endif  !  level_sumv(ka) > 0
      enddo !   k = start_level, end_level
      if( countu == 0 ) then
        if( me == 0 ) then
          write(6,*)' No large vertical u shear found above level '     &
     &          ,start_level - 1
        endif ! me == 0
      endif ! countu == 0
      if( countv == 0 ) then
        if( me == 0 ) then
          write(6,*)' No large vertical v shear found above level '     &
     &          ,start_level - 1
        endif ! me == 0
      endif ! countv == 0

      return     ! End of routine vert_diff_uv
      END SUBROUTINE vert_diff_uv

