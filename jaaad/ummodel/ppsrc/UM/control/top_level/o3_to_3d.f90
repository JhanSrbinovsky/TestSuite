
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to expand ozone from the ancillary to a full field.
!     Subroutine Interface:
      Subroutine O3_to_3D(lexpand_ozone, i_ozone_int,                   &
     &  rows, row_length, model_levels, ozone_levels,                   &
     &  halo_i, halo_j, off_x, off_y, at_extremity,                     &
     &  z_top_of_model,                                                 &
     &  theta, r_theta_levels, eta_theta_levels, exner_theta_levels,    &
     &  rho,   r_rho_levels,   eta_rho_levels,   exner_rho_levels,      &
     &  nd_o3, ozone_in,                                                &
     &  L_use_stochem_O3, nd_stochem, o3_stoch,                         &
     &  min_trop_level, max_trop_level,                                 &
     &  L_O3_trop_level,L_O3_trop_height,L_T_trop_level,L_T_trop_height,&
     &  O3_trop_level,O3_trop_height,T_trop_level,T_trop_height,        &
     &  proc_row_group,                                                 &
     &  global_row_length,                                              &
     &  ozone3D,                                                        &
     &  ErrorStatus, cmessage                                           &
     &  )
!
!
!
      Implicit None
!
! Description:
!   This routine takes the ozone field supplied in the ancillary
!   and converts it to a full 2-D (zonal mean) or 3-D field
!   as directed by the option selected.
!
! Method:
!   Essentially, the routine must carry out vertical interpolation
!   of the ozone field supplied, typically from an ancillary file,
!   on to the vertical levels at a grid-point. Various options are
!   permitted, the newer ones allowing interpolation in height to
!   match the new vertical structure of the model.
!      This code is largely developmental and will subsequently be
!   combined with changes and extensions to the ancillary system
!   for ozone.
!      The ozone concentrations are found on theta levels,
!   and the ozone tropopause is found on a rho level.
!   Rho level 2 is found between theta levels 1 and 2.
!
! Current Code Owner: J. M. Edwards
!
! History:
! Version  Date     Comment
! =======  ====     =======
! 5.2      28/11/00 Original Code.
!                   (J. M. Edwards)
!
! 5.3      28/09/01 Modified to calculate the ozone tropopause
!                   from the ozone concentrations in the ancillary file.
!                   The criteria in Bethan et al. 1996, QJR 122,
!                   pp. 929-944 are used.
!                                      E. Ostrom
!
! 5.5      24/02/03 Generate diagnostics and map ozone tropopause
!                   onto thermal tropopause keeping ozone column
!                   mass constant.                (J.-C. Thelen)
!  6.0   18/06/03 Fix control structure for ozone in
!                 in o3_to_3d         (J.-C. Thelen)
!  6.1   24/08/04  NEC supplied optimisations S. WIlson
!  6.1   20/08/03  Code for STOCHEM feedback.  C. Johnson
!  6.2   20/03/06  SX-6 optimisation by J-C Rioual. Lodged by M Saunby.
!  6.2   02/02/06 Modify argument list in call to tropin with SCM
!                 dummy variables to retain consistency in
!                 argument list.                       R. Wong

! Code description:
!   FORTRAN 90
!   This code is written to the programming standards of version 6
!   of UMDP3.
!
!     ----------------------------------------------------------------
!
!
!     Input arguments:
!
      Integer, intent(IN) :: proc_row_group, global_row_length
! Specify the total number of gridpoints in the east-west direction
! and which processor takes care of which row.
      Logical, intent(IN) :: lexpand_ozone
!       Flag for expansion of ozone from the ancillary
      Integer, intent(IN) :: i_ozone_int
!         Method of expanding the ozone in the dump to a full 3-D field
      Integer, intent(IN) :: rows
!         Number of EW rows
      Integer, intent(IN) :: row_length
!         Number of points on each row
!
      Integer, intent(IN) :: halo_i
!         Size of large EW Halo
      Integer, intent(IN) :: halo_j
!         Size of large NS Halo
      Integer, intent(IN) :: off_x
!         Size of EW Halo
      Integer, intent(IN) :: off_y
!         Size of NS Halo
      Logical, intent(IN) :: at_extremity(4)
!         Flags to indicate whether the processor is at the boundary
!         of the domain
!
      Logical, intent(IN) :: L_O3_trop_level, L_O3_trop_height
      Logical, intent(IN) :: L_T_trop_level, L_T_trop_height
      Logical, intent(IN) :: L_use_stochem_O3
! Flags indicating which diagnostics are on or off.
      Integer, intent(IN) :: model_levels
!         Number of vertical levels in the atmosphere
      Integer, intent(IN) :: ozone_levels
!         Number of levels on which ozone is specified: these are
!         contiguous levels reaching to the top of the model. If
!         used with the options IO3_3DSPEC or IO3_2DSPEC the value
!         on the lowest level of this array is used to set the mixing
!         ratio on lower layers of the model. Its interpretation under
!         other options is not yet finalized (but note that it is
!         used in radiation, so any changes here will need to be
!         accounted for there).
      Integer, intent(IN) :: nd_o3
!         Size of the array of ozone supplied from D1
      Integer, intent(IN) :: nd_stochem
!         Size of the array of ozone supplied from STOCHEM
!
      Integer, intent(IN) :: min_trop_level
!         Minimum permitted level of the tropopause
      Integer, intent(IN) :: max_trop_level
!         Maximum permitted level of the tropopause
!
      Real, Intent(IN)    :: z_top_of_model
!         Height of the top of the model
!
      Real, Intent(IN)    :: r_theta_levels(1-halo_i:row_length+halo_i, &
     &                                      1-halo_j:rows+halo_j,       &
     &                                      0: model_levels)
!         Radial coordinates on theta levels
      Real, Intent(IN)    :: exner_theta_levels(                        &
     &                         1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y,                      &
     &                         model_levels)
!         Exner function on theta levels
      Real, Intent(IN)    :: eta_theta_levels(0: model_levels)
!         Eta values on theta levels
      Real, Intent(IN)    :: theta(1-off_x:row_length+off_x,            &
     &                             1-off_y:rows+off_y,                  &
     &                             model_levels)
!         Potential temperatures on theta-levels
!
      Real, Intent(IN)    :: r_rho_levels(1-halo_i:row_length+halo_i,   &
     &                                    1-halo_j:rows+halo_j,         &
     &                                    model_levels)
!         Radial coordinates of rho levels
      Real, Intent(IN)    :: exner_rho_levels(                          &
     &                         1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y,                      &
     &                         model_levels)
!         Exner function on rho levels
      Real, Intent(IN)    :: eta_rho_levels(1: model_levels)
!         Eta values on rho levels
!
      Real, Intent(IN)    :: rho(1-off_x:row_length+off_x,              &
     &                           1-off_y:rows+off_y,                    &
     &                           model_levels)
!         Atmospheric densities on rho-levels
!
      Real, Intent(IN)    :: ozone_in(nd_o3)
!         Input ozone field in D1
      Real, Intent(IN)    :: o3_stoch(nd_stochem)
!         Input ozone field in D1 from STOCHEM
!
!
! SCM Dummy variables to keep call to tropin consistent.
       REAL                                                             &
     & scm_dummy_1d(1,1)                                                &
     &,scm_dummy_2d(1,1,0:model_levels)

!     Output arguments:

      Real, Intent(OUT)   :: ozone3D(row_length, rows, ozone_levels)
!         Expanded ozone field

!
      Real, intent(OUT) :: T_trop_level(row_length,rows)
      Real, intent(OUT) :: O3_trop_level(row_length,rows)
!Points to the lower boundary of the layer containing the thermal and
!ozone tropopause.
      Real, intent(OUT) :: O3_trop_height(row_length,rows)
      Real, intent(OUT) :: T_trop_height(row_length,rows)
!Height of the ozone and thermal tropopause.
!
!     Error Status:
      Integer, Intent(INOUT):: ErrorStatus
!         Error code
      Character*(*), Intent(INOUT) :: cmessage
!         Short error message
!
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
!+ ---------------------------------------------------------------------
!  Module to specify allowed methods of interpolating from the
!  ancillary file.
!
!  Current Code Owner: J. M. Edwards
!
!  History:
!
!  Version  Date      Comment.
!  5.2      14/11/00  Original code.
!                     (J. M. Edwards)
!
!- ---------------------------------------------------------------------
!
      Integer, parameter :: IO3_3DSPEC = 1
!       Ozone is provided as a full 3D field.
      Integer, parameter :: IO3_2DSPEC = 2
!       Ozone is expanded from a 2D field by direct copying.
      Integer, parameter :: IO3_2DMASSCON = 3
!       Ozone is expanded from a 2D field with conservation of mass.
      Integer, parameter :: IO3_TROP_MAP = 4
!       Ozone mixing ratios are set by mapping each height at each
!       grid-point in the real profile to a height in the ancillary
!       profile and using the mixing ratio there. The mapping is
!       set using the height of the tropopause
      Integer, parameter :: IO3_TROP_MAP_MASSCON = 5
!       Ozone mixing ratios are set as above, but scaled so as to
!       preserve the vertically integrated column ozone
!
! ----------------------------------------------------------------------
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! ----------------------- Header file O3CRITS  ------------------------
! Description: Parameters of ozone tropopause criteria.
!
! Current Code Owner: E. Ostrom
!
! History:
! Version  Date      Comment.
!  5.3     28/09/01  Original Code.   E. Ostrom
!
!----------------------------------------------------------------------
!     ! Criteria for establishing the level of the
!     ! ozone tropopause to enable mapping of ozone concentrations
!     ! in ozone ancillary file to model levels using the
!     ! thermal tropopause level in the model for scaling.

!     ! The gradient of O3 concentration has to exceed 60 ppbv per km
!     ! i.e. 99.423E-09 (= 60ppbv) / 1000 to convert to conc per meter
      Real, Parameter :: O3_grad_crit    = 99.423E-12

!     ! The concentration of O3 has to exceed 80 ppbv
      Real, Parameter :: O3_conc_crit    = 132.56E-09

!     ! The conc. of O3 has to exceed 110 ppbv in the stratosphere
      Real, Parameter :: O3_strat_crit   = 182.27E-09
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!     Local variables:
      Integer :: i
!         Loop variable for Latitude
      Integer :: j
!         Loop variable for Longitude
      Integer :: k
!         Loop variable for vertical
      Integer :: jk
!         Loop variable
      Integer :: jkp1
!         Loop variable
      Integer :: ijk
!         Loop variable
      Integer :: k_o3_trop(row_length)
!         Index of theta level just below ozone tropopause
      logical :: k_o3_trop_done(row_length)
      Integer :: k_anc(row_length)
!         Index of immediate theta level below the interpolated height
      Integer :: k_off
!         Offset of the ozone levels from the model levels
      Integer :: lvl, lsl, lso, nv, err_msg
! Variables needed for the call to GCG_RVECSUMR. For details see the
! documentation.
      Integer :: trindx(row_length, rows)
!         Points to the lower boundary of the layer containing the
!         thermal tropopause.
      Real :: height_above_surf
!         Height above the surface of the current grid-point
!         (on theta level).
      Real :: z_trop_O3(row_length)
!         Height of the ozone tropopause above the surface in the
!         ancillary profile.
      Real :: z_trop_pt(row_length)
!         Height of the thermal tropopause above the surface at the
!         current grid-point.
      Real :: z_int(row_length)
!         Interpolated height in the ancillary profile (on no level).
      Real :: T_filled_halo(1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y,                         &
     &                      model_levels)
!         Temperatures on theta levels with haloes added
      Real :: rho_anc(model_levels)
      Real :: layer_mass(row_length,rows)
!Mass on rho_levels
      Real :: avg_layer_mass(rows,model_levels)
!Mass on rho_levels averaged over longitude.
      Real :: ww1(rows)
!Working Array
!         Atmospheric densities on rho-levels in the ancillary profile
      Real :: o3_mass_anc(0: model_levels)
!         Column integrated mass of ozone in the ancillary profile
!         on rho levels, the unit is kg per kg per m^2
      Real :: o3_mass_cumul(row_length)
!         Column mass of ozone at the grid-point integrated from
!         the surface to the top of the current layer
      Real :: o3_mass_cumul_below(row_length)
!         Column mass of ozone at the grid-point integrated from
!         the surface to the top of the layer below the current one
      Real :: o3_gradient_anc(row_length, rows, model_levels)
!         Ozone concentration gradient on rho-levels in ancillary file.
!
!     Output arguments:
!- End of header
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!
! Additional loop variable for vectorization

      Integer :: kl
!NEC optimisation of the initialisation of k_anc array
!
!
!
!     Preliminary calculations:-
!
      If ( (i_ozone_int  ==  IO3_TROP_MAP).OR.                          &
     &     (i_ozone_int  ==  IO3_TROP_MAP_MASSCON) ) then
!
!       The (thermal) tropopause is required: this code is
!       copied from RAD_CTL, T_n being replaced by
!       its definition.
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              T_filled_halo(i,j,k)                                      &
     &          =theta(i,j,k)*exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
!
! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &    T_filled_halo, row_length, rows,                              &
     &    model_levels, off_x, off_y, fld_type_p, .false.)
! DEPENDS ON: fill_external_halos
        Call Fill_external_halos(T_filled_halo,row_length,rows,         &
     &    model_levels,off_x,off_y)
!
! DEPENDS ON: tropin
        Call tropin (T_filled_halo,                                     &
     &                 exner_rho_levels, exner_theta_levels,            &
     &                 row_length, rows, model_levels, off_x, off_y,    &
     &                 at_extremity,scm_dummy_1d,scm_dummy_2d,          &
     &                 min_trop_level, max_trop_level, trindx )
!

!     ! Provisional code: the ancillary density is set to a value of 1
!     ! until established.
      lvl = row_length
      lsl = row_length
      lso = 1
      nv  = rows

! Calculate the mass in each layer between rho-levels.

      Do k = 1,model_levels

         Do j = 1, rows
            Do i = 1,row_length
               if (k == 1) then
                  layer_mass(i,j)=rho(i,j,k+1)                          &
     &               *(r_rho_levels(i,j,k+1)-r_theta_levels(i,j,k-1))
               endif
               if ((k >  1).and.(k <  model_levels)) then
                  layer_mass(i,j)=0.5*(rho(i,j,k+1)+rho(i,j,k))         &
     &               *(r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
               endif
               if (k == model_levels) then
                  layer_mass(i,j)=rho(i,j,k)                            &
     &               *(r_theta_levels(i,j,k)-r_rho_levels(i,j,k))
               endif
            end do
         end do

! Average the mass per layer over longitude.

         call gcg_rvecsumr(lvl, lsl, lso, nv, layer_mass,               &
     &                  proc_row_group,err_msg,ww1)

         do j=1, rows
            avg_layer_mass(j,k)=ww1(j)/global_row_length
         end do
      end do
!
      Endif
!
!======================================================
! Set SCM dummy values to zero
       scm_dummy_1d(:,:)   = 0.0
       scm_dummy_2d(:,:,:) = 0.0
!======================================================
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!
      Select Case(i_ozone_int)
!
      Case(IO3_3DSPEC)
!
!       In this case we should have a 3D ancillary file so
!       we return an error message if not.
        If (lexpand_ozone) then
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("O3_to_3D", ErrorStatus,                         &
     &       "A 2D ozone ancillary has been specified with a" //        &
     &       "3D ozone option.")
        end if
!
!       The ozone field supplied is the intended 3-D field and is
!       simply copied across. (N.B. This is inefficient: in the future
!       we should be able to avoid this with a pointer and an
!       allocatable array.)
        Do k = 1, ozone_levels
          Do j = 1, rows
            Do i = 1, row_length
              ijk = i + (j-1)*row_length + (k-1)*row_length*rows
              ozone3D(i,j,k)=ozone_in(ijk)
            end do
          end do
        end do
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
!
      Case (IO3_2DSPEC)
!
!       In this case we should have a 2D ancillary file so
!       we return an error message if not.
        If (.NOT.lexpand_ozone) then
          ErrorStatus=132
! DEPENDS ON: ereport
          Call Ereport("O3_to_3D", ErrorStatus,                         &
     &       "A 3D ozone ancillary has been specified with a" //        &
     &       "2D ozone option.")
        end if
!
!       The ozone is copied round a latitude circle. This does not
!       conserve the total amount of ozone in the atmosphere because
!       it does not allow for orography.
        IF (L_use_stochem_O3) THEN
! Use STOCHEM O3, this is a 3-D field
          DO k = 1, ozone_levels
            DO j = 1, rows
              DO i = 1, row_length
                ijk = i + (j-1)*row_length + (k-1)*row_length*rows
                ozone3D(i,j,k)=o3_stoch(ijk)
              END DO
            END DO
          END DO
        ELSE
        Do k = 1, ozone_levels
          Do j = 1, rows
            Do i = 1, row_length
              jk = j + (k-1)*rows
              ozone3D(i,j,k)=ozone_in(jk)
            end do
          end do
        end do
        ENDIF
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
      Case (IO3_2DMASSCON)
!
! ! !   Since z_top_of_model * eta levels doesn't work in new dynamics
! ! !   This option is not possible and this option will be identical
! ! !   with the option above.
!
!       The ozone field is expanded to three dimensions in such a way
!       that mass-loading of each atmospheric layer is the same at
!       each grid-point on a latitude circle, regardless of the
!       orography.
!
!       The ozone levels do not start at the bottom of the model.

        k_off=model_levels-ozone_levels
!
!       In this provisional code the difference in density between
!       the ancillary profile and the profile at the grid-point at
!       the same level is ignored.
!
        Do k = 1, ozone_levels
          Do j = 1, rows
            Do i = 1, row_length
              jk = j + (k-1)*rows
              ozone3D(i,j,k)=ozone_in(jk)
!     &          * z_top_of_model
!     &          * (eta_theta_levels(k+k_off)
!     &           - eta_theta_levels(k-1+k_off))
!     &          / (r_theta_levels(i, j, k+k_off)
!     &           - r_theta_levels(i, j, k-1+k_off))

            end do
          end do
        end do
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
      Case (IO3_TROP_MAP)
!
!       Still unoptimizied code: this option has been
!       developed, but might not be the definite solution.
!
        Do j=1, rows
!

!          ! Initialize the pointer to theta level just below the ozone
!          ! tropopause level.
          k_o3_trop(:) = 0
          k_o3_trop_done(:) = .false.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
          Do k=min_trop_level, max_trop_level
!           ! Compute the ozone gradient in the ancillary file.
!           ! In order to access data on row 4, you need to index
!           ! j + 3 i.e. j + (k-1)*rows, hence the strange indexing
              jk   = j + (k-1)*rows
              jkp1 = j + (k)  *rows

            Do i=1, row_length
!           ! min_trop_level should be a rho-level, and rho levels
!           ! with the same index as theta-levels are higher up.

!           ! Set the pointer to the current rho-level. When the
!           ! tropopause criteria are met, this will be saved and will
!           ! indicate the rho-level, between two ozone levels that
!           ! contains the tropopause.
              if(.not.k_o3_trop_done(i)) then
                k_o3_trop(i) = k
!
                o3_gradient_anc(i,j,k) =                                &
     &              ( ozone_in(jkp1) - ozone_in(jk) ) /                 &
     &              ( r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
!
!           ! test whether we have found the ozone tropopause !
                If ( o3_gradient_anc(i,j,k) > O3_grad_crit .AND.        &
     &             ozone_in(jk)           > O3_conc_crit .AND.          &
     &             ozone_in(jkp1)         > O3_strat_crit     )         &
     &             k_o3_trop_done(i) = .true.
              endif
            end do
          end do


!          ! Calculate the height of the ozone tropopause above the
!          ! surface. This is assumed to lie on a rho-level, and
!          ! is used for the scaling of the ozone concentration
!          ! in the model.
!
          Do i=1, row_length
            z_trop_O3(i) = r_rho_levels(i, j, k_o3_trop(i))             &
     &                - r_theta_levels(i, j, 0)

!
             if (L_O3_trop_height) then
                O3_trop_height(i,j)=z_trop_O3(i)
             endif

             if (L_O3_trop_level) then
                O3_trop_level(i,j)=real(k_o3_trop(i))
             endif
!
!           Set the thermal tropopause at the grid-point. Note that as
!           TRINDX is counted upward from the surface, while the first
!           rho level is omitted from the physics, no offset to the
!           final index of r_rho_levels is required.
            z_trop_pt(i) = r_rho_levels(i, j, trindx(i, j))             &
     &                - r_theta_levels(i, j, 0)
!

             if (L_T_trop_height) then
                T_trop_height(i,j)=z_trop_pt(i)
             endif

             if (L_T_trop_level) then
                T_trop_level(i,j)=real(trindx(i,j))
             endif
          EndDo
!
!           Initialize the pointer that points to the lower boundary of
!           the layer containing the height we are currently inter-
!           polating ozone concentrations for.
!          k_anc(:) =0
!
          Do k=1, model_levels
            Do i=1, row_length
!
!           ! Find the height above the surface to the center of the
!           ! current model layer. I.E. find the height to the theta
!           ! level we want ozone on.
              height_above_surf                                         &
     &          = r_theta_levels(i, j, k) - r_theta_levels(i, j, 0)
!
              If (height_above_surf  <   z_trop_pt(i)) then
!               We are in the thermal troposphere.
!               Calculate the corresponding height in the ancillary
!               file using similarity.
                z_int(i) = height_above_surf                            &
     &            * z_trop_O3(i) / z_trop_pt(i)
              Else
!               We are in the thermal stratosphere.
                z_int(i) = z_trop_O3(i)                                 &
     &            + (height_above_surf - z_trop_pt(i))                  &
     &            * (z_top_of_model - z_trop_O3(i))                     &
     &            / (z_top_of_model - z_trop_pt(i))
              Endif
            EndDo
!
!           ! Increment k_anc until it points to the theta level in
!           ! the ancillary profile just below the interpolated height.
!
!          Do i=1, row_length
!              Do While ( (r_theta_levels(i,j,k_anc(i))
!     &                  - r_theta_levels(i, j, 0)  )   <  
!     &                   z_int(i) .AND. k_anc(i)  <=  model_levels )
!                k_anc(i) = k_anc(i) + 1
!              Enddo
!!             In order to point to theta level below int. height.
!              k_anc(i) = k_anc(i) - 1
!            EndDo

             k_anc(:)=-1
             do kl=0, model_levels
               Do i=1, row_length
                if ( (                                                  &
     &              ( r_theta_levels(i,j,kl)-r_theta_levels(i, j, 0))   &
     &               >=   z_int(i) )                                    &
     &              .and. ( k_anc(i)  ==  -1 ) ) then
                      k_anc(i)=kl-1
                end if
               end do
             end do

             where (k_anc==-1) k_anc=model_levels

          Do i=1, row_length
!
!             Interpolate within this layer.
              jk=j+(k_anc(i)-1)*rows

              IF      ( k_anc(i) == 0 ) THEN
                ozone3D(i, j, k) = ozone_in( j)
              ELSE IF ( k_anc(i) == model_levels ) THEN
                ozone3D(i, j, k) = ozone_in(jk)
              ELSE
              ozone3D(i, j, k) = ozone_in(jk)                           &
     &          + ( ozone_in(jk+rows) - ozone_in(jk) )                  &
     &          *  (z_int(i) - ( r_theta_levels(i,j,k_anc(i))           &
     &                      - r_theta_levels(i,j,0)     ) )             &
     &          / (r_theta_levels(i,j,k_anc(i)+1)                       &
     &               -  r_theta_levels(i,j,k_anc(i)) )
              ENDIF
            Enddo
!
          Enddo
!
        Enddo
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
      Case(IO3_TROP_MAP_MASSCON)
!
!       In this case we should have a 3D ancillary file so
!       we return an error message if not.
        If (.not.lexpand_ozone) then
          ErrorStatus=123
! DEPENDS ON: ereport
          Call Ereport("O3_to_3D", ErrorStatus,                         &
     &       "A 3D ozone ancillary has been specified with a" //        &
     &       "2D ozone option.")
        end if
!
!       Provisional unoptimizied code: these options are still being
!       developed.
!
!       Interpolation using cumulative ozone amounts is preferred to
!       conserved column masses.
!
        Do j=1, rows
!
!         Calculate cumulative amounts of ozone to the top of each
!         layer in the field supplied (i.e. rho level).
          o3_mass_anc(0)=0.0
!
!         When ozone_levels is less than model_levels,
!         the lowest levels all have the same ozone mass mixing ratio.
!

!
          Do k=1, model_levels-ozone_levels
            o3_mass_anc(k) = o3_mass_anc(k-1)                           &
     &        + ozone_in(j) * avg_layer_mass(j,k)
          Enddo
          Do k=model_levels-ozone_levels+1, model_levels
            jk=j+(k-1)*rows
            o3_mass_anc(k) = o3_mass_anc(k-1)                           &
     &        + ozone_in(jk) * avg_layer_mass(j,k)

          Enddo
!
!
!
!          ! Initialize the pointer to theta level just below the ozone
!          ! tropopause level.
          k_o3_trop(:) = 0
          k_o3_trop_done(:) = .false.
!
          Do k=min_trop_level+1, max_trop_level

!           ! Compute the ozone gradient in the ancillary file.
!           ! In order to access data on row 4, you need to index
!           ! j + 3 i.e. j + (k-1)*rows, hence the strange indexing

              jk   = j + (k-1)*rows
              jkp1 = j + (k)  *rows

!cdir novector
            Do i=1, row_length

!           ! min_trop_level should be a rho-level, and rho levels
!           ! with the same index as theta-levels are higher up.
!           !
!           ! Set the pointer to the current rho_level. When the
!           ! tropopause criteria are met, this will be saved and will
!           ! indicate the rho-level, between two ozone levels that
!           ! contains the tropopause.

              if(.not.k_o3_trop_done(i)) then
                k_o3_trop(i) = k
!
                o3_gradient_anc(i,j,k) =                                &
     &              ( ozone_in(jkp1) - ozone_in(jk) ) /                 &
     &              ( r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k))
!
!           ! test whether we have found the ozone tropopause !
                If (o3_gradient_anc(i,j,k) > O3_grad_crit .AND.         &
     &            ozone_in(jk)           > O3_conc_crit .AND.           &
     &            ozone_in(jkp1)         > O3_strat_crit     )          &
     &            k_o3_trop_done(i) = .true.
              endif
            end do
          end do

          Do i=1, row_length

!          ! Calculate the height of the ozone tropopause above the
!          ! surface. This is assumed to lie on a rho-level, and
!          ! is used for the scaling of the ozone concentration
!          ! in the model.
!
            z_trop_O3(i) = r_rho_levels(i, j, k_o3_trop(i))             &
     &                - r_theta_levels(i, j, 0)

            if (L_O3_trop_height) then
              O3_trop_height(i,j)=z_trop_O3(i)
            endif

            if (L_O3_trop_level) then
              O3_trop_level(i,j)=real(k_o3_trop(i))
            endif

!           Set the thermal tropopause at the grid-point. Note that as
!           TRINDX is counted upward from the surface, while the first
!           rho level is omitted from the physics, no offset to the
!           final index of r_rho_levels is required.
            z_trop_pt(i) = r_rho_levels(i, j, trindx(i, j))             &
     &                - r_theta_levels(i, j, 0)
!
!
            if (L_T_trop_height) then
              T_trop_height(i,j)=z_trop_pt(i)
            endif

            if (L_T_trop_level) then
              T_trop_level(i,j)=real(trindx(i,j))
            endif
          EndDo

!           Initialize the pointer, that points to the lower boundary
!           of the layer containing the height we are currently inter-
!           polating ozone concentrations for.
!
!          k_anc(:)=1
!
!           At the surface there is no ozone below the current level.
          o3_mass_cumul_below(:)=0.0
!
          Do k=1, model_levels

            if ((k >  1).and.(k <  model_levels)) then

            Do i=1, row_length
!
!           ! Find the height above the surface to the center of the
!           ! current model layer. I.E. find the height to the theta
!           ! level we want ozone on.

              height_above_surf                                         &
     &          = r_rho_levels(i, j, k+1) - r_theta_levels(i, j, 0)
              If (height_above_surf  <=  z_trop_pt(i)) then
!           !   We are in the thermal troposphere.
!           !   Calculate the corresponding height in the ancillary
!           !   file using similarity.
                z_int(i) = height_above_surf                            &
     &            * z_trop_O3(i) / z_trop_pt(i)
              Else
!           !   We are in the thermal stratosphere.
                z_int(i) = z_trop_O3(i)                                 &
     &            + (height_above_surf - z_trop_pt(i))                  &
     &            * (z_top_of_model - z_trop_O3(i))                     &
     &            / (z_top_of_model - z_trop_pt(i))
              Endif
            EndDo
!
!           ! Increment k_anc until it points to the theta level in
!           ! the ancillary profile just below the interpolated height.
!
!            Do i=1, row_length
!
!              Do While ( ( r_theta_levels(i,j,k_anc(i))
!     &                   - r_theta_levels(i,j,0)   )  <  
!     &                   z_int(i) )
!                if(k_anc(i) <  model_levels) then
!                  k_anc(i) = k_anc(i) + 1
!                else
!                  exit
!                endif
!              Enddo
!!!             In order to point to theta level below int. height.
!              k_anc(i) = k_anc(i) - 1
!            EndDo


             k_anc(:)=-1
             do kl=0, model_levels
               Do i=1, row_length
                if ( (                                                  &
     &              ( r_theta_levels(i,j,kl)-r_theta_levels(i, j, 0))   &
     &               >=   z_int(i) )                                    &
     &              .and. ( k_anc(i)  ==  -1 ) ) then
                      k_anc(i)=kl-1
                end if
               end do
             end do

             where (k_anc==-1) k_anc=model_levels




            endif
!
!             Interpolate within this layer.

            Do i=1, row_length

              if ((k >  1).and.(k <  model_levels)) then

!             Calculate the o3_mass_cumul at z_int using the o3_mass_anc
!             at the theta levels below and above z_int
!
                o3_mass_cumul(i) = o3_mass_anc(k_anc(i)-1)              &
     &          + ( o3_mass_anc(k_anc(i)) - o3_mass_anc(k_anc(i)-1) )   &
     &          *   (z_int(i) - ( r_rho_levels(i,j,k_anc(i))            &
     &                       - r_theta_levels(i,j,0)   ) )              &
     &            / ( r_rho_levels(i,j,k_anc(i)+1)                      &
     &              - r_rho_levels(i,j,k_anc(i))   )
!
              endif
              if (k == 1) then
                o3_mass_cumul(i)=o3_mass_anc(1)
              endif
              if (k == model_levels) then
                o3_mass_cumul(i)=o3_mass_anc(model_levels)
              endif
!
!           ! Calculate the ozone contentration on the theta level
!           ! at the model gridpoint, using
!           ! qO3 = delta_QO3 / (delta_rho * delta_z)
!
              if (k == 1) then
                ozone3D(i, j, k) =                                      &
     &                  (o3_mass_cumul(i) - o3_mass_cumul_below(i)) *   &
     &                  1. /( rho(i, j, k+1)                            &
     &                  * (r_rho_levels(i, j, k+1)                      &
     &                  - r_theta_levels(i, j, k-1) ) )
              endif
              if ((k >  1).and.(k <  model_levels)) then
                ozone3D(i, j, k) =                                      &
     &                  (o3_mass_cumul(i) - o3_mass_cumul_below(i)) *   &
     &                  2. / ( (rho(i, j, k) + rho(i, j, k+1) )         &
     &                  * (r_rho_levels(i, j, k+1)                      &
     &                  - r_rho_levels(i, j, k) ) )
!
              endif
              if (k == model_levels) then
                ozone3D(i, j, k) =                                      &
     &                  (o3_mass_cumul(i) - o3_mass_cumul_below(i)) *   &
     &                  1.0/( rho(i, j, k)                              &
     &                  * (r_theta_levels(i ,j ,k)                      &
     &                  - r_rho_levels(i, j, k) ) )
              endif

              o3_mass_cumul_below(i)=o3_mass_cumul(i)
!
            Enddo
!
          Enddo
!
        Enddo
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~72
      Case Default
!
        cmessage='*** Error: Unrecognized expansion of ozone.'
        ErrorStatus=123
        Return
!
      End Select
!

      Return
      END SUBROUTINE O3_to_3D
