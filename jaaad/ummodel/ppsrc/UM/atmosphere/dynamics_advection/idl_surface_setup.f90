
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_Surface_setup

      Subroutine IDL_Surface_setup(                                     &
     &                      Earth_radius, Pi                            &
     &,                     row_length, rows, model_levels              &
     &,                     global_row_length, global_rows              &
     &,                     halo_i, halo_j                              &
     &,                     me, n_proc, at_extremity, model_domain      &
     &,                     l_datastart, all_proc_group                 &
     &,                     delta_lambda, delta_phi, Base_phi           &
     &,                     n_rows, base_lambda                         &
!  VarRes Grid Spacing
     &,                     lambda_p, phi_p, lambda_u, phi_v            &
     &,                     lambda_p_end, phi_p_end                     &
     &,                     L_regular                                   &
     &,                     delta_x, delta_y                            &
     &,                     orography, orog_haloes                      &
     &,                     L_fix_orog_hgt_lbc, orog_hgt_lbc, rimwidth  &
     &,                     surface_type                                &
     &,                     h_o, lambda_fraction, phi_fraction          &
     &,                     half_width_x, half_width_y                  &
     &,                     xp, yp, Witch_power                         &
     &,                     L_code_test)

! Purpose:
!          Sets up surface
!          surface_number defines options
!
! Original Programmer:  T. Davies
! Current code owner: Andrew J. Malcolm
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   12/10/01   Original code.                      T.Davies
!   6.1   04/10/04   fixes for EW-cyclic LAM             A.Malcolm
!   6.2   31/01/06  Add option to set constant orog in lbc zone. YMTang
!  6.2   31/12/05  Variable resolution changes      Terry Davies
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      Implicit None

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a processor row
     &, rows                                                            &
                         ! number of rows in a processor theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, n_rows                                                          &
                         ! number of rows in a processor v field
     &, rimwidth         ! IN : Size of boundary region

      Real                                                              &
     &  orography(row_length,rows)                                      &
     &, orog_haloes(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
!    orog_haloes is r_theta_levels(0)
     &, h_o                                                             &
                              ! mountain height
     &, orog_hgt_lbc                                                    &
                              ! Height of orog in lbc zone
     &, lambda_fraction                                                 &
                              ! hill position as fraction of domain
     &, phi_fraction                                                    &
                              ! hill position as fraction of domain
     &, half_width_x                                                    &
                             ! hill half-width East-West
     &, half_width_y                                                    &
                             ! hill half-width North-South
     &, xp                                                              &
                         ! Plateau EW half-width  (surface_type 4)
     &, yp                                                              &
                         ! Plateau NS half-width  (surface_type 4)
     &, Witch_power                                                     &
                         ! Exponent power in Witch definition
     &, delta_x, delta_y

       Integer                                                          &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, global_row_length                                               &
                                ! number of points on a row
     &, global_rows                                                     &
                                ! number of rows in a theta field
     &, surface_type


      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, n_proc                                                          &
                   ! Total number of processors
     &, all_proc_group                                                  &
                       ! Group identifier for all processors.
     &, l_datastart(2)       ! First gridpoints held by this processor

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_code_test                                                     &
                         ! User switch
     &, L_fix_orog_hgt_lbc ! Switch to fix orog hgt in lbc zone

! Include physical constants
      Real                                                              &
           ! physical constants
     &  Earth_radius                                                    &
     &, Pi

      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, Base_phi                                                        &
     &, Base_lambda                                                     &
     &, lambda_p_end                                                    &
     &, phi_p_end

      Real                                                              &
           ! VarRes grid information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : n_rows+halo_j )

      Logical, Intent(In) :: L_regular


! local variables

      Integer                                                           &
     &  i, j, gi ,gj                                                    &
     &, haloi, haloj                                                    &
                      ! Local haloes - different for LAM/global
     &, info, power

      Real                                                              &
     &  x, y                                                            &
     &, xo, yo                                                          &
     &, width_x, width_y                                                &
     &, twopi                                                           &
     &, rad_to_deg                                                      &
     &, lambda_o, phi_o                                                 &
     &, lambda_deg, phi_deg                                             &
     &, h_agnesi                                                        &
     &, half_circumf                                                    &
     &, dist_x                                                          &
     &, dist_y                                                          &
     &, h_max

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
! Description: COMDECK containing surface types
!  for use in idealised problems
!
! Author : T. Davies
! History:
! Version  Date      Comment.
! 5.3      15/11/01  New code

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
      INTEGER, PARAMETER :: surface_dump=10

! No External routines

! ----------------------------------------------------------------------
! Section 1. Set domain information and check settings
! ----------------------------------------------------------------------
      rad_to_deg = 180.0 / Pi
      half_circumf = Earth_radius * Pi  !no need to vary with latitude
      If (Witch_power  <   1.1) Then
        power = nint(Witch_power)
      else
        power = 2
      endif !   Witch_power  <   1.1


      delta_x = Earth_radius * delta_lambda
      delta_y = Earth_radius * delta_phi

!  xo is equator value -
!  latitude accounted for in orog_haloes calculation
      If (model_domain  ==  mt_global) Then
!  For global domains, do not do haloes (apply swap bounds at end)
        haloi = 0
        haloj = 0
        lambda_o =  lambda_fraction * 2.0 * Pi
        phi_o = phi_fraction * Pi
        xo = Earth_radius * lambda_o
        yo = Earth_radius * phi_o
        phi_o = phi_o - 0.5 * Pi
        lambda_deg = lambda_o * rad_to_deg
        phi_deg = phi_o * rad_to_deg
      Else If (model_domain  ==  mt_cyclic_LAM .or.                     &
     &       model_domain  ==  mt_bi_cyclic_LAM) Then
!  For cyclic domains, do not do haloes (apply swap bounds at end)
        haloi = 0
        If(model_domain  ==  mt_cyclic_LAM)then
           haloj=halo_j
        else
           haloj=0
        endif
        If (L_regular) Then
        xo = lambda_fraction *  global_row_length * delta_x
        yo = phi_fraction * (global_rows - 1.0) * delta_y
        lambda_o = lambda_fraction * global_row_length * delta_lambda
        lambda_deg = lambda_o * rad_to_deg
        phi_o =  (phi_fraction - 0.5) * (global_rows - 1.0) * delta_phi
        phi_deg = phi_o * rad_to_deg
        Else
          lambda_o = lambda_fraction * (lambda_p_end-base_lambda)
          lambda_deg = lambda_o * rad_to_deg
          phi_o =  (phi_fraction - 0.5) * (phi_p_end-base_phi)
          phi_deg = phi_o * rad_to_deg
          xo = lambda_o * Earth_radius
          yo = phi_fraction *(phi_p_end-base_phi) * Earth_radius
        End If     ! L_regular
      Else If (model_domain  ==  mt_lam .and. L_regular .and.           &
     &              lambda_fraction  <   0.49)Then
!  For LAMs, fill haloes so extended halo regions have correct values
        haloi = halo_i
        haloj = halo_j
        xo = lambda_fraction *  global_row_length * delta_x
        yo = phi_fraction * (global_rows - 1.0) * delta_y
        lambda_o = lambda_fraction * global_row_length * delta_lambda
        lambda_deg = lambda_o * rad_to_deg
        phi_o =  (phi_fraction - 0.5) * (global_rows - 1.0) * delta_phi
        phi_deg = phi_o * rad_to_deg
      Else If (model_domain  ==  mt_lam .and. L_regular .and.           &
     &              lambda_fraction  >   0.49)Then
!  For LAMs, fill haloes so extended halo regions have correct values
        haloi = halo_i
        haloj = halo_j
        xo = lambda_fraction *  (global_row_length - 1.0) * delta_x
        yo = phi_fraction * (global_rows - 1.0) * delta_y
        lambda_o = lambda_fraction * (global_row_length - 1.0) *        &
     &                                        delta_lambda
        lambda_deg = lambda_o * rad_to_deg
        phi_o =  (phi_fraction - 0.5) * (global_rows - 1.0) * delta_phi
        phi_deg = phi_o * rad_to_deg
      Else If (model_domain  ==  mt_lam .and. .not. L_regular) Then
        haloi = halo_i
        haloj = halo_j
        lambda_o = lambda_fraction * (lambda_p_end-base_lambda)
        lambda_deg = lambda_o * rad_to_deg
        phi_o =  (phi_fraction - 0.5) * (phi_p_end-base_phi)
        phi_deg = phi_o * rad_to_deg
        xo = lambda_o * Earth_radius
        yo = phi_fraction *(phi_p_end-base_phi) * Earth_radius
      End If

      If(me  ==  0)then
        print*,' Model_domain = ',model_domain
        if ( model_domain  ==  mt_global ) then
          print*,' Global model_domain  =  ',model_domain
        elseif ( model_domain  ==  mt_LAM ) then
          print*,' Limited-area model_domain  =  ',model_domain
        elseif ( model_domain  ==  mt_cyclic_LAM ) then
          print*,' Cyclic (East-West) Limited-area model_domain  =  '   &
     &          ,model_domain
        elseif ( model_domain  ==  mt_bi_cyclic_LAM ) then
          print*,' Bi-Cyclic Limited-area model_domain  =  '            &
     &          ,model_domain
        else
          print*,' DANGER model domain ',model_domain,' NOT SUPPORTED'
          print*,' DANGER  MODEL WILL FAIL '
        endif ! model_domain  ==  mt_global
      End if ! (me == 0)

! ----------------------------------------------------------------------
! Section 2.  Insert mountain and set up orography arrays
! ----------------------------------------------------------------------

      If(surface_type  ==  surface_zero )then

        If(me  ==  0)then
          print*,' Flat surface_type = ',surface_type
          print*,' Orography set to 0.0 everywhere '
        End if ! (me == 0)

!  Make sure h_o = 0 for diagnostic printing
        h_o = 0.0

        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            orog_haloes(i,j) = 0.0
          end do
        end do

      Elseif(surface_type  ==  surface_ellipse)then
        If(me  ==  0)then
        print*,'surface_type = ',surface_type,' ** Witch of Agnesi **'
          print*,'Hill height = ',h_o,' metres. '
          print*,'Centre at ',lambda_deg,' degrees East or '            &
     &    ,lambda_o,' radians, ',xo,' metres from Greenwich meridian'
          if (phi_deg  >=  0.0) then
            print*,'Centre at ',phi_deg,' degrees North or '            &
     &            ,phi_o,' radians, ',yo,' metres from South Pole '
          else
            print*,'Centre at ',phi_deg,' degrees South or '            &
     &            ,phi_o,' radians, ',yo,' metres from South Pole '
          endif    !(phi_deg  >=  0.0)
          If(half_width_x  ==  half_width_y)then
            print*,' Circular hill half-width = ',half_width_x,' metres'
          Else
            print*,'Elliptical hill half_width_x = ',half_width_x,      &
     &             ' metres',' half_width_y = ',half_width_y,' metres'
          EndIf     !(half_width_x  ==  half_width_y)
        End if ! (me == 0)

        If (L_regular) then
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x = (gi-1) * delta_x
            y = (gj-1) * delta_y
            dist_x = ABS(x-xo)
!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x
            if(dist_x  >   half_circumf) then
              dist_x = 2.0 * half_circumf - x + xo
            endif     ! dist_x  >   half_circumf
            if(power  ==  1)then
              orog_haloes(i,j) = h_o /(1. +                             &
     &             (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
!   standard Witch of Agnesi includes ** 1.5 as below
!     &                         ( (y-yo)/half_width_y)**2.) ** 1.5
            else
              orog_haloes(i,j) = h_o /(1. +                             &
     &            (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)  &
     &             ** Witch_power
            endif   !(power  ==  1)
          end do
        end do
        Else
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            x = (lambda_p(i)-base_lambda) * Earth_radius
            y = (phi_p(i,j) - base_phi) * Earth_radius
            dist_x = ABS(x-xo)
!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x
            if(dist_x  >   half_circumf) then
              dist_x = 2.0 * half_circumf - x + xo
            endif     ! dist_x  >   half_circumf
            if(power  ==  1)then
              orog_haloes(i,j) = h_o /(1. +                             &
     &             (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
!   standard Witch of Agnesi includes ** 1.5 as below
!     &                         ( (y-yo)/half_width_y)**2.) ** 1.5
            else
              orog_haloes(i,j) = h_o /(1. +                             &
     &            (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)  &
     &             ** Witch_power
            endif   !(power  ==  1)
          end do
        end do
        End If   ! L_regular

      Elseif(surface_type  ==  surface_ridge) Then
!
! code to create a cos**2 ridge with amplitude 4*half_width_x
! Sets to 0 outside region
!
        width_x = 2.0 * half_width_x
        width_y = 2.0 * half_width_y

        If(me  ==  0)then
          print*,' surface_type = ',surface_type
          print*,' North-South ridge'
          print*,' Ridge height = ',h_o,' metres. Centre at '           &
     &          ,xo,' metres'
          print*,' Ridge half-width = ',half_width_x,' metres'
        End if ! (me == 0)

        If (L_regular) then
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x =  (gi-1)*delta_x
            if( abs(x-xo)  <=  width_x ) then
              orog_haloes(i,j) = h_o * cos(Pi*(x-xo)/(2.*width_x))**2.
            else
              orog_haloes(i,j) = 0.0
            endif
          end do
        end do
        Else
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            x = (lambda_p(i)-base_lambda) * Earth_radius
            if( abs(x-xo)  <=  width_x ) then
              orog_haloes(i,j) = h_o * cos(Pi*(x-xo)/(2.*width_x))**2.
            else
              orog_haloes(i,j) = 0.0
            endif
          end do
        end do
        End If   ! L_regular

      Elseif(surface_type  ==  surface_ridge_series) then

        width_x = 2.0 * half_width_x
        width_y = 2.0 * half_width_y
        twopi = 2.0 * Pi
        If(me  ==  0)then
          print*,' surface_type = ',surface_type
          print*,' North-South ridge series'
          print*,' Ridge height = ',h_o,' metres. Centre at '           &
     &          ,xo,' metres'
          print*,' Ridge half-width = ',half_width_x,' metres'
          gi = l_datastart(1)
          if (L_regular) Then
          x = (gi-1) * delta_x
          Else
            x = (lambda_p(1)-base_lambda) * Earth_radius
          End If
          orog_haloes(1,1) = h_o * cos(twopi * (x-xo) / width_x)
          print*,' First point gi value = ',gi,' First point x = ',x
          print*,' orog_haloes(1,1) = ',orog_haloes(1,1)
        End if ! (me == 0)

        if (L_regular) Then
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x =  (gi-1)*delta_x
            if( abs(x-xo)  <=  half_width_x ) then
              orog_haloes(i,j) = h_o * cos(twopi*(x-xo)/width_x)
            endif
          end do
        end do
        Else
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            x = (lambda_p(i)-base_lambda) * Earth_radius
            if( abs(x-xo)  <=  half_width_x ) then
              orog_haloes(i,j) = h_o * cos(twopi*(x-xo)/width_x)
            endif
          end do
        end do
        End If   ! L_regular

      Elseif(surface_type  ==  surface_plateau )then
        If(me  ==  0)then
          print*,'surface_type = ',surface_type,                        &
     &     ' ** FLATTENED  Witch of Agnesi **'
          print*,'** Constant height from centre to half-width '
          print*,'Plateau height = ',h_o,' metres. '
          print*,'Centre at ',lambda_deg,' degrees East or '            &
     &      ,lambda_o,' radians, ',xo,' metres from Greenwich meridian'
          if (phi_deg  >=  0.0) then
            print*,'Centre at ',phi_deg,' degrees North or '            &
     &      ,phi_o,' radians, ',yo,' metres from South Pole '
          else
            print*,'Centre at ',phi_deg,' degrees South or '            &
     &      ,phi_o,' radians, ',yo,' metres from South Pole '
          endif    !(phi_deg  >=  0.0)
          If(half_width_x  ==  half_width_y)then
          print*,' Circular hill half-width = ',half_width_x,' metres'
          Else
        print*,'Elliptical hill half_width_x = ',half_width_x,' metres' &
     &, ' half_width_y = ',half_width_y,' metres'
          EndIf     !(half_width_x  ==  half_width_y)
        End if ! (me == 0)

        if (L_regular) Then
        h_agnesi = 2.0 * h_o
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x = (gi-1) * delta_x
            y = (gj-1) * delta_y
            dist_x = ABS(x-xo)
!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x
            if(dist_x  >   half_circumf) then
              dist_x = 2.0 * half_circumf - x + xo
            endif     ! dist_x  >   half_circumf
            if(power  ==  1)then
              orog_haloes(i,j) = h_agnesi /(1. +                        &
     &             (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
!   standard Witch of Agnesi includes ** 1.5 as below
!     &                         ( (y-yo)/half_width_y)**2.) ** 1.5
            else
              orog_haloes(i,j) = h_agnesi /(1. +                        &
     &             (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2) &
     &             ** Witch_power
            endif   !(power  ==  1)

            if(orog_haloes(i,j)  >   h_o )then
              orog_haloes(i,j) = h_o
            endif         ! orography(i,j)  >   h_o
          end do
        end do
        Else
        h_agnesi = 2.0 * h_o
        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            x = (lambda_p(i)-base_lambda) * Earth_radius
            y = (phi_p(i,j) - base_phi) * Earth_radius
            dist_x = ABS(x-xo)
!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x
            if(dist_x  >   half_circumf) then
              dist_x = 2.0 * half_circumf - x + xo
            endif     ! dist_x  >   half_circumf
            if(power  ==  1)then
              orog_haloes(i,j) = h_agnesi /(1. +                        &
     &             (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
!   standard Witch of Agnesi includes ** 1.5 as below
!     &                         ( (y-yo)/half_width_y)**2.) ** 1.5
            else
              orog_haloes(i,j) = h_agnesi /(1. +                        &
     &             (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2) &
     &             ** Witch_power
            endif   !(power  ==  1)

            if(orog_haloes(i,j)  >   h_o )then
              orog_haloes(i,j) = h_o
            endif         ! orography(i,j)  >   h_o
          end do
        end do
        End If   ! L_regular

       Elseif(surface_type  ==  surface_massif)then
       If(me  ==  0)then
          print*,'surface_type = ',surface_type
          print*,' ** Massif plateau based upon 4 Witches of Agnesi **'
          print*,'Hill height = ',h_o,' metres. '
          print*,'Centre at ',lambda_deg,' degrees East or '            &
     &      ,lambda_o,' radians, ',xo,' metres from Greenwich meridian'
          if (phi_deg  >=  0.0) then
            print*,'Centre at ',phi_deg,' degrees North or '            &
     &            ,phi_o,' radians, ',yo,' metres from South Pole '
          else
            print*,'Centre at ',phi_deg,' degrees South or '            &
     &            ,phi_o,' radians, ',yo,' metres from South Pole '
          endif    !(phi_deg  >=  0.0)
          If(half_width_x  ==  half_width_y)then
            print*,' Circular hill half-width = ',half_width_x,' metres'
          Else
            print*,'Elliptical hill half_width_x = ',half_width_x,      &
     &             ' metres', ' half_width_y = ',half_width_y,' metres'
          EndIf     !(half_width_x  ==  half_width_y)
        End if ! (me == 0)

! There are nine sub-regions to set count 1-9 from bottom (S-Pole)
! and Eastwards
!               *            *
!           7   *     8      *    9
!               *            *
!     *************************************
!               * 5          * |
!           4   *   (xo,yo)  * 2yp   6
!               *            * |
!     *************************************
!               *<-- 2xp  -->*
!           1   *            *  3
!               *     2      *


        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            IF (L_regular) then
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x = (gi-1) * delta_x
            y = (gj-1) * delta_y
            Else
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            End If
            if(x  <   xo - xp) then
              if( y  <   yo - yp)then       ! region 1
                dist_x =  ABS(x-xo+xp)
                dist_y =  ABS(y-yo+yp)
              elseif( y  >   yo + yp)then   ! region 7
                dist_x =  ABS(x-xo+xp)
                dist_y =  ABS(y-yo-yp)
              else                         ! region 4
                dist_x =  ABS(x-xo+xp)
                dist_y =  0.0
              endif  !  y  <   yo- yp
            elseif(x  >   xo + xp) then
              if( y  <   yo - yp)then       ! region 3
                dist_x =  ABS(x-xo-xp)
                dist_y =  ABS(y-yo+yp)
              elseif( y  >   yo + yp)then   ! region 9
                dist_x =  ABS(x-xo-xp)
                dist_y =  ABS(y-yo-yp)
              else                         ! region 6
                dist_x =  ABS(x-xo-xp)
                dist_y =  0.0
              endif  !  y  <   yo- yp
            else
              if( y  <   yo - yp)then       ! region 2
                dist_x =  0.0
                dist_y =  ABS(y-yo+yp)
              elseif( y  >   yo + yp)then   ! region 8
                dist_x =  0.0
                dist_y =  ABS(y-yo-yp)
              else                         ! region 5
                dist_x =  0.0
                dist_y =  0.0
              endif  !  y  <   y- yp
            endif   ! x position

!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x

            if(dist_x  >   half_circumf ) then
              if(x  <   xo - xp) then
                dist_x = 2.0 * half_circumf - x + xo - xp
              elseif(x  >   xo + xp) then
                dist_x = 2.0 * half_circumf - x + xo + xp
              endif !x  <   xo - xp
            endif     ! dist_x  >   half_circumf
            if(power  ==  1)then
              orog_haloes(i,j) = h_o /(1. +                             &
     &            (dist_x/half_width_x)**2 + (dist_y/half_width_y)**2)
!   standard Witch of Agnesi includes ** 1.5 as below
!     &                         ( (y-yo)/half_width_y)**2.) ** 1.5
            else
              orog_haloes(i,j) = h_o /(1. +                             &
     &            (dist_x/half_width_x)**2 + (dist_y/half_width_y)**2)  &
     &             ** Witch_power
            endif   !(power  ==  1)
          end do
        end do

      Elseif(surface_type  ==  surface_mask)then

       If(me  ==  0)then
        print*,'surface_type = ',surface_type,' ** Real orography **'
        print*,'within box East-West half-width ',half_width_x,' metres'&
     &, ' North-South half-width ',half_width_y,' metres'
        print*,'Centre at ',lambda_deg,' degrees East or '              &
     & ,lambda_o,' radians, '
        print*,xo,' metres from Greenwich meridian'
         if (phi_deg  >=  0.0) then
           print*,'Centre at ',phi_deg,' degrees North or '             &
     &                        ,phi_o,' radians '
           print*,yo,' metres from South Pole '
         else
           print*,'Centre at ',phi_deg,' degrees South or '             &
     &                        ,phi_o,' radians, '
           print*,yo,' metres from South Pole '
         endif    !(phi_deg  >=  0.0)
       End if ! (me == 0)

!  Zero r_theta_levels(i,j,0)
      do j = 1-halo_j, rows+halo_j
        do i = 1-halo_i, row_length+halo_i
          orog_haloes(i,j)  = 0.0
        end do
      end do
!  Copy orography into r_theta_levels(i,j,0)
      do j = 1, rows
        do i = 1, row_length
          orog_haloes(i,j)  = orography(i,j)
        end do
      end do

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(orog_haloes, row_length, rows, 1,                &
     &                 halo_i, halo_j, fld_type_p,.FALSE.)



!   Reset h_o and use to store max orography
        h_max = 0.0

        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            IF (L_regular) then
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x = (gi-1) * delta_x
            y = (gj-1) * delta_y
            Else
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            End If
            dist_x = ABS(x-xo)
            dist_y = ABS(y-yo)
            if(dist_x  >   half_width_x .or.                            &
     &         dist_y  >   half_width_y) then
              orog_haloes(i,j) = 0.0
            endif     ! test if outside ellipse of real orography
!  Set h_o to max orography found
            if (orog_haloes(i,j)  >   h_max)then
              h_max = orog_haloes(i,j)
            endif
          end do
        end do
! find max  over all processors
        call gc_rmax(1, n_proc, info, h_max)
! Put max in h_o so that z_orog_print can calculate grid over max h_o
        h_o = h_max

      Elseif(surface_type  ==  surface_gauss)then

        If(me  ==  0)then
          print*,' surface_type = ',surface_type,' ** Gaussian **'
          print*,' Hill height = ',h_o,' metres. East-West centre at '  &
     &          ,xo,' metres, North-South centre at ',yo,' metres'
          If(half_width_x  ==  half_width_y)then
            print*,' Circular hill half-width = ',half_width_x,' metres'
          Else
        print*,'Elliptical hill half_width_x = ',half_width_x,' metres' &
     &    , ' half_width_y = ',half_width_y,' metres'
          EndIf     !(half_width_x  ==  half_width_y)
        End if ! (me == 0)

        do j = 1-haloj, rows+haloj
          do i = 1-haloi, row_length+haloi
            IF (L_regular) then
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x = (gi-1) * delta_x
            y = (gj-1) * delta_y
            Else
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            End If
            dist_x = ABS(x-xo)
            if(dist_x  >   half_circumf) then
              dist_x = 2.0 * half_circumf - x + xo
            endif     ! dist_x  >   half_circumf
            orog_haloes(i,j) = h_o * exp( -1.0 * (                      &
     &                   ( dist_x/half_width_x)**2. +                   &
     &                         ((y-yo)/half_width_y)**2.))
          end do
        end do

      Elseif(surface_type  ==  surface_dump)then

!  Zero r_theta_levels(i,j,0)
        do j = 1-halo_j, rows+halo_j
          do i = 1-halo_i, row_length+halo_i
            orog_haloes(i,j)  = 0.0
          end do
        end do
!  Copy orography into r_theta_levels(i,j,0)
        do j = 1, rows
          do i = 1, row_length
            orog_haloes(i,j)  = orography(i,j)
          end do
        end do

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(orog_haloes, row_length, rows, 1,              &
     &                 halo_i, halo_j, fld_type_p,.FALSE.)

        !--------------------------------------------------------
        !
        !          Set orography in the lbc zone
        !
        !--------------------------------------------------------

        ! If (L_fix_orog_hgt_lbc) is true and this is a limited
        ! area model domain, then set the height of
        ! the orography in the lateral boundary zone to a constant
        ! (external halo and rimwidth) and linearly interpolate
        ! in a zone in the interior to the input orography
        ! for a smooth transition (currently hardwired to 4 points).

        If (L_fix_orog_hgt_lbc .and.                                    &
     &      model_domain == mt_cyclic_LAM) Then

          If (me == 0) Then
            Write(6,*) ' '
            Write(6,*) ' L_fix_orog_hgt_lbc = .True.'
            Write(6,*) ' Setting constant height orography'
            Write(6,*) ' in the lbc zone.',orog_hgt_lbc,' (m)'
          End If

          If (at_extremity(PSouth)) Then
            Do j = 1-halo_j,  rimwidth
              Do i = 1-halo_i, row_length+halo_i
                orog_haloes(i,j) = orog_hgt_lbc
              End Do
            End Do
          End If

          If (at_extremity(PNorth)) Then
            Do j = 1-halo_j,  rimwidth
              Do i = 1-halo_i, row_length+halo_i
                orog_haloes(i,rows-j+1) = orog_hgt_lbc
              End Do
            End Do
          End If

          If (at_extremity(PWest)) Then
            Do j = 1, rows
              Do i = 1-halo_i, rimwidth
                orog_haloes(i,j) = orog_hgt_lbc
              End Do
            End Do
          End If

          If (at_extremity(PEast)) Then
            Do j = 1, rows
              Do i = 1-halo_i, rimwidth
                orog_haloes(rows-i+1, j) = orog_hgt_lbc
              End Do
            End Do
          End If

          If (at_extremity(PSouth)) Then
            Do j = rimwidth + 1,  rimwidth + 3
              Do i =  rimwidth + 1, row_length - rimwidth
                dist_x = j - rimwidth
                h_max = orography(i, rimwidth + 4) - orog_hgt_lbc
                orog_haloes(i,j) = orog_hgt_lbc + 0.25 * dist_x * h_max
              End Do
            End Do
          End If

          If (at_extremity(PNorth)) Then
            Do j = rimwidth + 1,  rimwidth + 3
              Do i =  rimwidth + 1, row_length - rimwidth
                dist_x = j - rimwidth
                gj = rows-j+1
                y = orography(i, rows - rimwidth - 3) - orog_hgt_lbc
                orog_haloes(i,gj) = orog_hgt_lbc + 0.25 * dist_x * y
              End Do
            End Do
          End If

          If (at_extremity(PWest)) Then
            Do j = rimwidth + 4, rows - rimwidth - 3
              Do i =  rimwidth + 1,  rimwidth + 3
                dist_x = i - rimwidth
                h_max = orography(rimwidth + 4 , j) - orog_hgt_lbc
                orog_haloes(i,j) = orog_hgt_lbc + 0.25 * dist_x * h_max
              End Do
            End Do
          End If

          If (at_extremity(PEast)) Then
            Do j = rimwidth + 4, rows - rimwidth - 3
              Do i =  rimwidth + 1,  rimwidth + 3
                dist_x = i - rimwidth
                gi = row_length-i+1
                y = orography(row_length-rimwidth-3, j) - orog_hgt_lbc
                orog_haloes(gi, j)= orog_hgt_lbc + 0.25 * dist_x * y
              End Do
            End Do
          End If

          If (at_extremity(PSouth) .and. at_extremity(PEast)) Then
            h_max = orography(rimwidth+4, rimwidth+4) - orog_hgt_lbc
            Do i = 1, 3
              gi = rimwidth+i
              dist_x = i
              Do j = gi, rimwidth+4, 1
                orog_haloes(j, gi) = orog_hgt_lbc + 0.25 * dist_x*h_max
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
            End Do
          End If

          If (at_extremity(PSouth) .and. at_extremity(PWest)) Then
            h_max = orography(row_length-rimwidth-3, rimwidth + 4)      &
     &                      - orog_hgt_lbc
            Do i = 1, 3
              gi = row_length-rimwidth+1-i
              gj = rimwidth+i
              dist_x = i
              Do j = gi, row_length-rimwidth-1,-1
                orog_haloes(j, gj) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
              Do j = gj, rimwidth+4 ,1
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
            End Do
          End If

          If (at_extremity(PNorth) .and. at_extremity(PWest)) Then
            h_max = orography(row_length-rimwidth-3, rows-rimwidth-3)   &
     &                      - orog_hgt_lbc
            Do i = 1, 3
              gi = row_length-rimwidth+1-i
              gj = rows-rimwidth+1-i
              dist_x = i
              Do j = gi, row_length-rimwidth-1,-1
                orog_haloes(j, gj) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
              Do j = gj, rows-rimwidth-1,-1
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
            End Do
          End If

          If (at_extremity(PNorth) .and. at_extremity(PEast)) Then
            h_max = orography(rimwidth+4,rows-rimwidth-3)-orog_hgt_lbc
            Do i = 1, 3
              gi = rimwidth +i
              gj = rows-rimwidth+1-i
              dist_x = i
              Do j = gi, rimwidth+4, 1
                orog_haloes(j, gj) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
              Do j = gj, rows-rimwidth-1,-1
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              End Do
            End Do
          End If

        End If ! on (L_fix_orog_hgt_lbc .and. mt_lam)
        !----------------------------------------------------------

!   Reset h_o and use to store max orography
        h_max = 0.0

        do j=1,rows
          do i=1,row_length
            IF (L_regular) then
            gi = l_datastart(1) + i - 1
            gj = l_datastart(2) + j - 1
            x = (gi-1) * delta_x
            y = (gj-1) * delta_y
            Else
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            End If
            dist_x = ABS(x-xo)
            dist_y = ABS(y-yo)
!  Set h_o to max orography found
            if (orography(i,j)  >   h_max)then
              h_max = orography(i,j)
            endif
          end do
        end do
! find max  over all processors
        call gc_rmax(1, n_proc, info, h_max)
! Put max in h_o so that z_orog_print can calculate grid over max h_o
        h_o = h_max
        If(me  ==  0)then
          print*,'surface_type = ',surface_type
          print*,' ** Orography in input dump used everywhere **'
          print*,'Maximum orographic height = ',h_o,' metres'
        End if ! (me == 0)

      Else

        print*,' surface_type = ',surface_type,' option unavailable '

      Endif        ! on surface_type

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(orog_haloes, row_length, rows, 1,              &
     &                 halo_i, halo_j, fld_type_p,.FALSE.)
! Next series of loops adjust the mountain profile so as to ensure that
! it is zero along Northern and Southern LAM boundaries.
!  Not needed for ridge (type 2) nor input orography (type 10)

      If(   surface_type  /=  surface_ridge                             &
     & .and. surface_type  /=  surface_dump )then

        print*,' Northern and Southern boundaries set to zero'

        If (at_extremity(PSouth)) Then
          Do j = 1, 2
            do i = 1-haloi, row_length+haloi
              orog_haloes(i,j) = 0.0
            End Do
          End Do
        End If

        If (at_extremity(PNorth)) Then
          Do j = rows-1, rows
            do i = 1-haloi, row_length+haloi
              orog_haloes(i,j) = 0.0
            End Do
          End Do
        End If

      EndIf  !        surface_type  /=  surface_ridge
             ! .and.  surface_type  /=  surface_dump

!  Copy orog_haloes into D1 array
      do j = 1, rows
        do i = 1, row_length
          orography(i,j) = orog_haloes(i,j)
        end do
      end do

      Return
      END SUBROUTINE IDL_Surface_setup

! End IDL_Surface_setup
