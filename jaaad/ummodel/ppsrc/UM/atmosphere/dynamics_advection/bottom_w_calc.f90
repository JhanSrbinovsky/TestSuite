
! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************
!
! Subroutine Bottom_w_Calc.

      Subroutine Bottom_w_Calc(                                         &
     &                         r_theta_levels, u, v,                    &
     &                         w, w_adv, sec_theta_latitude,            &
     &                         delta_lambda, delta_phi,                 &
     &                         lambda_p, phi_p, lambda_u, phi_v,        &
     &                         rows, n_rows, row_length, model_levels,  &
     &                         model_domain, L_regular,                 &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         proc_row_group, at_extremity,            &
     &                         g_row_length)

! Purpose:
!          Calculates etadot from given u, v and w fields.
!
! Method:
!          Is described in ;
!
! Original Progammer: C. J. Smith/ Terry davies
! Current code owner: Andy Malcolm
!
! History:
! Version   Date       Comment
! -------  -------     -------
!  5.3     10/10/01    Deck introduced                   Terry Davies
!  5.4     28/08/02    Bug Fix (Bi-cyclic LAM)           Carol Roadnight
!  6.2    25/12/05  Variable resolution changes            Yongming Tang
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, g_row_length                                                    &
                         ! global number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels
     &, off_x                                                           &
                   ! Size of small halo in i direction
     &, off_y                                                           &
                   ! Size of small halo in j direction
     &, halo_i                                                          &
                ! Size of halo in i direction
     &, halo_j                                                          &
                ! Size of halo in j direction
     &, proc_row_group ! Group id for processors on the same row

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south,east or west of the processor grid
     &, L_regular

      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )                           &
     &, phi_v    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : n_rows+halo_j )

      Real                                                              &
           ! trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

! Primary Arrays

      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k                                                         &
                          ! Loop indices
     &, j0                                                              &
     &, j1

      Real                                                              &
     &  Recip_delta_lambda                                              &
     &, Recip_delta_phi                                                 &
     &, weight1                                                         &
     &, weight2                                                         &
     &, sum_row

      Integer info

! arrays

      Real                                                              &
     &  u_on_theta(0:row_length, rows)                                  &
     &, v_on_theta(row_length, 0:n_rows)                                &
     &, term(0:row_length, 0:rows)

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

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Interpolate u and v to w points.
!              Not done at top level.
! ----------------------------------------------------------------------

      j0 = 1
      j1 = rows
      If (model_domain  /=  mt_bi_cyclic_LAM) then
      If (at_extremity(PSouth)) j0 = 2
      If (at_extremity(PNorth)) j1 = rows-1
      Endif
      recip_delta_lambda = 1./ delta_lambda
      recip_delta_phi = 1./ delta_phi

      k = 0

!  Level 1 u assumed to be same as surface

      Do j = j0, j1
        Do i = 0, row_length
          u_on_theta (i,j) = u(i,j,k+1)
        End Do
      End Do

!  Level 1 v assumed to be same as surface
      Do j = j0-1, n_rows
        Do i = 1, row_length
          v_on_theta (i,j) =  v(i,j,k+1)
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2.   Calculate u.grad(r)
! ----------------------------------------------------------------------

      IF (L_regular) then

! Calculate u dr/d lambda
      Do j = j0, j1
        Do i = 0, row_length
          term(i,j) = u_on_theta(i,j) *                                 &
     &               (r_theta_levels(i+1,j,k) -                         &
     &                r_theta_levels(i,j,k) ) * recip_delta_lambda      &
     &                * 2.0 / (r_theta_levels(i+1,j,k) +                &
     &                         r_theta_levels(i,j,k) )
        End Do
      End Do

! average quantity and store in w(i,j,0)
      Do j = j0, j1
        Do i = 1, row_length
          w(i,j,k) = (term(i,j) + term(i-1,j) ) * 0.5 *                 &
     &                        sec_theta_latitude(i,j)
        End Do
      End Do

! Calculate v dr/d phi
      Do j = j0-1, n_rows
        Do i = 1, row_length
          term(i,j) = v_on_theta(i,j) *                                 &
     &               (r_theta_levels(i,j+1,k) -                         &
     &                r_theta_levels(i,j,k) ) * recip_delta_phi         &
     &                  * 2.0 / (r_theta_levels(i,j+1,k) +              &
     &                           r_theta_levels(i,j,k) )
        End Do
      End Do

! average quantity and add on to u term
      Do j = j0, j1
        Do i = 1, row_length
          w(i,j,k) = (term(i,j) + term(i,j-1)) * 0.5 +                  &
     &                        w(i,j,k)
          w_adv(i,j,k) = w(i,j,k)
        End Do
      End Do

      ELSE

! Calculate u dr/d lambda
      Do j = j0, j1
        Do i = 0, row_length
        weight1 =  lambda_p(i+1)-lambda_u(i)
        weight2 =  lambda_u(i) - lambda_p(i)
        term(i,j) = u_on_theta(i,j)                                     &
     &            * (r_theta_levels(i+1,j,k) - r_theta_levels(i,j,k))   &
     &            / ( weight1 * r_theta_levels(i,j,k) +                 &
     &                weight2 * r_theta_levels(i+1,j,k) )
        End Do
      End Do

! average quantity and store in w(i,j,0)
      Do j = j0, j1
        Do i = 1, row_length
        recip_delta_lambda = 1.0/(lambda_u(i)-lambda_u(i-1))
        weight1 = (lambda_u(i)-lambda_p(i)) * recip_delta_lambda
        weight2 = 1.0-weight1
        w(i,j,k) = ( weight1 * term(i-1,j) + weight2 * term(i,j) )      &
     &                   *  sec_theta_latitude(i,j)
        End Do
      End Do

! Calculate v dr/d phi
      Do j = j0-1, n_rows
        Do i = 1, row_length
          weight1 =  phi_p(i,j+1) - phi_v(i,j)
          weight2 =  phi_v(i,j) - phi_p(i,j)
          term(i,j) = v_on_theta(i,j) *                                 &
     &               (r_theta_levels(i,j+1,k) - r_theta_levels(i,j,k))/ &
     &                              ( weight1 * r_theta_levels(i,j,k) + &
     &                                weight2 * r_theta_levels(i,j+1,k))
        End Do
      End Do

! average quantity and add on to u term
      Do j = j0, j1
        Do i = 1, row_length
          recip_delta_phi  =  1.0 / (phi_v(i,j) - phi_v(i,j-1))
          weight1 =  (phi_v(i,j) - phi_p(i,j)) * recip_delta_phi
          w(i,j,k) = ( weight1 * term(i,j-1) +                          &
     &                (1.0 - weight1) * term(i,j) ) +  w(i,j,k)
          w_adv(i,j,k) = w(i,j,k)
        End Do
      End Do

      End If   ! L_regular


! northern and southern boundary updates.
! save values for summing.

      If (model_domain  ==  mt_global) Then    ! global model domain

        If(at_extremity(PSouth))then

          Call gcg_rvecsumr(row_length, row_length, 1, 1,               &
     &                     term(1,1), proc_row_group, info, sum_row)

          sum_row = sum_row / g_row_length
          Do i = 1, row_length
            w(i,1,k) = sum_row
          End Do
        End If

        If(at_extremity(PNorth))then
          Call gcg_rvecsumr(row_length, row_length, 1, 1,               &
     &                   term(1,n_rows), proc_row_group, info, sum_row)

          sum_row = sum_row / g_row_length
          Do i = 1, row_length
            w(i,rows,k) = sum_row
          End Do
        End If

      Else     ! limited area model.
!  Assume values outside boundary are zero.

       If(at_extremity(PSouth) .and.                                    &
     &                 model_domain  /=  mt_bi_cyclic_LAM) then
          Do i= 1,row_length
            w(i,1,k) = 0.
          End Do
        End If

       If(at_extremity(PNorth) .and.                                    &
     &                 model_domain  /=  mt_bi_cyclic_LAM) then
          j = rows
          Do i=1,row_length
            w(i,j,k) = 0.
          End Do
        End If

      End If   !( model_domain  ==  mt_global )

! End of routine  Bottom_w_Calc.
      return
      END SUBROUTINE Bottom_w_Calc

