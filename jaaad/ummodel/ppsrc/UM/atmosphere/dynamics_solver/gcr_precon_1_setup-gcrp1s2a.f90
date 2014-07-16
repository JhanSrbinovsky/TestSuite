
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine GCR_precon_1_setup

      Subroutine GCR_precon_1_setup(                                    &
     &                     HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,          &
     &                     HM_Czz, HM_Cz, HM_C3, HM_C4,                 &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     r_theta_levels, r_rho_levels,                &
     &                     row_length, rows, model_levels,              &
     &                     FV_sec_theta_latitude,                       &
     &                     model_domain,                                &
     &                     weight_upper, weight_lower,                  &
     &                     a0, a1, factor,                              &
     &                     offx,offy,n_proc, global_row_length,         &
     &                     at_extremity,                                &
     &                     n_procx,n_procy,proc_row_group,              &
     &                     halo_i,halo_j                                &
     &                     )

! Purpose:
!          Calculates constant terms used inpre-conditioning operator
!          applied to field.
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Original Progammer: Mark H. Mawson
! Current code owner: Andrew J. Malcolm
!
! History:
! Date     Version     Comment
! ----     -------     -------
! 22/03/00 5.1        Changed j_start/stop for cyclic LAM  Andy Malcolm
!LL   5.1   11/02/00  Use DOMTYP parameters                    P.Burton
!   5.3     15/09/01  add mt_bi_cyclic_LAM code            A. Malcolm
!     5.3   19/10/01  Use appropriate gcg routines.   S. Cusack
!     6.2   21/10/05  Remove unused external. P.Selwood.
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit none

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, model_domain

      Integer                                                           &
     &  offx,offy,n_proc, global_row_length,n_procx,n_procy             &
     &  ,halo_i,halo_j

! interpolation weights for moving between theta levels and rho levels
      Real                                                              &
     &  weight_upper(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          &
     &, weight_lower(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)

      Real                                                              &
     &  HM_Cxx1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cxx2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Cyy2(1-offx:row_length+offx,1-offy:rows+offy,model_levels)   &
     &, HM_Czz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)    &
     &, HM_Cz(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C3(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, HM_C4(1-offx:row_length+offx,1-offy:rows+offy,model_levels)     &
     &, FV_sec_theta_latitude(1-offx:row_length+offx,1-offy:rows+offy)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i                      &
     &                 ,1-halo_j:rows+halo_j,0:model_levels)            &
     &, r_rho_levels (1-halo_i:row_length+halo_i                        &
     &                 ,1-halo_j:rows+halo_j,model_levels)              &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! Arguments with Intent OUT. ie: Output variables.
      Real                                                              &
     &  a0(row_length,rows,model_levels)                                &
     &, a1(row_length,rows,model_levels)                                &
     &, factor(row_length,rows,model_levels)

! Local Variables.

      Integer                                                           &
     &  i,j,k                                                           &
     &, istat                                                           &
     &, j0                                                              &
     &, j1                                                              &
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop

      Real                                                              &
     &  factor_1                                                        &
     &, factor_2                                                        &
     &, factor_3

!  parallel variables
      Integer                                                           &
     &  proc_row_group

      Real                                                              &
     &  sum_n(model_levels)                                             &
     &, sum_s(model_levels)                                             &
     &, sum_n_component(row_length,model_levels)                        &
     &, sum_s_component(row_length,model_levels)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Local arrays.

      Real                                                              &
     &  a2(row_length,rows,model_levels)

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
!   External routines.
      External                                                          &
     &   gcg_rvecsumr

!-----------------------------------------------------------------------
!     Section 1.
!-----------------------------------------------------------------------

        if(at_extremity(PSouth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
        j0= 2
      Else
        j0 = 1
      End If
        if(at_extremity(PNorth) .and.                                   &
     &     model_domain  /=  mt_bi_cyclic_LAM)then
        j1 = rows - 1
      Else
        j1 = rows
      End If

      If (model_domain  ==  mt_global .or.                              &
     &    model_domain  ==  mt_bi_cyclic_LAM) Then
! Solve over full domain
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
      Else If(model_domain  ==  mt_lam) Then
! Solve over interior points only.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
        if(at_extremity(PEast)) i_stop = row_length-1
        if(at_extremity(PWest)) i_start = 2
      Elseif (model_domain  ==  mt_cyclic_LAM) then
! Solve over interior points only periodic in x => i_start=1,
! i_stop=row_length.
        i_start = 1
        i_stop = row_length
        j_start = 1
        j_stop = rows
        if(at_extremity(PSouth)) j_start = 2
        if(at_extremity(PNorth)) j_stop = rows-1
      End If

      Do k = 1, model_levels
        If ( k  ==  1) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +        &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k) )
              a0(i,j,k) = - a1(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do

        Else if (k  ==  model_levels) Then
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a2(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -      &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0(i,j,k) = - a2(i,j,k) - HM_C4(i,j,k)
            End Do
          End Do
        Else
          factor_1 = 1. / (eta_theta_levels(k) -                        &
     &                     eta_theta_levels(k-1))
          factor_2 = 1. / (eta_rho_levels(k+1) -                        &
     &                     eta_rho_levels(k))
          factor_3 = 1. / (eta_rho_levels(k) -                          &
     &                     eta_rho_levels(k-1))

          Do j = j_start, j_stop
            Do i = i_start, i_stop
              a1(i,j,k) = factor_2 * (factor_1 * HM_Czz(i,j,k) +        &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k) *                 &
     &                                weight_upper(i,j,k) )
              a2(i,j,k) = factor_3 * (factor_1 * HM_Czz(i,j,k-1) -      &
     &                    HM_C3(i,j,k) * HM_Cz(i,j,k-1) *               &
     &                                weight_lower(i,j,k) )
              a0(i,j,k) = - a2(i,j,k) - a1(i,j,k)  - HM_C4(i,j,k)
            End Do
          End Do
        End If

      End Do

! Include horizontal second derivative terms in matrix.

! Interior points - Common to all model_domains.

      Do k = 1, model_levels
        Do j = j0, j1
          Do i = i_start, i_stop
            a0(i,j,k) = a0(i,j,k) - FV_sec_theta_latitude(i,j) *        &
     &                  (HM_Cxx1(i-1,j,k)*HM_Cxx2(i-1,j,k) +            &
     &                   HM_Cxx1(i,j,k)*HM_Cxx2(i,j,k) +                &
     &                   HM_Cyy1(i,j,k)*HM_Cyy2(i,j,k) +                &
     &                   HM_Cyy1(i,j-1,k)*HM_Cyy2(i,j-1,k) )
          End Do
        End Do
      End Do

! Treatment of North/South boundary.

      If (Model_domain  ==  mt_Global )Then

        If(at_extremity(PSouth))then
          Do k = 1, model_levels
            Do i = 1,row_length
              sum_s_component(i,k)= - HM_Cyy1(i,1,k)*HM_Cyy2(i,1,k)
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Do k = 1, model_levels
            Do i = 1,row_length
              sum_n_component(i,k)= - HM_Cyy1(i,rows-1,k)*              &
     &                                HM_Cyy2(i,rows-1,k)
            End Do
          End Do
        End If

        If(at_extremity(PSouth))then
          Call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      model_levels,sum_s_component,               &
     &                      proc_row_group,istat,sum_s)
          Do k = 1, model_levels
            sum_s(k) = sum_s(k) * FV_sec_theta_latitude(1,1)            &
     &                   /global_row_length
            Do i = 1,row_length
              a0(i,1,k) = a0(i,1,k) + sum_s(k)
            End Do
          End Do
        End If
        If(at_extremity(PNorth))then
          Call gcg_rvecsumr(row_length,row_length,1,                    &
     &                      model_levels,sum_n_component,               &
     &                      proc_row_group,istat,sum_n)
          Do k = 1, model_levels
            sum_n(k) = sum_n(k) * FV_sec_theta_latitude(1,rows)         &
     &                   /global_row_length
            Do i = 1,row_length
              a0(i,rows,k) = a0(i,rows,k) + sum_n(k)
            End Do
          End Do
        End If

      End If

! Solve tridiagonal system.
! solution is in Soln
! reduce matrix to upper diagonal form.

! for global parallel code a0 only valid on interior points , NOT halos
      Do j = j_start, j_stop
        Do i = i_start, i_stop
          a0(i,j,1) = 1./a0(i,j,1)
        End Do
      End Do

      Do k= 2, model_levels
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            factor(i,j,k) = a2(i,j,k) * a0(i,j,k-1)
            a0(i,j,k) = 1./(a0(i,j,k) - factor(i,j,k)*a1(i,j,k-1))
          End Do
        End Do
      End Do

!     end of routine GCR_precon_1_setup

      return
      END SUBROUTINE GCR_precon_1_setup

