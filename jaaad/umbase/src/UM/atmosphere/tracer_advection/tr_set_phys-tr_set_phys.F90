#if defined(A11_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE TR_Set_Phys(                                           &
     &                      super_array_size, super_tracer_phys,        &
     &                      L_CO2_interactive, CO2,                     &
     &                      L_Murk_advect, murk,                        &
     &                      L_Soot, soot_new, soot_agd, soot_cld,       &
     &                      L_SULPC_SO2, SO2, SO4_aitken, so4_accu,     &
     &                                   so4_diss,                      &
     &                      L_sulpc_nh3, nh3,                           &
     &                      L_sulpc_dms, dms,                           &
     &                      L_dust, DUST_DIV1, DUST_DIV2, DUST_DIV3,    &
     &                              DUST_DIV4, DUST_DIV5, DUST_DIV6,    &
     &                      L_biomass, bmass_new, bmass_agd, bmass_cld, &
     &                      L_ocff, ocff_new, ocff_agd, ocff_cld,       &
     &                      L_USE_CARIOLLE, OZONE_TRACER,               &
     &                      tracer_phys1, tracer, tracer_ukca,          &
     &                      row_length, rows,                           &
     &                      model_levels, tr_levels, tr_vars, tr_ukca,  &
     &                      offx, offy, model_domain,                   &
     &                      L_init, halo_i_in,halo_j_in                 &
     &                      )

! Purpose: Interface routine to initialise tracer fields
!
! Method:
!
! Original Programmer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version  Date        Comment
! -------  -------     -----
! 5.5      17/12/02    Original Deck introduced     A. Malcolm
! 6.4      17/12/06    Major rewrite                A. Malcolm
! 6.4      16/02/07    Fix for tracers in LAM's     A. Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, tr_levels                                                       &
                         ! number of tracer levels
     &, tr_vars                                                         &
                         ! number of tracers
     &, tr_ukca                                                         &
                         ! number of ukca tracers
     &, super_array_size                                                &
     &, offx                                                            &
     &, offy
      Integer, Intent(In) :: model_domain
      Integer, Intent(In) :: halo_i_in
      Integer, Intent(In) :: halo_j_in

      Real, Intent(In) ::                                               &
     &  CO2      (1-offx:row_length+offx, 1-offy:rows+offy,             &
     &            model_levels)                                         &
     &, murk      (1-offx:row_length+offx, 1-offy:rows+offy,            &
     &             model_levels)                                        &
     &, soot_new      (1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, soot_agd      (1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, soot_cld      (1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, SO2             (1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels)                &
     &, SO4_aitken      (1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels)                &
     &, so4_accu        (1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels)                &
     &, so4_diss        (1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels)                &
     &, nh3      (1-offx:row_length+offx,                               &
     &            1-offy:rows+offy, model_levels)                       &
     &, dms      (1-offx:row_length+offx,                               &
     &            1-offy:rows+offy, model_levels)                       &
     &, dust_div1      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, dust_div2      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, dust_div3      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, dust_div4      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, dust_div5      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, dust_div6      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, bmass_new      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, bmass_agd      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, bmass_cld      (1-offx:row_length+offx,                         &
     &            1-offy:rows+offy, model_levels)                       &
     &, ocff_new      (1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, ocff_agd      (1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, ocff_cld      (1-offx:row_length+offx,                          &
     &                 1-offy:rows+offy, model_levels)                  &
     &, tracer      (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &               tr_levels,1:tr_vars)                               &
     &, tracer_ukca (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &               tr_levels,1:tr_ukca)                               &
! Add cariolle specific parameters for ozone tracer     
     &, OZONE_TRACER(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)                                          

      Real, Intent(Out) :: tracer_phys1(1-halo_i_in:row_length+halo_i_in, &
                                        1-halo_j_in:rows+halo_j_in,     &
                                        tr_levels,tr_vars)     
      Real, Intent(Out) :: super_tracer_phys(1-halo_i_in:row_length+halo_i_in, &
                                             1-halo_j_in:rows+halo_j_in, &
                                             model_levels,super_array_size) 

      Logical                                                           &
     &  L_CO2_interactive                                               &
     &, L_Murk_advect                                                   &
     &, L_soot                                                          &
     &, L_SULPC_SO2                                                     &
     &, L_sulpc_nh3                                                     &
     &, L_sulpc_dms                                                     &
     &, L_biomass                                                       &
     &, L_dust                                                          &
     &, L_ocff                                                          &
     &, L_USE_CARIOLLE

      Logical :: L_init   ! flag for setting halo values 

!      local variables
      Integer                                                           &
     &  i,j,k,l,count        !loop variables

      Integer :: array_size_count

#include "parparm.h"
#include "domtyp.h"

      array_size_count=0

! ----------------------------------------------------------------------
! Section 1.1  carbon cycle.
! ----------------------------------------------------------------------
      If (L_CO2_interactive) then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) = CO2(i,j,k)
            End Do
          End Do
        End Do
      end if   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 1.2  Soot cycle.
! ----------------------------------------------------------------------
      If (L_soot) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  soot_new(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  soot_agd(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  soot_cld(i,j,k)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      end if    ! L_soot

! ----------------------------------------------------------------------
! Section 1.3  Biomass aerosol.
! ----------------------------------------------------------------------
      If (l_Biomass) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  bmass_new(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  bmass_agd(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  bmass_cld(i,j,k)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      endif  ! l_biomass

! ----------------------------------------------------------------------
! Section 1.4  sulphur cycle.
! ----------------------------------------------------------------------
      If (L_SULPC_SO2) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  so4_aitken(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  so4_accu(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  so4_diss(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+3) =  so2(i,j,k)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +3

        If(L_sulpc_nh3) then
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =  nh3(i,j,k)
              End Do
            End Do
          End Do
        end if  ! L_sulpc_nh3

        If(L_sulpc_dms) then
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =  dms(i,j,k)
              End Do
            End Do
          End Do
        end if  ! L_sulpc_dms
      end if  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 1.5  mineral dust.
! ----------------------------------------------------------------------
      If (L_dust) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  dust_div1(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  dust_div2(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  dust_div3(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+3) =  dust_div4(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+4) =  dust_div5(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+5) =  dust_div6(i,j,k)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +5
      end if    ! L_dust

! ----------------------------------------------------------------------
! New addition. Fossil-fuel organic carbon aerosol
! ----------------------------------------------------------------------
      If (L_ocff) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  ocff_new(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+1) =  ocff_agd(i,j,k)
              super_tracer_phys(i,j,k,array_size_count+2) =  ocff_cld(i,j,k)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      end if    ! L_ocff

! ----------------------------------------------------------------------
! Section 1.5.1  cariolle ozone tracer.
! ----------------------------------------------------------------------
      If (L_USE_CARIOLLE) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  OZONE_TRACER(i,j,k)
            End Do
          End Do
        End Do
      end if    ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Section 1.6.a  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF (model_levels==tr_levels.and. tr_vars>0) THEN
        do count=1,tr_vars
          array_size_count=array_size_count +1
          DO K = 1, MODEL_LEVELS
            DO J = 1, ROWS
              DO I = 1, ROW_LENGTH
                super_tracer_phys(i,j,k,array_size_count) = tracer(i,j,k,count)
              End Do
            End Do
          End Do
        End Do

!       IF(L_init) then
!       ! Call SWAPBOUNDS and set external halos for free tracers
! DEPENDS ON: swap_bounds
!         Call Swap_Bounds(                                             &
!    &     super_tracer_phys(1-halo_i_in,1-halo_j_in,1,                 &
!    &                       super_array_size-tr_vars+1),               &
!    &       row_length, rows,                                          &
!    &       model_levels*tr_vars,                                      &
!    &       halo_i_in, halo_j_in, fld_type_p,  .false. )
!         If (model_domain == mt_lam) Then
! DEPENDS ON: set_external_halos
!           Call SET_EXTERNAL_HALOS(                                    &
!    &      super_tracer_phys(1-halo_i_in,1-halo_j_in,1,                &
!!   &                        super_array_size-tr_vars+1),              &
!    &      row_length, rows,                                           &
!    &      tr_vars*model_levels, halo_i_in, halo_j_in, 0.0)
!         End If
!       End If ! L_INIT
      End IF  ! tr_vars>0

! ----------------------------------------------------------------------
! Section 1.7  UKCA tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      IF (model_levels==tr_levels .and. tr_ukca > 0) THEN
        DO count=1,tr_ukca
          array_size_count=array_size_count +1
          DO K = 1, model_levels
            DO J = 1, rows
              DO I = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =             &
     &               tracer_ukca(i,j,k,count)
              END DO
            END DO
          END DO
        END DO
      END IF     ! tr_ukca>0

! ----------------------------------------------------------------------
! Section 1.99  Murk cycle.  This must be the last Full level field in
!                            the super_array
! ----------------------------------------------------------------------
      If (L_Murk_advect) then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) =  murk(i,j,k)
            End Do
          End Do
        End Do
      end if  ! L_Murk_advect

      if(L_init)then
! ----------------------------------------------------------------------
! Call SWAPBOUNDS and set external halos for all arrays if required
! ----------------------------------------------------------------------

! DEPENDS ON: swap_bounds
        Call Swap_Bounds(                                               &
     &       super_tracer_phys,                                         &
     &       row_length, rows,                                          &
     &       model_levels*array_size_count,                             &
     &       halo_i_in, halo_j_in, fld_type_p,  .false. )

        If (model_domain == mt_lam) Then
! DEPENDS ON: set_external_halos
          Call SET_EXTERNAL_HALOS(super_tracer_phys, row_length, rows,  &
     &        array_size_count*model_levels, halo_i_in, halo_j_in, 0.0)
        End If

      Endif  !L_init
        IF (L_Murk_advect) then

! DEPENDS ON: fill_external_halos
          CALL FILL_EXTERNAL_HALOS(                                     &
     &  super_tracer_phys(1-halo_i_in,1-halo_j_in,1,array_size_count)   &
     &     , row_length, rows,                                          &
     &                     model_levels, halo_i_in, halo_j_in)

        End IF  ! L_Murk_advect
      
!ajm     Endif  !L_init
! ----------------------------------------------------------------------
! Section 1.6.b  free tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      If (tr_vars > 0.and. tr_levels /= model_levels) Then
        do l = 1 , tr_vars
          Do k = 1, tr_levels
            Do j = 1, rows
              Do i = 1, row_length
                tracer_phys1(i,j,k,l) = TRACER(i,j,k,l)
              End Do
            End Do
          End Do
        End Do

        If(L_INIT)then
! Call SWAPBOUNDS and set external halos for free tracers if necessary
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &         tracer_phys1,                                            &
     &         row_length, rows, tr_levels*tr_vars,                     &
     &         halo_i_in, halo_j_in, fld_type_p,  .false. )
          If (model_domain == mt_lam) Then
! DEPENDS ON: set_external_halos
            Call SET_EXTERNAL_HALOS(tracer_phys1, row_length,           &
     &        rows, tr_vars*tr_levels, halo_i_in, halo_j_in, 0.0)
          End If
        End If !L_init
      End If  ! tr_vars > 0

      Return
      END SUBROUTINE TR_Set_Phys
#endif
