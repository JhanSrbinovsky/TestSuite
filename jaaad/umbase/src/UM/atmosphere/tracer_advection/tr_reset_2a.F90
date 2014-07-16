#if defined(A11_2A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
        SUBROUTINE TR_Reset(                                            &
                            super_array_size, super_tracer_phys,        &
                            L_CO2_interactive, CO2,                     &
                            L_Murk_advect, murk,                        &
                            L_Soot, soot_new, soot_agd, soot_cld,       &
                            L_SULPC_SO2, SO2, SO4_aitken, so4_accu,     &
                                         so4_diss,                      &
                            L_sulpc_nh3, nh3,                           &
                            L_sulpc_dms, dms,                           &
                            L_dust, DUST_DIV1, DUST_DIV2, DUST_DIV3,    &
                                    DUST_DIV4, DUST_DIV5, DUST_DIV6,    &
                            L_biomass, bmass_new, bmass_agd, bmass_cld, &
                            L_ocff, ocff_new, ocff_agd, ocff_cld,       &
                            L_USE_CARIOLLE, OZONE_TRACER,               &
                            tracer_phys2, tracer, tracer_ukca,          &
                            row_length, rows,                           &
                            model_levels, tr_levels, tr_vars, tr_ukca,  &
                            offx, offy                                  &
                            )

! Purpose: Interface routine to calculate increase in tracer fields
!          during the call to atmos_physics2
!
! Method:
!
! Original Programmer: Andrew J. Malcolm
! Current code owner: Andrew J. Malcolm
!
! History:
! Version  Date        Comment
! -------  -------     -----
! 6.4      17/12/06    Original Deck introduced     Andy Malcolm
! 6.4      16/02/07    fix for tracers in LAM's     Andy Malcolm
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      Implicit None

! Arguments with Intent IN. ie: Input variables.

! Arguments with Intent IN. ie: Input variables.

      Integer, Intent(In) :: row_length  ! number of points on a row
      Integer, Intent(In) :: rows        ! number of rows in a theta field
      Integer, Intent(In) :: model_levels   ! number of model levels         
      Integer, Intent(In) :: tr_levels   ! number of tracer levels         
      Integer, Intent(In) :: tr_vars   ! number of tracers         
      Integer, Intent(In) :: tr_ukca   ! number of ukca tracers         
      Integer, Intent(In) :: super_array_size  
      Integer, Intent(In) :: offx      ! halo size in x direction
      Integer, Intent(In) :: offy      ! halo size in y direction

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
     &            1-offy:rows+offy, model_levels)                       &
     &, ocff_agd      (1-offx:row_length+offx,                          &
     &            1-offy:rows+offy, model_levels)                       &
     &, ocff_cld      (1-offx:row_length+offx,                          &
     &            1-offy:rows+offy, model_levels)                       &
     &, tracer      (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &               tr_levels,1:tr_vars)                               &
     &, tracer_ukca (1-offx:row_length+offx, 1-offy:rows+offy,          &
     &               tr_levels,1:tr_ukca)                               &
! Add cariolle specific parameters for ozone tracer     
     &, OZONE_TRACER(1-offx:row_length+offx, 1-offy:rows+offy,          &
     &           model_levels)              

      Real, Intent(InOut) ::                                            &
     &  super_tracer_phys(row_length, rows,                             &
     &               model_levels,super_array_size)                     &
     &, tracer_phys2(row_length, rows, tr_levels,tr_vars)                

      Logical, Intent(In):: L_CO2_interactive 
      Logical, Intent(In):: L_Murk_advect     
      Logical, Intent(In):: L_soot
      Logical, Intent(In):: L_SULPC_SO2
      Logical, Intent(In):: L_SULPC_nh3
      Logical, Intent(In):: L_SULPC_dms
      Logical, Intent(In):: L_biomass  
      Logical, Intent(In):: L_dust
      Logical, Intent(In):: L_ocff
      Logical, Intent(In):: L_USE_CARIOLLE

!      local variables
      Integer                                                           &
     &  i,j,k,l,count        !loop variables

      Integer :: array_size_count

      array_size_count=0

! ----------------------------------------------------------------------
! Section 1.1  carbon cycle.
! ----------------------------------------------------------------------
      If (L_CO2_interactive) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              super_tracer_phys(i,j,k,array_size_count) = CO2(i,j,k)    &
     &             -super_tracer_phys(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      End If   ! L_CO2_interactive

! ----------------------------------------------------------------------
! Section 1.2  Soot cycle.
! ----------------------------------------------------------------------
      If (L_soot) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
         super_tracer_phys(i,j,k,array_size_count) =  soot_new(i,j,k)   &
     &            -super_tracer_phys(i,j,k,array_size_count)
         super_tracer_phys(i,j,k,array_size_count+1) =  soot_agd(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+1)
         super_tracer_phys(i,j,k,array_size_count+2) =  soot_cld(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      End If    ! L_soot

! ----------------------------------------------------------------------
! Section 1.3  Biomass aerosol.
! ----------------------------------------------------------------------
      If (l_Biomass) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  bmass_new(i,j,k)   &
     &            -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  bmass_agd(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  bmass_cld(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      End If  ! l_biomass

! ----------------------------------------------------------------------
! Section 1.4  sulphur cycle.
! ----------------------------------------------------------------------
      If (L_SULPC_SO2) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  so4_aitken(i,j,k)  &
     &             -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  so4_accu(i,j,k)  &
     &             -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  so4_diss(i,j,k)  &
     &             -super_tracer_phys(i,j,k,array_size_count+2)
        super_tracer_phys(i,j,k,array_size_count+3) =  so2(i,j,k)       &
     &             -super_tracer_phys(i,j,k,array_size_count+3)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +3

        If (L_sulpc_nh3) Then
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
          super_tracer_phys(i,j,k,array_size_count) =  nh3(i,j,k)       &
     &            -super_tracer_phys(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End If  ! L_sulpc_nh3

        If (L_sulpc_dms) Then
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  dms(i,j,k)         &
     &            -super_tracer_phys(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End If  ! L_sulpc_dms
      End If  ! L_sulpc_SO2

! ----------------------------------------------------------------------
! Section 1.5  mineral dust.
! ----------------------------------------------------------------------
      If (L_dust) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  dust_div1(i,j,k)   &
     &            -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  dust_div2(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  dust_div3(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+2)
        super_tracer_phys(i,j,k,array_size_count+3) =  dust_div4(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+3)
        super_tracer_phys(i,j,k,array_size_count+4) =  dust_div5(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+4)
        super_tracer_phys(i,j,k,array_size_count+5) =  dust_div6(i,j,k) &
     &            -super_tracer_phys(i,j,k,array_size_count+5)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +5
      End If    ! L_dust

! ----------------------------------------------------------------------
! New addition. Fossil-fuel organic carbon aerosol.
! ----------------------------------------------------------------------
      If (l_ocff) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  ocff_new(i,j,k)    &
     &            -super_tracer_phys(i,j,k,array_size_count)
        super_tracer_phys(i,j,k,array_size_count+1) =  ocff_agd(i,j,k)  &
     &            -super_tracer_phys(i,j,k,array_size_count+1)
        super_tracer_phys(i,j,k,array_size_count+2) =  ocff_cld(i,j,k)  &
     &            -super_tracer_phys(i,j,k,array_size_count+2)
            End Do
          End Do
        End Do
        array_size_count=array_size_count +2
      End If  ! l_ocff

! ----------------------------------------------------------------------
! Section 1.5.1  cariolle ozone tracer.
! ----------------------------------------------------------------------
      If (L_USE_CARIOLLE) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  OZONE_TRACER(i,j,k)&
                 -super_tracer_phys(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      End If    ! L_USE_CARIOLLE

! ----------------------------------------------------------------------
!  ANY NEW NAMED TRACER SPECIES SHOULD BE ADDED HERE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.6.a  free tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      If (model_levels == tr_levels .and. tr_vars > 0) Then
        Do count=1,tr_vars
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) = tracer(i,j,k,count) &
     &            -super_tracer_phys(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End Do
      End If  ! tr_vars > 0

! ----------------------------------------------------------------------
! Section 1.7  ukca tracers  (model_levels=tr_levels)
! ----------------------------------------------------------------------
      If (model_levels == tr_levels .and. tr_ukca > 0) Then
        Do count=1,tr_ukca
          array_size_count=array_size_count +1
          Do k = 1, model_levels
            Do j = 1, rows
              Do i = 1, row_length
                super_tracer_phys(i,j,k,array_size_count) =             &
     &             tracer_ukca(i,j,k,count)                             &
     &            -super_tracer_phys(i,j,k,array_size_count)
              End Do
            End Do
          End Do
        End Do
      End If  ! tr_ukca > 0

! ----------------------------------------------------------------------
! Section 1.99  Murk cycle.  This must be the last Full level field in
!                           the super_array
! ----------------------------------------------------------------------
      If (L_Murk_advect) Then
        array_size_count=array_size_count +1
        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
        super_tracer_phys(i,j,k,array_size_count) =  murk(i,j,k)        &
     &            -super_tracer_phys(i,j,k,array_size_count)
            End Do
          End Do
        End Do
      End If  ! L_Murk_advect

! ----------------------------------------------------------------------
! Section 1.7.b  free tracers  (model_levels/=tr_levels)
! ----------------------------------------------------------------------
      If (tr_vars > 0 .and. tr_levels /= model_levels) Then
        Do l = 1 , tr_vars
          Do k = 1, tr_levels
            Do j = 1, rows
              Do i = 1, row_length
                tracer_phys2(i,j,k,l) = tracer(i,j,k,l)                 &
     &            -tracer_phys2(i,j,k,l)
              End Do
            End Do
          End Do
        End Do
      End If  ! tr_vars > 0

      Return
      END SUBROUTINE TR_Reset
#endif
