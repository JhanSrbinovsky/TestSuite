#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the output dump lookup tables

Module Rcf_Setup_Lookup_Mod

!  Subroutine Rcf_Setup_Lookup - sets up lookups for output dump
!
! Description:
!   The lookup headers are filled in - with calculated addressing
!   and other namelist derived variables
!
! Method:
!    UMDP F3 defines the lookups.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_Lookup( Hdr_In, Hdr_Out )

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Exppx_mod, Only : &
    Rcf_Exppx

Use Rcf_Ppx_Info_mod, Only : &
    STM_Record_Type,         &
    NSectP

Use Rcf_Submodel_Mod, Only : &
    Internal_model_list, &
    N_Internal_model

Use Rcf_NRecon_Mod, Only : &
    Recondat_Node,         &
    RecondatList,          &
    DumpProgLevs

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Model_Mod, Only : &
    Necf

Use Rcf_HeadAddress_Mod, Only :&
    FH_HorizGrid,         RC_LatSpacing,        RC_LongSpacing, &
    RC_PoleLat,           RC_PoleLong,          FH_VertCoord,   &
    FH_VertCoord_CP,      RC_FirstLat,          LDC_ZseaTheta,  &
    LDC_CkTheta,          LDC_ZseaRho,          LDC_CkRho

Use Rcf_generate_heights_mod, Only : &
    height_gen_original,             &
    height_gen_smooth

Use Rcf_Readnl_horizont_Mod, Only : &
    Iproj

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Items_Mod, Only : &
    num_items,          area_array,             item_array

Use Rcf_Recon_Mod, Only : &
    Lcal360,          &
    Lozone_zonal,     &
    Dump_Pack,        &
    Rimwidtha,        &
    Var_Recon

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec,        &
    stashcode_ozone,           &
    stashcode_ozone_tracer,    &
    stashcode_o3_prod_loss,    &
    stashcode_o3_p_l_vmr,      &
    stashcode_o3_vmr,          &
    stashcode_o3_p_l_temp,     &
    stashcode_o3_temp,         &
    stashcode_o3_p_l_colo3,    &
    stashcode_o3_colo3,        &
    stashcode_u

Implicit None

! Arguments
Type (Um_Header_type), Intent(In) :: Hdr_In
Type (Um_Header_type), Target     :: Hdr_Out

! Comdecks
#include "c_mdi.h"
#include "cppxref.h"
#include "clookadd.h"
#include "c_model_id.h"

! Local vars.
Integer                :: int_val = 1      ! used for transfers
Integer, Pointer       :: Lookup(:,:)
Type (STM_Record_Type) :: STM_Record
Integer                :: i, j, jj, k ! Loopers
Integer                :: icount     ! Counter for calc. Lookup(40)
Integer                :: sec_item
Integer                :: item
Integer                :: section
Integer                :: model
Integer                :: k_out
Integer                :: whole      ! Whole packing indicator -
                                     ! Lookup(21) - see UMDP F3
Integer                :: n1, n2, n3 ! Packing and compression parts
                                     ! of above
Integer                :: length
Integer                :: start_address
Integer                :: n_levels
Integer                :: n_plevels
Integer                :: lev          ! level
Integer                :: bot_lev      ! bottom level for field
Integer                :: lblev_val    ! value for lblev
Integer                :: area
Integer                :: ppxref_grid_type
Integer                :: ErrorStatus
Integer                :: Area_Expand( Hdr_Out % Len2Lookup )
Character (Len=*), Parameter :: RoutineName='Rcf_Setup_Lookup'
Character (Len=80)     :: Cmessage
Real                   :: depth
Real                   :: level( Hdr_Out % Len1LevDepC )
Real                   :: zsea  !} height defining constants
Real                   :: Ck    !}
Real                   :: riv_bzx
Real                   :: riv_bzy

! Function & Subroutine calls:
Integer            :: get_um_version_id

!------------------------------------------------------------------
! First set all elements of lookup to 0 - except dates which
! come from the input lookup (there must be 1).  If reconfiguring
! for VAR then need to make sure date and time is from an 
! instantaneous field (problem when time accumulated fields are at
! the beginning) therefore find where u velocity is.
!------------------------------------------------------------------
Lookup => Hdr_Out % Lookup
k = 1
If (Var_Recon) Then
  Do i = 1, Hdr_In % Len2Lookup
    If (Hdr_In % Lookup(item_code, i) == stashcode_u) Then
      k = i
      Exit
    End If
  End Do
End If

Do i = lbyr, lbdayd
  Lookup( i, : ) = Hdr_In % Lookup( i, k )
End Do
Lookup( 13 : 45, : )  = 0
Lookup( 46 : 64, : ) = Transfer( 0.0, int_val)

K_OUT=1

! 5: Loop through prognostic items and initialise Lookup
Do J=1,N_INTERNAL_MODEL
  Do JJ=0,NSectP
    recondat_node => RecondatList(J,JJ)
    Do While (Associated(recondat_node % recondat_info)) 
      ! 5.1: Extract addressing and number of levels from linked list
      SEC_ITEM      = recondat_node % recondat_info % sec_item
      N_LEVELS      = recondat_node % recondat_info % rlevs
      LENGTH        = recondat_node % recondat_info % len
      START_ADDRESS = recondat_node % recondat_info % raddress
      N_PLEVELS     = recondat_node % recondat_info % rplevs

      area = 1
      Do k = 1, num_items
        If ( sec_item == item_array(k) ) Then
          area = area_array(k)
        End If
      End Do

      ICOUNT=0
      Do  K=K_OUT,K_OUT+N_LEVELS*N_PLEVELS-1

        area_expand(k) = area

        ! 5.3.1: Initialise Lookup
        Lookup( item_code, K ) = SEC_ITEM
        Lookup( model_code,K )=INTERNAL_MODEL_LIST(J)

        item    = Mod(sec_item,1000)
        section = (Lookup(item_code,K) - item)/1000
        model   = internal_model_list(j)
        STM_record = Rcf_Exppx( model, section, item )

        ! 5.3.2: Set addressing information
        Lookup( lblrec,K)=LENGTH/(N_LEVELS*N_PLEVELS)

        Lookup( naddr, K)=START_ADDRESS+ICOUNT
        ICOUNT=ICOUNT+Lookup( lblrec,K)

        !-----------------------------------------------------------
        !        Calculate levels from level dependent constants
        !        Set levels for multi-level fields only
        !        Levels not set for single level fields
        !-----------------------------------------------------------
        If (N_LEVELS > 1) Then

          ! Set the level
          Call Rcf_Level_Code( STM_record % lb_code, bot_lev,          &
          Output_Grid)

          ! Set level number. For ozone we only use top of the
          ! atmosphere levels. Otherwise we count from the surface.
          lev = Mod(K - K_OUT, N_levels) + bot_lev

          If (lev == 0) Then
            lblev_val = 9999                      ! surface
          Else
            lblev_val = lev
          End If

          If ( section == stashcode_prog_sec .AND.                     & 
               item    == stashcode_ozone) Then
            lev = lev + Output_Grid % model_levels -                   &
            Output_Grid % ozone_levels
          End If
          
          LOOKUP( lblev,K ) = lblev_val

          If ( Output_Grid % height_gen_method == height_gen_smooth) Then

            If (STM_record % lbvc_code == 9 .OR.                  &
            STM_record % lbvc_code ==65 ) Then  ! Hybrid/Eta levels

              If (STM_record % lv_code == ppx_theta_level) Then

                If (lev == Output_Grid % model_levels ) Then ! Top level

                  ! Level is current theta level
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Lower boundary is rho level below
                  zsea = Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                  Ck   = Hdr_Out % LevDepC( lev, LDC_CkRho )
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )

                  ! Upper boundary is rho level above - calculated
                  zsea  = 2.0 * Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta) &
                  - Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                  Ck    = 0.0    ! above 1st constant rho level
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                Else If ( lev > Output_Grid % model_levels ) Then

                  Lookup( bulev, lev )  = Transfer( 0., int_val )
                  Lookup( blev, lev )   = Transfer( 0., int_val )
                  Lookup( brlev, lev )  = Transfer( 0., int_val )
                  Lookup( bhulev, lev ) = Transfer( 0., int_val )
                  Lookup( bhlev, lev )  = Transfer( 0., int_val )
                  Lookup( bhrlev, lev ) = Transfer( 0., int_val )

                Else

                  ! Level is current theta level
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaTheta )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Upper boundary is rho level above
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZseaRho )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkRho )
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                  ! Set the lower boundary
                  If (lev <= 1) Then
                    ! Lowest theta level has a lower boundary of the
                    ! physical surface, as defined in the physics.
                    zsea = 0.0
                    Ck   = 1.0
                  Else
                    ! Lower boundary is rho level below
                    zsea = Hdr_Out % LevDepC( lev, LDC_ZseaRho )
                    Ck   = Hdr_Out % LevDepC( lev, LDC_CkRho )
                  End If
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )
                End If

              Else If (STM_record % lv_code == ppx_rho_level ) Then
                If ( lev == Output_Grid % model_levels + 1) Then

                  ! Level is top rho level - calculate
                  zsea = 2.0 * Hdr_Out % LevDepC( lev, LDC_ZseaTheta ) - &
                  Hdr_Out % LevDepC( lev - 1, LDC_ZseaRho )
                  Ck   = 0.0      ! Above 1st const rho by definition
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Upper level boundary same as level
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                  ! Lower boundary is top theta level
                  zsea = Hdr_Out % LevDepC( lev, LDC_ZseaTheta )
                  Ck   = Hdr_Out % LevDepC( lev, LDC_CkTheta )
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )

                Else

                  ! Level above is Theta level
                  zsea = Hdr_Out % LevDepC( lev+1, LDC_ZSeaTheta )
                  Ck   = Hdr_Out % LevDepC( lev+1, LDC_CkTheta )
                  Lookup( bulev, K )  = Transfer( zsea, int_val )
                  Lookup( bhulev, K ) = Transfer( Ck, int_val )

                  ! Level is Rho level
                  zsea = Hdr_Out % LevDepC( lev, LDC_ZSeaRho )
                  Ck   = Hdr_Out % LevDepC( lev, LDC_CkRho )
                  Lookup( blev, K )  = Transfer( zsea, int_val )
                  Lookup( bhlev, K ) = Transfer( Ck, int_val )

                  ! Level below is Theta level
                  If (lev <=0) Then  ! orography
                    zsea = 0.0
                    Ck   = 1.0
                  Else
                    zsea = Hdr_Out % LevDepC( lev, LDC_ZSeaTheta )
                    Ck   = Hdr_Out % LevDepC( lev, LDC_CkTheta )
                  End If
                  Lookup( brlev, K )  = Transfer( zsea, int_val )
                  Lookup( bhrlev, K ) = Transfer( Ck, int_val )

                End If    ! number of levels
              End If      ! rho levels
            End If        ! Hybrid Levels

          Else If (Output_Grid % height_gen_method ==                  &
          height_gen_original) Then

            ! Alternative lbvc code=65 enabled in case adopted in future
            If (STM_record % lbvc_code == 9 .or. &
            STM_record % lbvc_code ==65 )  Then   ! Hybrid/Eta levels
              Lookup( bhlev, K )  = Transfer( RMDI, int_val )
              Lookup( bhulev, K ) = Transfer( RMDI, int_val )
              Lookup( bhrlev, K ) = Transfer( RMDI, int_val )

              If (STM_record % lv_code == ppx_theta_level) Then

                If ( lev == Output_Grid % model_levels ) Then ! Top level
                  Lookup( bulev, K ) = Transfer(                           &
                  2.0 * Output_Grid % eta_theta_levels( lev ) -    &
                  Output_Grid % eta_rho_levels( lev ), int_val )
                  Lookup( blev, K ) = Transfer(                            &
                  Output_Grid % eta_theta_levels( Lev), int_val )
                  Lookup( brlev, K ) = Transfer(                           &
                  Output_Grid % eta_rho_levels  ( Lev), int_val )
                Else If ( lev > Output_Grid % model_levels ) Then
                  Lookup( bulev, K ) = Transfer( 0., int_val )
                  Lookup( blev, K  ) = Transfer( 0., int_val )
                  Lookup( brlev,K )  = Transfer( 0., int_val )
                Else
                  Lookup( bulev, K ) = Transfer(                           &
                  Output_Grid % eta_rho_levels( Lev + 1), int_val )
                  Lookup( blev, K ) = Transfer(                            &
                  Output_Grid % eta_theta_levels( Lev), int_val )

                  If (lev <= 0) Then
                    Lookup( brlev, K ) = Transfer( 0., int_val )
                  Else
                    Lookup( brlev, K ) = Transfer(                         &
                    Output_Grid % eta_rho_levels  ( Lev), int_val )
                  End If
                End If

              Else If ( STM_record % lv_code == ppx_rho_level ) Then

                If ( lev == Output_Grid % model_levels + 1) Then
                  Lookup( bulev, K ) = Transfer(                           &
                  2.0 * Output_Grid % eta_theta_levels(lev - 1) -  &
                  Output_Grid % eta_rho_levels( lev - 1), int_val )
                  Lookup( blev, K  ) = Lookup( bulev, K )
                  Lookup( brlev,K )  = Transfer(                           &
                  Output_Grid % eta_theta_levels( lev - 1), int_val)
                Else
                  Lookup( bulev, K ) = Transfer(                           &
                  Output_Grid % eta_theta_levels( Lev), int_val )
                  Lookup( blev, K ) = Transfer(                            &
                  Output_Grid % eta_rho_levels  ( Lev), int_val )

                  If ( lev <= 0 ) Then   ! bottom level
                    Lookup( brlev, K ) = Transfer( 1.0, int_val )
                  Else
                    Lookup( brlev, K ) = Transfer(                         &
                    Output_Grid % eta_theta_levels( Lev - 1), int_val)
                  End If
                End If
              Else
                Write (6,*) 'ERROR:- lv_code=',STM_record % lv_code
                Write (6,*) 'model, section, item = ', model, section, item
                ErrorStatus = 10
                Cmessage = 'lv_code not right from STASHmaster'
                Call Ereport(  RoutineName, ErrorStatus, Cmessage )

              End If

            End If   ! Hybrid levels
          End If   ! Height Gen methods
        End If    ! N_LEVELS > 1

        !--------------------------------------------------------
        ! Set pseudo level number if variable is on pseudo levels
        !--------------------------------------------------------
        If (STM_Record % pt_code > 0) Then
          Lookup(lbplev,K) = k - k_out + 1
        Endif

      End Do ! K=K_OUT...
      K_OUT=K_OUT+(N_LEVELS*N_PLEVELS)
      recondat_node => recondat_node % next
    End Do ! While associated
  End Do ! All Sections
End Do ! All internal/sub models
!-------------------------------------------------------------------
! Initialise LOOKUP fields from PPXREF
!-------------------------------------------------------------------

Do K=1, Hdr_Out % Len2Lookup
  ITEM=MOD(Lookup(item_code,K),1000)
  SECTION=(Lookup(item_code,K)-ITEM)/1000
  MODEL=Lookup(model_code,K)
  STM_Record = Rcf_Exppx( model, section, item)

  If ( Hdr_Out % FIXHD( FH_HorizGrid ) < 100) Then
    Lookup( lbcode,K)= 1
    Lookup( lbhem, K)= Hdr_Out % FixHd( FH_HorizGrid )
  Else
    Lookup( lbcode,K)= 101 !100 added for non-standard polar axis
    Lookup( lbhem, K)= Hdr_Out % FixHd( FH_HorizGrid ) - 100
  End If

  ! LBCs have an lbhem value of 99
  If (STM_Record % grid_type == ppx_atm_lbc_theta .OR.          &
  STM_Record % grid_type == ppx_atm_lbc_u      .OR.         &
  STM_Record % grid_type == ppx_atm_lbc_v) Then
    Lookup( lbhem, K)= 99
  End If

  Lookup( lbext,K )  = 0 ! No extra data
  Lookup( lbrel,K )  = 2 ! Header release number currently 2
  Lookup( lbfc,K )   = STM_Record % field_code
  Lookup( lbvc,K )   = STM_Record % lbvc_code
  Lookup( lbegin,K ) = 0
  Lookup( lbnrec,K ) = 0
  Lookup( lbproj,K ) = IPROJ
  Lookup( lbtyp,K )  = STM_Record % cf_fieldcode

  If (Lookup( lblev,K )  ==  0) Then
    Lookup( lblev,K ) = STM_Record % cf_levelcode
  End If

  ! DEPENDS ON: get_um_version_id
  Lookup( lbsrce,K )=get_um_version_id(model_id)

  If (Lookup( data_type,K )  ==  0) Then
    Lookup( data_type,K ) = STM_Record % data_type
  End If

  If (Lookup( lbpack,K )  ==  0) Then
    Lookup( lbpack,K ) = STM_Record % dump_packing
    If (DUMP_PACK.eq.2 .or. DUMP_PACK.eq.3 ) Then
      ! Do not pack data ; Override packing indicator from PPXREF
      N1 = 0   !   No packing
      Lookup( lbpack,K ) = (Lookup( lbpack,K )/10)*10 + N1
    End If
  End If

  Lookup( bmdi,K ) = Transfer( rmdi, int_val )
  Lookup( bmks,K ) = Transfer( 1.0,  int_val )
End Do

!-------------------------------------------------------------------
! Change LOOKUP to allow for change in horizontal dimensions
!-------------------------------------------------------------------

Do K=1, Hdr_Out % Len2Lookup

  ITEM=MOD(Lookup( item_code,K ),1000)
  SECTION=(Lookup( item_code,K )-ITEM)/1000
  MODEL=Lookup( model_code,K )
  STM_Record = Rcf_Exppx( model, section, item)

  If (AREA_Expand(K) == 1) Then

    ! Get N2 and N3 from whole value of LBPACK
    WHOLE=Lookup( lbpack,K )
    N2=MOD(INT(WHOLE/10),10)
    N3=MOD(INT(WHOLE/100),10)
    
    ! Find grid type to calculate grid info
    PPXREF_GRID_TYPE = STM_Record % grid_type

    If (N2 == 2 .AND. N3 == 1) Then
      Lookup( lbrow,K ) = 0
      Lookup( lbnpt,K ) = 0
    ElseIf (N2 == 1 .AND. N3 == 1) Then
      Lookup( lbrow,K ) = 0
      Lookup( lbnpt,K ) = 0
    Else
      ! Only deal with C-grid
      Select Case( ppxref_grid_type )
      Case (ppx_atm_cuall)
        Lookup( lbrow,K ) = Output_Grid % Glob_u_rows
        Lookup( lbnpt,K ) = Output_Grid % Glob_u_row_length

      Case (ppx_atm_cvall)
        Lookup( lbrow,K ) = Output_Grid % Glob_v_rows
        Lookup( lbnpt,K ) = Output_Grid % Glob_v_row_length

      Case (ppx_atm_lbc_theta) ! LBC theta points
        Lookup( lbrow,K ) = rimwidtha
        Lookup( lbnpt,K ) = Output_Grid % Glob_p_row_length

      Case (ppx_atm_lbc_u) ! LBC U points
        Lookup( lbrow,K ) = rimwidtha
        Lookup( lbnpt,K ) = Output_Grid % Glob_u_row_length - 1
                            ! Note this is actual rather than
                            ! oversized and only exists for LAM

      Case (ppx_atm_lbc_v) ! LBC V points
        Lookup( lbrow,K ) = rimwidtha
        Lookup( lbnpt,K ) = Output_Grid % Glob_v_row_length

      Case (ppx_atm_river)
        Lookup( lbrow,K ) = Output_Grid % Glob_r_rows
        Lookup( lbnpt,K ) = Output_Grid % Glob_r_row_length

      Case Default
        Lookup( lbrow,K ) = Output_Grid % Glob_p_rows
        If (LOZONE_ZONAL .AND. SECTION == 0    &     
                         .AND. (  ITEM == stashcode_ozone        &   !STASH = 60 
                         .OR.     ITEM == stashcode_o3_prod_loss &   !STASH = 481
                         .OR.     ITEM == stashcode_o3_p_l_vmr   &   !STASH = 482
                         .OR.     ITEM == stashcode_o3_vmr       &   !STASH = 483
                         .OR.     ITEM == stashcode_o3_p_l_temp  &   !STASH = 484
                         .OR.     ITEM == stashcode_o3_temp      &   !STASH = 485
                         .OR.     ITEM == stashcode_o3_p_l_colo3 &   !STASH = 486
                         .OR.     ITEM == stashcode_o3_colo3)    &   !STASH = 487
                                             ) Then

          Lookup( lbnpt,K ) = 1
        Else
          Lookup( lbnpt,K ) = Output_Grid % Glob_p_row_length
        End If
      End Select
    End If

    Lookup( bplat,K ) = Transfer( Hdr_Out % RealC( RC_PoleLat  ),    &
                                  int_val)
    Lookup( bplon,K ) = Transfer( Hdr_Out % RealC( RC_PoleLong ),    &
                                  int_val)
    Lookup( bgor, K ) = Transfer( 0.                            ,    &
                                  int_val)

    Lookup( bdy,K ) = Transfer( Hdr_Out % RealC(2), int_val )
    Lookup( bzy,K ) = Transfer( Hdr_Out % RealC(3) -               &
                                 Hdr_Out % RealC(2), int_val )

    Lookup( bzx,K ) = Transfer( Hdr_Out % RealC(4) -                  &
                             Hdr_Out % RealC(1), int_val)
    If (Lookup( lbnpt,K ) ==  1) Then
      Lookup( bdx,K ) = Transfer( 360., int_val )
    Else
      Lookup( bdx,K ) = Transfer( Hdr_Out % RealC( RC_LongSpacing),   &
                                  int_val )
    End If

    If (PPXREF_GRID_TYPE == ppx_atm_cuall) Then

      Lookup( bzy,K ) = Transfer( Hdr_Out % RealC(3) -           &
                                 Hdr_Out % RealC(2), int_val )
      Lookup( bzx,K ) = Transfer( Hdr_Out % RealC(4) -           &
                                 Hdr_Out % RealC(1) * .5, int_val )
    Else If (PPXREF_GRID_TYPE == ppx_atm_cvall) Then

      Lookup( bzy,K ) = Transfer( Hdr_Out % RealC(3) -           &
                                 Hdr_Out % RealC(2) *.5, int_val )
      Lookup( bzx,K ) = Transfer( Hdr_Out % RealC(4) -           &
                                 Hdr_Out % RealC(1), int_val )
    Else If (PPXREF_GRID_TYPE == ppx_atm_river) Then
      Lookup( bdy,K ) = Transfer( - 180.0/Output_Grid % Glob_r_rows, &
                            int_val)
      riv_bzy = - Hdr_Out % RealC(3) &
                  - 0.5*(- 180.0/Output_Grid % Glob_r_rows)
      Lookup( bzy,K ) = Transfer( riv_bzy,int_val)
      Lookup( bdx,K ) =                                          &
                  Transfer( 360.0/Output_Grid % Glob_r_row_length, &
                       int_val)
      riv_bzx = Hdr_Out % RealC(4) &
                    - 0.5*(360.0/Output_Grid % Glob_r_row_length)
      Lookup( bzx,K ) = Transfer( riv_bzx, int_val )
    End If
  End If

  ! Set lookup 13 from LCAL360. prefromed in UI prior to vn 3.5
  If (LCAL360) Then
    Lookup( lbtim,K ) = 2
  Else
    Lookup( lbtim,K ) = 1
  End If

End Do

! Set BZY and BZX to RMDI in variable resolution header 
  If (Hdr_Out %RealC(1)< 0) Then
    Do K=1, Hdr_Out % Len2Lookup 
      Lookup( bzy,K ) = Transfer( rmdi, int_val )
      Lookup( bzx,K ) = Transfer( rmdi, int_val )
    End Do
  End If

Return
End Subroutine Rcf_Setup_Lookup
End Module Rcf_Setup_Lookup_Mod

#endif
