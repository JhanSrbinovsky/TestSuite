#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ sets up the output grid land-sea mask

Module Rcf_Setup_LSM_out_mod

!  Subroutine Rcf_Setup_LSM_Out - output land-sea mask setup
!
! Description:
! This subroutine sets up the Land-Sea Mask for the output grid
! (if not ancillary) and computes gather indexes etc for coastal
! adjustment should these be required. This needs to be done for
! each PE (with whole grid) as for the weights for horizontal
! interpolation as this is where the coastal adjustment occurs.
!
! Method:
!   All pes work out the global mask (if not ancillary).
!   This is scattered to create the local mask.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains


Subroutine Rcf_Setup_LSM_Out( fields_in, field_count_in, fields_out,  &
                              field_count_out )

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_lsm,             &
    stashcode_prog_sec

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_DecompTP_Mod, Only : &
    decomp_rcf_output

Use Rcf_FreeUMhdr_Mod, Only : &
    Rcf_FreeUMhdr

Use Rcf_Items_Mod, Only : &
    Num_items,        &
    Item_Array,       &
    Source_Array

Use Rcf_Parvars_mod, Only : &
    mype,               &
    nproc,              &
    current_decomp_type,&
    gc_all_proc_group

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_ReadUMhdr_mod, Only : &
    Rcf_ReadUMhdr

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Lsm_Mod, Only : &
    glob_lsm_in,    &
    glob_lsm_out,   &
    local_lsm_out,  &
    glob_land_out,  &
    local_land_out

Use Rcf_Recon_Mod, Only : &
    Lspiral_S

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    Rcf_Free_Unit

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Change_Decomposition_Mod, Only : &
    Rcf_Change_Decomposition

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Interp_Weights_Mod, Only : &
    bl_index_b_l,               bl_index_b_r,              &
    weight_b_l,                 weight_t_r,                &
    weight_t_l,                 weight_b_r

Use Rcf_Lsm_Mod, Only : &
    N_Coastal_Points,           coast_index_in,            &
    coast_index_out,            index_targ_land,           &
    index_targ_sea,             land_unres_index,          &
    sea_unres_index,            n_land_points_unres,       &
    n_sea_points_unres

Implicit None

! Arguments
Integer, Intent(In)              :: field_count_in
Integer, Intent(In)              :: field_count_out
Type (field_type), Pointer       :: fields_in( : )
Type (field_type), Pointer       :: fields_out( : )

! Note that the two fields for lsm contain sizes but *NOT* data!!!
! Data is held in the lsm module.

! Comdecks
#include "cppxref.h"

! Local Data
Integer                :: i           ! Looper
Integer                :: ij          ! Looper
Integer                :: k           ! Looper
Integer                :: kk          ! Looper
Integer                :: ijmin       ! used for output regulation
Integer                :: pos         ! position in fields array
Integer                :: msg         ! tag for comms
Integer                :: dump_pos_tmp
Integer                :: ErrorStatus
Integer                :: Land_Points ! local counter
Integer                :: Orig_Decomp ! Temporary for decomp change
Integer, Allocatable   :: int_lsm_in(:)  ! integer version
Integer, Allocatable   :: int_lsm_out(:) ! integer version
Integer, Allocatable   :: land_sea_index(:) ! temp indexing
Logical                :: LSM         ! is output lsm ancillary?
Logical, Allocatable   :: Land_sea_Mask_temp(:) ! temp space
Type (Um_header_type)  :: hdr_anc     ! Header for ancillary
Character (Len=*), Parameter :: RoutineName = 'Setup_LSM_Out'
Character (Len=80)           :: Cmessage
Type (field_type), Pointer   :: lsm_in      ! input lsm field
Type (field_type), Pointer   :: lsm_out     ! output lsm field


! External f77 routines
External Coast_AJ, GC_IBcast, File_Open, File_Close

!----------------------------------------------------------------
! Locate the lsm fields from fields arrays
!----------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields_in, field_count_in, pos, .TRUE. )

If (pos == 0) Then
  ErrorStatus = -10
  Cmessage = 'Land-Sea Mask is not in input file'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

  Nullify( lsm_in )
Else
  lsm_in => fields_in( pos )
End If

Call Rcf_Locate( stashcode_prog_sec, stashcode_lsm,                  &
                 fields_out, field_count_out, pos, .TRUE.)

If (pos == 0 ) Then
  Nullify( lsm_out )
Else
  lsm_out => fields_out( pos )
End If

!------------------------------------------------------------------
! Set decomposition to rcf_output as all computation will be for
! output lsm
!------------------------------------------------------------------
orig_decomp = current_decomp_type
If ( orig_decomp /= decomp_rcf_output ) Then
  Call Rcf_Change_Decomposition( decomp_rcf_output )
End If

!------------------------------------------------------------------
! Only need to set the output LSM et al if it is required in the
! output dump!
!-----------------------------------------------------------------
If ( Associated( lsm_out ) ) Then

!------------------------------------------------------------------
! Check that output LSM memory is allocated - rather belt and braces
!-----------------------------------------------------------------
  If (.NOT. Allocated( Local_LSM_Out ) ) Then
    ErrorStatus = 30
    Cmessage = 'Local output LSM space not allocated for output!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  If (.NOT. Allocated( Glob_LSM_Out ) ) Then
    ErrorStatus = 40
    Cmessage = 'Global output LSM space not allocated for output!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

!----------------------------------------------------------------
! Establish if the LSM has to be read in from  ancillary -
! set the LSM flag
!----------------------------------------------------------------

  LSM = .FALSE.
  Do i = 1, Num_Items
    If ( Item_Array(i) == stashcode_lsm .AND. Source_Array(i) == 2) Then
      LSM = .TRUE.
    End If
  End Do

!---------------------------------------------------------------
! Do we have to read in ancillary LSM? If so, do it here
!---------------------------------------------------------------
  If (LSM) Then

    If ( PrintStatus >= PrStatus_Normal .AND. mype == 0 ) Then
      Write (6,*) 'Reading in Ancillary Land-Sea Mask'
    End If

    Call Rcf_Get_Unit( Hdr_Anc % UnitNum )

! DEPENDS ON: file_open
    Call File_Open( Hdr_Anc % UnitNum,'MASK',4,0,0,ErrorStatus)
    If ( ErrorStatus /= 0 ) Then
      Cmessage = 'Error Opening Land/Sea Mask Ancillary'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    Call Rcf_ReadUMhdr( Hdr_Anc )

! Read Local LSM
! Note that a temporary overwrite of the dump position is done
! as the lsm_out field doesn't correspond to the ancillary dump...
    dump_pos_tmp = lsm_out % dump_pos
    lsm_out % dump_pos = 1
    Call Rcf_Alloc_Field( lsm_out )

    Call Rcf_Read_Field( lsm_out, Hdr_Anc,  decomp_rcf_output )
    Local_Lsm_Out(:) = lsm_out % Data_Log(:,1)

    Call Rcf_Dealloc_Field( lsm_out )
    lsm_out % dump_pos = dump_pos_tmp

! Do the comms to get it as a global LSM on all PEs
    Call Rcf_Gather_Field( local_lsm_out, glob_lsm_out,               &
                           lsm_out % row_len, lsm_out % rows,         &
                           lsm_out % glob_row_len,                    &
                           lsm_out % glob_rows, 0,                    &
                           gc_all_proc_group )

  ! Note that this communication should ideally use GC_BBcast and
  ! use Kind to determine number of bytes, but we'll not bother.
  ! After all Rcf_Gather_Field (above) and the I/O assume that a
  ! Logical is the same size as a Real as is an Integer. We will do
  ! likewise.

    msg = 801
    Call GC_IBcast( msg, lsm_out % glob_level_size, 0, nproc,  &
                    ErrorStatus, glob_lsm_out )

    If ( ErrorStatus /= 0 ) Then
      Cmessage = 'Problem broadcasting global land-sea mask from PE0'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

! DEPENDS ON: file_close
    Call File_Close( Hdr_Anc % UnitNum,'MASK',4,0,0,ErrorStatus )
    If ( ErrorStatus /= 0 ) Then
      Cmessage = 'Problem closing Land-Sea ancillary'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    Call Rcf_Free_Unit( Hdr_Anc % UnitNum )
    Call Rcf_FreeUMhdr( Hdr_Anc )

  End If

!----------------------------------------------------------------
! Have we changed resolution or LSM? If so, turn on horizontal
! interpolation.
!----------------------------------------------------------------
  If (LSM .AND. (.NOT. h_int_active) ) Then

    If ( Associated( lsm_in ) ) Then
      If ( lsm_out % glob_level_size /= lsm_in % glob_level_size ) Then

        h_int_active = .TRUE.
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) Then
          Write(6,*) 'Horizontal interpolation is switched on because &
                    &of a change in resolution'
        End If

      Else

       Do i = 1, lsm_out % glob_level_size
         If ( glob_lsm_out( i ) .neqv. glob_lsm_in( i ) ) Then
           If (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) Then
             Write(6,*) 'Horizontal interpolation is switched on &
                        &because of a change in Land-Sea Mask'
           End If
           h_int_active = .TRUE.
           Exit
         End If
       End Do

      End If
    End If
  End If


!-------------------------------------------------------------
! Temporary space for integer version of lsms
!-------------------------------------------------------------
  If ( Associated( lsm_in ) ) Then
    Allocate( int_lsm_in( lsm_in % glob_level_size ) )
  Else
    Allocate( int_lsm_in( lsm_out % glob_level_size ) )
  End If
  Allocate( int_lsm_out( lsm_out % glob_level_size ) )
  Allocate( land_sea_mask_temp( lsm_out % glob_level_size ) )
  Allocate( land_sea_index( lsm_out % glob_level_size ) )

!------------------------------------------------------------
! Set up integer versions of LSM (in and out)
!------------------------------------------------------------

  If ( Associated (lsm_in ) ) Then
    Do i = 1, lsm_in % glob_level_size
      If ( glob_lsm_in( i ) ) Then
        int_lsm_in( i ) = 1
      Else
        int_lsm_in( i ) = 0
      End If
    End Do
  Else
    int_lsm_in( : ) = 0
  End If

  Do i = 1, lsm_out % glob_level_size
    If ( glob_lsm_out( i ) ) Then
      int_lsm_out( i ) = 1
    Else
      int_lsm_out( i ) = 0
    End If
  End Do

!--------------------------------------------------------------
! If interpolation is turned on we need to calculate gather
! indexes etc. for coastal adjustment
!--------------------------------------------------------------
  If ( h_int_active .AND. Associated( lsm_in) ) Then

!-------------------------------------------------------------
! Allocate space required for gather indexes
!-------------------------------------------------------------
    ! Global
    Allocate( Coast_Index_In( lsm_out % glob_level_size ) )
    Allocate( Coast_Index_Out( lsm_out % glob_level_size ) )
    Allocate( Index_Targ_Land( lsm_out % glob_level_size ) )
    Allocate( Index_Targ_Sea( lsm_out % glob_level_size ) )

! Note that the following 2 arrays are only used for non-spiral
! adjustment - thus spiral saves memory...
    If ( .NOT. LSpiral_S ) Then
      Allocate( Land_Unres_Index( lsm_out % glob_level_size ) )
      Allocate( Sea_Unres_Index( lsm_out % glob_level_size ) )
    End If

! DEPENDS ON: coast_aj
    Call Coast_AJ( bl_index_b_l, bl_index_b_r, weight_t_r, weight_b_r, &
                   weight_t_l, weight_b_l, lsm_in % glob_row_len,      &
                   lsm_in % glob_rows, lsm_out % glob_level_size,      &
                   int_lsm_in, int_lsm_out, coast_index_out,           &
                   coast_index_in, n_coastal_points, lsm,              &
                   index_targ_sea, n_sea_points_unres, index_targ_land,&
                   n_land_points_unres )

    If ( .NOT. lsm ) Then    ! Coast_aj will have estimated a lsm for us
      Do i = 1, lsm_out % glob_level_size
        If ( int_lsm_out( i ) == 0 ) Then
          glob_lsm_out( i ) = .FALSE.
        Else
          glob_lsm_out( i ) = .TRUE.
        End If
      End Do
    End If

    If (PrintStatus >= PrStatus_Normal .AND. mype == 0 ) Then
      WRITE(6,'('' COASTAL PTS ='',I10)')N_COASTAL_POINTS
      WRITE(6,'('' UNRES SEA PTS ='',I10)')N_SEA_POINTS_UNRES
      WRITE(6,'('' UNRES LAND PTS ='',I10)')N_LAND_POINTS_UNRES
    End If

    ! Print out the input land-sea mask if appropriate.
    If ( PrintStatus >= PrStatus_Diag .AND.  mype == 0) Then
      Write (6,'(/,'' Input Land Sea Mask.'')')
      ij = lsm_in % glob_row_len * (lsm_in % glob_rows - 1 ) + 1
      Do k = lsm_in % glob_rows, 1, -1
        ijmin = min( ij+149, ij+ lsm_in % glob_row_len - 1)
        Write (6,'('' '',150I1)')( int_lsm_in(i),i = IJ,IJMIN)
        ij = ij - lsm_in % glob_row_len
      End Do
    End If

!------------------------------------------------------------------
! Set up gather indices to satify unresolved land and sea points
!------------------------------------------------------------------

! Compute gather index for sea points minus unresolved points

    If (.NOT.LSPIRAL_S) Then
      Do I=1,lsm_out % glob_level_size
        LAND_SEA_MASK_TEMP(I)=.NOT.glob_lsm_out(I)
      End Do
      Do I=1,N_SEA_POINTS_UNRES
        If(.NOT.glob_lsm_out(INDEX_TARG_SEA(I))) Then
          LAND_SEA_MASK_TEMP(INDEX_TARG_SEA(I))=.FALSE.
        End If
      End Do

      LAND_POINTS = 0
      Do I=1,lsm_out % glob_level_size
        If(LAND_SEA_MASK_TEMP(I))Then
          LAND_POINTS=LAND_POINTS + 1
          LAND_SEA_INDEX(LAND_POINTS) = I
        End If
      End Do

      If (land_points < 0 .OR. land_points > lsm_out % glob_level_size)&
                                                     Then
        Write (6,*) 'land_points = ', land_points
        Cmessage = 'Value for land_points is not valid'
        ErrorStatus =  45
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

! Assign each unresolved sea pt to nearest non-unresolved sea pt

      Do I=1,N_SEA_POINTS_UNRES

        If(INDEX_TARG_SEA(I).LE.LAND_SEA_INDEX(1))Then
          SEA_UNRES_INDEX(I)=LAND_SEA_INDEX(1)

        Else If(INDEX_TARG_SEA(I).GT.LAND_SEA_INDEX(LAND_POINTS))Then

          SEA_UNRES_INDEX(I)=LAND_SEA_INDEX(LAND_POINTS)

        ELSE

          Do KK=1,LAND_POINTS-1
            If (INDEX_TARG_SEA(I) .GE. LAND_SEA_INDEX(KK) .AND.        &
                INDEX_TARG_SEA(I).LT.LAND_SEA_INDEX(KK+1)) Then
              SEA_UNRES_INDEX(I)=LAND_SEA_INDEX(KK)
            End If
          End Do

        End If
      End Do

! Compute gather index for land points minus unresolved points

      Do I=1,lsm_out % glob_level_size
        LAND_SEA_MASK_TEMP(I)=glob_lsm_out(I)
      End Do
      Do I=1,N_LAND_POINTS_UNRES
        If (glob_lsm_out(INDEX_TARG_LAND(I))) Then
          LAND_SEA_MASK_TEMP(INDEX_TARG_LAND(I))=.FALSE.
        End If
      End Do

      LAND_POINTS = 0
      Do I=1,lsm_out % glob_level_size
        If(LAND_SEA_MASK_TEMP(I))Then
          LAND_POINTS=LAND_POINTS + 1
          LAND_SEA_INDEX(LAND_POINTS) = I
        End If
      End Do

! Assign each unresolved land pt to nearest non-unresolved land pt

      Do I=1,N_LAND_POINTS_UNRES

        If(INDEX_TARG_LAND(I).LE.LAND_SEA_INDEX(1))Then
          LAND_UNRES_INDEX(I)=LAND_SEA_INDEX(1)
        Else If(INDEX_TARG_LAND(I).GT.LAND_SEA_INDEX(LAND_POINTS))Then
         LAND_UNRES_INDEX(I)=LAND_SEA_INDEX(LAND_POINTS)
        ELSE

          Do KK=1,LAND_POINTS-1
            If(INDEX_TARG_LAND(I).GE.LAND_SEA_INDEX(KK).AND.        &
              INDEX_TARG_LAND(I).LT.LAND_SEA_INDEX(KK+1))Then
              LAND_UNRES_INDEX(I)=LAND_SEA_INDEX(KK)
            End If
          End Do

        End If
      End Do

    End If

  Else If (.NOT. LSM) Then  ! h_int_active (or no input)

! Reuse the input Land-Sea Mask
    Do i = 1, lsm_out % glob_level_size
      If ( int_lsm_in( i ) == 0 ) Then
        glob_lsm_out( i ) = .FALSE.
      Else
        glob_lsm_out( i ) = .TRUE.
      End If
    End Do


  End If     ! h_int_active


!--------------------------------------------------------------
! Check number of land points against that specified in umui
!--------------------------------------------------------------
  glob_land_out = 0
  Do i = 1, lsm_out % glob_level_size
    If ( glob_lsm_out( i ) ) Then
      glob_land_out = glob_land_out + 1
    End If
  End Do

  If ( glob_land_out /= Output_Grid % glob_land_field) Then

    Write (6,*) 'Reconfiguration Error'
    Write (6,*) 'No of land points in output land_sea mask     = ', &
                glob_land_out
    Write (6,*) 'No of land points specified in namelist RECON = ', &
                Output_Grid % glob_land_field
    Write (6,*) 'Please reprocess the job with the correct number of &
                &land points in the UMUI panal'

    ErrorStatus = 50
    Cmessage='Number of land points does not agree with input namelist!'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

  End If

!----------------------------------------------------------------
! Print out the output land-sea mask if appropriate.
!----------------------------------------------------------------
  If ( PrintStatus >= PrStatus_Diag .AND. Associated( lsm_in) &
       .AND. mype == 0) Then
    WRITE (6,'(/,'' Output Land Sea Mask.'')')
    ij = lsm_out % glob_row_len * (lsm_out % glob_rows - 1 ) + 1
    Do k = lsm_out % glob_rows, 1, -1
      ijmin = min( ij+149, ij+ lsm_out % glob_row_len - 1)
      If ( h_int_active ) Then
        WRITE(6,'('' '',150I1)')( int_lsm_out(i),i = IJ,IJMIN)
      Else
        WRITE(6,'('' '',150I1)')( int_lsm_in(i),i = IJ,IJMIN)
      End If
      ij = ij - lsm_out % glob_row_len
    End Do
  End If

!-------------------------------------------------------------
! Scatter the output grid from PE 0 (global should be same on
! all PEs) across LPG and setup local size.
! Only need to do this if NOT ancillary lsm.
!-------------------------------------------------------------
  If ( .NOT. LSM ) Then
    ! Strictly need a real field to do comms, so convert LSM to real
    Call Rcf_Scatter_Field( local_lsm_out, glob_lsm_out,              &
                            lsm_out % row_len, lsm_out % rows,        &
                            lsm_out % glob_row_len,                   &
                            lsm_out % glob_rows, 0, gc_all_proc_group )

  End If

!--------------------------------------------------------------
! Clear up allocated space etc
!--------------------------------------------------------------
  Deallocate( int_lsm_in )
  Deallocate( int_lsm_out )
  Deallocate( land_sea_mask_temp )
  Deallocate( land_sea_index )

End If

!--------------------------------------------------------------
! Count local land-field size
!--------------------------------------------------------------
local_land_out = 0
If ( Associated( lsm_out ) ) Then
  Do i = 1, lsm_out % level_size
    If ( local_lsm_out( i ) ) Then
      local_land_out = local_land_out + 1
    End If
  End Do
End If

If ( PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
  Write (6,*) 'No. of land points is ', glob_land_out
  Write (6,*) 'Local no. of land points is ',local_land_out,&
                ' on PE ', mype
End If

Output_Grid % loc_land_field = local_land_out

!---------------------------------------------------------------
! Set local sizes for land-only fields in fields_out
!---------------------------------------------------------------
Do i = 1, field_count_out
  If ( fields_out( i ) % stashmaster % grid_type ==                 &
                                       ppx_atm_compressed) Then
    If ( .NOT. Associated(lsm_out) ) Then
      ErrorStatus = 60
      Cmessage = 'Output field requested is land-packed but there is&
                 & no output Land-Sea Mask'
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    End If

    fields_out( i ) % level_size = local_land_out
  End If
End Do

!----------------------------------------------------------------
! Change decomposition back to original one
!----------------------------------------------------------------
If ( orig_decomp /= current_decomp_type ) Then
  Call Rcf_Change_Decomposition( orig_decomp )
End If

Return
End Subroutine Rcf_Setup_LSM_Out
End Module Rcf_Setup_LSM_Out_Mod
#endif
