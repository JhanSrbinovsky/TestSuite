
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Checks that cloud fraction fields are consistent with q fields.

Module Rcf_Cloud_Frac_Chk_Mod

!  Subroutine Rcf_Cloud_Frac_Chk
!
! Description:
!   Ensures that cloud fraction fields are consistent with
!   the humidity (qcf and qcl) fields.
!
! Method:
!   Checks carried out are :-
!   1. Ensure that frozen cloud fraction is zero if QCF is zero.
!   2. Ensure that liquid cloud fraction is zero if QCL is zero.
!   3. Ensure that both area & bulk cloud fractions are zero if both
!      frozen and liquid cloud fractions are zero.
!
! Current Code Owner: D. Robinson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   12/10/00   Original code.  D. Robinson
!   5.3   24/09/01   Correct logic for deallocation of qcl/qcf.
!                    P.Selwood
!   5.5   06/03/03   Correct bug, swapped pos_qcl and pos_qcf. R.Forbes
!   5.5   28/02/03   Include qcf2 (second ice variable) in cloud
!                    fraction consisteny check.      R.M.Forbes
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Cloud_Frac_Chk( fields, field_count, grid, decomp, hdr )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Write_Field_Mod, Only : &
    Rcf_Write_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_done

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_qcf,             &
    stashcode_qcl,             &
    stashcode_qcf2,            &
    stashcode_area_cf,         &
    stashcode_bulk_cf,         &
    stashcode_frozen_cf,       &
    stashcode_liquid_cf,       &
    stashcode_prog_sec

Use Rcf_Parvars_Mod, Only : &
    nproc,                  &
    mype

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,                &
    PrStatus_Diag

Implicit None

! Arguments
Type( field_type ), Pointer        :: fields(:)
Type( grid_type ), Intent(In)      :: grid
Type( um_header_type ), Intent(In) :: hdr
Integer, Intent(In)                :: field_count
Integer, Intent(In)                :: decomp

! Local variables
Integer                            :: pos_qcf
Integer                            :: pos_qcl
Integer                            :: pos_qcf2
Integer                            :: pos_cf_a
Integer                            :: pos_cf_b
Integer                            :: pos_cf_f
Integer                            :: pos_cf_l
Integer                            :: cldfr_a_changed
Integer                            :: cldfr_b_changed
Integer                            :: cldfr_f_changed
Integer                            :: cldfr_l_changed

Integer                            :: i
Integer                            :: lev
Integer                            :: istat
Integer                            :: count
Type( field_type ), Pointer        :: qcf
Type( field_type ), Pointer        :: qcf2
Type( field_type ), Pointer        :: qcl
Type( field_type ), Pointer        :: cldfr_a
Type( field_type ), Pointer        :: cldfr_b
Type( field_type ), Pointer        :: cldfr_f
Type( field_type ), Pointer        :: cldfr_l

! Note that any negative qcf, qcl and cloud fractions are reset
! to zero earlier in Rcf_post_process_atmos

!-------------------------------------------------------------------
! Loacte qcl, qcf and cloud fractions in output dump
!-------------------------------------------------------------------

Call Rcf_Locate ( stashcode_prog_sec, stashcode_qcl,       &
                  fields, field_count, pos_qcl, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_qcf,       &
                  fields, field_count, pos_qcf, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_qcf2,      &
                  fields, field_count, pos_qcf2, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_area_cf,   &
                  fields, field_count, pos_cf_a, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_bulk_cf,   &
                  fields, field_count, pos_cf_b, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_frozen_cf, &
                  fields, field_count, pos_cf_f, .TRUE. )
Call Rcf_Locate ( stashcode_prog_sec, stashcode_liquid_cf, &
                  fields, field_count, pos_cf_l, .TRUE. )

!--------------------------------------------------------------------
! Check that frozen cloud fraction is zero if qcf is zero.
!--------------------------------------------------------------------

! Need to check that level size == for qcl and clf_l ?

cldfr_f_changed = 0

If (pos_qcf /= 0 .AND. pos_cf_f /= 0 ) Then

  qcf => fields( pos_qcf )
  Call Rcf_Alloc_Field( qcf )
  Call Rcf_Read_Field( qcf, hdr, decomp )

  cldfr_f => fields( pos_cf_f )
  Call Rcf_Alloc_Field( cldfr_f )
  Call Rcf_Read_Field( cldfr_f, hdr, decomp )

  If (pos_qcf2 /= 0) Then
    qcf2 => fields( pos_qcf2 )
    Call Rcf_Alloc_Field( qcf2 )
    Call Rcf_Read_Field( qcf2, hdr, decomp )
  End If

  If (qcf     % interp == interp_done .or.   &
      cldfr_f % interp == interp_done ) Then

    Do lev = 1, grid % wet_levels
    count = 0
      Do i = 1, qcf % level_size

        If (pos_qcf2 /= 0) Then ! qcf2 field is present
                                ! so reset cloud fraction only
                                ! if both qcf and qcf2 are zero
          If ( qcf     % data (i,lev) == 0.0  .and.   &
               qcf2    % data (i,lev) == 0.0  .and.   &
               cldfr_f % data (i,lev) >  0.0 ) Then

            cldfr_f % data (i,lev) = 0.0
            cldfr_f_changed = 1
            count = count + 1

          End If

        Else ! qcf2 field not present, only use qcf

          If ( qcf     % data (i,lev) == 0.0  .and.   &
               cldfr_f % data (i,lev) >  0.0 ) Then

            cldfr_f % data (i,lev) = 0.0
            cldfr_f_changed = 1
            count = count + 1

          End If

        End If ! on pos_qcf2


      End Do

      Call gc_isum (1, nproc, istat, count)
      If (PrintStatus >= PrStatus_Diag .AND. mype == 0) Then
        write (6,*) ' Level ',lev, &
        ' Cloud fracs (frozen) reset to zero ',count
      End If

    End Do

  End If

  If (pos_qcf2 /= 0 ) Then
    Call Rcf_Dealloc_Field( qcf2 )   ! qcf2 no longer required
  End If
  
  Call Rcf_Dealloc_Field( qcf )   ! qcf no longer required

End If


!--------------------------------------------------------------------
! Check that liquid cloud fraction is zero if qcl is zero.
!--------------------------------------------------------------------

! Need to check that level size == for qcl and clf_l ?

cldfr_l_changed = 0

If (pos_qcl /= 0 .AND. pos_cf_l /= 0 ) Then

  qcl => fields( pos_qcl )
  Call Rcf_Alloc_Field( qcl )
  Call Rcf_Read_Field( qcl, hdr, decomp )

  cldfr_l => fields( pos_cf_l )
  Call Rcf_Alloc_Field( cldfr_l )
  Call Rcf_Read_Field( cldfr_l, hdr, decomp )

  If (qcl     % interp == interp_done .or.   &
      cldfr_l % interp == interp_done ) Then

    Do lev = 1, grid % wet_levels
      count = 0
      Do i = 1, qcl % level_size

        If ( qcl     % data (i,lev) == 0.0  .and.   &
             cldfr_l % data (i,lev) >  0.0 ) Then

          cldfr_l % data (i,lev) = 0.0
          cldfr_l_changed = 1
          count = count + 1

        End If

      End Do

      Call gc_isum (1, nproc, istat, count)
      If (PrintStatus >= PrStatus_Diag .AND. mype == 0) Then
        write (6,*) ' Level ',lev, &
        ' Cloud fracs (liquid) reset to zero ',count
      End If

    End Do

  End If
  
  Call Rcf_Dealloc_Field( qcl )   ! qcl no longer required

End If


!--------------------------------------------------------------------
! Check that area and bulk cloud fractions are set to zero
! if both frozen and liquid cloud fractions are zero.
!--------------------------------------------------------------------

! need to check that cld_f and cld_l have already been read in ?
! need to check that level_size == for all 4 cf fields ?

cldfr_a_changed = 0
cldfr_b_changed = 0

If (pos_cf_f /= 0 .AND. pos_cf_l /= 0 ) Then
  If (pos_cf_a /= 0 .AND. pos_cf_b /= 0 ) Then

    cldfr_a => fields( pos_cf_a )
    Call Rcf_Alloc_Field( cldfr_a )
    Call Rcf_Read_Field( cldfr_a, hdr, decomp )

    cldfr_b => fields( pos_cf_b )
    Call Rcf_Alloc_Field( cldfr_b )
    Call Rcf_Read_Field( cldfr_b, hdr, decomp )

    If (qcf     % interp == interp_done .or.   &
        qcl     % interp == interp_done .or.   &
        cldfr_a % interp == interp_done .or.   &
        cldfr_b % interp == interp_done .or.   &
        cldfr_f % interp == interp_done .or.   &
        cldfr_l % interp == interp_done ) Then

      Do lev = 1, grid % wet_levels
        count = 0
        Do i = 1, qcl % level_size

          If ( cldfr_f % data (i,lev) == 0.0  .and.   &
               cldfr_l % data (i,lev) == 0.0  .and.   &
              (cldfr_a % data (i,lev) >  0.0  .or.    &
               cldfr_b % data (i,lev) >  0.0) ) Then

            cldfr_a % data (i,lev) = 0.0
            cldfr_b % data (i,lev) = 0.0
            cldfr_a_changed = 1
            cldfr_b_changed = 1
            count = count + 1

          End If

        Enddo

        Call gc_isum (1, nproc, istat, count)
        If (PrintStatus >= PrStatus_Diag .AND. mype == 0) Then
          write (6,*) ' Level ',lev, &
          ' Cloud fracs (area & bulk) reset to zero ',count
        End If

      End Do

    End If

  End If
End If

!---------------------------------------------------------------------
! Synchronise `changed' flags
!---------------------------------------------------------------------
Call GC_Imax( 1, nproc, istat, cldfr_a_changed )
Call GC_Imax( 1, nproc, istat, cldfr_b_changed )
Call GC_Imax( 1, nproc, istat, cldfr_f_changed )
Call GC_Imax( 1, nproc, istat, cldfr_l_changed )

!---------------------------------------------------------------------
! Write out changed fields
!---------------------------------------------------------------------

If (cldfr_a_changed == 1) Then
  Call Rcf_Write_Field( cldfr_a, hdr, decomp )
End If
If (cldfr_b_changed == 1) Then
  Call Rcf_Write_Field( cldfr_b, hdr, decomp )
End If
If (cldfr_f_changed == 1) Then
  Call Rcf_Write_Field( cldfr_f, hdr, decomp )
End If
If (cldfr_l_changed == 1) Then
  Call Rcf_Write_Field( cldfr_l, hdr, decomp )
End If

If (pos_cf_f /= 0 .AND. pos_cf_l /= 0 ) Then
  If (pos_cf_a /= 0 .AND. pos_cf_b /= 0 ) Then
    Call Rcf_Dealloc_Field( cldfr_a )
    Call Rcf_Dealloc_Field( cldfr_b )
  End If
End If

If (pos_qcf /= 0 .AND. pos_cf_f /= 0 ) Then
  Call Rcf_Dealloc_Field( cldfr_f )
End If

If (pos_qcl /= 0 .AND. pos_cf_l /= 0 ) Then
  Call Rcf_Dealloc_Field( cldfr_l )
End If

Return
End Subroutine Rcf_Cloud_Frac_Chk
End Module Rcf_Cloud_Frac_Chk_Mod
