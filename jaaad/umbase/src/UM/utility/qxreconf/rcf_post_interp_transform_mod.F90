#if defined(RECON)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs transforms after a field has been interpolated.

Module Rcf_Post_Interp_Transform_Mod

!  Subroutine Rcf_Post_Interp_Transform
!
! Description:
!   Wrapper to perform tranformations/processing on a field after
!   it has been interpolated.
!
! Method:
!   Choice of transform/processing is based on stashcode.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.2   15/11/00   Add checks for negative Qs. P.Selwood
!   5.2   16/11/00   Add checks for negative QCF, QCL and
!                    cloud fractions. D.Robinson
!   5.3   05/09/01   Set range of w levels to zero. D.Robinson
!   5.5   07/02/03   Add checks for negative canopy water.
!                    Also for negative qcf2, qrain, and qgraup.
!                    Refactor checks for minimum. T. White
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Use Rcf_Parvars_mod, Only : &
    mype,                   &
    nproc

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_PrintStatus_Mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

Implicit None

Private field_set_to_min

Contains


Subroutine Rcf_Post_Interp_Transform( output_field, fields_out, &
                                      field_count_out )

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_P_Exner

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_done

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_cca,                 stashcode_cc_lwp,        &
    stashcode_w,                   stashcode_w_adv,         &
    stashcode_exner,               stashcode_q,             &
    stashcode_qcf,                 stashcode_qcl,           &
    stashcode_area_cf,             stashcode_bulk_cf,       &
    stashcode_liquid_cf,           stashcode_frozen_cf,     &
    stashcode_mean_canopyw,        stashcode_can_water_tile,&
    stashcode_qcf2,                stashcode_qrain,         &
    stashcode_qgraup,              stashcode_prog_sec,      &
    stashcode_tracer_sec,          stashcode_ukca_sec      

Use Rcf_Recon_Mod, Only : &
    q_min,                &
    model_levels,         &
    w_zero_start,         &
    w_zero_end

Use Ereport_Mod, Only : &
    Ereport

! Arguments
Type( field_type ), Intent( InOut ) :: output_field
Type( field_type ), Pointer         :: fields_out(:)
Integer,            Intent(In)      :: field_count_out

! Local variables
Integer                             :: level
Character (Len=*), Parameter        :: RoutineName='Rcf_Post_Interp_Transform'
Character (Len=80)                  :: Cmessage
Integer                             :: ErrorStatus


!---------------------------------------------------------------
! Only do transforms if interpolation is switched on
!---------------------------------------------------------------
If ( output_field % interp == interp_done ) Then

  Select Case( output_field % stashmaster % section )
  
  Case( stashcode_prog_sec )

    ! Which fields do we wish to apply transforms to?
    Select Case( output_field % stashmaster % item )

    Case( stashcode_exner )
      ! convert interpolated P back to Exner
      If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
        Write (6,'(a25)') 'Converting P to exner'
      End If
      Call Rcf_Conv_P_Exner( output_field )

    Case( stashcode_w, stashcode_w_adv )
      ! Zero first level (surface)
      output_field % Data( :,1 ) = 0.0
      If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
        Write (6,*) ' Setting w to zero, level 0'
      End If

      Do level = w_zero_start+1, w_zero_end+1
        output_field % Data( :, level ) = 0.0
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
          Write (6,*) ' Setting w to zero, level ',level-1
        End If
      End Do

      If (w_zero_end /= model_levels) Then
      output_field % Data( :, output_field % levels ) = 0.0
        If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
          Write (6,*) ' Setting w to zero, level ',model_levels
        End If
      End If

    Case( stashcode_cca )
      Call field_set_to_min( "Convective Cloud Amount", &
                              output_field )

    Case( stashcode_cc_lwp )
      Call field_set_to_min( "Liquid Water Path", &
                              output_field )

    Case( stashcode_q )
      Call field_set_to_min( "Specific Humidity", &
                              output_field,       &
                              field_min=q_min )

    Case( stashcode_qcf )
      Call field_set_to_min( "QCF", &
                              output_field )

    Case( stashcode_qcl )
      Call field_set_to_min( "QCL", &
                              output_field )

    Case( stashcode_qcf2)
      Call field_set_to_min( "QCF2", &
                              output_field )

    Case( stashcode_qrain)
      Call field_set_to_min( "QRain", &
                              output_field )

    Case( stashcode_qgraup)
      Call field_set_to_min( "Qgraup", &
                              output_field )

    Case( stashcode_area_cf )
      Call field_set_to_min( "Area Cloud Fraction", &
                              output_field )

    Case( stashcode_liquid_cf )
      Call field_set_to_min( "Liquid Cloud Fraction", &
                              output_field )

    Case( stashcode_bulk_cf )
      Call field_set_to_min( "Bulk Cloud Fraction", &
                              output_field )

    Case( stashcode_frozen_cf )
      Call field_set_to_min( "Frozen Cloud Fraction", &
                              output_field )

    Case( stashcode_mean_canopyw )
      Call field_set_to_min( "Canopy Water", &
                              output_field )

    Case( stashcode_can_water_tile )
      Call field_set_to_min( "Canopy Water on Tiles", &
                              output_field )

    End Select
  End Select
End If

Return
End Subroutine Rcf_Post_Interp_Transform


Subroutine field_set_to_min( field_name, output_field, &
                             field_min)
! Iterates across whole of output_field, resetting any
! values below the minimum to be at the minimum.
! Default minimum is zero.

! Arguments
Character( * ),     Intent( In )      :: field_name
Type( field_type ), Intent( InOut )   :: output_field

!Optional Argument
Real, Optional,     Intent( In )      :: field_min

! Local variables
Integer                               :: i
Integer                               :: j
Integer                               :: count
Integer                               :: istat

Real                                  :: min

! Default minimum is zero
If( Present( field_min ) ) Then
   min = field_min
Else
   min = 0.0
Endif

! Check for fields too low and reset to minimum.
Do j = 1, output_field % levels
   count = 0
   Do i = 1, output_field % level_size
      If ( output_field % Data(i,j) < min ) Then
         output_field % Data(i,j) = min
         count = count + 1
      End If
   End Do

   Call gc_isum (1, nproc, istat, count)
   If (PrintStatus >= PrStatus_Normal .AND. mype == 0) Then
      write (6,*) ' Level ',j, &
           trim(field_name)//' reset to minimum ',count
   End If

End Do

Return
End Subroutine field_set_to_min
End Module Rcf_Post_Interp_Transform_Mod
#endif
