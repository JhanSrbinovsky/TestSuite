
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+Determine STASH input length per vertical level for prog var

Module Rcf_Address_Length_Mod

!  Subroutine Rcf_Address_Length - determines field size
!
! Description:
!    Calculates size of field for output dump addressing (single level).
!
! Method:
!    Calculates field sizes based on Grid_Type from stashmaster
!    Generally single level, but old LBC needs level info.
!
!    Based on UM 4.5/5.0 code.
!
! Current Code Owner: P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Contains

SUBROUTINE Rcf_Address_Length( IGPL, halo_type, LEN )

Use Rcf_Submodel_Mod

Use Rcf_CntlAtm_Mod

Use Rcf_Lsm_Mod, Only : &
    glob_land_out

Use Rcf_Model_Mod, Only :  &
    ZonAvOzone

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_Parvars_Mod

Use Rcf_Recon_Mod, Only : &
    Tr_Vars,              &
    Rimwidtha,            &
    Rimwidtho

Use Rcf_Level_Code_Mod, Only : &
    Rcf_Level_Code

Use Rcf_Readnl_Horizont_Mod, Only : &
    Extended_Halo_Size_EW,          &
    Extended_Halo_Size_NS

IMPLICIT NONE

! Subroutine arguments:
Integer, Intent(In)   :: IGPL      ! Grid code
Integer, Intent(In)   :: halo_type ! code from stashmaster
Integer, Intent(Out)  :: LEN       ! Length

! Local scalars:
INTEGER               :: IP   !  pressure levels
INTEGER               :: IQ   !  wet levels
INTEGER               :: IT   !  tracer levels
INTEGER               :: IX1  !  leftmost point
INTEGER               :: IX2  !  rightmost point
INTEGER               :: IY1  !  lower point
INTEGER               :: IY2  ! upper point

Integer               :: halo_i  ! EW halo size for LBC addressing
Integer               :: halo_j  ! NS halo size for LBC addressing

! Function and subroutine calls:
EXTERNAL LLTORC

!- End of Header ---------------------------------------------------

! Determine row/column nos. for global domain on output grid
! DEPENDS ON: lltorc
CALL LLTORC(IGPL,90,-90,0,360,IY1,IY2,IX1,IX2, Output_Grid)

! What sizes of halo will we have (used for LBC fields)
Select Case( halo_type )
  Case( halo_type_no_halo )
    halo_i = 0
    halo_j = 0

  Case( halo_type_single )
    halo_i = 1
    halo_j = 1

  Case( halo_type_extended )
    halo_i = Extended_Halo_Size_EW
    halo_j = Extended_Halo_Size_NS

End Select

Select Case (IGPL)
  Case ( 21 )
    LEN= glob_land_out     ! Land compressed

  Case ( 22 )
    IF (ZonAvOzone) THEN   !  Zonal
      LEN = IY2-IY1+1
    ELSE                   !  Full fields
      LEN = (IX2-IX1+1)*(IY2-IY1+1)
    ENDIF

  Case (23)
    LEN = Output_Grid % glob_r_field

  Case ( 25 )              ! Old LBC
    CALL Rcf_Level_Code( 2, IP, Output_Grid)   ! pressure levels
    CALL Rcf_Level_Code( 3, IQ, Output_Grid)   ! wet levels
    CALL Rcf_Level_Code(11, IT, Output_Grid)   ! tracer levels
    LEN=( Output_Grid % GLOB_P_ROWS +                               &
          Output_Grid % GLOB_P_ROW_LENGTH-2*RIMWIDTHA)*2*RIMWIDTHA* &
        (1+3*IP+2*IQ+TR_VARS*IT)-2*IP*4*RIMWIDTHA


  Case ( 26 )              ! LBC pressure grid
    Len = 2 * (RimwidthA + halo_j) *                                 &
          (Output_Grid % glob_p_row_length + (2 * halo_i)) +         &
          2 * (RimwidthA + halo_i) *                                 &
          (Output_Grid % glob_p_rows - (2 * RimwidthA))

  Case ( 27 )              ! LBC u  C-grid
    Len = 2 * (RimwidthA + halo_j) *                                 &
          (Output_Grid % glob_p_row_length - 1 + (2 * halo_i)) +     &
          2 * (RimwidthA + halo_i) *                                 &
          (Output_Grid % glob_p_rows - (2 * RimwidthA))

  Case ( 28 )              ! LBC v C-grid
    Len = 2 * (RimwidthA + halo_j) *                                 &
          (Output_Grid % glob_p_row_length + (2 * halo_i)) +         &
          2 * (RimwidthA + halo_i) *                                 &
          (Output_Grid % glob_p_rows - 1 - (2 * RimwidthA))

  Case Default
    LEN =(IX2-IX1+1)*(IY2-IY1+1)

End Select

RETURN
END Subroutine Rcf_Address_Length

End Module Rcf_Address_Length_Mod
