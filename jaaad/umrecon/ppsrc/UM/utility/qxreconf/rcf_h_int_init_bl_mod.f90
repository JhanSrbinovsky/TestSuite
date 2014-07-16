
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialises bilinear interpolation weights

Module Rcf_H_Int_Init_BL_Mod

!  Subroutine Rcf_H_Int_Init_BL
!
! Description:
!   Sets up the interpolation and rotation weights for bilinear
!   horizontal interpolation.
!
! Method:
!   Weights stored in the interp_weights_mod module
!   Seperate calculations for P, U, V and P zonal points.
!
!   Based (with many changes) on UM4.5 code.
!
! Current Code Owner: P.Selwood
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   20/09/01   Make sure comments are not considered as cpp
!                    continuations.  P.Selwood
!   5.5   12/02/03   Correct nested LAM check and clarified
!                    comments. T.White
!   5.5   27/02/03   Corrected interpolation of ozone data from
!                    zonal to full field.  P.Dando
!   6.0   19/06/03   Correction for LAM to LAM when crossing Greenwich
!                    meridian. Roddy Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_H_Int_Init_BL( grid_in, grid_out, hdr_in, hdr_out, &
                              grib )

Use Rcf_Interp_Weights_Mod            ! All of it (well almost)

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Headaddress_Mod, Only :  &
    FH_RowDepCStart

Use RCF_Recon_Mod, Only : &
    LOZONE_ZONAL

Implicit None

! Arguments
Type( grid_type ),      Intent( In )   :: grid_in
Type( grid_type ),      Intent( In )   :: grid_out
Type( um_header_type ), Intent( In )   :: hdr_in
Type( um_header_type ), Intent( In )   :: hdr_out
Logical,                Intent( In )   :: Grib

! Comdecks
!
! Description:
!  Comdeck contains parameters for horizontal interpolation routines
!
! Current Code Owner: D.M.Goddard
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   4.0  18/04/95   Original code. D.M. Goddard
!
! Declarations:

! Global parameters:
      INTEGER      BGRID            !I_GRID_IN/OUT value for B_GRID
      INTEGER      CGRID            !I_GRID_IN/OUT value for C_GRID
      INTEGER      IDLAT            !Posn of lat increment in REALHD
      INTEGER      IDLON            !Posn of lon increment in REALHD
      INTEGER      IPLAT            !Posn of N pole latitide in REALHD
      INTEGER      IPLON            !Posn of N pole longitude in REALHD
      INTEGER      ISLAT            !Posn of start latitide in REALHD
      INTEGER      ISLON            !Posn of start longitude in REALHD
      INTEGER      ISTAG            !Posn of grid staggering in FIXHD
      INTEGER      ITYPE            !Posn of domain type in FIXHD
      REAL         SMALL            !Used in IF test to check that
                                    !  target is within source area.

      PARAMETER(BGRID=1)
      PARAMETER(CGRID=2)
      PARAMETER(IDLAT=2)
      PARAMETER(IDLON=1)
      PARAMETER(IPLAT=5)
      PARAMETER(IPLON=6)
      PARAMETER(ISLAT=3)
      PARAMETER(ISLON=4)
      PARAMETER(ISTAG=9)
      PARAMETER(ITYPE=4)
      PARAMETER(SMALL=0.001)


!- End of COMDECK declaration

! Local variables
Integer            :: i                ! looper
Integer            :: j                ! looper
Integer            :: ij               ! looper
Integer            :: p_points_out     ! Total points in output grid
Integer            :: ErrorStatus
Real               :: DELTA_LAT_SOURCE !\                              .
Real               :: DELTA_LAT_TARGET ! \Grid spacing
Real               :: DELTA_LON_SOURCE ! /
Real               :: DELTA_LON_TARGET !/
Real               :: MAX_LAMBDA_IN    !\                              .
Real               :: MAX_LAMBDA_OUT   ! \                             .
Real               :: MAX_PHI_IN       !  \  Max and min of lat and long
Real               :: MAX_PHI_OUT      !   \ of source and target grids
Real               :: MIN_LAMBDA_IN    !   / when converted to standard
Real               :: MIN_LAMBDA_OUT   !  /  coordinates. Used to check
Real               :: MIN_PHI_IN       ! /   that target area is within
Real               :: MIN_PHI_OUT      !/source area.
Real               :: NPOLE_LAT_SOURCE !\                              .
Real               :: NPOLE_LAT_TARGET ! \North pole coordinatest
Real               :: NPOLE_LON_SOURCE ! /
Real               :: NPOLE_LON_TARGET !/
Real               :: START_LAT_SOURCE !\                              .
Real               :: START_LAT_TARGET ! \Coords of first data point
Real               :: START_LON_SOURCE ! /
Real               :: START_LON_TARGET !/
Real               :: D_LAMBDA_IN( 2 )    !\                           .
Real               :: D_LAMBDA_OUT( 2 )   ! \ Offsets for wind grids
Real               :: D_PHI_IN( 2 )       ! /
Real               :: D_PHI_OUT( 2 )      !/
Real               :: LAMBDA_IN( grid_in % glob_p_row_length )
                                    ! Longitude coords of source p-grid
Real               :: LAMBDA_INN( grid_in % glob_p_field )
                                    ! Longitude coords of source p-grid
Real               :: LAMBDA_OUT( grid_out % glob_p_field )
                                    ! Longitude coords of target grid
Real               :: LAMBDA_ROT( grid_in % glob_p_field )
                                    ! Standard lon coords of source grid
Real               :: LAMBDA_TMP( grid_out % glob_p_field )
                                    ! Standard lon coords of target grid
Real               :: PHI_IN( grid_in % glob_p_rows )
                                    ! Latitude coords of source p-grid
Real               :: PHI_INN( grid_in % glob_p_field )
                                    ! Latitude coords of source p-grid
Real               :: PHI_OUT( grid_out % glob_p_field )
                                    ! Latitude coords of target grid
Real               :: PHI_TMP( grid_out % glob_p_field )
                                    ! Standard lat coords of target grid
Logical            :: Rot_In
Logical            :: Rot_Out
Logical            :: Cyclic
Logical            :: Translate ! Do we need to translate longitudinal
                                ! coords to be contiguous?
Character (Len=*), Parameter     :: RoutineName = 'rcf_h_int_init_bl'
Character (Len=80)               :: Cmessage

External NEAR_PT, H_INT_CO, EQTOLL, LLTOEQ

!--------------------------------------------------------------
! 1: Test if requested horizontal interpolation is sensible
!--------------------------------------------------------------

! 1.1: Hemispheric or LAM -> global not allowed
If ( Hdr_Out % FixHd(ITYPE) == 0 .AND. Hdr_In % FixHd(ITYPE) > 0) Then
  Write(6,'('' *ERROR* Trying to interpolate from a hemispheric &
 & or LAM to a global domain'')')
  Cmessage = 'Cannot interpolation from hemisphere/LAM to global'
  ErrorStatus = 2
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! 1.2: LAM -> hemispheric not allowed
If ( Hdr_Out % FixHd(ITYPE) < 3 .AND. Hdr_In % FixHd(ITYPE) > 2) Then
  Write(6,'('' *ERROR* Trying to interpolate from a limited area &
    & domain to a global or hemispheric domain'')')
  Cmessage = 'Cannot interpolate from LAM to global or hemisphere'
  ErrorStatus = 3
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

!-----------------------------------------------------------------
! 2: Initialise local constants
!-----------------------------------------------------------------

! 2.1: Grid spacing
delta_lat_source = Hdr_In % RealC( idlat )
delta_lat_target = Hdr_Out % RealC( idlat )
delta_lon_source = Hdr_In % RealC( idlon )
delta_lon_target = Hdr_Out % RealC( idlon )

! 2.2: Coordinates of north pole on grid
npole_lat_source = Hdr_In % RealC( iplat )
npole_lat_target = Hdr_Out % RealC( iplat )
npole_lon_source = Hdr_In % RealC( iplon )
npole_lon_target = Hdr_Out % RealC( iplon )

! 2.3: Coordinates of top left hand  p-point on grid
start_lat_source = Hdr_In % RealC( islat )
start_lat_target = Hdr_Out % RealC( islat )
start_lon_source = Hdr_In % RealC( islon )
start_lon_target = Hdr_Out % RealC( islon )

! 2.4: Logical to indicate output grid is rotated
ROT_out = npole_lat_target /= 90. .OR. npole_lon_target /= 0.

! 2.5: Logical to indicate input grid is rotated
ROT_in = npole_lat_source /= 90. .OR. npole_lon_source /= 0.

! 2.6: Logical to indicate if input data cyclic
CYCLIC = Hdr_In % FixHd( ITYPE ) < 3

!-------------------------------------------------------------------
! 3. Weights and indices for P points
!-------------------------------------------------------------------

! 3.1: Lat and lon of new grid
!If output grid is variable resolution
If (hdr_out % FIXHD(FH_RowDepCStart) > 1) Then
  ij=0
    Do j=1,grid_out % glob_p_rows
      Do i=1,grid_out % glob_p_row_length
        ij=ij+1
        lambda_out(ij)=hdr_out % ColDepC(i,1)
        phi_out(ij)   =hdr_out % RowDepC(j,1)
      End Do
     End Do
Else
  ij = 0
  Do j = 1, grid_out % glob_p_rows
    Do i = 1, grid_out % glob_p_row_length
      ij = ij + 1
      lambda_out(ij) = start_lon_target + delta_lon_target * (i-1)
      phi_out(ij)    = start_lat_target + delta_lat_target * (j-1)
    End Do
  End Do
End If 
       
! 3.2: Convert output grid to standard lat-lon if output grid rotated
If (ROT_OUT) Then
! DEPENDS ON: eqtoll
  Call EQTOLL(phi_out, lambda_out, phi_tmp, lambda_tmp,             &
              npole_lat_target, npole_lon_target,                   &
              grid_out % glob_p_field )

  ! 3.2.1: Calculate coefficients for rotating wind
! DEPENDS ON: w_coeff
  Call W_COEFF(COEFF1, COEFF2, lambda_tmp, lambda_out,               &
               npole_lat_target, npole_lon_target,                   &
               grid_out % glob_p_field )

  ! 3.2.2: Copy across latitude and longitude of target grid in standard
  !        coords for next step
  Do i = 1, grid_out % glob_p_field
    lambda_out(i) = lambda_tmp(i)
    phi_out(i)    = phi_tmp(i)
  End Do
End If

! 3.3: Convert output grid to input lat-lon if input grid rotated
If (ROT_IN) Then
! DEPENDS ON: lltoeq
  Call LLTOEQ(phi_out, lambda_out, phi_out, lambda_out,            &
              npole_lat_source, npole_lon_source,                  &
              grid_out % glob_p_field )
End If

! 3.4: LLTOEQ ensures that lambda_out is between 0 and 360

! 3.5: Lat and lon of old grid
Do j = 1, grid_in % glob_p_rows
  phi_in(j)    = start_lat_source + delta_lat_source * (J-1)
End Do
Do I = 1, grid_in % glob_p_row_length
  lambda_in(I) = start_lon_source + delta_lon_source * (I-1)
End Do

! 3.6: Ensure that 0<=lambda_in<=360
If (ROT_IN) Then
  lambda_in = Modulo(lambda_in, 360.0)
End If

! 3.7: Check if target area contained within source area when
!      LAM->LAM.
Translate = .False. !Default is no translation required.
If (ROT_IN) Then

  max_lambda_in = lambda_in( grid_in % glob_p_row_length )
  min_lambda_in = lambda_in( 1 )
  max_phi_in = phi_in( grid_in % glob_p_rows )
  min_phi_in = phi_in( 1 )

! If max < min, our input grid range lies across zero. We
! must rescale so the range is
! max_lambda_in-360<=lambda<=max_lambda_in
  Translate = (max_lambda_in < min_lambda_in)
  If (Translate) Then
    !max_lambda unchanged.
    min_lambda_in = min_lambda_in - 360.0
    Where (lambda_out > max_lambda_in + SMALL) &
      lambda_out = lambda_out - 360.0
    Where (lambda_in  > max_lambda_in + SMALL) &
      lambda_in  = lambda_in  - 360.0
  End If

! At this point, if any of lambda_out lies outside of lambda_in,
! then lambda_out is not contained within lambda_in. However, we
! cannot guarantee that the max/min taken next are true (if the
! output grid spans max_lambda_in, and we have rescaled, they
! are inaccurate. If necessary for the error message, we
! recalculate below)
  max_lambda_out = Maxval(lambda_out)
  min_lambda_out = Minval(lambda_out)
  max_phi_out    = Maxval(phi_out)
  min_phi_out    = Minval(phi_out)

  If (max_phi_in < max_phi_out - SMALL       .OR. &
      min_phi_in > min_phi_out + SMALL       .OR. &
      max_lambda_in < max_lambda_out - SMALL .OR. &
      min_lambda_in > min_lambda_out + SMALL) Then
! Calculate correct values for max/min output grid.
! Put everything on a 0-360 scale.
    max_lambda_in  = Modulo(max_lambda_in,  360.0) ! Noop
    min_lambda_in  = Modulo(min_lambda_in,  360.0)
    max_lambda_out = Modulo(max_lambda_out, 360.0)
    min_lambda_out = Modulo(min_lambda_out, 360.0)
! If lambda_out crosses zero, wrap it left in thirty degree
! increments until it is contiguously within range.
    If ( Abs(max_lambda_out - min_lambda_out) < 2.0*SMALL ) Then
      lambda_out = Modulo(lambda_out, 360.0)
      Do i = 330, 0, -30 !This assumes that the output grid
                         !spans less than 330 degrees longitude.
         Where (lambda_out > Real(i) ) &
           lambda_out = lambda_out - 360.0
         max_lambda_out = Modulo(Maxval(lambda_out), 360.0)
         min_lambda_out = Modulo(Minval(lambda_out), 360.0)
         If ( Abs(max_lambda_out - min_lambda_out) > 2.0*SMALL ) &
           Exit
      End Do
    End If

    WRITE(6,'('' *ERROR* Target LAM not contained within &
             &source LAM'')')
    WRITE(6,'(''max_phi_in,max_phi_out'',2F8.2)')              &
                max_phi_in,max_phi_out
    WRITE(6,'(''min_phi_in,min_phi_out'',2F8.2)')              &
                min_phi_in,min_phi_out
    WRITE(6,'(''max_lambda_in,max_lambda_out'',2F8.2)')        &
                max_lambda_in,max_lambda_out
    WRITE(6,'(''min_lambda_in,min_lambda_out'',2F8.2)')        &
                min_lambda_in,min_lambda_out
    Cmessage = 'Target LAM not within source LAM'
    ErrorStatus=4
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If
EndIf

! 3.7.1: Calculate wind rotation coefficients from source to standard
!      lat_long grid. Used to rotate winds in CONTROL when they are
!      first read in. A grid supported only

If (ROT_in) Then
  IJ=0
  Do J=1, grid_in % glob_p_rows
    Do I=1, grid_in % glob_p_row_length
      IJ=IJ+1
      lambda_inn(IJ) = start_lon_source + delta_lon_source * (I-1)
      phi_inn(IJ)    = start_lat_source + delta_lat_source * (J-1)
    End Do
  End Do

! DEPENDS ON: eqtoll
  Call EQTOLL(phi_inn, lambda_inn, phi_inn,lambda_rot,           &
              npole_lat_source, npole_lon_source,                &
              grid_in % glob_p_field)

! Calculate coefficients for rotating wind
! DEPENDS ON: w_coeff
  Call W_COEFF(COEFF3, COEFF4, lambda_rot, lambda_inn,           &
               npole_lat_source, npole_lon_source,               &
               grid_in % glob_p_field)

End If


! 3.8: Initialise Indices and weights for Bi-linear interpolation
!        As these are required for coastal adjustment, this must be
!        called whether bi-linear interpolation is requested or not.

! 3.8.1: Call H_INT_CO



! DEPENDS ON: h_int_co
Call H_INT_CO(BL_INDEX_B_L(1,1), BL_INDEX_B_R(1,1),           &
              WEIGHT_T_R(1,1), WEIGHT_B_R(1,1),               &
              WEIGHT_T_L(1,1), WEIGHT_B_L(1,1),               &
              lambda_in,phi_in, lambda_out, phi_out,          &
              grid_in % glob_p_row_length,                    &
              grid_in % glob_p_rows, grid_out % glob_p_field, CYCLIC)




! 3.8.2: Initialise index of nearest P-points on source grid
! DEPENDS ON: near_pt
Call NEAR_PT (BL_INDEX_B_L(1,1), BL_INDEX_B_R(1,1),                 &
              WEIGHT_T_R(1,1), WEIGHT_B_R(1,1),                     &
              WEIGHT_T_L(1,1), WEIGHT_B_L(1,1),                     &
              grid_out % glob_p_field, grid_in % glob_p_row_length, &
              BL_INDEX_NEAREST)

!--------------------------------------------------------------------
! 4: Weights and indices for U/V points
!--------------------------------------------------------------------
! 4.0.1: Offsets

If ( Hdr_In % FixHd(ISTAG) == 3) Then
  d_lambda_in(1) = 0.5
  d_phi_in(1)    = 0.0
  d_lambda_in(2) = 0.0
  d_phi_in(2)    = 0.5
Else
  d_lambda_in(1) = 0.5
  d_phi_in(1)    = 0.5
  d_lambda_in(2) = 0.5
  d_phi_in(2)    = 0.5
End If

If ( Hdr_Out % FixHd(ISTAG) == 3) Then
  d_lambda_out(1) = 0.5
  d_phi_out(1)    = 0.0
  d_lambda_out(2) = 0.0
  d_phi_out(2)    = 0.5
Else
  d_lambda_out(1) = 0.5
  d_phi_out(1)    = 0.5
  d_lambda_out(2) = 0.5
  d_phi_out(2)    = 0.5
End If

! A grid - note that only interpolating *from* A grid (not to it)
If (grib) Then
  d_lambda_in(:)  = 0.
  d_phi_in(:)     = 0.
End If

!-------------------------------------------------------------------
! U first
!-------------------------------------------------------------------
!If output grid is variable resolution
If (hdr_out % FIXHD(FH_RowDepCStart) > 1) Then
  ij=0
  Do j=1,grid_out % glob_u_rows
    Do i=1,grid_out % glob_u_row_length
      ij=ij+1
      lambda_out(ij)=hdr_out % ColDepC(i,2)
      phi_out(ij)   =hdr_out % RowDepC(j,1)
    End Do
  End Do
Else
  IJ=0
  Do J = 1, grid_out % glob_u_rows
    Do I = 1, grid_out % glob_u_row_length
      IJ=IJ+1
      lambda_out(IJ) = start_lon_target + delta_lon_target  &
                * (I - 1 + d_lambda_out(1) )
      phi_out(IJ) = start_lat_target + delta_lat_target     &
                * (J - 1 + d_phi_out(1))
    End Do
  End Do 
End If
 
! 4.2: Lat and lon of source grid
Do J = 1, grid_in % glob_u_rows
  phi_in(J) = start_lat_source + delta_lat_source              &
            * (J - 1 + d_phi_in(1))
End Do
Do I = 1, grid_in % glob_u_row_length
  lambda_in(I) = start_lon_source + delta_lon_source           &
               * (I - 1 + d_lambda_in(1))
End Do

! 4.3: Convert to standard lat-lon if new grid rotated
If (ROT_OUT) Then
! DEPENDS ON: eqtoll
  Call EQTOLL(phi_out, lambda_out, phi_out, lambda_out,        &
              npole_lat_target, npole_lon_target,              &
              grid_out % glob_u_field )
End If

! 4.6: Convert target grid to source lat-lon if source grid rotated
If (ROT_IN) Then
! DEPENDS ON: lltoeq
  Call LLTOEQ(phi_out, lambda_out, phi_out, lambda_out,        &
              npole_lat_source, npole_lon_source,              &
              grid_out % glob_u_field )
End If

! 4.7: Scale longitude if LAM target grid
!        to make it monotonically increasing
If (Translate) Then
  Where (lambda_out > max_lambda_out + SMALL) &
    lambda_out = lambda_out - 360.0
End If

! 4.8: Scale longitude if LAM source grid
!     to make it monotonically increasing
If (Translate) Then
  Where (lambda_in > max_lambda_out + SMALL) &
    lambda_in = lambda_in - 360.0
End If

! 4.9: Indices and weights for horizontal interpolation
! DEPENDS ON: h_int_co
CALL H_INT_CO( BL_INDEX_B_L(1,2), BL_INDEX_B_R(1,2),      &
               WEIGHT_T_R(1,2), WEIGHT_B_R(1,2),          &
               WEIGHT_T_L(1,2), WEIGHT_B_L(1,2),          &
               lambda_in, phi_in, lambda_out, phi_out,    &
               grid_in % glob_u_row_length,               &
               grid_in % glob_u_rows, grid_out % glob_u_field, CYCLIC )

!-------------------------------------------------------------------
!  Then V
!-------------------------------------------------------------------
!If output grid is variable resolution
If (hdr_out % FIXHD(FH_RowDepCStart) > 1) Then
  ij=0
  Do j=1,grid_out % glob_V_rows
    Do i=1,grid_out % glob_V_row_length
      ij=ij+1
      lambda_out(ij)=hdr_out % ColDepC(i,1)
      phi_out(ij)   =hdr_out % RowDepC(j,2)
    End Do
  End Do
Else
  IJ=0
  Do J = 1, grid_out % glob_v_rows
    Do I = 1, grid_out % glob_v_row_length
      IJ=IJ+1
      lambda_out(IJ) = start_lon_target + delta_lon_target     &
                   * (I - 1 + d_lambda_out(2) )
      phi_out(IJ) = start_lat_target + delta_lat_target        &
                   * (J - 1 + d_phi_out(2))
    End Do
  End Do
End If
 
! 4.2: Lat and lon of source grid
Do J = 1, grid_in % glob_v_rows
  phi_in(J) = start_lat_source + delta_lat_source              &
            * (J - 1 + d_phi_in(2))
End Do
Do I = 1, grid_in % glob_u_row_length
  lambda_in(I) = start_lon_source + delta_lon_source           &
               * (I - 1 + d_lambda_in(2))
End Do

! 4.3: Convert to standard lat-lon if new grid rotated
If (ROT_OUT) Then
! DEPENDS ON: eqtoll
  Call EQTOLL(phi_out, lambda_out, phi_out, lambda_out,        &
              npole_lat_target, npole_lon_target,              &
              grid_out % glob_v_field )
End If

! 4.6: Convert target grid to source lat-lon if source grid rotated
If (ROT_IN) Then
! DEPENDS ON: lltoeq
  Call LLTOEQ(phi_out, lambda_out, phi_out, lambda_out,        &
              npole_lat_source, npole_lon_source,              &
              grid_out % glob_v_field )
End If

! 4.7: Scale longitude if LAM target grid
!        to make it monotonically increasing
If (Translate) Then
  Where (lambda_out > max_lambda_out + SMALL) &
    lambda_out = lambda_out - 360.0
End If

! 4.8: Scale longitude if LAM source grid
!     to make it monotonically increasing
If (Translate) Then
  Where (lambda_in > max_lambda_out + SMALL) &
    lambda_in = lambda_in - 360.0
End If

! 4.9: Indices and weights for horizontal interpolation
! DEPENDS ON: h_int_co
CALL H_INT_CO( BL_INDEX_B_L(1,3), BL_INDEX_B_R(1,3),      &
               WEIGHT_T_R(1,3), WEIGHT_B_R(1,3),          &
               WEIGHT_T_L(1,3), WEIGHT_B_L(1,3),          &
               lambda_in, phi_in, lambda_out, phi_out,    &
               grid_in % glob_v_row_length,               &
               grid_in % glob_v_rows, grid_out % glob_v_field, CYCLIC )

!--------------------------------------------------------------
! 5: Weights and indices for zonal mean P points:
!--------------------------------------------------------------
! 5.1: Lat and lon of target grid

If ( LOZONE_ZONAL ) Then
!If output grid is variable resolution 
!for phi_out for OZONE_ZONAL. C.Wang 17/05/07
  If (hdr_out % FIXHD(FH_RowDepCStart) > 1) Then
    lambda_out(1) = hdr_out % ColDepC(1,1)
    Do j=1, grid_out % glob_p_rows
      phi_out(j)   =hdr_out % RowDepC(j,1)
    End Do
  Else
    lambda_out(1) = start_lon_target
    Do j=1, grid_out % glob_p_rows
      phi_out(j) = start_lat_target + delta_lat_target * (j - 1)
    End Do
  End If
    
  p_points_out = grid_out % glob_p_rows
Else

! Interpolating Zonal ozone to full field. Set up output grid as
! for any other p field but warn user that the zonal data has been
! copied to full field.
! If output grid is variable resolution
  If (hdr_out % FIXHD(FH_RowDepCStart) > 1) Then
    ij=0
    Do j=1,grid_out % glob_p_rows
      Do i=1,grid_out % glob_p_row_length
          ij=ij+1
          lambda_out(ij)=hdr_out % ColDepC(i,1)
          phi_out(ij)   =hdr_out % RowDepC(j,1)
      End Do
    End Do
  Else
    IJ=0
    Do j = 1, grid_out % glob_p_rows
      Do i = 1, grid_out % glob_p_row_length
        IJ = IJ + 1
        lambda_out(IJ) = start_lon_target + delta_lon_target * (i - 1)
        phi_out(IJ)    = start_lat_target + delta_lat_target * (j - 1)
      End Do
    End Do
  End If
  p_points_out = grid_out % glob_p_field

! Write warning message that ozone has been interpolated from zonal to
! full field.

  Cmessage = 'Interpolating ozone from zonal to full field'
  ErrorStatus = -10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

End If


! 5.2: Lat and lon of source grid
lambda_in(1) = start_lon_source
Do J = 1, grid_in % glob_p_rows
  phi_in(J) = start_lat_source + delta_lat_source * (J-1)
End Do

! 5.3: Convert output grid to standard lat-lon if target grid rotated
If (ROT_OUT) Then
! DEPENDS ON: eqtoll
  Call EQTOLL(phi_out, lambda_out, phi_out, lambda_out,      &
              npole_lat_target, npole_lon_target,            &
              p_points_out )
End If

! 5.4: Convert targget grid to source lat-lon if source grid rotated
If (ROT_IN) Then
! DEPENDS ON: lltoeq
  Call LLTOEQ(phi_out, lambda_out, phi_out, lambda_out,       &
              npole_lat_source, npole_lon_source,             &
              p_points_out )
End If

! 5.5: Initialise indices and weights for horizontal interpolation
! DEPENDS ON: h_int_co
Call H_INT_CO( BL_INDEX_B_L(1,4), BL_INDEX_B_R(1,4),            &
               WEIGHT_T_R(1,4), WEIGHT_B_R(1,4),                &
               WEIGHT_T_L(1,4), WEIGHT_B_L(1,4),                &
               lambda_in, phi_in, lambda_out, phi_out,          &
               1, grid_in % glob_p_rows, p_points_out,&
               cyclic )




Return
End Subroutine Rcf_H_Int_Init_BL
End Module Rcf_H_Int_Init_BL_Mod
