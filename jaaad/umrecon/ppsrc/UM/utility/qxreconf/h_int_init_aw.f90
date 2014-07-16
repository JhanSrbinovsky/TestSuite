
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Initialises arrays used in area weighed horizontal interpolation.
      SUBROUTINE H_INT_INIT_AW(ICOF,IDIM,P_FIELD_OUT                    &
     &,                        P_ROWS_IN,P_ROWS_OUT                     &
     &,                        ROW_LENGTH_IN,ROW_LENGTH_OUT             &
     &,                        U_FIELD_IN,U_FIELD_OUT                   &
     &,                        U_ROWS_IN,U_ROWS_OUT                     &
     &,                        GLOBAL,GRIB,FIXHD_IN,FIXHD_OUT           &
     &,                        REALHD_IN,REALHD_OUT,AW_AREA_BOX         &
     &,                        AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP      &
     &,                       BL_INDEX_B_L,BL_INDEX_B_R,BL_INDEX_NEAREST&
     &,                        AW_COLAT_T,AW_LONG_L                     &
     &,                        WEIGHT_T_R,WEIGHT_B_R                    &
     &,                        WEIGHT_T_L,WEIGHT_B_L)
!
! Subroutine Interface:

      Use Ereport_Mod, Only :                                           &
     &    Ereport

      Use Rcf_Grid_Type_Mod, Only :                                     &
     &    Input_Grid,                                                   &
     &    Output_Grid

      IMPLICIT NONE
!
! Description:
!   Initialises arrays used in horizontal interpolation.
!   This replaces routine SETWTS1 (A Dickinson) whose function it
!   incorporates.
!
! Method:
!   Sets up gather index and weight arrays for later call to H_INT_BL h
!   Also sets up rotation coefficients for use in ROTATE with rotated
!   grids.
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.0      12/04/95  Original code. D.M. Goddard
! 4.1      12/06/96  Extended to cope with C_grid. D.M. Goddard
!LL  5.1   10/04/00   New reconfiguration &  S->N ordering support.
!LL                   P.Selwood.
!LL  5.3  24/09/01  Portability changes.    Z. Gardner
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v7 programming standards.
!
! System component covered: S121
! System Task:              S1
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
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

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER      ICOF             !Second dimension of coefficents ary
      INTEGER      IDIM             !Second dimension of index arrays
      INTEGER      P_FIELD_OUT      !No of P pts on target grid
      INTEGER      P_ROWS_IN        !No of P rows on source grid
      INTEGER      P_ROWS_OUT       !No of P rows on target grid
      INTEGER      ROW_LENGTH_IN    !No of pts per row on source grid
      INTEGER      ROW_LENGTH_OUT   !No of pts per row on target grid
                                    ! grid
      INTEGER      U_FIELD_OUT      !No of U pts on target grid
      INTEGER      U_FIELD_IN       !No of U pts on source grid
      INTEGER      U_ROWS_IN        !No of U rows on source grid
      INTEGER      U_ROWS_OUT       !No of U rows on target grid
      LOGICAL      GLOBAL           !T= Global; F= LAM.
      LOGICAL      GRIB             !=T if winds imported on A-grid

!   Array  arguments with intent(in):
      INTEGER      FIXHD_IN(*)      !Fixed length header for source grid
      INTEGER      FIXHD_OUT(*)     !Fixed length header for target grid
      REAL         REALHD_IN(*)     !Real constants from source grid
      REAL         REALHD_OUT(*)    !Real constants from target grid

!   Array  arguments with intent(Out):

      INTEGER      AW_INDEX_TARG_LHS(ROW_LENGTH_OUT+1,IDIM)
                                    !Index of source box overlapping
                                    !lhs of target grid-box
      INTEGER      AW_INDEX_TARG_TOP(P_ROWS_OUT+1,IDIM)
                                    !Index of source box overlapping
                                    !top of target grid-box
      INTEGER      BL_INDEX_B_L(P_FIELD_OUT,IDIM)
                                    !Gather index for bottom l.h.c of
                                    !source grid box. 1=P-pts; 2=UV-pts
      INTEGER      BL_INDEX_B_R(P_FIELD_OUT,IDIM)
                                    !Gather index for bottom r.h.c of
                                    !source grid box. 1=P-pts; 2=UV-pts
      INTEGER      BL_INDEX_NEAREST(P_FIELD_OUT)
                                    !Gather index for nearest point on
                                    !source grid for each target P-pt
      REAL         AW_AREA_BOX(IDIM) !area of grid box in sq units of
                                    !  source grid
      REAL         AW_COLAT_T(P_ROWS_OUT+1,IDIM)
                                    !Colatitude of top of target grd-box
                                    ! (in units of DELTA_LAT_SRCE)
      REAL         AW_LONG_L(ROW_LENGTH_OUT+1,IDIM)
                                    !Left longitude of target grid-box
                                    ! (in units of DELTA_LONG_SRCE)
      REAL         WEIGHT_T_R(P_FIELD_OUT,IDIM) ! Weights for bilinear
      REAL         WEIGHT_B_R(P_FIELD_OUT,IDIM) !\horizontal interpolatn
      REAL         WEIGHT_T_L(P_FIELD_OUT,IDIM) !/ 1=P-pts; 2=U-pts;
      REAL         WEIGHT_B_L(P_FIELD_OUT,IDIM) ! 3=V-pts;4=zonal mean


      INTEGER      ErrorStatus          ! Error flag (0 = OK)

! Local parameters:
      Character (Len=*), Parameter :: RoutineName='HINTIAW'

! Local scalars:
      INTEGER      I                !\                                 .
      INTEGER      IJ               ! \ Loop
      INTEGER      J                ! /variables
      INTEGER      K                !/
      REAL         DELTA_LAT_SOURCE !\                                 .
      REAL         DELTA_LAT_TARGET ! \Grid spacing
      REAL         DELTA_LON_SOURCE ! /
      REAL         DELTA_LON_TARGET !/
      REAL         NPOLE_LAT_SOURCE !\                                 .
      REAL         NPOLE_LAT_TARGET ! \North pole coordinatest
      REAL         NPOLE_LON_SOURCE ! /
      REAL         NPOLE_LON_TARGET !/
      REAL         START_LAT_SOURCE !\                                 .
      REAL         START_LAT_TARGET ! \Coordinates of first data point
      REAL         START_LON_SOURCE ! /
      REAL         START_LON_TARGET !/
      LOGICAL      CYCLIC           !T= Data cyclic
      CHARACTER (LEN=80) :: Cmessage

! Local dynamic arrays:
      INTEGER      I_GRID_IN(ICOF)
      INTEGER      I_GRID_OUT(ICOF)
      REAL         D_LAT_IN(ICOF)
      REAL         D_LON_IN(ICOF)
      REAL         D_LAT_OUT(ICOF)
      REAL         D_LON_OUT(ICOF)
      REAL         LAMBDA_IN(ROW_LENGTH_IN)
                                  ! Latitude coords of source p-grid
      REAL         LAMBDA_OUT(P_FIELD_OUT)
                                    ! Latitude coords of target grid
      REAL         PHI_IN(P_ROWS_IN)
                                    ! Longitude coords of source p-grid
      REAL         PHI_OUT(P_FIELD_OUT)
                                    ! Longitude coords of target grid
! Function & Subroutine calls:
      External H_INT_CO,NEAR_PT,BOX_BND

!- End of header

! 1: Test if requested horizontal interpolation is sensible

! 1.1: Hemispheric or LAM -> global not allowed
      IF(FIXHD_OUT(ITYPE) == 0.AND.FIXHD_IN(ITYPE) >  0)THEN
        WRITE(6,'('' *ERROR* Trying to interpolate from a hemispheric'' &
     &,   '' or LAM to a global domain'')')
        Cmessage = 'Cannot interpolate from hemisphere or LAM to global'
        ErrorStatus = 2
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

! 1.2: LAM -> hemispheric not allowed
      IF(FIXHD_OUT(ITYPE) <  3.AND.FIXHD_IN(ITYPE) >  2)THEN
        WRITE(6,'('' *ERROR* Trying to interpolate from a limited area''&
     &,           '' domain to a global or hemispheric domain'')')
        Cmessage = 'Cannot interpolate from LAM to global or hemisphere'
        ErrorStatus = 3
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      END IF

! 2: Initialise local constants

! 2.1: Grid spacing
      DELTA_LAT_SOURCE=REALHD_IN(IDLAT)
      DELTA_LAT_TARGET=REALHD_OUT(IDLAT)
      DELTA_LON_SOURCE=REALHD_IN(IDLON)
      DELTA_LON_TARGET=REALHD_OUT(IDLON)

! 2.2: Coordinates of north pole on grid
      NPOLE_LAT_SOURCE=REALHD_IN(IPLAT)
      NPOLE_LAT_TARGET=REALHD_OUT(IPLAT)
      NPOLE_LON_SOURCE=REALHD_IN(IPLON)
      NPOLE_LON_TARGET=REALHD_OUT(IPLON)

! 2.3: Coordinates of top left hand  p-point on grid
      START_LAT_SOURCE=REALHD_IN(ISLAT)
      START_LAT_TARGET=REALHD_OUT(ISLAT)
      START_LON_SOURCE=REALHD_IN(ISLON)
      START_LON_TARGET=REALHD_OUT(ISLON)

  ! 2.4: Logical to indicate if input data cyclic
      CYCLIC=FIXHD_IN(ITYPE) <  3

  ! 2.5: Initialise I_GRID_IN and IGRID_OUT to 1 ie p-grid
      I_GRID_IN(1)=1
      I_GRID_IN(2)=1
      I_GRID_OUT(1)=1
      I_GRID_OUT(2)=1

! 3: Weights and indices for P points:

  ! 3.1: If source or target grids have different poles
  !      abort with error message

      IF(NPOLE_LAT_SOURCE /= NPOLE_LAT_TARGET.AND.                      &
     &   NPOLE_LON_SOURCE /= NPOLE_LON_TARGET)THEN

        WRITE(6,*) 'Source and target grids have different poles'
        WRITE(6,*) 'Reconfigure onto a grid with the same pole as'      &
     &,            'target grid'
        WRITE(6,*) 'before attempting area weighted interpolation'

        Cmessage = 'Grids have different poles!'
        ErrorStatus = 4
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      ENDIF

  ! 3.2: Calculate area weighted indices for
  !      interpolating from the source grid onto the target grid
! DEPENDS ON: box_bnd
      CALL BOX_BND(AW_INDEX_TARG_LHS(1,1),AW_LONG_L(1,1)                &
     &,            AW_INDEX_TARG_TOP(1,1),AW_COLAT_T(1,1)               &
     &,            AW_AREA_BOX(1)                                       &
     &,            ROW_LENGTH_OUT,P_ROWS_OUT                            &
     &,            ROW_LENGTH_IN,P_ROWS_IN                              &
     &,            DELTA_LON_TARGET,DELTA_LAT_TARGET                    &
     &,            START_LON_TARGET,START_LAT_TARGET                    &
     &,            DELTA_LON_SOURCE,DELTA_LAT_SOURCE                    &
     &,            START_LON_SOURCE,START_LAT_SOURCE                    &
     &,              I_GRID_OUT(1),I_GRID_IN(1),GLOBAL)


! 4: Weights and indices for U and V points:

  ! 4.1: Calculate offsets for source winds
      IF(GRIB) THEN

    ! 4.1.1: Source winds on A grid
        D_LAT_IN(1)=0.0
        D_LON_IN(1)=0.0
        D_LAT_IN(2)=0.0
        D_LON_IN(2)=0.0

      ELSEIF(FIXHD_IN(ISTAG) == 3)THEN

    ! 4.1.2: Source winds on C grid
        D_LON_IN(1)=0.5
        D_LAT_IN(1)=0.0
        D_LON_IN(2)=0.0
        D_LAT_IN(2)=0.5
        I_GRID_IN(1)=1
        I_GRID_IN(2)=2

      ELSE

    ! 4.1.3: Source winds on B grid
        D_LAT_IN(1)=0.5
        D_LON_IN(1)=0.5
        D_LAT_IN(2)=0.5
        D_LON_IN(2)=0.5
        I_GRID_IN(1)=2
        I_GRID_IN(2)=2
      END IF

  ! 4.2: Calculate offsets for target winds

      IF(FIXHD_OUT(ISTAG) == 3)THEN

    ! 4.2.1: Target winds on C grid
        D_LON_OUT(1)=0.5
        D_LAT_OUT(1)=0.0
        D_LON_OUT(2)=0.0
        D_LAT_OUT(2)=0.5
        I_GRID_OUT(1)=1
        I_GRID_OUT(2)=2

      ELSE

    ! 4.2.2: Target winds on B grid
        D_LAT_OUT(1)=0.5
        D_LON_OUT(1)=0.5
        D_LAT_OUT(2)=0.5
        D_LON_OUT(2)=0.5
        I_GRID_OUT(1)=2
        I_GRID_OUT(2)=2

      ENDIF

! Loop over u, then v
      DO K=1,2

  ! 4.3: Calculate area weighted indices for

  !      interpolating from the source grid onto the target grid

  ! This should be fixed for all C grid stuff now.

        IF (K  ==  1) THEN            ! U grid
! DEPENDS ON: box_bnd
        CALL BOX_BND(AW_INDEX_TARG_LHS(1,K+1),AW_LONG_L(1,K+1)          &
     &,              AW_INDEX_TARG_TOP(1,K+1),AW_COLAT_T(1,K+1)         &
     &,              AW_AREA_BOX(K+1)                                   &
     &,              Output_Grid % Glob_u_row_length                    &
     &,              Output_Grid % Glob_u_rows                          &
     &,              Input_Grid % Glob_u_row_length                     &
     &,              Input_Grid % Glob_u_rows                           &
     &,              DELTA_LON_TARGET,DELTA_LAT_TARGET                  &
     &,              START_LON_TARGET+D_LON_OUT(K)*DELTA_LON_TARGET     &
     &,              START_LAT_TARGET+D_LAT_OUT(K)*DELTA_LAT_TARGET     &
     &,              DELTA_LON_SOURCE,DELTA_LAT_SOURCE                  &
     &,              START_LON_SOURCE+D_LON_IN(K)*DELTA_LON_SOURCE      &
     &,              START_LAT_SOURCE+D_LAT_IN(K)*DELTA_LAT_SOURCE      &
     &,              I_GRID_OUT(K),I_GRID_IN(K),GLOBAL)

        ELSE                    ! V Grid
! DEPENDS ON: box_bnd
        CALL BOX_BND(AW_INDEX_TARG_LHS(1,K+1),AW_LONG_L(1,K+1)          &
     &,              AW_INDEX_TARG_TOP(1,K+1),AW_COLAT_T(1,K+1)         &
     &,              AW_AREA_BOX(K+1)                                   &
     &,              Output_Grid % Glob_v_row_length                    &
     &,              Output_Grid % Glob_v_rows                          &
     &,              Input_Grid % Glob_v_row_length                     &
     &,              Input_Grid % Glob_v_rows                           &
     &,              DELTA_LON_TARGET,DELTA_LAT_TARGET                  &
     &,              START_LON_TARGET+D_LON_OUT(K)*DELTA_LON_TARGET     &
     &,              START_LAT_TARGET-D_LAT_OUT(K)*DELTA_LAT_TARGET     &
     &,              DELTA_LON_SOURCE,DELTA_LAT_SOURCE                  &
     &,              START_LON_SOURCE+D_LON_IN(K)*DELTA_LON_SOURCE      &
     &,              START_LAT_SOURCE-D_LAT_IN(K)*DELTA_LAT_SOURCE      &
     &,              I_GRID_OUT(K),I_GRID_IN(K),GLOBAL)

        END IF

      END DO

! 5: Weights and indices for zonal mean P points:

! 5.1: Reset I_GRID_IN and I_GRID_OUT for p grid
      I_GRID_IN(1)=1
      I_GRID_OUT(1)=1
! 5.2: Calculate area weighted indices for
  !      interpolating from the source grid onto the target grid
! DEPENDS ON: box_bnd
      CALL BOX_BND(AW_INDEX_TARG_LHS(1,4),AW_LONG_L(1,4)                &
     &,            AW_INDEX_TARG_TOP(1,4),AW_COLAT_T(1,4)               &
     &,            AW_AREA_BOX(4)                                       &
     &,            1,P_ROWS_OUT                                         &
     &,            1,P_ROWS_IN                                          &
     &,            DELTA_LON_TARGET,DELTA_LAT_TARGET                    &
     &,            START_LON_TARGET,START_LAT_TARGET                    &
     &,            DELTA_LON_SOURCE,DELTA_LAT_SOURCE                    &
     &,            START_LON_SOURCE,START_LAT_SOURCE                    &
     &,              I_GRID_OUT(1),I_GRID_IN(1),GLOBAL)


! 6: Weights and indices for Coastal adjustment and Integer fields

  ! 6.1: Lat and lon of target grid
      IJ=0
      DO J=1,P_ROWS_OUT
        DO I=1,ROW_LENGTH_OUT
          IJ=IJ+1
          LAMBDA_OUT(IJ)=START_LON_TARGET+DELTA_LON_TARGET*(I-1)
          PHI_OUT(IJ)=START_LAT_TARGET+DELTA_LAT_TARGET*(J-1)
        END DO
      END DO

  ! 6.2: Lat and lon of source grid
      DO J=1,P_ROWS_IN
        PHI_IN(J)=START_LAT_SOURCE+DELTA_LAT_SOURCE*(J-1)
      END DO
      DO I=1,ROW_LENGTH_IN
        LAMBDA_IN(I)=START_LON_SOURCE+DELTA_LON_SOURCE*(I-1)
      END DO

  ! 6.3: Initialise Indices and weights for Bi-linear interpolation

! DEPENDS ON: h_int_co
      CALL H_INT_CO(BL_INDEX_B_L(1,1),BL_INDEX_B_R(1,1)                 &
     &,             WEIGHT_T_R(1,1),WEIGHT_B_R(1,1)                     &
     &,             WEIGHT_T_L(1,1),WEIGHT_B_L(1,1)                     &
     &,             LAMBDA_IN,PHI_IN,LAMBDA_OUT,PHI_OUT                 &
     &,             ROW_LENGTH_IN,P_ROWS_IN,P_FIELD_OUT,CYCLIC)


  ! 6.8: Initialise index of nearest P-points on source grid
! DEPENDS ON: near_pt
      CALL NEAR_PT(BL_INDEX_B_L(1,1),BL_INDEX_B_R(1,1)                  &
     &,            WEIGHT_T_R(1,1),WEIGHT_B_R(1,1)                      &
     &,            WEIGHT_T_L(1,1),WEIGHT_B_L(1,1)                      &
     &,            P_FIELD_OUT,ROW_LENGTH_IN,BL_INDEX_NEAREST)


      RETURN
      END SUBROUTINE H_INT_INIT_AW
