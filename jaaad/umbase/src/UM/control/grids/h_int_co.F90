#if defined(C92_2A) || defined(RECON) || defined(VAROPSVER) \
 || defined(MAKEBC) || defined(FRAMES)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE H_INT_CO----------------------------------------------
!LL
!LL  Purpose:  Calculates bi-linear horizontal interpolation
!LL            coefficients and gather indices for interpolating
!LL            between generalised latitude-longitude grids (eg
!LL            global, regional or rotated lat-lon grid) in which the
!LL            gridlength may vary with latitude and/or longitude. The
!LL            interpolation is carried out by subroutine
!LL            H_INT. Gather indices point to bottom left hand
!LL            corner and bottom right hand corner of each grid box on
!LL            source grid enclosing a target point. Two indices are
!LL            needed to cater for east-west (lambda direction) cyclic
!LL            boundaries when the source data is global. If a target po
!LL            falls outside the domain of the source data, one sided
!LL            differencing is used. The source latitude coordinates
!LL            must be supplied in decreasing order. The source long-
!LL            itude coordinates must be supplied in increasing order,
!LL            starting at any value, but not wrapping round. The
!LL            target points may be specified in any order.
!LL
!LL A.Dickinson <- programmer of some or all of previous code or changes
!LL J.Gregory   <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.4   24/06/94  New checks to keep within array bounds for
!LL                   non_cyclic cases. D. Robinson
!     4.5   29/07/98  Optimisation changes for T3E
!                     Author D.M. Goddard
!LL   5.1   10/04/00 South->North ordering fix. P.Selwood.
!LL   5.2   24/07/00  Adjust target longitudes to ensure they are
!LL                   larger than the minimum source longitude
!LL                   even if the source grid is not cyclic. M J Bell
!     5.3   24/09/01 Change parenthesis to get around optimiser
!                    bug. P.Selwood.
!     5.5   17/02/03 Allow Wave model to use boundary code.
!                                                         D.Holmes-Bell
!     6.0   05/09/03 Added new def for use with makebc. R.Sempers
!     6.2   23/11/05 Removed all references to the wavemodel.
!                    T.Edwards
!     6.2   06/12/05 Remove goto statements from loops to
!                    allow vectorisation
!                    Jean-Christophe Rioual/R Sempers
!LL
!LL  Programming standard:
!LL           Unified Model Documentation Paper No 3
!LL           Version No 1 15/1/90
!LL
!LL  System component: S121,S122
!LL
!LL  System task: S1
!LL
!LL  Documentation:
!LL            The interpolation formulae are described in
!LL            unified model on-line documentation paper S1.
!LL
!LL  -----------------------------------------------------------------
!
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE H_INT_CO                                               &
     &(INDEX_B_L,INDEX_B_R,WEIGHT_T_R,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L  &
     &,LAMBDA_SRCE,PHI_SRCE,LAMBDA_TARG,PHI_TARG                        &
     &,POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE,POINTS,CYCLIC)

      IMPLICIT NONE

      INTEGER                                                           &
     & POINTS_LAMBDA_SRCE                                               &
                          !IN Number of lambda points on source grid
     &,POINTS_PHI_SRCE                                                  &
                          !IN Number of phi points on source grid
     &,POINTS                                                           &
                          !IN Total number of points on target grid
     &,INDEX_B_L(POINTS)                                                &
                          !OUT Index of bottom lefthand corner
                          !    of source gridbox
     &,INDEX_B_R(POINTS)  !OUT Index of bottom righthand corner
                          !    of source gridbox

      REAL                                                              &
     & LAMBDA_TARG(POINTS)                                              &
                           !IN Lambda coords of target grid in degrees
                           !   using same rotation as source grid
     &,PHI_TARG(POINTS)                                                 &
                           !IN Phi coords of target grid in degrees
                           !   using same rotation as source grid
     &,WEIGHT_T_R(POINTS)                                               &
                           !OUT Weight applied to value at top right
                           !    hand corner of source gridbox
     &,WEIGHT_B_L(POINTS)                                               &
                           !OUT Weight applied to value at bottom left
                           !    hand corner of source gridbox
     &,WEIGHT_B_R(POINTS)                                               &
                           !OUT Weight applied to value at bottom right
                           !    hand corner of source gridbox
     &,WEIGHT_T_L(POINTS)                                               &
                           !OUT Weight applied to value at top left
                           !    hand corner of source gridbox
     &,LAMBDA_SRCE(POINTS_LAMBDA_SRCE)                                  &
                                       !IN Lambda coords of source grid
                           !    in degrees
     &,PHI_SRCE(POINTS_PHI_SRCE) !IN Phi coords of target grid in degree

      LOGICAL                                                           &
     & CYCLIC              !IN =T, then source data is cyclic
                           !   =F, then source data is non-cyclic

! Local arrays:---------------------------------------------------------
      REAL                                                              &
     & T_LAMBDA(POINTS) !Local value of target longitude

      INTEGER                                                           &
     & IXP1(POINTS)                                                     &
                       !Longitudinal index plus 1
     &,IX(POINTS)                                                       &
                       !Longitudinal index
     &,IY(POINTS)      !Latitudinal index

! new logical arrays used for vectorizing the set up of
! longitudinal and latitudinal indexes

      LOGICAL :: SETX(POINTS), SETY(POINTS)
! External subroutines called:------------------------------------------
! None
!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------
      REAL                                                              &
     & A                                                                &
                      !Longitudinal weight
     &,B              !Latitudinal weight

      INTEGER                                                           &
     & I,J            !Loop indices
! ----------------------------------------------------------------------

!L    1. Initialise arrays

!         1.1 Scale target longitude so that it falls between
!             LAMBDA_SRCE(1) and LAMBDA_SRCE(1)+360

      DO I=1,POINTS
          T_LAMBDA(I)=MOD(((LAMBDA_TARG(I)-LAMBDA_SRCE(1))+720.),360.)  &
     &    +LAMBDA_SRCE(1)
      END DO

      IF(CYCLIC)THEN
        DO I=1,POINTS
          IX(I)=0
          IY(I)=1
          SETX(I)=.FALSE.
          SETY(I)=.FALSE.
        END DO
      ELSE
        DO I=1,POINTS
          IX(I)=1
          IY(I)=1
          SETX(I)=.FALSE.
          SETY(I)=.FALSE.
        END DO
      END IF

!
!L 2. Calculate lat and lon index of bottom left hand corner of
!L    source grid box enclosing each target point.

! Longitude

      DO J=1,POINTS_LAMBDA_SRCE
         DO I=1,POINTS
           IF ((.NOT. SETX(I)) .AND.                                    &
     &             (LAMBDA_SRCE(J)  <=  T_LAMBDA(I))) THEN
                  IX(I)=J
           ELSE
                  SETX(I)=.TRUE.
           END IF
         END DO
      END DO

       DO J=1,POINTS_PHI_SRCE
          DO I=1,POINTS
             IF ((.NOT. SETY(I)) .AND.                                  &
     &               ((PHI_SRCE(J)  <=  PHI_TARG(I)))) THEN
               IY(I)=J
             ELSE
               SETY(I)=.TRUE.
             END IF
          END DO
       END DO


!L  3. Correct 1-D indices for wrap around etc and then calculate
!L     2-D indices of bottom left and bottom right hand corner
!L     of each grid box.

      IF(CYCLIC)THEN
!     3.1 Cyclic case

      DO I=1,POINTS

! Set index for cyclic wrap around
        IF(IX(I) <  1)THEN
          IX(I)=POINTS_LAMBDA_SRCE
          T_LAMBDA(I)=T_LAMBDA(I)+360.
        ENDIF

! Set index for one sided difference if target point to north or
! south of source area.
        IY(I)=MAX(IY(I),1)
        IY(I)=MIN(IY(I),POINTS_PHI_SRCE-1)

! 2-D indices
        INDEX_B_L(I)=IX(I)+(IY(I)-1)*POINTS_LAMBDA_SRCE
        INDEX_B_R(I)=INDEX_B_L(I)+1

! Correct for cyclic boundaries if target point outside source grid.

        IXP1(I)=IX(I)+1
        IF(IX(I) == POINTS_LAMBDA_SRCE)THEN
          INDEX_B_R(I)=INDEX_B_R(I)-POINTS_LAMBDA_SRCE
          IXP1(I)=IXP1(I)-POINTS_LAMBDA_SRCE
        ENDIF

      ENDDO

      ELSE

!     3.2 Non cyclic case

      DO I=1,POINTS

! Set index for one sided difference if outside source area
        IX(I)=MAX(IX(I),1)
        IX(I)=MIN(IX(I),POINTS_LAMBDA_SRCE-1)
        IF (IX(I) <  1) THEN ! IX(I) < 1 if POINTS_LAMBDA_SRCE = 1
          IX(I)=1
        ENDIF

        IXP1(I)=IX(I)+1
        IXP1(I)=MIN(IXP1(I),POINTS_LAMBDA_SRCE)

! Set index for one sided difference if outside source area
        IY(I)=MAX(IY(I),1)
        IY(I)=MIN(IY(I),POINTS_PHI_SRCE-1)
        IF (IY(I) <  1) THEN ! IY(I) < 1 if POINTS_PHI_SRCE = 1
          IY(I)=1
        ENDIF


! 2-D indices
        INDEX_B_L(I)=IX(I)+(IY(I)-1)*POINTS_LAMBDA_SRCE
        INDEX_B_R(I)=INDEX_B_L(I)+1

      ENDDO

      ENDIF

!L 4. Compute interpolation weights

      DO I=1,POINTS

! Calculate basic weights (equation 2.2)
        A=(AMOD(360.+LAMBDA_SRCE(IXP1(I))-LAMBDA_SRCE(IX(I)),360.))
        IF(A /= 0.)THEN
          A=(T_LAMBDA(I)-LAMBDA_SRCE(IX(I)))/A
        ELSE
          A=0.
        ENDIF

        B=(PHI_TARG(I)-PHI_SRCE(IY(I)))/                                &
     &  ABS(PHI_SRCE(IY(I))-PHI_SRCE(IY(I)+1))

! Calculate bi-linear interpolation weights as per equation (2.1)

        WEIGHT_T_R(I)=A*B
        WEIGHT_B_L(I)=(1.-A)*(1.-B)
        WEIGHT_T_L(I)=(1.-A)*B
        WEIGHT_B_R(I)=A*(1.-B)

      ENDDO

      RETURN
      END SUBROUTINE H_INT_CO
#endif
