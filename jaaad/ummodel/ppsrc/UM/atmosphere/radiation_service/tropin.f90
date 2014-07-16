



! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  SUBROUTINE TROPIN------------------------------------------------
!LL
!LL  Purpose:  Finds the tropopause & returns index of where it is.
!LL      Not suitable for single-column model use, as it does
!LL             horizontal filling-in, & includes dynamical allocation.
!LL
!LL    Based on routine TROP, but
!LL      1) taking temperature rather than theta as input
!LL      2) rather than returning the temperature, pressure and height
!LL      of a continuously-varying tropopause found by extrapolating
!LL      lapse rates, it returns the index of the adjacent layer
!LL      boundary (N meaning the bottom of the Nth model layer)
!LL      3) "filling in" where no tropopause is found.
!LL
!LL        Author:  William Ingram
!LL
!LL  Model                Modification history:
!LL version  Date
!LL   4.2   23/9/96       New code, based on TROP
!LL   6.2   27/01/06    Code added to allow SCM to call TROPIN.
!LL                     Default SCM tropopause model level based
!LL                     on mean tropopause height varying by lat
!LL                     -itude.                           R.Wong
!       6.2             21/02/06   Updefs Added for version
!                                  control of radiation code
!                                            (J.-C. Thelen)
!LL
!LL    Note that the definition of "tropopause" matches TROP, which
!LL      is not quite the WMO definition, though the same critical
!LL      lapse rate is used (unless comdeck C_LAPSE is altered).
!LL      For details of the interpolation assumptions, see UMDP S1
!LL      section 3.2.2 or Swinbank & Wilson (1990: SRFRTN 48).
!LL      Any physical changes to one routine should be considered for
!LL                                         mirroring in the other.
!LLEND----------------------------------------------------------------
!
!*L  Arguments:-------------------------------------------------------
      SUBROUTINE TROPIN(T, Exner_rho_levels, Exner_theta_levels,        &
     &                  ROW_LENGTH, Rows, P_LEVELS, off_x, off_y,       &
     &                  at_extremity,true_latitude,height_theta,        &
     &                  MIN_TROP_LEVEL, MAX_TROP_LEVEL, IT)

      IMPLICIT NONE

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north, sout
                         ! east or west of the processor grid

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

      INTEGER                                                           &
             !, INTENT (IN)
     &     ROW_LENGTH,                                                  &
                                        !   Number of points per row
     &     ROWS,                                                        &
                                        !   Number of rows
     &     P_LEVELS,                                                    &
                                        !   Number of model levels
     &     MIN_TROP_LEVEL,                                              &
     &     MAX_TROP_LEVEL,                                              &
     &     off_x, off_y                 ! halo sizes
!     ! Limits on where the tropopause may be deemed to be - between
!     !  the MIN_TROP_LEVELth and MAX_TROP_LEVELth layers (with the
!     !  convention used here for layer boundaries, the actual index
!     !  IT returned has MIN_TROP_LEVEL < IT =< MAX_TROP_LEVEL.)

      REAL                                                              &
          !, INTENT(IN)
 !  Temperature at layer centres
     &     T(1-off_x:row_length+off_x, 1-off_y:rows+off_y, P_LEVELS)    &
     &,    EXNER_rho_levels(1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y,P_LEVELS)                &
     &,    EXNER_theta_levels(1-off_x:row_length+off_x,                 &
     &                      1-off_y:rows+off_y,P_LEVELS)                &
     &,    true_latitude(1,1)                                           &
                                           ! Latitude in radians
                                           ! for SCM
     &,    height_theta(1,1,0:p_levels)    ! For SCM


      INTEGER                                                           &
             !, INTENT (OUT)
     &     IT(row_length, rows)
!     ! Integer indexing the tropopause, taken to be @ a layer boundary
!     !   with the convention that N means the bottom of layer N.

! Workspace usage:-----------------------------------------------------

      INTEGER                                                           &
             !, INTENT (OUT)
     &     IT_work(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

      REAL LAPSE_RATE(1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &                MIN_TROP_LEVEL+1:MAX_TROP_LEVEL+1)
!     ! Lapse rate between layer centres

      LOGICAL LTROP(1-off_x:row_length+off_x, 1-off_y:rows+off_y)
!     !  Logical array which indicates whether we are still seeking a
!     !    tropopause (if not, it has already been found lower down)

!*---------------------------------------------------------------------
! Define local variables:----------------------------------------------

      INTEGER I, J, K, j0, j1, j0f, j1f,                                &
                                          ! Loopers over level & point
     &     POINT,                                                       &
                            ! Point counter at ends of rows
     &     KP1,                                                         &
!     !  K+1, except where this would cause out-of-bounds reference
     &     NNEIGH,                                                      &
                            ! Number of well-defined tropopauses among
!     ! the 8 nearest neighbours of a point without one of its own
     &     FILLIN,                                                      &
!     !  Used to fill in points without a clearly-defined tropopause
     &     DTI              ! Default tropopause index for points where
!     ! not even one nearest neighbour has a well-defined tropopause

      REAL                                                              &
     &     DEL_EXNER_J,                                                 &
                            ! Differences of Exner pressure across
     &     DEL_EXNER_JM1,                                               &
                            !                           half layers
     &     DENOM,                                                       &
                            ! Denominator in lapse rate expression
     &     Z_TROP           ! Height in metres of default tropopause
                            ! in SCUM code.
!*L------------------COMDECK C_PI---------------------------------------
!LL
!LL 4.0 19/09/95  New value for PI. Old value incorrect
!LL               from 12th decimal place. D. Robinson
!LL 5.1 7/03/00   Fixed/Free format P.Selwood
!LL

      ! Pi
      Real, Parameter :: Pi                 = 3.14159265358979323846

      ! Conversion factor degrees to radians
      Real, Parameter :: Pi_Over_180        = Pi/180.0

      ! Conversion factor radians to degrees
      Real, Parameter :: Recip_Pi_Over_180  = 180.0/Pi

!*----------------------------------------------------------------------

!*L------------------COMDECK C_G----------------------------------------
! G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE

      Real, Parameter :: G = 9.80665

!*----------------------------------------------------------------------
!*L------------------COMDECK C_R_CP-------------------------------------
! History:
! Version  Date      Comment.
!  5.0  07/05/99  Add variable P_zero for consistency with
!                 conversion to C-P 'C' dynamics grid. R. Rawlins
!  5.1  07/03/00  Fixed/Free format conversion   P. Selwood

! R IS GAS CONSTANT FOR DRY AIR
! CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
! PREF IS REFERENCE SURFACE PRESSURE

      Real, Parameter  :: R      = 287.05
      Real, Parameter  :: CP     = 1005.
      Real, Parameter  :: Kappa  = R/CP
      Real, Parameter  :: Pref   = 100000.

      ! Reference surface pressure = PREF
      Real, Parameter  :: P_zero = Pref
!*----------------------------------------------------------------------
!*L------------------COMDECK C_LAPSE ----------------------------------
      Real, Parameter :: lapse      = 0.0065  ! Near surface lapse rate
      Real, Parameter :: lapse_trop = 0.002   ! Tropopause lapse rate
!*---------------------------------------------------------------------

      REAL CP_OVER_G, P_EXNER_500, P_EXNER_50
      PARAMETER ( CP_OVER_G = CP / G )

!----------------------------------------------------------------------

!     ! 1.  Set up local constants and initialise arrays

! set up j0, j1 loop limits. AT north and south boundary ignore halo
! rows

      j0 = 1 - off_y
      j1 = rows + off_y
      j0f = 1
      j1f = rows
      If (at_extremity(PSouth)) Then
      j0 = 1
      j0f = 2
      End If
      If (at_extremity(PNorth)) Then
      j1 = rows
      j1f = rows - 1
      End If
      P_EXNER_500 = (500./1000.)**KAPPA
      P_EXNER_50  =  (50./1000.)**KAPPA

      DTI = ( MIN_TROP_LEVEL + MAX_TROP_LEVEL ) / 2

      DO J = j0, j1
        DO I=1-off_x, row_length+off_x
          LTROP(I,J) = .TRUE.
        END DO
      ENDDO

!L    ! Compute lapse rate between full levels: equation 3.16, UMDP S1

      DO K=MIN_TROP_LEVEL+1, MIN(MAX_TROP_LEVEL+1,P_LEVELS)
        DO J = j0, j1
          DO I=1-off_x, row_length+off_x
!         ! Exner pressure difference across half layers
            DEL_EXNER_J = exner_rho_levels(i,j,k) /                     &
     &                    exner_theta_levels(i,j,k) - 1
            DEL_EXNER_JM1 = 1 - exner_rho_levels(i,j,k)                 &
     &                          / exner_theta_levels(i,j,k-1)
!         ! Denominator
            DENOM = T(I,j,k-1) * DEL_EXNER_JM1                          &
     &            + T(I,J,k) * DEL_EXNER_J
!         !  Lapse rate between level k-1 and k
            LAPSE_RATE(I,J,k) = ( T(I,J,k-1) - T(I,J,k) ) /             &
     &                          ( CP_OVER_G * DENOM )
          END DO
        END DO
      END DO

!L    ! 2.  Find level of tropopause, where it is well defined

      DO K=MIN_TROP_LEVEL+1, MAX_TROP_LEVEL

! 'K+1' level for lapse rate test; allows K iteration up to P_LEVELS
        KP1=MIN(K+1,P_LEVELS)

        DO J = j0, j1
          DO I=1-off_x, row_length+off_x

!         ! Not-quite-WMO criteria for interval containing tropopause
!         ! (where 'interval' stretches between layer centres k and k-1)

            IF ( exner_theta_levels(i,j,k-1)  >   P_EXNER_50 .AND.      &
     &           exner_theta_levels(i,j,k)  <   P_EXNER_500 .AND.       &
     &           LAPSE_RATE(I,J,k)  <   LAPSE_TROP .AND.                &
     &           LAPSE_RATE(I,J,kp1)  <   LAPSE_TROP .AND. LTROP(I,J) ) &
     &      THEN
              LTROP(I,J)=.FALSE.
              IT_work(I,J) = K
            ENDIF
          END DO
        END DO
      END DO

!L    ! 4.  Fill in it array from it_work, where set, else
!L          where the above criteria did not find a tropopause

!     !  Run through all internal points and where no tropopause was
!     !   found, set the level to the average of those found at the 8
!     !   surrounding points.  If none of these did find one, then DTI,
!     !   the middle of the permitted range, is set.
!     !   (Using integers means rounding down - probably the best thing
!     !   overall as the lack of vertical resolution to pick out a
!     !   tropopause precisely is likely to mean they'll be diagnosed
!     !   too high, if anything, at neighbouring points.  Similarly
!     !   it's not worth worrying about level number being non-linear
!     !   in height or pressure at those points where ipso facto the
!     !   result is so arbitrary.)





! Loop over non-halo points
! Don't do poles as search over near neighbours is expensive.
! If poles not set then set to dti

      DO J= j0f, j1f
        DO I= 1, ROW_LENGTH

          IF ( LTROP(I,J) ) THEN

            FILLIN = 0
            NNEIGH = 0
            IF ( .NOT. LTROP(i-1,j-1) ) THEN
              FILLIN = FILLIN + IT_work(i-1,j-1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i,j-1) ) THEN
              FILLIN = FILLIN + IT_work(i,j-1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i+1,j-1) ) THEN
              FILLIN = FILLIN + IT_work(i+1,j-1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i-1,j) ) THEN
              FILLIN = FILLIN + IT_work(i-1,j)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i+1,j) ) THEN
              FILLIN = FILLIN + IT_work(i+1,j)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i-1,j+1) ) THEN
              FILLIN = FILLIN + IT_work(i-1,J+1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i,j+1) ) THEN
              FILLIN = FILLIN + IT_work(i,j+1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(i+1,j+1) ) THEN
              FILLIN = FILLIN + IT_work(i+1,j+1)
              NNEIGH = NNEIGH + 1
            ENDIF

            IF ( NNEIGH  ==  0 ) THEN
               IT(i,j) = DTI
             ELSE
               IT(i,j) = FILLIN / NNEIGH
            END IF

          Else
            IT(i,j) = it_work(i,j)
          END IF

        END DO
      END DO

! South pole loop: only executed if j0f > 1
      DO J= 1, j0f - 1
        DO I= 1, ROW_LENGTH
          IF ( LTROP(I,J) ) THEN
            IT(i,j) = DTI
          Else
            IT(i,j) = it_work(i,j)
          END IF
        END DO
      END DO

! North pole loop: only executed if j1f < rows
      DO J= j1f+1, rows
        DO I= 1, ROW_LENGTH
          IF ( LTROP(I,J) ) THEN
            IT(i,j) = DTI
          Else
            IT(i,j) = it_work(i,j)
          END IF
        END DO
      END DO

      RETURN



      END SUBROUTINE TROPIN
