
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculates level dependent constants from ETA half levels.
!
! Subroutine Interface:
      SUBROUTINE ABCALC (MODE_L,MODE_H,MODE_C,LEVELS,ETA_P,ETA_S,ETAH   &
     &  ,AK,BK,AKH,BKH,IERR)

      IMPLICIT NONE
!
! Description:
!   Calculates a set of A and B values to define a
!   set of hybrid model levels from ETA half levels.
!   Original code written for user interface by
!   Richard Swinbank 15/06/90.
!
! Method:
!   <Say how it does it: include references to external documentation>
!   <If this routine is very complex, then include a "pseudo code"
!    description of it to make its structure and method clear>
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.5       24/02/95 Original code. (D.M. Goddard)
! 6.0       05/09/03 Added extra def for use with makebc. R.Sempers
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v7 programming standards.
!
! System component covered: None
! System Task:              None
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
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

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER      MODE_L               ! Method of calculation of
                                        ! eta level (ETAK)
                                        ! from layers (ETAKH)
                                        ! 1 - LOGARTHMIC  MEAN
                                        ! 2 - OLD MET O 20 METHOD
                                        ! 3 - ARITHMETIC MEAN
                                        ! 4 - arithmetic mean of half
                                        !     level exner
                                        ! 5 - pressure weighted mean of
                                        !     half level exner

      INTEGER      MODE_H               ! Method of calculating AKH and
                                        ! BKH from ETAKH, ETA_S and
                                        ! ETA_P
                                        ! 1 - Cubic variation of A, B
                                        !      with ETA

      INTEGER      MODE_C               ! Method of calculating AK and
                                        ! BK from AKH, BKH AND ETAK
                                        ! 1 - Linear interpolation
                                        ! 2 - Half-height method

      INTEGER      LEVELS               ! Number of levels

      REAL         ETA_S                ! Eta values at which levels
                                        !become sigma surf

      REAL         ETA_P                ! Eta values at which levels
                                        ! become P surfaces

!   Array  arguments with intent(in):
      REAL         ETAH(LEVELS+1)       ! Eta values at model layer
                                        ! boundaries ETAKH

!   Scalar arguments with intent(out):
      INTEGER      IERR                 ! Error code (>0 if error)

!   Array  arguments with intent(out):
      REAL         AK(LEVELS)           ! Value to define hybrid levels
      REAL         BK(LEVELS)           ! Value to define hybrid levels
      REAL         AKH(LEVELS+1)        ! Value to defn layer boundaries
      REAL         BKH(LEVELS+1)        ! Value to defn layer boundaries

! Local parameters:
      INTEGER      ILEVP                ! Maximum number of levels
      INTEGER      ILEVP1               ! Maximum number of half levels
      REAL         TINY                 ! Smallest allowed A value
      PARAMETER (ILEVP=75, ILEVP1=76)
      PARAMETER (TINY=1.0E-8)

! Local scalars:
      INTEGER      JL                   ! Loop index
      REAL         Z                    !!
      REAL         Z1                   !! Working parameters used to
      REAL         Z2                   !!calculate ETAL

! Local dynamic arrays:
      REAL        AH(ILEVP1)            ! Hybrid coordinate parameter
      REAL        AL(ILEVP)             ! A values on model levels
      REAL        ETAD(ILEVP)           ! Thickness of model layers
      REAL        ETAL(ILEVP)           ! Eta values on model levels

!- End of header


!-------------------------------------------------------------------
! 1 initialisation
!-------------------------------------------------------------------


      IF(LEVELS >  ILEVP) THEN
        WRITE(6,'('' NO. LEVELS ('',I3,'') > MAXIMUM ('',I3,'')'')')    &
     &    LEVELS,ILEVP
        IERR=1
        RETURN
      END IF
      IF(ETA_S <  ETA_P) THEN
        WRITE(6,'('' ETA_S < ETA_P '',2F12.8)')  ETA_S,ETA_P
        IERR=2
        RETURN
      ELSE
        IERR=0
      END IF

!-------------------------------------------------------------------
! 2 Calculate derived values
!-------------------------------------------------------------------

! Calculate layer thickness
      DO JL=1,LEVELS
        ETAD(JL)=ETAH(JL)-ETAH(JL+1)
      END DO

! Calculate (Layer centre) value of ETA
      IF(MODE_L == 1) THEN
        DO JL=1,LEVELS
          IF(ETAH(JL+1) <= 0.0) THEN
            ETAL(JL)=ETAH(JL)*EXP(-1.0)
          ELSE
            ETAL(JL)=SQRT(ETAH(JL)*ETAH(JL+1))
          END IF
        END DO
      ELSE IF(MODE_L == 2) THEN
        DO JL=1,LEVELS
          IF(ETAH(JL+1) <= 0.0) THEN
            ETAL(JL)=ETAH(JL)*EXP(-1.0)
          ELSE
            ETAL(JL)=(ETAH(JL)-ETAH(JL+1))/LOG(ETAH(JL)/ETAH(JL+1))
          END IF
        END DO
      ELSE IF(MODE_L == 3) THEN
        DO JL=1,LEVELS
          ETAL(JL) = 0.5 * (ETAH(JL) + ETAH(JL+1))
        END DO
      ELSE IF(MODE_L == 4) THEN
! Method 4 - unified model - pre vn2.6
        DO  JL=1,LEVELS
          Z1=ETAH(JL)**KAPPA
          Z2=ETAH(JL+1)**KAPPA
          Z=0.5*(Z1+Z2)
          ETAL(JL) = Z**(1.0/KAPPA)
        END DO
      ELSE IF(MODE_L == 5) THEN
! Method 5 - unified model -  vn2.6 onwards
        DO JL=1,LEVELS
          Z1=ETAH(JL)**(KAPPA + 1.0)
          Z2=ETAH(JL+1)**(KAPPA + 1.0)
          Z=(Z2-Z1)/(( ETAH(JL+1) - ETAH(JL) )*(KAPPA + 1.0) )
          ETAL(JL) = Z**(1.0/KAPPA)
        END DO
      END IF

! Calculate hybrid coordinate parameter A
      Z1=(ETA_S-ETA_P)*(ETA_S-ETA_P)*(ETA_S-ETA_P)
      DO JL=1,LEVELS+1
        IF(ETAH(JL) >= ETA_S) THEN
          AH(JL)=0.0
        ELSE IF(ETAH(JL) <= ETA_P) THEN
          AH(JL)=ETAH(JL)
        ELSE
          Z2=ETAH(JL)*(ETA_S+ETA_P)-2.0*ETA_P*ETA_P
          AH(JL)=(ETA_S-ETAH(JL))*(ETA_S-ETAH(JL))*Z2/Z1
          IF(ABS(AH(JL)) <  TINY) AH(JL)=0.0
        END IF
      END DO

! Calculate A values at LEVELS
      IF(MODE_C == 1) THEN
        DO JL=1,LEVELS
          AL(JL)=AH(JL)+(AH(JL+1)-AH(JL))*                              &
     &     (ETAL(JL)-ETAH(JL))/(ETAH(JL+1)-ETAH(JL))
          IF(ABS(AL(JL)) <  TINY) AL(JL)=0.0
        END DO
      ELSE IF(MODE_C == 2) THEN
        DO JL=1,LEVELS
          Z1=(ETAH(JL  )-AH(JL  )) * ETAH(JL  )**(KAPPA-1.0)
          Z2=(ETAH(JL+1)-AH(JL+1)) * ETAH(JL+1)**(KAPPA-1.0)
          Z=(0.5*(Z1+Z2)) * ETAL(JL)**(-(KAPPA-1.0))
          AL(JL)=ETAL(JL)-Z
          IF(ABS(AL(JL)) <  TINY) AL(JL)=0.0
        END DO
      END IF

! Set A and B Arrays
      DO JL=1,LEVELS
        AK(JL)=PREF*AL(JL)
        BK(JL)=ETAL(JL)-AL(JL)
        IF(ABS(BK(JL)) <  TINY) BK(JL)=0.0
      END DO
      DO JL=1,LEVELS+1
        AKH(JL)=PREF*AH(JL)
        BKH(JL)=ETAH(JL)-AH(JL)
        IF(ABS(BKH(JL)) <  TINY) BKH(JL)=0.0
      END DO

      RETURN
      END SUBROUTINE ABCALC
