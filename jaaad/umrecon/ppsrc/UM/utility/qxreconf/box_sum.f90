
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE BOX_SUM
!LL
!LL  NOT SUITABLE FOR SINGLE COLUMN USE
!LL
!LL
!LL  SUITABLE FOR ROTATED GRIDS
!LL
!LL  ORIGINAL VERSION FOR CRAY Y-MP/IBM
!LL  WRITTEN 12/07/91 BY C. WILSON
!LL
!LL  CODE REVIEWED BY R.SMITH ??/??/??
!LL
!LL  VERSION NO. 2 DATED 16/09/91
!LL         COSMOS DSN MS15.CWUM.JOBS(BOXSUM2)
!LL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
!LL  VERSION 1, DATED 12/09/89
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.0      12/04/95  Imported into Unified model. D.M. Goddard
! 4.1      12/06/96  Corrections for zonal means. D.M. Goddard
! 4.4      30/09/97  Corrections for portable model. D.M. Goddard
!LL
!LL  SYSTEM TASK:  S1 (part,extension for area mean interpolation)
!LL
!LL  PURPOSE:
!LL  Routine sums contributions from gridboxes for source data on a
!LL  regular lat-long grid to form means for gridboxes of a regular
!LL  lat-long grid specified as target.
!LL  Both grids are defined with the same pole and orientation;
!LL  the original data must be interpolated onto a rotated
!LL  grid ,if the target grid is a rotated grid, BEFORE calling this
!LL  routine.
!LL  The algorithms are general and will cope with either a finer
!LL  or coarser resolution source grid.
!LL
!LL  DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1
!LL                  BY A.DICKINSON/C WILSON VERSION ??DATED ??/??/91
!LL
!LL
!LLEND-------------------------------------------------------------
!
!*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE BOX_SUM                                                &
     &  (SOURCE_ROW_LENGTH,SOURCE_ROWS,ROW_LENGTH,ROWS,                 &
     &   LONG_L,COLAT_T,I_L,J_T,GLOBAL,BOXSUM,SOURCE)

      IMPLICIT NONE

      INTEGER                                                           &
     &  SOURCE_ROW_LENGTH                                               &
                           !IN    Number of points per row (source data)
                           !      on rotated grid if necessary
     &, SOURCE_ROWS                                                     &
                           !IN    Number of rows of source data
                           !      on rotated grid if necessary
     &, ROW_LENGTH                                                      &
                           !IN    Number of points per row target area
     &, ROWS               !IN    Number of rows of target area
!
      INTEGER                                                           &
     & I_L(ROW_LENGTH+1)                                                &
                        !IN Index of first source gridbox to overlap
                        !   with left hand side of target gridbox
     &,J_T(ROWS+1)      !IN Index of first source gridbox to overlap
                        !   top of target gridbox
!
!L N.B.I_L(I) is the first source gridbox to overlap LH side of target
!L            box  I of a row
!L     I_L(I+1) is the last source gridbox to overlap RH side of target
!L            box  I of a row
!L     J_T(J) is the first source gridbox to overlap top of target
!L            box on row J
!L     J_T(J+1) is the last source gridbox to overlap bottom of target
!L            box on row J
!L
!L REAL value of:-
!L     I_L(I) is also used to measure the 'longitude' of the RHS of the
!L            source gridbox
!L     J_T(J) is also used to measure the 'colatitude' of the bottom
!L            of the source gridbox
!

      REAL                                                              &
     & SOURCE(SOURCE_ROW_LENGTH,SOURCE_ROWS)!IN  source data
      REAL BOXSUM(ROW_LENGTH,ROWS)                                      &
                                   !OUT Sum of data on target grid
     &,    LONG_L(ROW_LENGTH +1)                                        &
                                    !IN Left longitude of gridbox (in
                                    ! units of souce gridbox EW length)
     &,    COLAT_T(ROWS +1)         !IN Colatitude of top of gridbox (in
                                    ! units of source gridbox NS length)

      LOGICAL GLOBAL       !IN    true if global area required
!*---------------------------------------------------------------------

!*L  WORKSPACE USAGE:-------------------------------------------------
!   Define local workspace arrays:
!    1 real of length row_length
      REAL                                                              &
     &     EW_SUM(ROW_LENGTH)                                           &
                                  ! summed WE source data
     &,    EW_WEIGHT(ROW_LENGTH)                                        &
                                  ! summed WE weights for source data
     &,    BOX_WEIGHT(ROW_LENGTH)  ! summed weights for target boxes
!
!*---------------------------------------------------------------------
!
!*L EXTERNAL SUBROUTINES CALLED---------------------------------------
! None
!*------------------------------------------------------------------

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!    DEFINE LOCAL VARIABLES

      INTEGER I,J,I1,I2,IT,J1,J2,JT,K ! loop counters

      REAL RH_BOX

!L *********************************************************************
!L 1.0 Sum source boxes (whole and partial) contributions to target box
!L *********************************************************************

      DO 150 J=1,ROWS
        J1 = J_T(J)
        J2 = J_T(J+1)

        DO I=1,ROW_LENGTH
          BOX_WEIGHT(I)=0.0
        ENDDO

        DO 140 JT=J1,J2

!L *********************************************************************
!L 1.1 Sum  EW (whole and partial) contributions to target grid boxes
!L *********************************************************************

          DO I=1,ROW_LENGTH
            EW_SUM(I)=0.0
            EW_WEIGHT(I)=0.0
          ENDDO

          DO 120 I=1,ROW_LENGTH
            I1 = I_L(I)
            I2 = I_L(I+1)
            IF(I1 >  I2.AND.GLOBAL) THEN
!   If grid box spans zero longitude need to split summation

              DO 101 IT=I1,SOURCE_ROW_LENGTH
                IF(IT == I1) THEN
!   Left side partial contribution
                  RH_BOX = LONG_L(I+1)
                  IF(RH_BOX <  LONG_L(I)) RH_BOX=RH_BOX                 &
     &                                        +SOURCE_ROW_LENGTH
!XX               IF(NINT(SOURCE(I1,JT)) /= IMDI) THEN
                  IF(NINT(SOURCE(I1,JT)) /= NINT(RMDI)) THEN
!XX               IF(SOURCE(I1,JT) >= 0.0) THEN
                    EW_WEIGHT(I) =                                      &
     &                EW_WEIGHT(I) + (MIN(REAL(I1),RH_BOX) - LONG_L(I))
                    EW_SUM(I) = (MIN(REAL(I1),RH_BOX) - LONG_L(I))      &
     &                      *SOURCE(I1,JT) + EW_SUM(I)
                  ELSEIF(NINT(SOURCE(I1,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 1 ? I1 JT data',I1,JT,      &
     &              source(i1,jt)
                  ENDIF

                ELSE

!   Whole contributions
!XX               IF(NINT(SOURCE(IT,JT)) /= IMDI) THEN
                  IF(NINT(SOURCE(IT,JT)) /= NINT(RMDI)) THEN
!XX               IF(SOURCE(IT,JT) >= 0.0) THEN
                    EW_WEIGHT(I) = EW_WEIGHT(I)+ 1.0
                    EW_SUM(I) = EW_SUM(I) + SOURCE(IT,JT)
                  ELSEIF(NINT(SOURCE(IT,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 2 ? IT JT data',IT,JT,      &
     &              source(it,jt)
                  ENDIF

                ENDIF

101           CONTINUE

              DO 102 IT=1,I2
                IF(IT == I2) THEN
!   Right side partial contribution
!XX               IF(NINT(SOURCE(I2,JT)) /= IMDI) THEN
                  IF(NINT(SOURCE(I2,JT)) /= NINT(RMDI)) THEN
!XX               IF(SOURCE(I2,JT) >= 0.0) THEN
                    EW_WEIGHT(I) = EW_WEIGHT(I)+(LONG_L(I+1)-(I2-1))
                    EW_SUM(I)=EW_SUM(I) +                               &
     &                               (LONG_L(I+1)-(I2-1))*SOURCE(I2,JT)
                  ELSEIF(NINT(SOURCE(I2,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 3 ? I2 JT data',I2,JT,      &
     &              source(i2,jt)
                  ENDIF
                ELSE

!   Whole contributions
!XX               IF(NINT(SOURCE(IT,JT)) /= IMDI) THEN
                  IF(NINT(SOURCE(IT,JT)) /= NINT(RMDI)) THEN
!XX               IF(SOURCE(IT,JT) >= 0.0) THEN
                    EW_WEIGHT(I) = EW_WEIGHT(I)+ 1.0
                    EW_SUM(I) = EW_SUM(I) + SOURCE(IT,JT)
                  ELSEIF(NINT(SOURCE(IT,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 4 ? IT JT data',IT,JT,      &
     &              source(it,jt)
                  ENDIF

                ENDIF

102           CONTINUE

          ELSE IF(I1 <  I2)THEN ! no zero meridian crossing
            DO 110 IT=I1,I2
              IF(IT == I1) THEN
!   Left side partial contribution
!XX             IF(NINT(SOURCE(I1,JT)) /= IMDI) THEN
                IF(NINT(SOURCE(I1,JT)) /= NINT(RMDI)) THEN
!XX             IF(SOURCE(I1,JT) >= 0.0) THEN
                  EW_WEIGHT(I) = EW_WEIGHT(I) +                         &
     &                      (MIN(REAL(I1),LONG_L(I+1)) - LONG_L(I))
                  EW_SUM(I) = (MIN(REAL(I1),LONG_L(I+1)) - LONG_L(I))   &
     &                      *SOURCE(I1,JT) + EW_SUM(I)
                  ELSEIF(NINT(SOURCE(I1,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 5 ? I1 JT data',I1,JT,      &
     &              source(i1,jt)
                ENDIF

              ELSE IF(IT == I2) THEN
!   Right side partial contribution
!XX             IF(NINT(SOURCE(I2,JT)) /= IMDI) THEN
                IF(NINT(SOURCE(I2,JT)) /= NINT(RMDI)) THEN
!XX             IF(SOURCE(I2,JT) >= 0.0) THEN
                  EW_WEIGHT(I) = EW_WEIGHT(I)+ (LONG_L(I+1)-(I2-1))
                  EW_SUM(I)=EW_SUM(I)+(LONG_L(I+1)-(I2-1))*SOURCE(I2,JT)
                  ELSEIF(NINT(SOURCE(I2,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 6 ? I2 JT data',I2,JT,      &
     &              source(i2,jt)
                ENDIF

              ELSE

!   Whole contributions

!XX             IF(NINT(SOURCE(IT,JT)) /= IMDI) THEN
                IF(NINT(SOURCE(IT,JT)) /= NINT(RMDI)) THEN
!XX             IF(SOURCE(IT,JT) >= 0.0) THEN
                  EW_WEIGHT(I) = EW_WEIGHT(I)+ 1.0
                  EW_SUM(I) = EW_SUM(I) + SOURCE(IT,JT)
                  ELSEIF(NINT(SOURCE(IT,JT)) == NINT(RMDI)) THEN
!                   missing data
                  ELSE
                    write (6,*) '-ve source 7 ? IT JT data',IT,JT,      &
     &              source(it,jt)
                ENDIF

              ENDIF
110           CONTINUE

          ELSE
!Zonal mean no need to average in EW direction
            DO K=1,ROW_LENGTH
              EW_WEIGHT(K)=1.0
              EW_SUM(K)=SOURCE(1,JT)
            END DO



          ENDIF

120       CONTINUE

!L *********************************************************************
!L 1.3 Add summed EW  box contributions into rows J  target grid boxes
!L *********************************************************************

          IF(JT == J1) THEN
!   Top row
            DO 130 I=1,ROW_LENGTH
              BOXSUM(I,J) = (MIN(REAL(J1),COLAT_T(J+1)) - COLAT_T(J))   &
     &                       *EW_SUM(I)
              BOX_WEIGHT(I) = (MIN(REAL(J1),COLAT_T(J+1)) - COLAT_T(J)) &
     &                       *EW_WEIGHT(I) + BOX_WEIGHT(I)
130         CONTINUE

          ELSE IF(JT == J2) THEN
!   Bottom of row J
            DO 131 I=1,ROW_LENGTH
              BOXSUM(I,J) = BOXSUM(I,J) +                               &
     &         (1-(J2 - COLAT_T(J+1)))*EW_SUM(I)
              BOX_WEIGHT(I) = BOX_WEIGHT(I) +                           &
     &         (1-(J2 - COLAT_T(J+1)))*EW_WEIGHT(I)
131         CONTINUE

          ELSE
!   Whole contributions to row J
            DO 132 I=1,ROW_LENGTH
              BOXSUM(I,J) = BOXSUM(I,J) + EW_SUM(I)
              BOX_WEIGHT(I) = BOX_WEIGHT(I) + EW_WEIGHT(I)
132         CONTINUE

          ENDIF

140     CONTINUE

          DO 142 I=1,ROW_LENGTH
            IF(BOX_WEIGHT(I) /= 0.0) THEN
              BOXSUM(I,J) = BOXSUM(I,J) / BOX_WEIGHT(I)
            ELSE
              BOXSUM(I,J) = RMDI
            ENDIF
142       CONTINUE

150   CONTINUE

      RETURN
      END SUBROUTINE BOX_SUM
