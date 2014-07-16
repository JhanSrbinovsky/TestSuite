

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE PRE_AREAVER(GAPS_LAMBDA_SRCE,LAMBDA_SRCE               &
     &,GAPS_PHI_SRCE,PHI_SRCE,CYCLIC_SRCE,LROW_SRCE,WANT,MASK_SRCE      &
     &,GAPS_LAMBDA_TARG,LAMBDA_TARG,GAPS_PHI_TARG,PHI_TARG              &
     &,CYCLIC_TARG,SPHERICAL                                            &
     &,MAXL,COUNT_TARG,BASE_TARG,INDEX_SRCE,WEIGHT,WEIGHT2,ICODE,CMESSAGE)
!LL
!LL   Subroutine PRE_AREAVER
!LL
!LL   Calculate weights for area-averaging data on the source grid to
!LL   data on the target grid.
!LL
!LL  J.Gregory  <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
!LL                   portability.  Author Tracey Smith.
!LL   3.3   3.9.93    Correct out-of-bounds errors for IXL and IYT
!LL   4.1   22.5.96   J.M.Gregory  Add argument SPHERICAL to produce
!LL                   correct weights for a spherical surface.
!LL   5.3   22.10.01  Change for south -> north ordering. C.F.Durman.
!  5.5  28.02.03  Add defined RIVERS. C. Bunton
!  6.0  12.09.03  Change DEF from A20 to A26. D. Robinson.
!
!      ****ANY CHANGES MADE TO THIS ROUTINE SHOULD BE MIRRORED IN
!          SUBROUTINE PRE_AREAVG_ICE******
!
!LL
!LL   Standard, paper 4 version 4 (14.12.90)
!LL
!L
!L    The source grid and target grid are each specified by a lambda set
!L    and a phi set of coordinates delimiting the boxes. These sets are
!L    supplied in 1D arrays aa_bb for coordinate aa=LAMBDA,PHI on grid
!L    bb=SRCE,TARG.  The number of gaps is specified by GAPS_aa_bb,
!L    which is equal to the number of lines IF (CYCLIC_bb) and aa is
!L    LAMBDA, otherwise one less. (By "gap" we mean the interval
!L    spanning a box in one coordinate only. The total number of boxes,
!L    or grid points, is the product of the numbers of gaps in the two
!L    coordinates.) Whether the axes are cyclic is not known until
!L    run-time, so the dimensions of the arrays LAMBDA_bb are not known
!L    at compile-time, and they are dimensioned assumed-size. There are
!L    no restrictions on the base meridian of lambda, and it does not
!L    have to be the same for source and target grids. The lambda
!L    coordinates should be in increasing order (running from left to
!L    right), the phi increasing (top to bottom). The coordinates must
!L    be given in degrees for cyclic axes, because a range of 360 is
!L    assumed, or IF (SPHERICAL), when trigonometric functions are
!L    used. IF (SPHERICAL), the weights computed are for a spherical
!L    surface, assuming that LAMBDA is longitude and PHI latitude.
!L    Otherwise, LAMBDA and PHI are treated as Cartesian axes on a
!L    plane.
!L
!L    The array MASK_SRCE is the land/sea mask for the source grid. The
!L    logical value indicating points to be used should be supplied in
!L    WANT. The first dimension of MASK_SRCE should be supplied in
!L    LROW_SRCE. Specifying this separately allows for the possibility
!L    of extra rows and columns in MASK_SRCE which are to be ignored.
!L
!L    The arrays COUNT_TARG and BASE_TARG should be dimensioned in the
!L    calling program to the number of boxes in the target array.
!L
!L    The arrays INDEX_SRCE and WEIGHT are returned in the form which
!L    the area-averaging routine DO_AREAVER expects. They are continuous
!L    lists comprising consecutive groups of entries. There is a group
!L    for each target point, for which the number of entries is spec-
!L    ified by COUNT_TARG, and the groups appear in the normal order of
!L    grid points. The size required for INDEX_SRCE and WEIGHT depends
!L    on how many source boxes go into each target box, on average, and
!L    is not known at compile-time. The maximum that could be needed is
!L    (GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)*(GAPS_PHI_SRCE+GAPS_PHI_TARG)
!L    and the size to which the arrays are actually dimensioned should
!L    be supplied in MAXL. The size used is returned in MAXL. It is the
!L    responsibility of the calling routine to provide enough space.
!L
      IMPLICIT NONE
!*L
      INTEGER                                                           &
     & GAPS_LAMBDA_SRCE                                                 &
                               !IN number of lambda gaps in source grid
     &,GAPS_PHI_SRCE                                                    &
                               !IN number of phi gaps in source grid
     &,GAPS_LAMBDA_TARG                                                 &
                               !IN number of lambda gaps in target grid
     &,GAPS_PHI_TARG                                                    &
                               !IN number of phi gaps in target grid
     &,LROW_SRCE                                                        &
                               !IN first dimension of MASK_SRCE
     &,MAXL                                                             &
                               !INOUT maximum entries in output lists
     &,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                       &
!                              !OUT no. of source boxes in target box
     &,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)                        &
!                              !OUT first index in list for target box
     &,INDEX_SRCE(MAXL)                                                 &
                               !OUT list of source box indices
     &,ICODE                   !OUT return code
      LOGICAL                                                           &
     & CYCLIC_SRCE                                                      &
                               !IN source grid is cyclic
     &,CYCLIC_TARG                                                      &
                               !IN target grid is cyclic
     &,WANT                                                             &
                               !IN indicator of wanted points in mask
     &,MASK_SRCE(LROW_SRCE,*)                                           &
                               !IN land/sea mask for source grid
     &,SPHERICAL               !IN calculate weights for a sphere
      REAL                                                              &
     & LAMBDA_SRCE(*)                                                   &
                               !IN source lambda line coordinates
     &,PHI_SRCE(*)                                                      &
                               !IN source phi line coordinates
     &,LAMBDA_TARG(*)                                                   &
                               !IN target lambda line coordinates
     &,PHI_TARG(*)                                                      &
                               !IN target phi line coordinates
     &,WEIGHT(MAXL)                                                     &
                               !OUT list of weights for source boxes    
     &,WEIGHT2(MAXL)
                               !OUT weights o2a coastal flux regridding
      CHARACTER*(80) CMESSAGE  !OUT error message
!*
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
!
      INTEGER                                                           &
     & LINES_LAMBDA_SRCE                                                &
                               ! number of source lambda lines
     &,LINES_PHI_SRCE                                                   &
                               ! number of source phi lines
     &,LINES_LAMBDA_TARG                                                &
                               ! number of target lambda lines
     &,LINES_PHI_TARG                                                   &
                               ! number of target phi lines
     &,COUNT_LAMBDA(GAPS_LAMBDA_TARG)                                   &
!                              ! number of source lambda gaps per target
     &,COUNT_PHI(GAPS_PHI_TARG)                                         &
!                              ! number of source phi gaps per target
     &,INDEX_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)                  &
!                              ! source lambda gap indices
     &,INDEX_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)                           &
!                              ! source phi gap indices
     &,IX1,IY1,IX2,IY2                                                  &
                               ! working SRCE/TARG LAMBDA/PHI indices
     &,IX1N,IX1W                                                        &
                               ! working indices
     &,IXL(GAPS_LAMBDA_TARG+1)                                          &
                               ! source line on the left of target line
     &,IX2N                                                             &
                               ! target gap to the right of IX2
     &,IYT(GAPS_PHI_TARG+1)                                             &
                               ! source line above target line
     &,IXP,IYP,IP                                                       &
                               ! pointers into lists
     &,IX,IY,I                 ! loop indices
      REAL                                                              &
     & BASLAM                                                           &
                               ! minimum lambda for TEMP coordinates
     &,BTARG                                                            &
                               ! width of target gap
     &,BLO,BHI                                                          &
                               ! limits of gap
     &,TEMP_SRCE(GAPS_LAMBDA_SRCE+1)                                    &
!                              ! adjusted version of LAMBDA_SRCE
     &,TEMP_TARG(GAPS_LAMBDA_TARG+1)                                    &
!                              ! adjusted version of LAMBDA_TARG
     &,FRAC_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)                   &
!                              ! fractions of target lambda gaps
     &,FRAC_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)                            &
!                              ! fractions of target phi gaps
     &,SUM                                                              &
!                              ! sum of weights
     &,SUM2                    ! sum of weights in coastal zone
!
!L    1  Set up the lambda coordinates to make them easier to handle.
!
!L    1.1  Produce in TEMP_SRCE a monotonically increasing set of angles
!L    equivalent to LAMBDA_SRCE i.e. equal under modulo 360.
!
      IF (CYCLIC_SRCE) THEN
        LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE
      ELSE
        LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE+1
      ENDIF
      BASLAM=LAMBDA_SRCE(1)
      DO 10 IX1=1,LINES_LAMBDA_SRCE
        IF (LAMBDA_SRCE(IX1) <  BASLAM) THEN
          TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)+360.
        ELSE
          TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)
        ENDIF
   10 CONTINUE ! Labelled DO for the sake of fpp
!
!L    1.2  Produce in TEMP_TARG a set of angles equivalent to
!L    LAMBDA_TARG i.e. equal under modulo 360, but all in the range
!L    BASLAM to BASLAM+360, where BASLAM=min(TEMP_LAMBDA).
!
      IF (CYCLIC_TARG) THEN
        LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG
      ELSE
        LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG+1
      ENDIF
      DO 20 IX2=1,LINES_LAMBDA_TARG
        TEMP_TARG(IX2)=MOD(LAMBDA_TARG(IX2)-BASLAM,360.)
        IF (TEMP_TARG(IX2) <  0.) TEMP_TARG(IX2)=TEMP_TARG(IX2)+360.
        TEMP_TARG(IX2)=TEMP_TARG(IX2)+BASLAM
   20 CONTINUE ! Labelled DO for the sake of fpp
!
!L    2  For each target lambda line, find the index of the next source
!L    lambda line to the left.
!
      DO IX2=1,LINES_LAMBDA_TARG
        DO IX1=1,LINES_LAMBDA_SRCE
          IF (TEMP_TARG(IX2) >= TEMP_SRCE(IX1)) IXL(IX2)=IX1
        ENDDO
      ENDDO
!
!L    3  Find which source lambda gaps cover each target gap and the
!L    fractions they contribute.
!
!     At this point IXL(target_line) gives the index of the next source
!     lambda line to the left of the target lambda line, wrapping round
!     if the source grid is cyclic. This is also the index of the source
!     gap in which the target line falls. Similarly, the index of the
!     target line is also that of the target gap of which it is the
!     left-hand limit. Therefore also IXL(target_gap+1, wrapping round
!     if reqd.), is the index of the source gap which contains the
!     right-hand limit of the target gap. For each target gap, we loop
!     over all source gaps and find the fraction covered by each. Record
!     the fraction and the source index in cumulative lists. If the
!     source grid is not cyclic, parts of the target gap lying outside
!     the source grid are neglected.
!
      IXP=0
      DO IX2=1,GAPS_LAMBDA_TARG
        IX=0
        IX2N=MOD(IX2,LINES_LAMBDA_TARG)+1
        BTARG=TEMP_TARG(IX2N)-TEMP_TARG(IX2)
        IF (BTARG <  0.) THEN
          BTARG=BTARG+360.
          IX1N=IXL(IX2N)+LINES_LAMBDA_SRCE
        ELSE
          IX1N=IXL(IX2N)
        ENDIF
        DO IX1W=IXL(IX2),IX1N
          IX1=MOD(IX1W-1,LINES_LAMBDA_SRCE)+1
          IF (CYCLIC_SRCE.OR.IX1 /= LINES_LAMBDA_SRCE) THEN
            IF (IX1W == IXL(IX2)) THEN
              BLO=TEMP_TARG(IX2)
            ELSE
              BLO=TEMP_SRCE(IX1)
            ENDIF
            IF (IX1W == IX1N) THEN
              BHI=TEMP_TARG(IX2N)
            ELSE
              BHI=TEMP_SRCE(MOD(IX1,LINES_LAMBDA_SRCE)+1)
            ENDIF
            IF (BHI <  BLO) BHI=BHI+360.
            IF (ABS(BHI-BLO) >  1E-7*ABS(BTARG)) THEN
              IX=IX+1
              INDEX_LAMBDA(IXP+IX)=IX1
              FRAC_LAMBDA(IXP+IX)=(BHI-BLO)/BTARG
            ENDIF
          ENDIF
        ENDDO
        COUNT_LAMBDA(IX2)=IX
        IXP=IXP+COUNT_LAMBDA(IX2)
      ENDDO
!
!L    4  For each target phi line, find the index of the next source phi
!L    line above. Comments as for section 2, without wrap-round. Note
!L    that this assumes that the atmosphere gridbox ordering is from
!L    south to north; a north to south atmosphere would require LE
!L    instead of GE in the following IF test.
!
      LINES_PHI_SRCE=GAPS_PHI_SRCE+1
      LINES_PHI_TARG=GAPS_PHI_TARG+1
      DO IY2=1,LINES_PHI_TARG
        IYT(IY2)=0
        DO IY1=1,LINES_PHI_SRCE
          IF (PHI_TARG(IY2) >= PHI_SRCE(IY1)) IYT(IY2)=IY1
        ENDDO
      ENDDO
!
!L    5  Find which source phi gaps cover each target gap and the
!L    fractions they contribute. Comments as for section 3, without
!L    wrap-round.
!
      IYP=0
      DO IY2=1,GAPS_PHI_TARG
        IY=0
        IF (SPHERICAL) THEN
!     Contain angle between +-90. There is no real area outside
!     these limits on a sphere.
          BTARG=SIN(MAX(MIN(PHI_TARG(IY2+1),90.),-90.)*PI_OVER_180)     &
     &    -SIN(MAX(MIN(PHI_TARG(IY2),90.),-90.)*PI_OVER_180)
        ELSE
           BTARG=PHI_TARG(IY2+1)-PHI_TARG(IY2)
        ENDIF
        DO IY1=IYT(IY2),IYT(IY2+1)
          IF (.NOT.(IY1 == 0.OR.IY1 == LINES_PHI_SRCE)) THEN
            IF (IY1 == IYT(IY2)) THEN
              BLO=PHI_TARG(IY2)
            ELSE
              BLO=PHI_SRCE(IY1)
            ENDIF
            IF (IY1 == IYT(IY2+1)) THEN
              BHI=PHI_TARG(IY2+1)
            ELSE
              BHI=PHI_SRCE(IY1+1)
            ENDIF
            IF (SPHERICAL) THEN
              BLO=MAX(MIN(BLO,90.),-90.)
              BHI=MAX(MIN(BHI,90.),-90.)
            ENDIF
            IF (ABS(BHI-BLO) >  1E-7*ABS(BTARG)) THEN
              IY=IY+1
              INDEX_PHI(IYP+IY)=IY1
!             Both numerator and denominator in the following are -ve.
              IF (SPHERICAL) THEN
                FRAC_PHI(IYP+IY)                                        &
     &          =(SIN(BHI*PI_OVER_180)-SIN(BLO*PI_OVER_180))/BTARG
              ELSE
                FRAC_PHI(IYP+IY)=(BHI-BLO)/BTARG
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        COUNT_PHI(IY2)=IY
        IYP=IYP+COUNT_PHI(IY2)
      ENDDO
!
!L    6  For each target box, loop over contributing source boxes and
!L    calculate the weights for each one, ignoring land boxes.
!
!     After the first pass for each target box, go back and normalise
!     the weights to compensate for land source boxes and any outside
!     the source grid. Record the source box index and the weight in
!     cumulative lists.
!
      IP=0
      IYP=0
      DO IY2=1,GAPS_PHI_TARG
        IXP=0
        DO IX2=1,GAPS_LAMBDA_TARG
          I=0
          SUM=0
          SUM2=0
          DO IY=IYP+1,IYP+COUNT_PHI(IY2)
            DO IX=IXP+1,IXP+COUNT_LAMBDA(IX2)
              IF (MASK_SRCE(INDEX_LAMBDA(IX),INDEX_PHI(IY))             &
     &        .EQV.WANT) THEN
                I=I+1
                INDEX_SRCE(IP+I)=INDEX_LAMBDA(IX)                       &
     &          +(INDEX_PHI(IY)-1)*GAPS_LAMBDA_SRCE
                WEIGHT(IP+I)=FRAC_LAMBDA(IX)*FRAC_PHI(IY)
                SUM=SUM+WEIGHT(IP+I)
              ENDIF
              SUM2=SUM2+(FRAC_LAMBDA(IX)*FRAC_PHI(IY))
            ENDDO
          ENDDO
          COUNT_TARG(IX2,IY2)=I
          BASE_TARG(IX2,IY2)=IP
          DO I=1,COUNT_TARG(IX2,IY2)
            WEIGHT2(IP+I)=WEIGHT(IP+I)/SUM2
            WEIGHT(IP+I)=WEIGHT(IP+I)/SUM
          ENDDO
          IP=IP+COUNT_TARG(IX2,IY2)
          IXP=IXP+COUNT_LAMBDA(IX2)
        ENDDO
        IYP=IYP+COUNT_PHI(IY2)
      ENDDO
      MAXL=IP
!
      ICODE=0
      CMESSAGE=' '
      RETURN
      END SUBROUTINE PRE_AREAVER
!
!------------------------------------------------------------
!
