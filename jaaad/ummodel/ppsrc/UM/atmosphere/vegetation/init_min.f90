
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Ensures that PFT fractions are greater than non-zero minimum fraction
!
! Subroutine Interface:
      SUBROUTINE INIT_MIN(LAND_PTS,FRAC,CS)

      IMPLICIT NONE
!
! Description:
!   If fractions of any PFTs are less than a non-zero minimum fraction
!   on land points that are not entirely (or mostly) covered by ice,
!   water or urban, initialise the PFT fractions to the minimum fraction
!   and take the excess proportionally from other PFTs and soil to
!   ensure that the total fractional cover of all PFTs + soil remains
!   unchanged.
!
! Method:
!   For PFTs with fraction < minimum fraction, reset fraction to minimum
!   fraction and find the total increase for all PFTs.  For all PFTS,
!   define "available fraction" as the difference between fraction
!   and minimum fraction, and find "available fraction" from sum
!   of all PFT available fractions plus fraction of soil (this is the
!   "available fraction" for soil; the minimum fraction for soil is
!   zero).  Reduce fractions of all PFTs and soil by amounts weighted
!   by "available fraction" / "total available fraction" such that the
!   sum of the reductions equals the total increase made earlier.  On
!   points with insufficent veg or soil to do this, take no action as
!   vegetation will not be modelled on these points.
!
!
! Current Code Owner: Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.5    8/5/98    Original code.  Richard Betts
!   5.2   15/11/00   Re-Written for New Dynamics    M. Best
!  6.2  01/03/06  initialise all 4 soil carbon pools.       C.D. Jones
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------ nstypes.h ----------------------------------
!jhan:further renovation of ths file may be necessary params are dependent on dataset
!jhan: ALSO nstypes_cable.h should be unecessary nsoil/soil is only difference
      !--- Number of non-vegetation surface types
      Integer, Parameter :: NNVG  = 4

      !--- Number of plant functional types.
      Integer, Parameter :: NPFT  = 13
      
      !--- Number of surface types.
      Integer, Parameter :: NTYPE =17 
      
      !--- Index of the surface type 'Soil'
      !Integer, Parameter :: SOIL  = 16 
      !dhb599, 20110615: change made as per Peter Vohralik, item 1:
      Integer, Parameter :: SOIL  = 14

!--- Land surface types :
!--- original veg. tiles 
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!--- for testing these tiles are set = 1:5 
!     6 - Broadleaf Tree
!     7 - Needleleaf Tree
!     8 - C3 Grass
!     9 - C4 Grass
!    10 - Shrub
!--- for testing these tiles are set = 0
!    11 - 0 
!    11 - 0
!    11 - 0
!--- original non-veg tiles moved to these indices
!     14 - Urban
!     15 - Water
!     16 - Soil
!     17 - Ice



! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER                                                           &
     & LAND_PTS            ! IN Number of land point to be processed.

      REAL                                                              &
     & CS(LAND_PTS,4)                                                   &
                           ! INOUT Soil carbon content (kg C/m2).
     &,FRAC(LAND_PTS,NTYPE) ! INOUT Fractions of surface types.

      INTEGER                                                           &
     & L                                                                &
                           ! Loop counter for land points
     &,N                   ! Loop counter for surface types

      REAL                                                              &
     & FRAC_AVAIL(LAND_PTS,NTYPE)                                       &
                                   ! LOCAL The part of FRAC that is
!                                  !       available for "donation"
     &,TOT_FRAC_NEED(LAND_PTS)                                          &
                                   ! LOCAL Total fraction needed to make
!                                  !       PFT fractions up to minimum
     &,TOT_FRAC_AVAIL(LAND_PTS)    ! LOCAL Total fractional area
!                                  !       available to give to PFTs
!                                  !       with less than minimum frac.

!----------------------------------------------------------------------
! Local parameters
!----------------------------------------------------------------------
! CSMIN start

      ! Minimum soil carbon (kg C/m2).
      REAL,PARAMETER:: CS_MIN = 1.0E-6

! CSMIN end
! Start seed
! Description:
!   This file sets the values of the variables FRAC_MIN and FRAC_SEED
!
! Current Code Owner:
!
! History:
! Version  Date     Comment
! -------  ----     -------
!   5.3   25/09/01  Portability changes.  Z. Gardner
!   5.5   17/04/03  Remove reference to obsolete section
!                   A03_7A. T.White
!
      ! Minimum areal fraction for PFTs.
      REAL, PARAMETER:: FRAC_MIN  = 1.0E-6

      ! "Seed" fraction for PFTs.
      REAL, PARAMETER:: FRAC_SEED = 0.01
! End Seed

      DO L=1,LAND_PTS
        TOT_FRAC_NEED(L) = 0.0
        TOT_FRAC_AVAIL(L) = 0.0

!-----------------------------------------------------------------------
! Find total fraction available for donation to PFTs with less than
! the minimum coverage
!-----------------------------------------------------------------------
        DO N=1,NPFT
          IF (FRAC(L,N) <  FRAC_MIN) THEN
              TOT_FRAC_NEED(L) = TOT_FRAC_NEED(L) +                     &
     &        (FRAC_MIN - FRAC(L,N))
          ELSE IF (FRAC(L,N) >= FRAC_MIN) THEN
            FRAC_AVAIL(L,N) = FRAC(L,N)-FRAC_MIN
            TOT_FRAC_AVAIL(L) = TOT_FRAC_AVAIL(L) + FRAC_AVAIL(L,N)
          ENDIF
        ENDDO
        N=SOIL
        FRAC_AVAIL(L,N) = FRAC(L,N)
        TOT_FRAC_AVAIL(L) = TOT_FRAC_AVAIL(L) + FRAC(L,N)

!-----------------------------------------------------------------------
! If sufficient total fraction is available, modify fractions of veg and
! soil and also modify soil carbon.  If insufficient fraction available,
! do neither of these as TRIFFID will not operate on such points.
!-----------------------------------------------------------------------
        IF (TOT_FRAC_AVAIL(L) >= TOT_FRAC_NEED(L)) THEN

!-----------------------------------------------------------------------
! i)  If PFT fraction is less than the minimum fraction, increase it
!     to the minimum fraction.
!-----------------------------------------------------------------------
          DO N=1,NPFT
            IF (FRAC(L,N) <  FRAC_MIN) THEN
              FRAC(L,N) = FRAC_MIN
              FRAC_AVAIL(L,N) = 0.0
            ELSEIF (FRAC(L,N) == FRAC_MIN) THEN
              FRAC_AVAIL(L,N) = 0.0
            ENDIF
          ENDDO

!-----------------------------------------------------------------------
! ii) Scale other PFTs and soil to keep total coverage of veg+soil
!     unchanged.  The relative proportions of the soil fraction and the
!     PFT fractions greater than the minimum fraction remain constant.
!-----------------------------------------------------------------------
          DO N=1,NPFT
            FRAC(L,N) = FRAC(L,N) -                                     &
     &      ( (FRAC_AVAIL(L,N)/TOT_FRAC_AVAIL(L)) * TOT_FRAC_NEED(L) )
          ENDDO

          N=SOIL
          FRAC(L,N) = FRAC(L,N) -                                       &
     &    ( (FRAC_AVAIL(L,N)/TOT_FRAC_AVAIL(L)) * TOT_FRAC_NEED(L) )

!-----------------------------------------------------------------------
! iii) If soil carbon content is less than minimum allowed, increase
!      it to the minimum.
!-----------------------------------------------------------------------
          IF ((CS(L,1)+CS(L,2)+CS(L,3)+CS(L,4)) <  CS_MIN) THEN
            CS(L,1) = CS_MIN*0.1
            CS(L,2) = CS_MIN*0.2
            CS(L,3) = CS_MIN*0.3
            CS(L,4) = CS_MIN*0.4
          ENDIF

        ENDIF

      ENDDO

      RETURN
      END SUBROUTINE INIT_MIN
