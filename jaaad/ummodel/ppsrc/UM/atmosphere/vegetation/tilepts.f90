
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Counts the number of points containing each surface type and creates
! a TILE_INDEX array specifying the location of these points on the land
! grid.
!
! Subroutine Interface:

      SUBROUTINE TILEPTS(LAND_PTS,FRAC,TILE_PTS,TILE_INDEX)

! module for land-surface namelist
      USE LAND_SURF_MOD, ONLY :                                         &
     & ALL_TILES

      IMPLICIT NONE
!
! Description:
!
! Method:
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    16/10/97   Original code.  Peter Cox
!   5.2    15/11/00   Re-Written for New Dynamics      M. Best
!
      INTEGER, INTENT(IN) :: LAND_PTS  !num land points to process

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



      REAL, INTENT(IN) :: FRAC(LAND_PTS,NTYPE)  !fractions of surface types
      
      INTEGER, INTENT(OUT) ::  TILE_PTS(NTYPE)                      &
                                                ! Number of land points which
                                                ! include the nth surface type
     &,                        TILE_INDEX(LAND_PTS,NTYPE)
                                                ! Indices of land points which
                                                ! include the nth surface type
      
      INTEGER :: N,L,C   !local counters: type, land pts, count
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
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

!-----------------------------------------------------------------------
! Create the TILE_INDEX array of land points with each surface type
!-----------------------------------------------------------------------
      TILE_PTS (:) = 0
      TILE_INDEX(:,:) = 0
      
      DO N=1,NTYPE
        C=0
!CDIR NODEP        
        DO L=1,LAND_PTS
          IF ( ((ALL_TILES == 0).AND.                                  &
     &          (FRAC(L,N) >  0.0))                                    &
     &        .OR.                                                     &
     &         ((ALL_TILES == 1).AND.                                  &
     &          ( (N .LT. NTYPE .AND. FRAC(L,NTYPE) .LT. 0.5)          &
     &           .OR.                                                  &
     &            (N .EQ. NTYPE .AND. FRAC(L,NTYPE) .GT. 0.5))) )THEN
            C = C + 1
            TILE_INDEX(C,N) = L
          ENDIF
        ENDDO
        TILE_PTS(N) = C
      ENDDO

      RETURN
      END SUBROUTINE TILEPTS
