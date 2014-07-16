

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set fields of climatological aerosols in HADCM3.
!
! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.
!
! Method:
!   Straightforward.
!
! Current Owner of Code: James Manners
!
!- ---------------------------------------------------------------------
      SUBROUTINE R2_SET_AERO_CLIM_HADCM3(N_PROFILE, NLEVS, N_LAYER      &
     &   , I_GATHER, L_EXTRA_TOP                                        &
     &   , L_CLIM_AERO_HGT, BL_DEPTH, T, N_LEVELS_BL                    &
     &   , LAND, LYING_SNOW, PSTAR, P_LAYER_BOUNDARIES, TRINDX          &
     &   , L_VOLCTS, VOLCMASS							&
     &   , AEROSOL_MIX_RATIO_CLIM                                       &
     &   , NPD_FIELD, NPD_PROFILE, NPD_LAYER                            &
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     COMDECKS INCLUDED.
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
!
!     DUMMY ARGUMENTS.
!
!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
     &     NPD_FIELD                                                    &
!             FIELD SIZE IN CALLING PROGRAM
     &   , NPD_PROFILE                                                  &
!             SIZE OF ARRAY OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
     &     N_PROFILE                                                    &
!             NUMBER OF PROFILES
     &   , NLEVS                                                        &
!             Number of atmospheric layers used outside the radiation
!             scheme
     &   , N_LAYER
!             Number of atmospheric layers seen in radiation
!
!     Variables related to options for setting the field
      LOGICAL, INTENT(IN) ::                                            &
     &     L_CLIM_AERO_HGT                                              &
!             Flag to use the depth of the boundary layer to set
!             the layers in which a boundary-layer aerosol is used
     &   , L_EXTRA_TOP
!             Flag to use an extra top layer in radiative calculations
      INTEGER, INTENT(IN) ::                                            &
     &     N_LEVELS_BL
!             Common number of layers taken to be occupied by the
!             boundary-layer aerosol if this the boundary layer
!             depth is not used to determine the number separately
!             at each grid-point
      REAL, INTENT(IN) ::                                               &
     &     BL_DEPTH(NPD_FIELD)                                          &
!             Depth of the boundary layer
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             Temperatures of atmospheric layers
!
!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
     &     I_GATHER(NPD_FIELD)
!             LIST OF POINTS TO GATHER
!
!     GENERAL ATMOSPHERIC PROPERTIES:
      INTEGER                                                           &
                !, INTENT(IN)
     &     TRINDX(NPD_FIELD)
!             LAYER BOUNDARY OF TROPOPAUSE
      REAL                                                              &
                !, INTENT(IN)
     &     PSTAR(NPD_FIELD)                                             &
!             SURFACE PRESSURES
     &,    P_LAYER_BOUNDARIES(NPD_FIELD,0:NLEVS)			&
!             PRESSURE AT BOUNDARIES OF LAYERS
     &   , VOLCMASS(NPD_FIELD)
!             Mass of stratospheric volcanic aerosol at each point

      LOGICAL  L_VOLCTS

!
!     SURFACE FIELDS
      LOGICAL                                                           &
                !, INTENT(IN)
     &     LAND(NPD_FIELD)
!             LAND-SEA MASK
      REAL                                                              &
                !, INTENT(IN)
     &     LYING_SNOW(NPD_FIELD)
!             DEPTH OF LYING SNOW
!
      REAL                                                              &
                !, INTENT(OUT)
     &     AEROSOL_MIX_RATIO_CLIM(NPD_PROFILE, 0: NPD_LAYER, 5)
!             MIXING RATIOS OF CLIMATOLOGICAL AEROSOLS
!
!
!
!     LOCAL VARIABLES:
      INTEGER                                                           &
     &     I                                                            &
!             LOOP VARIABLE
     &   , I_UM                                                         &
!             Index of a level in the UM's upward convention
     &   , J                                                            &
!             LOOP VARIABLE
     &   , L                                                            &
!             LOOP VARIABLE
     &   , LG                                                           &
!             INDEX FOR GATHERING
     &   , BL_TOP_LYR(NPD_PROFILE)                                      &
!             Topmost layer occupied by boundary layer aerosol,
!             indexed using the top-down convention of radiation.
     &   , N_POINTS_SET_Z
!             Number of points where the top of the boundary layer
!             has not yet been reached
      LOGICAL                                                           &
     &     L_SET_Z(NPD_PROFILE)
!             Array to flag points where the top of the boundary layer
!             has not yet been reached
      REAL                                                              &
     &     PRESSURE_WT(NPD_FIELD)                                       &
!             ARRAY FOR SCALING AEROSOL AMOUNTS FOR DIFFERENT SURFACE
!             PRESSURES
     &   , P_TOA(NPD_PROFILE)                                           &
!             Pressure at the top of the atmosphere seen by radiation
     &   , Z_BOTTOM(NPD_PROFILE)                                        &
!             Height of the bottom of the current layer above
!             the surface
     &   , DZ
!             Depth of the current layer
!
!     TOTAL COLUMN MASS (KG M-2) OF EACH AEROSOL SPECIES IN
!     THE BOUNDARY LAYER, THE FREE TROPOSPHERE AND THE STRATOSPHERE
!     RESPECTIVELY. THIS MODEL ASSUMES THAT THERE ARE FIVE AEROSOLS.
      REAL                                                              &
     &     BL_OCEANMASS(5)                                              &
     &   , BL_LANDMASS(5)                                               &
     &   , FREETROP_MASS(5)                                             &
     &   , STRAT_MASS(5)
!
!     INITIALIZATION FOR THE CLIMATOLOGICAL AEROSOL MODEL
      DATA BL_LANDMASS/2.77579E-5, 6.70018E-5, 0.0, 9.57169E-7, 0.0/
      DATA BL_OCEANMASS/1.07535E-5, 0.0, 2.043167E-4, 0.0, 0.0/
      DATA FREETROP_MASS/3.46974E-6, 8.37523E-6, 0.0, 1.19646E-7, 0.0/
      DATA STRAT_MASS/0.0, 0.0, 0.0, 0.0, 1.86604E-6/
!
!
!
!     TROPOSPHERIC AEROSOL LOADING IS A SIMPLE FUNCTION OF SURFACE
!     PRESSURE: HALVING PSTAR HALVES THE TROPOSPHERIC AEROSOL BURDEN.
!     THE STRATOSPHERIC BURDEN IS INDEPENDENT OF PSTAR.  NOTE THE
!     FACTOR MULTIPLING AEROSOL AMOUNTS USES A REFERENCE PRESSURE
!     OF 1013 mbars.
      DO L=1, N_PROFILE
        PRESSURE_WT(L)=PSTAR(I_GATHER(L))*(1.0/1.013E5)
      END DO
!
!     For each of the 5 aerosol species, the column amount in the
!     boundary layer, free troposphere and stratosphere are known for
!     a standard atmosphere over ocean and land. These can be used
!     to find mixing ratios for the UM by dividing total aerosol by
!     total air mass (and using pressure weighting in the
!     troposphere).
!
!     Firstly, mixing ratios are set for the 5 aerosol species in the
!     stratosphere.
      IF (L_EXTRA_TOP) THEN
!       With an extra layer the radiative atmosphere is notionally
!       extended to zero pressure.
        DO L=1, N_PROFILE
          P_TOA(L)=0.0E+00
        ENDDO
      ELSE
!       Otherwise the top of the atmosphere seen elsewhere in the
!       model is used.
        DO L=1, N_PROFILE
          P_TOA(L)=P_LAYER_BOUNDARIES(I_GATHER(L), NLEVS)
        ENDDO
      ENDIF


! route through here depends on L_VOLCTS to control
! time-varying volcanic forcing.
      if (L_VOLCTS) then
       DO I=1,5
        DO L=1, N_PROFILE
         LG=I_GATHER(L)
         IF ( I .LT. 5 ) THEN 
           AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-TRINDX(LG),I) = 0.
         ELSE
           AEROSOL_MIX_RATIO_CLIM(L,NLEVS+1-TRINDX(LG),I) =             &
     &                                                VOLCMASS(LG) * G /&
     &        (P_LAYER_BOUNDARIES(LG,TRINDX(LG)-1)                      &
     &       - P_TOA(L))
         ENDIF
        END DO
       END DO
      else   ! L_VOLCTS is F
       DO I=1,5
        DO L=1, N_PROFILE
          LG=I_GATHER(L)
          AEROSOL_MIX_RATIO_CLIM(L, N_LAYER+1-TRINDX(LG), I)            &
     &      =STRAT_MASS(I)*G/                                           &
     &        (P_LAYER_BOUNDARIES(LG,TRINDX(LG)-1)                      &
     &       - P_TOA(L))
        END DO
       END DO
      endif   ! L_VOLCTS

      DO I=1,5
        DO L=1, N_PROFILE
          LG=I_GATHER(L)
            DO J=(TRINDX(LG)+1),N_LAYER
              AEROSOL_MIX_RATIO_CLIM(L, N_LAYER+1-J, I)=                &
     &          AEROSOL_MIX_RATIO_CLIM(L, N_LAYER+1-TRINDX(LG), I)
            END DO
         END DO
       END DO
!
!      At each point the number of layers considered to contain
!      boundary layer aerosol must be determined. In the original
!      form of the scheme this number is fixed, but in the newer
!      form it is determined from the boundary layer depth.
       IF (L_CLIM_AERO_HGT) THEN
!
!        Initialize:
         DO L=1, N_PROFILE
           BL_TOP_LYR(L)=N_LAYER
           L_SET_Z(L)=.TRUE.
           Z_BOTTOM(L)=0.0E+00
         ENDDO
         N_POINTS_SET_Z=N_PROFILE
           I=N_LAYER
!
         DO WHILE (N_POINTS_SET_Z >  0)
!
!          Assign the UM's indexing over layers: the UM now indexes the
!          boundaries of layers from 0 at the bottom to NLEVS at the
!          top, while the radiation code uses the same range but starts
!          at the top (possibly with an extra layer at the top of
!          the model). I and I_UM refer to layers, not to the edges of
!          layers. BL_TOP_LYR now holds the topmost layer containing
!          boundary layer aerosol indexed as in the radiation code.
           I_UM=N_LAYER+1-I
           DO L=1, N_PROFILE
             IF (L_SET_Z(L)) THEN
               LG=I_GATHER(L)
               DZ=R*T(L, I)*LOG(P_LAYER_BOUNDARIES(LG, I_UM-1)          &
     &           /P_LAYER_BOUNDARIES(LG, I_UM))/G
               IF ( (MAX(BL_DEPTH(LG), 5.0E+02) >=                      &
     &               Z_BOTTOM(L)+0.5E+00*DZ).AND.                       &
     &              (I_UM <= TRINDX(LG)-2) ) THEN
!                The top of the boundary layer is above the middle
!                of this layer, so we take the layer to contain
!                boundary layer aerosol, provisionally setting
!                the index of the top of the boundary layer to this
!                level: it will be overwritten if higher layers are
!                also in the boundary layer. A upper limit is applied
!                to the number of the layers that can be filled with
!                the boundary layer aerosol so that there is at least
!                one layer in the free troposphere.
                 BL_TOP_LYR(L)=I
!                Increment the height of the bottom of the layer
!                for the next pass.
                 Z_BOTTOM(L)=Z_BOTTOM(L)+DZ
               ELSE
!                The top of the boundary layer has been passed at this
!                point: it does not need to be considered further.
                 L_SET_Z(L)=.FALSE.
                 N_POINTS_SET_Z=N_POINTS_SET_Z-1
               ENDIF
             ENDIF
           ENDDO
           I=I-1
         ENDDO
       ELSE
         DO L=1, N_PROFILE
!          We do not allow the stratosphere to meet the boundary
!          layer, so we ensure that, counting up from the surface,
!          the highest layer allowed to be filled with boundary-
!          layer aerosol is the TRINDX-2nd: as a result, there must
!          be at least one layer in the free troposphere.
           BL_TOP_LYR(L)=N_LAYER+1                                      &
     &       -MIN(N_LEVELS_BL, TRINDX(I_GATHER(L))-2)
         ENDDO
       ENDIF
!
!      Now, the mixing ratios are set for the 5 aerosol species
!      in the free troposphere. Initially set the mixing ratio
!      in the lowest layer of the free troposphere.
       DO I=1,5
         DO L=1, N_PROFILE
           LG=I_GATHER(L)
           AEROSOL_MIX_RATIO_CLIM(L, BL_TOP_LYR(L)-1, I)                &
     &       =FREETROP_MASS(I)*G*                                       &
     &       PRESSURE_WT(L)/                                            &
     &       (P_LAYER_BOUNDARIES(LG, N_LAYER+1-BL_TOP_LYR(L))-          &
     &        P_LAYER_BOUNDARIES(LG,TRINDX(LG)-1) )
         END DO
       END DO
!      Fill in the remaining levels where necessary.
       DO L=1, N_PROFILE
         LG=I_GATHER(L)
         DO I=1,5
           DO J=N_LAYER+2-TRINDX(LG), BL_TOP_LYR(L)-2
             AEROSOL_MIX_RATIO_CLIM(L, J, I)=                           &
     &         AEROSOL_MIX_RATIO_CLIM(L, BL_TOP_LYR(L)-1, I)
           END DO
         END DO
       END DO
!
!      Now, the boundary layer mixing ratios are set for the
!      5 aerosol species. A continental aerosol is used over most land
!      areas, but not over ice sheets, which are identified by the
!      criterion used in the boundary layer scheme that the mass of
!      lying snow exceeds 5000 kgm-2. Over ice sheets a maritime
!      aerosol is used.
       DO I=1,5
         DO L=1, N_PROFILE
           LG=I_GATHER(L)
           IF ( LAND(LG).AND.(LYING_SNOW(LG) <  5.0E+03) ) THEN
             AEROSOL_MIX_RATIO_CLIM(L, BL_TOP_LYR(L), I)                &
     &         =BL_LANDMASS(I)*G*PRESSURE_WT(L)                         &
     &           /(PSTAR(LG)                                            &
     &           -P_LAYER_BOUNDARIES(LG, N_LAYER+1-BL_TOP_LYR(L)))
           ELSE
             AEROSOL_MIX_RATIO_CLIM(L, BL_TOP_LYR(L), I)                &
     &         =BL_OCEANMASS(I)*G*PRESSURE_WT(L)                        &
     &           /(PSTAR(LG)                                            &
     &           -P_LAYER_BOUNDARIES(LG, N_LAYER+1-BL_TOP_LYR(L)))
           END IF
         END DO
       END DO
       DO I=1,5
         DO L=1, N_PROFILE
           DO J=BL_TOP_LYR(L)+1, N_LAYER
            AEROSOL_MIX_RATIO_CLIM(L,J,I)=                              &
     &         AEROSOL_MIX_RATIO_CLIM(L, BL_TOP_LYR(L), I)
           END DO
         END DO
       END DO
!
!
!
      RETURN
      END SUBROUTINE R2_SET_AERO_CLIM_HADCM3
