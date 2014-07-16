#if defined(A03_8A) || defined(A03_8B) || defined(A03_8C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!**********************************************************************
! Routine to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes
!
!
!    Programming standard:
!
!**********************************************************************
      SUBROUTINE SF_STOM  (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX          &
     &,                    VEG_PTS,VEG_INDEX                            &
     &,                    FT,CO2,CO2_3D,CO2_DIM_LEN                    &
     &,                    CO2_DIM_ROW,L_CO2_INTERACTIVE                &
     &,                    FSMC,HT,IPAR,LAI,PSTAR                       &
     &,                    Q1,RA,TSTAR                                  &
     &,                    CAN_RAD_MOD,ILAYERS,FAPARV                   &
     &,                    GPP,NPP,RESP_P,RESP_W,GC)


      IMPLICIT NONE

      INTEGER                                                           &
     & ROW_LENGTH                                                       &
                                  ! IN Number of points on a row
     &,ROWS                                                             &
                                  ! IN Number of rows in a theta field
     &,LAND_PTS                                                         &
                                  ! IN Number of land points to be
!                                 !    processed.
     &,LAND_INDEX(LAND_PTS)                                             &
                                  ! IN Index of land points on the
!                                 !    P-grid.
     &,VEG_PTS                                                          &
                                  ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_PTS)                                              &
                                  ! IN Index of vegetated points
!                                 !    on the land grid.
     &,CO2_DIM_LEN                                                      &
                                  ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                ! IN Number of CO2 field rows.

      INTEGER                                                           &
     & FT                         ! IN Plant functional type.
     
      LOGICAL L_CO2_INTERACTIVE   ! switch for 3D CO2 field
      
      INTEGER                                                           &
     &  CAN_RAD_MOD                                                     &
!                                  !Switch for canopy radiation model
     & ,ILAYERS
!                                  !No of layers in canopy radiation model

      REAL                                                              &
     & CO2                                                              &
                                  ! IN Atmospheric CO2 concentration
     &,CO2_3D(CO2_DIM_LEN,CO2_DIM_ROW)                                  &
!                                 ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
     &,FSMC(LAND_PTS)                                                   &
                                  ! IN Soil water factor.
     &,HT(LAND_PTS)                                                     &
                                  ! IN Canopy height (m).
     &,IPAR(ROW_LENGTH,ROWS)                                            &
                                  ! IN Incident PAR (W/m2).
     &,LAI(LAND_PTS)                                                    &
                                  ! IN Leaf area index.
     &,PSTAR(LAND_PTS)                                                  &
                                  ! IN Surface pressure (Pa).
     &,FAPARV(LAND_PTS,ILAYERS)                                         &
                                  ! IN Profile of absorbed PAR.
     &,Q1(ROW_LENGTH,ROWS)                                              &
                                  ! IN Specific humidity at level 1
     &,RA(LAND_PTS)                                                     &
                                  ! IN Aerodynamic resistance (s/m).
     &,TSTAR(LAND_PTS)                                                  &
                                  ! IN Surface temperature (K).
     &,GPP(LAND_PTS)                                                    &
                                  ! OUT Gross Primary Productivity
!                                 !     (kg C/m2/s).
     &,NPP(LAND_PTS)                                                    &
                                  ! OUT Net Primary Productivity
!                                 !     (kg C/m2/s).
     &,RESP_P(LAND_PTS)                                                 &
                                  ! OUT Plant respiration rate
!                                 !     (kg C/m2/sec).
     &,RESP_W(LAND_PTS)                                                 &
                                  ! OUT Wood respiration rate
!                                 !     (kg C/m2/sec).
     &,GC(LAND_PTS)                                                     &
                                  ! INOUT Canopy resistance to H2O
!                                 !       (m/s).


     &,VCMV(LAND_PTS,ILAYERS)                                           &
                                  ! Vcm-layer
     &,FAPARV_layer(LAND_PTS,ILAYERS)                                   &
                                  ! work/out to leaf ..absorbed par(layers)
     &,CI_CA(LAND_PTS,ILAYERS)                                          &
                                  ! work to store ci_Ca per layer
     &,ALAYER(LAND_PTS,ILAYERS)   ! work to store ALper layer


!  External routines called :-
      EXTERNAL QSAT,CANOPY,LEAF_C3,LEAF_C4,LEAF,LEAF_LIMITS



      REAL                                                              &
     & ANETC(LAND_PTS)                                                  &
                                  ! WORK Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
     &,CO2C(LAND_PTS)                                                   &
                                  ! WORK Canopy level CO2 concentration
!                                 !      (kg CO2/kg air).
     &,CI(LAND_PTS)                                                     &
                                  ! WORK Internal CO2 pressure (Pa).
     &,DQ(LAND_PTS)                                                     &
                                  ! WORK Specific humidity deficit
!                                 !      (kg H2O/kg air).
     &,DQC(LAND_PTS)                                                    &
                                  ! WORK Canopy level specific humidity
!                                 !      deficit (kg H2O/kg air).
     &,FPAR(LAND_PTS)                                                   &
                                  ! WORK PAR absorption factor.
     &,LAI_BAL(LAND_PTS)                                                &
                                  ! WORK Leaf area index in balanced
!                                 !      growth state.
     &,NL(LAND_PTS)                                                     &
                                  ! WORK Mean leaf nitrogen
!                                 !      concentration (kg N/kg C).
     &,NL_BAL(LAND_PTS)                                                 &
                                  ! WORK Mean leaf nitrogen
!                                 !      concentration in balanced
!                                 !      growth state (kg N/kg C).
     &,N_LEAF(LAND_PTS)                                                 &
                                  ! WORK Nitrogen contents of the leaf,
     &,N_ROOT(LAND_PTS)                                                 &
                                  !      root,
     &,N_STEM(LAND_PTS)                                                 &
                                  !      and stem (kg N/m2).
     &,QS(LAND_PTS)                                                     &
                                  ! WORK Saturated specific humidity
!                                 !      (kg H2O/kg air).
     &,RA_RC(LAND_PTS)                                                  &
                                  ! WORK Ratio of aerodynamic resistance
!                                 !      to canopy resistance.
     &,RDC(LAND_PTS)                                                    &
                                  ! WORK Canopy dark respiration,
!                                 !      without soil water dependence
!                                 !      (mol CO2/m2/s).
     &,RESP_P_G(LAND_PTS)                                               &
                                  ! WORK Plant growth respiration rate
!                                 !      (kg C/m2/sec).
     &,RESP_P_M(LAND_PTS)                                               &
                                  ! WORK Plant maintenance respiration
!                                 !      rate (kg C/m2/sec).
     &,ROOT(LAND_PTS)             ! WORK Root carbon (kg C/m2).

      INTEGER                                                           &
     & I,J,K,L,M,N                  ! WORK Loop counters.


!-----------------------------------------------------------------------
! Parameters
!-----------------------------------------------------------------------
      REAL, PARAMETER      :: MERCADO_KN = 0.78
                                  ! Nitrogen Allocation coefficient:
                                  ! Tellus 59B, 553-565 (2007), page 555.
      REAL, PARAMETER      :: CCONU = 12.0E-3
                                  ! kg C in 1 mol CO2
! Gas constant for dry air is defined 
      REAL, PARAMETER      :: O2 = 0.23
                                  ! Atmospheric concentration of
!                                 ! oxygen (kg O2/kg air).

      INTEGER, PARAMETER   :: ITER = 3
                                  ! Number of iterations to
!                                 ! determine the canopy climate.

#include "nstypes.h"
#include "trif.h"
#include "ccarbon.h"
#include "c_r_cp.h"          
                          ! For R = 287.05, Gas constant for dry air


!-----------------------------------------------------------------------
! Error reporting
!-----------------------------------------------------------------------
      CHARACTER(LEN=256)      :: RoutineName = "sf_stom.F90"     ! Name of the routine
      INTEGER                 :: ErrorStatus                     ! Error code
      CHARACTER(LEN=256)      :: Message                         ! Text for output.


      REAL                                                              &
     & ANETL(LAND_PTS)                                                  &
                                  ! WORK Net leaf photosynthesis
!                                 !      (mol CO2/m2/s/LAI).
     &,APAR(LAND_PTS)                                                   &
                                  ! WORK PAR absorbed by the top leaf
!                                 !      (W/m2).
     &,CA(LAND_PTS)                                                     &
                                  ! WORK Canopy level CO2 pressure
!                                 !      (Pa).
     &,DQM(LAND_PTS)                                                    &
                                  ! WORK Canopy level humidity
!                                 !      deficit (mol H2O/m3).
     &,GL(LAND_PTS)                                                     &
                                  ! WORK Leaf conductance for H2O
!                                 !      (m/s).
     &,OA(LAND_PTS)                                                     &
                                  ! WORK Atmospheric O2 pressure
!                                 !      (Pa).
     &,RD(LAND_PTS)                                                     &
                                ! WORK Dark respiration of top leaf
!                               !      (mol CO2/m2/s).
     &,WCARB(LAND_PTS)                                                  &
                                ! WORK Carboxylation, ...
     &,WLITE(LAND_PTS)                                                  &
                                !      ... Light, and ...
     &,WEXPT(LAND_PTS)                                                  &
                                !      ... export limited gross ...
!                               !      ... photosynthetic rates ...
!                               !      ... (mol CO2/m2/s).
     &,WLITEV(LAND_PTS)         ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for each layer
!                                 !      (mol CO2/m2/s).
      REAL                                                              &
     & DLAI(LAND_PTS)           ! WORK LAI Increment.

      INTEGER                                                           &

     & CLOS_INDEX(LAND_PTS)                                             &
                                  ! WORK Index of land points
!                                 !      with closed stomata.
     &,CLOS_PTS                                                         &
                                  ! WORK Number of land points
!                                 !      with closed stomata.
     &,OPEN_INDEX(LAND_PTS)                                             &
                                  ! WORK Index of land points
!                                 !      with open stomata.
     &,OPEN_PTS                   ! WORK Number of land points
!                                 !      with open stomata.

! Initialisations
      ANETL(:)=0.0
      NL(:)=0.0
!-----------------------------------------------------------------------
! Calculate the surface to level 1 humidity deficit and the surface
! density of the air
!-----------------------------------------------------------------------
! DEPENDS ON: qsat
      CALL QSAT(QS,TSTAR,PSTAR,LAND_PTS)
      DO M=1,VEG_PTS
        L = VEG_INDEX(M)
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        DQ(L) = MAX(0.0,(QS(L) - Q1(I,J)))
        CI(L)=0.0
      ENDDO

!-----------------------------------------------------------------------
! Calculate the PAR absorption factor
!-----------------------------------------------------------------------
      IF (CAN_RAD_MOD == 1) THEN
        DO M=1,VEG_PTS
          L = VEG_INDEX(M)
  
          FPAR(L) = (1 - EXP(-KPAR(FT)*LAI(L))) / KPAR(FT)
  
        ENDDO


!-----------------------------------------------------------------------
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration. Intially set the canopy humidity
! deficit using the previous value of GC.
!-----------------------------------------------------------------------
        DO K=1,ITER

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit and CO2 concentration
!-----------------------------------------------------------------------
          DO M=1,VEG_PTS
            L = VEG_INDEX(M)
            RA_RC(L) = RA(L) * GC(L)
            DQC(L) = DQ(L) / (1 + RA_RC(L))
          ENDDO
          IF (L_CO2_INTERACTIVE) THEN
!  use full 3D CO2 field
            DO M=1,VEG_PTS
              L = VEG_INDEX(M)
              J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
              I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
              CO2C(L) = CO2_3D(I,J)
            ENDDO
          ELSE
!  just use single CO2_MMR value
            DO M=1,VEG_PTS
              L = VEG_INDEX(M)
              CO2C(L) = CO2
            ENDDO
          ENDIF

!-----------------------------------------------------------------------
! Call CANOPY to calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
! DEPENDS ON: canopy
          CALL CANOPY (ROW_LENGTH,ROWS,LAND_PTS,LAND_INDEX              &
     &,              VEG_PTS,VEG_INDEX                                  &
     &,              FT,DQC,IPAR,TSTAR,CO2C,O2,PSTAR                    &
     &,              FPAR,FSMC,LAI                                      &
     &,              GC,ANETC,CI,RDC)

        ENDDO

        DO M=1,VEG_PTS
          L = VEG_INDEX(M)

!-----------------------------------------------------------------------
! Assume that root biomass is equal to balanced growth leaf biomass
!-----------------------------------------------------------------------
          LAI_BAL(L) = (A_WS(FT)*ETA_SL(FT)*HT(L)/A_WL(FT))             &
     &             **(1.0/(B_WL(FT)-1))
          ROOT(L) = SIGL(FT) * LAI_BAL(L)

!-----------------------------------------------------------------------
! Calculate the actual and balanced mean leaf nitrogen concentration
! assuming perfect light acclimation
!-----------------------------------------------------------------------
          NL(L) = (FPAR(L) / LAI(L)) * NL0(FT)
          NL_BAL(L) = (1 - EXP(-KPAR(FT)*LAI_BAL(L)))                   &
     &            / (KPAR(FT)*LAI_BAL(L)) * NL0(FT)

!-----------------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------
          N_LEAF(L) = NL(L) * SIGL(FT) * LAI(L)
          N_ROOT(L) = NR_NL(FT) * NL_BAL(L) * ROOT(L)
          N_STEM(L) = NS_NL(FT) * NL_BAL(L) * ETA_SL(FT) * HT(L) *      &
     &          LAI(L)

!-----------------------------------------------------------------------
! Calculate the Gross Primary Productivity, the plant maintenance
! respiration rate, and the wood maintenance respiration rate
! in kg C/m2/sec
!-----------------------------------------------------------------------
          GPP(L) = CCONU * (ANETC(L) + RDC(L)*FSMC(L))
          RESP_P_M(L) = CCONU * RDC(L)                                &
     &       * (N_LEAF(L)*FSMC(L) + N_STEM(L) + N_ROOT(L)) / N_LEAF(L)
          RESP_W(L) = CCONU * RDC(L) * N_STEM(L) / N_LEAF(L)

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------
          RESP_P_G(L) = R_GROW(FT) * (GPP(L) - RESP_P_M(L))
          RESP_P(L) = RESP_P_M(L) + RESP_P_G(L)
          NPP(L) = GPP(L) - RESP_P(L)

        ENDDO

        ELSEIF (CAN_RAD_MOD == 2) THEN

          DO M=1,VEG_PTS
            L = VEG_INDEX(M)
            ANETC(L) = 0.0
            GC(L) = 0.0
            GL(L)=0.0                      !lina
            RDC(L) = 0.0
            WLITEV(L) = WLITE(L)
          ENDDO

          DO N=1,ILAYERS
            DO M=1,VEG_PTS
              L = VEG_INDEX(M)
              GL(L)=0.0
            ENDDO

            DO K=1,ITER

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit and CO2 concentration
!-----------------------------------------------------------------------
            DO M=1,VEG_PTS
              L = VEG_INDEX(M)
              RA_RC(L) = RA(L) * GL(L)
              DQC(L) = DQ(L) / (1 + RA_RC(L))
            ENDDO
            IF (L_CO2_INTERACTIVE) THEN
!  use full 3D CO2 field
              DO M=1,VEG_PTS
                L = VEG_INDEX(M)
                J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
                I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
                CO2C(L) = CO2_3D(I,J)
              ENDDO
            ELSE
!  just use single CO2_MMR value
              DO M=1,VEG_PTS
                L = VEG_INDEX(M)
                CO2C(L) = CO2
              ENDDO
            ENDIF

!-----------------------------------------------------------------------
! Call CANOPY to calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
            DO M =1,VEG_PTS
              L = VEG_INDEX(M)
              J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
              I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
              CA(L) = CO2C(L) / EPCO2 * PSTAR(L)
              OA(L) = O2 / EPO2 * PSTAR(L)
              DQM(L) = DQC(L) / EPSILON * PSTAR(L) / (R * TSTAR(L))
!-----------------------------------------------------------------------
!     Calculate the PAR absorbed by the top leaf
!-----------------------------------------------------------------------
              APAR(L) = (1 - OMEGA(FT)) * IPAR(I,J)
              NL(L)       = NL0(FT)
!----------------------------------------------------------------------
!     Calculate the PAR absorbed at each layer
!-----------------------------------------------------------------------
              DLAI(L) = LAI(L)/FLOAT(ILAYERS)    ! if lai profile unknown
!             DLAI(L)=LAI(L)*LAI_STEP(N)         ! if lai profile known
              FAPARV_layer(L,N)=FAPARV(L,N)*DLAI(L) !lina
            ENDDO

!     lina trying   VN ie. max along the canopy
!-----------------------------------------------------------------------
! Calculate the limiting factors for leaf photosynthesis
!-----------------------------------------------------------------------
! The NL*exp((N-1)/FLOAT(ILAYERS)*(-MERCADO_KN)) term is the non-uniform 
! distribution.


! DEPENDS ON: leaf_limits
            CALL LEAF_LIMITS (LAND_PTS,VEG_PTS,VEG_INDEX,FT             &
     &        ,  NL*exp((N-1)/FLOAT(ILAYERS)*(-MERCADO_KN))             &
     &        ,  DQc,APAR,TSTAR,CA,OA                                   &
     &        ,PSTAR,FSMC                                               &
!    &  ,     NL,                    DQC,APAR,TSTAR,CA,OA,PSTAR,FSMC
     &        ,                   CLOS_PTS,OPEN_PTS,CLOS_INDEX          &
     &        ,OPEN_INDEX,     CI,RD,WCARB,WEXPT,WLITE                  &
     &        ,FAPARV_LAYER(:,N),VCMv(:,N))


            DO M=1,OPEN_PTS
              L = VEG_INDEX(OPEN_INDEX(M))
              J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
              I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
              WLITEV(L)=WLITE(L)/APAR(L)*FAPARV(L,N)*IPAR(I,J)
            ENDDO

! DEPENDS ON: leaf
            CALL LEAF (LAND_PTS,VEG_PTS,VEG_INDEX,FT                    &
     &        ,                CLOS_PTS,OPEN_PTS,CLOS_INDEX             &
     &        ,OPEN_INDEX,FSMC,TSTAR,CA,CI,RD,WCARB                     &
     &        ,WEXPT,WLITEV,GL,ANETL)

          ENDDO                 ! K-ITER

          DO M=1,VEG_PTS
            L = VEG_INDEX(M)
            J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
            I = LAND_INDEX(L) - (J-1)*ROW_LENGTH

            CI_CA(L,N)=CI(L)/CA(L)
            ANETC(L) = ANETC(L) + ANETL(L) * LAI(L)/FLOAT(ILAYERS)
            ALAYER(L,N)=ANETL(L)*LAI(L)/FLOAT(ILAYERS)
            GC(L) = GC(L) + GL(L) * LAI(L)/FLOAT(ILAYERS)
            RDC(L) = RDC(L)+ RD(L) * LAI(L)/FLOAT(ILAYERS)
          ENDDO

        ENDDO                   ! N LAYERS

!-----------------------------------------------------------------------
!     Calculate plant level respiration, NPP and GPP
!-----------------------------------------------------------------------

        DO M=1,VEG_PTS
          L = VEG_INDEX(M)

!-----------------------------------------------------------------------
! Assume that root biomass is equal to balanced growth leaf biomass
!-----------------------------------------------------------------------
          LAI_BAL(L) = (A_WS(FT)*ETA_SL(FT)*HT(L)/A_WL(FT))             &
     &               **(1.0/(B_WL(FT)-1))
          ROOT(L) = SIGL(FT) * LAI_BAL(L)

!-----------------------------------------------------------------------
! Calculate the actual and balanced mean leaf nitrogen concentration
! assuming perfect light acclimation
!-----------------------------------------------------------------------
          NL(L) = NL0(FT)
          NL_BAL(L) = NL0(FT)
          RDC(L) = RD(L) * LAI(L)
!---------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------
          N_LEAF(L) = NL(L) * SIGL(FT) * LAI(L)
          N_ROOT(L) = NR_NL(FT) * NL_BAL(L) * ROOT(L)
          N_STEM(L) = NS_NL(FT) * NL_BAL(L) * ETA_SL(FT) * HT(L) *      &
     &                  LAI(L)

!-----------------------------------------------------------------------
! Calculate the Gross Primary Productivity, the plant maintenance
! respiration rate, and the wood maintenance respiration rate
! in kg C/m2/sec
!-----------------------------------------------------------------------
          GPP(L) = CCONU * (ANETC(L) + RDC(L)*FSMC(L))
          RESP_P_M(L) = CCONU * RDC(L)                                &
     &       * (N_LEAF(L)*FSMC(L) + N_STEM(L) + N_ROOT(L)) / N_LEAF(L)
          RESP_W(L) = CCONU * RDC(L) * N_STEM(L) / N_LEAF(L)

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------
          RESP_P_G(L) = R_GROW(FT) * (GPP(L) - RESP_P_M(L))
          RESP_P(L) = RESP_P_M(L) + RESP_P_G(L)
          NPP(L) = GPP(L) - RESP_P(L)

        ENDDO

      ELSE

        IF (CAN_RAD_MOD == 3) THEN

          write(6,*)'CAN_RAD_MOD ERROR:'
          write(6,*)'YOU ARE RUNNING WITH CAN_RAD_MOD=3'
          write(6,*)'THIS VALUE IS NOT SUPPORTED IN THE UM'
          write(6,*)'SINCE THE COMPUTATIONAL COST OF RUNNING'
          write(6,*)'WITH THE FULL COMPLEXITY OPTION (=2) IS '
          write(6,*)'MINIMAL COMPARED TO ATM CALCULATIONS'
          write(6,*)'USE CAN_RAD_MOD=1(OLD) OR =2(NEW LIGHT MOD)'
          write(6,*)''
          write(6,*)'THE MODEL RUN HAS BEEN STOPPED IN sf_stom'
          write(Message,*) 'CAN_RAD_MOD value disallowed in sf_stom'
          ErrorStatus = 1
          CALL EREPORT (RoutineName, ErrorStatus,Message)


        ELSE

          write(6,*)'CAN_RAD_MOD ERROR:'
          write(6,*)'CAN_RAD_MOD = ',CAN_RAD_MOD
          write(6,*)'THIS VALUE IS NOT SUPPORTED - CHECK SETUP'
          write(6,*)'THE VARIABLE NITROGEN MOD IS ONLY DESIGNED FOR'
          write(6,*)'USE HERE WITH CAN_RAD_MOD=1,2'
          write(6,*)''
          write(6,*)'THE MODEL RUN HAS BEEN STOPPED IN sf_stom'
          ErrorStatus = 2
          write(Message,*) 'CAN_RAD_MOD value disallowed in sf_stom'
          CALL EREPORT (RoutineName, ErrorStatus,Message)

        ENDIF
      ENDIF


      RETURN
      END SUBROUTINE SF_STOM

!***********************************************************************
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)
!
! Written by Peter Cox (May 1995)
!***********************************************************************
#endif
