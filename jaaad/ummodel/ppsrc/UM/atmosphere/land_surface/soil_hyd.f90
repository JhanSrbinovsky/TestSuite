
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOIL_HYD-----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SOIL_HYD (NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX              &
     &,                    BEXP,DZ,EXT,FW,KS,KSZ,SATHH,TIMESTEP,V_SAT   &
     &,                    SLOW_RUNOFF,SMCL,STHU,SURF_ROFF,W_FLUX       &
     &,                    STF_SLOW_RUNOFF                              &
     &,                    ZW,STHZW,ZDEPTH,QBASE,QBASE_L                &
     &,                    DUN_ROFF,DRAIN,L_TOP,L_SOIL_SAT_DOWN         &
     &,                    LTIMER                                       &
     &)

! Include L_VG from land surface module
      USE land_surf_mod, ONLY :                                         &
     & L_VG_SOIL

      IMPLICIT NONE
!
! Description:
!     Increments the layer soil moisture contents and calculates
!     calculates gravitational runoff.
!
!
! Documentation : UM Documentation Paper 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.2    15/11/00   New Deck         M. Best
!  5.5  12/02/03   Code added for large-scale hydrology.
!                                          Nic Gedney.
!  6.0    02/09/03   Removed a multiple declaration. Luke Jones.
!  7.1    23/05/08  Add switch L_VG_SOIL to use Van Genuchten hydrology
!                   instead of Clapp Hornberger hydrology
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
! C_DENSTY for subroutine SF_EXCH
      REAL,PARAMETER:: RHOSEA = 1026.0 ! density of sea water (kg/m3)
      REAL,PARAMETER:: RHO_WATER = 1000.0! density of pure water (kg/m3)
! C_DENSTY end
! C_POND start

      ! Maximum amount of ponded water (kg/m2).
      REAL,PARAMETER:: POND_MAX = 0.0

! Maximum ponded fraction of the gridbox.
      REAL,PARAMETER:: POND_FRAC_MAX = 0.5

! C_POND end
! C_SOILH start
      ! No. of soil layers (must = NSOIL).
      REAL,PARAMETER:: PSOIL=4

      ! Tunable characteristic freq (rad/s)
      REAL,PARAMETER:: OMEGA1=3.55088E-4

      ! Density of lying snow (kg per m**3)
      REAL,PARAMETER:: RHO_SNOW=250.0

      ! Depth of `effective' snow surface layer (m)
      REAL,PARAMETER:: DEFF_SNOW=0.1

      ! Thermal conductivity of lying snow (Watts per m per K).
      REAL,PARAMETER:: SNOW_HCON=0.265

      ! Thermal capacity of lying snow (J/K/m3)
      REAL,PARAMETER:: SNOW_HCAP=0.63E6

! C_SOILH end
! C_TOPOG start
! 5.5 17/02/03    Required for large-scale hydrology L_TOP code.
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 4.0
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 5.5
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                                                            &
                            ! IN Number of soil moisture levels.
     &,SOIL_PTS             ! IN Number of soil points.

      REAL                                                              &
     & TIMESTEP             ! IN Model timestep (s).

      LOGICAL                                                           &
     & STF_SLOW_RUNOFF                                                  &
                            ! IN Stash flag for sub-surface runoff.
     &,L_TOP                                                            &
                            ! IN Flag for TOPMODEL-based hydrology.
     &,L_SOIL_SAT_DOWN      ! IN Direction of super-sat soil moisture



!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL                                                              &
     & BEXP(NPNTS)                                                      &
                            ! IN Clapp-Hornberger exponent.
     &,DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,EXT(NPNTS,NSHYD)                                                 &
                            ! IN Extraction of water from each soil
!                           !    layer (kg/m2/s).
     &,FW(NPNTS)                                                        &
                            ! IN Throughfall from canopy plus snowmelt
!                           !    minus surface runoff (kg/m2/s).
     &,SATHH(NPNTS)                                                     &
                            ! IN Saturated soil water pressure (m).
     &,V_SAT(NPNTS)                                                     &
                            ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
     &,KS(NPNTS)                                                        &
                            ! IN Saturated hydraulic
!                           !    conductivity at surface (kg/m2/s).
     &,KSZ(NPNTS,NSHYD)                                                 &
                            ! IN Saturated hydraulic conductivity
!                           !    in each soil layer (kg/m2/s).
     &,QBASE(NPNTS)                                                     &
                            ! IN Base flow (kg/m2/s).
     &,QBASE_L(NPNTS,NSHYD+1)                                           &
!                           ! IN Base flow from each level (kg/m2/s).
!                           !    (kg/m2/s).
     &,ZDEPTH(0:NSHYD)                                                  &
                            ! IN Soil layer depth at lower boundary (m).
     &,ZW(NPNTS)            ! IN Mean water table depth (m).

!
      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & SLOW_RUNOFF(NPNTS)                                               &
                            ! OUT Drainage from the base of the
!                           !     soil profile (kg/m2/s).
     &,SMCLSAT(NPNTS,NSHYD)                                             &
                            ! OUT The saturation moisture content of
!                           !     each layer (kg/m2).
     &,W_FLUX(NPNTS,0:NSHYD)                                            &
                            ! OUT The fluxes of water between layers
!                           !     (kg/m2/s).
     &,DUN_ROFF(NPNTS)                                                  &
                            ! OUT Cumm. Dunne surface runoff (kg/m2/s).
     &,DRAIN(NPNTS)         ! OUT Drainage out of nshyd'th level (kg/m2/s).

!   Array arguments with intent(INOUT) :
      REAL                                                              &
     & SMCL(NPNTS,NSHYD)                                                &
                            ! INOUT Total soil moisture contents
!                           !       of each layer (kg/m2).
     &,STHU(NPNTS,NSHYD)                                                &
                            ! INOUT Unfrozen soil moisture content of ea
!                           !       layer as a fraction of saturation.
     &,SURF_ROFF(NPNTS)                                                 &
                            ! INOUT Surface runoff (kg/m2/s).
     &,STHZW(NPNTS)         ! INOUT soil moist fraction in deep layer.


!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL HYD_CON,DARCY,GAUSS,CALC_ZW,HYD_CON_VG,DARCY_VG

      EXTERNAL TIMER

!-----------------------------------------------------------------------

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.

      REAL                                                              &
     & GAMCON               ! WORK Constant (s/mm).

! Local arrays:
      REAL                                                              &
     & A(NPNTS,NSHYD)                                                   &
                            ! WORK Matrix elements corresponding
!                           !      to the coefficients of DSTHU(n-1).
     &,B(NPNTS,NSHYD)                                                   &
                            ! WORK Matrix elements corresponding
!                           !      to the coefficients of DSTHU(n).
     &,C(NPNTS,NSHYD)                                                   &
                            ! WORK Matrix elements corresponding
!                           !      to the coefficients of DSTHU(n+1).
     &,D(NPNTS,NSHYD)                                                   &
                            ! WORK Matrix elements corresponding
!                           !      to the RHS of the equation.
     &,DSMCL(NPNTS,NSHYD)                                               &
                            ! WORK Soil moisture increment
!                           !      (kg/m2/timestep).
     &,DSTHU(NPNTS,NSHYD)                                               &
                            ! WORK Increment to STHU (/timestep).
     &,DSTHUMIN(NPNTS,NSHYD)                                            &
                            ! WORK Minimum value of DSTHU.
     &,DSTHUMAX(NPNTS,NSHYD)                                            &
                            ! WORK Maximum value of DSTHU.
     &,DWFLUX_DSTHU1(NPNTS,NSHYD)                                       &
                                  ! WORK The rate of change of the expli
!                           !      flux with STHU1 (kg/m2/s).
     &,DWFLUX_DSTHU2(NPNTS,NSHYD)                                       &
                                  ! WORK The rate of change of the expli
!                           !      flux with STHU2 (kg/m2/s).
     &,SMCLU(NPNTS,NSHYD)                                               &
                            ! WORK Unfrozen soil moisture contents
!                           !      of each layer (kg/m2).
     &,SMCLZW(NPNTS)                                                    &
                            ! WORK moisture content in deep layer.
     &,SMCLSATZW(NPNTS)                                                 &
                            ! WORK moisture content in deep layer
!                           !      at saturation.
     &,DW                                                               &
     &,DWZW


! Local parameters:
      REAL                                                              &
     & GAMMA                ! Forward timestep weighting.
      PARAMETER (GAMMA=1.0)


      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SOILHYD ',103)
      ENDIF
!-----------------------------------------------------------------------
! Calculate the unfrozen soil moisture contents and the saturation
! total soil moisture for each layer.
!-----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I)
          SMCLU(I,N)=STHU(I,N)*SMCLSAT(I,N)
          DSTHUMIN(I,N)=-STHU(I,N)
          DSTHUMAX(I,N)=1.0-SMCL(I,N)/SMCLSAT(I,N)
          DWFLUX_DSTHU1(I,N)=0.0
          DWFLUX_DSTHU2(I,N)=0.0
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Top boundary condition and moisture in deep layer:

!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        W_FLUX(I,0)=FW(I)
        SMCLSATZW(I)=RHO_WATER*V_SAT(I)                                 &
     &       *(ZW_MAX-ZDEPTH(NSHYD))
        SMCLZW(I)=STHZW(I)*SMCLSATZW(I)

      ENDDO

!-----------------------------------------------------------------------
! Calculate the Darcian fluxes and their dependencies on the soil
! moisture contents.
!-----------------------------------------------------------------------

! If L_VG_SOIL is T then Van Genuchten formulation is used otherwise
! Clapp Hornberger using Cosby parameters is used

      IF(L_VG_SOIL) THEN

! DEPENDS ON: hyd_con_vg
        CALL HYD_CON_VG (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP                 &
     &,               KSZ(1,NSHYD),STHU(1,NSHYD)                        &
     &,               W_FLUX(1,NSHYD),DWFLUX_DSTHU1(1,NSHYD)            &
     &,               LTIMER)

        DO N=2,NSHYD
! DEPENDS ON: darcy_vg
          CALL DARCY_VG (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KSZ(1,N-1)      &
     &,               SATHH,STHU(1,N-1),DZ(N-1),STHU(1,N),DZ(N)         &
     &,               W_FLUX(1,N-1),DWFLUX_DSTHU1(1,N-1)                &
     &,               DWFLUX_DSTHU2(1,N-1),LTIMER)

        ENDDO

      ELSE


! DEPENDS ON: hyd_con
        CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP                    &
     &,               KSZ(1,NSHYD),STHU(1,NSHYD)                        &
     &,               W_FLUX(1,NSHYD),DWFLUX_DSTHU1(1,NSHYD)            &
     &,               LTIMER)

        DO N=2,NSHYD
! DEPENDS ON: darcy
          CALL DARCY (NPNTS,SOIL_PTS,SOIL_INDEX,BEXP,KSZ(1,N-1),SATHH   &
     &,               STHU(1,N-1),DZ(N-1),STHU(1,N),DZ(N),W_FLUX(1,N-1) &
     &,               DWFLUX_DSTHU1(1,N-1),DWFLUX_DSTHU2(1,N-1)         &
     &,               LTIMER)
         ENDDO

      ENDIF

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation in deep layer:
!-----------------------------------------------------------------------
      IF(L_TOP)THEN
         DO J=1,SOIL_PTS
            I=SOIL_INDEX(J)
            DWZW=(SMCLZW(I)-SMCLSATZW(I))/TIMESTEP+                     &
     &           (W_FLUX(I,NSHYD)-QBASE_L(I,NSHYD+1))
            IF(DWZW >  0.0 ) THEN
               W_FLUX(I,NSHYD)=W_FLUX(I,NSHYD)-DWZW
            ENDIF
         ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Calculate the explicit increments.
! The calculation of the increments depends on the direction of
! super-saturated soil moisture (down if L_SOIL_SAT_DOWN, else up)
!-----------------------------------------------------------------------
      IF (L_SOIL_SAT_DOWN) THEN
        ! super-saturated soil moisture moves down
        DO N=1,NSHYD
          DO J=1,SOIL_PTS
            I=SOIL_INDEX(J)
            DSMCL(I,N)=(W_FLUX(I,N-1)-W_FLUX(I,N)-EXT(I,N))*TIMESTEP
            IF(L_TOP)DSMCL(I,N)=DSMCL(I,N)-QBASE_L(I,N)*TIMESTEP

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation.
!-----------------------------------------------------------------------
            IF (DSMCL(I,N) >  (SMCLSAT(I,N)-SMCL(I,N))) THEN
              DSMCL(I,N)=SMCLSAT(I,N)-SMCL(I,N)
              W_FLUX(I,N)=W_FLUX(I,N-1)-DSMCL(I,N)/TIMESTEP-EXT(I,N)
              IF(L_TOP)W_FLUX(I,N)=W_FLUX(I,N)-QBASE_L(I,N)
            ENDIF
!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent negative soil moisture.
! This MAY no longer be required!!!!
!-----------------------------------------------------------------------
            IF(L_TOP)THEN
              IF(SMCL(I,N)+DSMCL(I,N) <  0.0.AND.N >  1)THEN
                DSMCL(I,N)=SMCL(I,N)
                W_FLUX(I,N)=W_FLUX(I,N-1)-EXT(I,N-1)-QBASE_L(I,N)         &
     &          -DSMCL(I,N)/TIMESTEP
              ENDIF
            ENDIF


          ENDDO
        ENDDO
      ELSE
        ! super-saturated soil moisture moves up
        DO N=NSHYD,1,-1
          DO J=1,SOIL_PTS
            I=SOIL_INDEX(J)
            DSMCL(I,N)=(W_FLUX(I,N-1)-W_FLUX(I,N)-EXT(I,N))*TIMESTEP
            IF(L_TOP)DSMCL(I,N)=DSMCL(I,N)-QBASE_L(I,N)*TIMESTEP

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation.
!-----------------------------------------------------------------------
            IF (DSMCL(I,N) >  (SMCLSAT(I,N)-SMCL(I,N))) THEN
              DSMCL(I,N)=SMCLSAT(I,N)-SMCL(I,N)
              W_FLUX(I,N-1)=DSMCL(I,N)/TIMESTEP+W_FLUX(I,N)+EXT(I,N)
              IF(L_TOP)W_FLUX(I,N-1)=W_FLUX(I,N-1)+QBASE_L(I,N)
            ENDIF
!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent negative soil moisture.
! This MAY no longer be required!!!!
!-----------------------------------------------------------------------
            IF(L_TOP)THEN
              IF(SMCL(I,N)+DSMCL(I,N) <  0.0.AND.N >  1)THEN
                DSMCL(I,N)=SMCL(I,N)
                W_FLUX(I,N)=W_FLUX(I,N-1)-EXT(I,N-1)-QBASE_L(I,N)         &
     &          -DSMCL(I,N)/TIMESTEP
              ENDIF
            ENDIF


          ENDDO
        ENDDO

      END IF

!-----------------------------------------------------------------------
! Calculate the matrix elements required for the implicit update.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        GAMCON=GAMMA*TIMESTEP/SMCLSAT(I,1)
        A(I,1)=0.0
        B(I,1)=1.0+GAMCON*DWFLUX_DSTHU1(I,1)
        C(I,1)=GAMCON*DWFLUX_DSTHU2(I,1)
        D(I,1)=DSMCL(I,1)/SMCLSAT(I,1)
      ENDDO

      DO N=2,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          GAMCON=GAMMA*TIMESTEP/SMCLSAT(I,N)
          A(I,N)=-GAMCON*DWFLUX_DSTHU1(I,N-1)
          B(I,N)=1.0-GAMCON*(DWFLUX_DSTHU2(I,N-1)-DWFLUX_DSTHU1(I,N))
          C(I,N)=GAMCON*DWFLUX_DSTHU2(I,N)
          D(I,N)=DSMCL(I,N)/SMCLSAT(I,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Solve the triadiagonal matrix equation.
!-----------------------------------------------------------------------
! DEPENDS ON: gauss
      CALL GAUSS(NSHYD,NPNTS,SOIL_PTS,SOIL_INDEX,A,B,C,D                &
     &,          DSTHUMIN,DSTHUMAX,DSTHU)

!-----------------------------------------------------------------------
! Diagnose the implicit fluxes.
!-----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          DSMCL(I,N)=DSTHU(I,N)*SMCLSAT(I,N)
          W_FLUX(I,N)=W_FLUX(I,N-1)-EXT(I,N)-DSMCL(I,N)/TIMESTEP
          IF(L_TOP)W_FLUX(I,N)=W_FLUX(I,N)-QBASE_L(I,N)
        ENDDO
      ENDDO
      IF(L_TOP)THEN
         DO J=1,SOIL_PTS
            I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Limit explicit fluxes to prevent negative soil moisture in deep layer:
!-----------------------------------------------------------------------
            IF(SMCL(I,NSHYD)+DSMCL(I,NSHYD) <  0.0)THEN
               DSMCL(I,NSHYD)=SMCL(I,NSHYD)
               W_FLUX(I,NSHYD)=W_FLUX(I,NSHYD-1)                        &
     &              -EXT(I,NSHYD-1)-QBASE_L(I,NSHYD)                    &
     &              -DSMCL(I,NSHYD)/TIMESTEP
            ENDIF

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation in the deep layer.
! Adjust drainage flux out of layer above:
!-----------------------------------------------------------------------
            DWZW=(SMCLZW(I)-SMCLSATZW(I))/TIMESTEP+                     &
     &           (W_FLUX(I,NSHYD)-QBASE_L(I,NSHYD+1))
            IF(DWZW >  0.0 ) THEN
               W_FLUX(I,NSHYD)=W_FLUX(I,NSHYD)-DWZW
               DSMCL(I,NSHYD)=DSMCL(I,NSHYD)+                           &
     &              DWZW*TIMESTEP
            ENDIF
         ENDDO

!-----------------------------------------------------------------------
! Limit the explicit fluxes to prevent supersaturation in soil layers:
!-----------------------------------------------------------------------
         DO N=NSHYD,1,-1
            DO J=1,SOIL_PTS
               I=SOIL_INDEX(J)
               DW=(SMCL(I,N)+DSMCL(I,N)-SMCLSAT(I,N))/TIMESTEP
               IF (DW >= 0.0) THEN
                  DSMCL(I,N)=SMCLSAT(I,N)-SMCL(I,N)
                  W_FLUX(I,N-1)=W_FLUX(I,N-1)-DW
                  IF(N /= 1)DSMCL(I,N-1)=DSMCL(I,N-1)+DW*TIMESTEP
               ENDIF
            ENDDO
         ENDDO
      ENDIF


!-----------------------------------------------------------------------
! Update the prognostic variables.
!-----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCLU(I,N)=SMCLU(I,N)+DSMCL(I,N)
          SMCL(I,N)=SMCL(I,N)+DSMCL(I,N)
          STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)
        ENDDO
      ENDDO

      IF(L_TOP)THEN

!-----------------------------------------------------------------------
! Diagnose the new water table depth.
! Assume local equilibrium psi profile
!-----------------------------------------------------------------------

        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCLZW(I)=(SMCLZW(I)-(QBASE_L(I,NSHYD+1)-W_FLUX(I,NSHYD))*    &
     &         TIMESTEP)
! Update prognostic deep layer soil moisture fraction:
          STHZW(I)=SMCLZW(I)/SMCLSATZW(i)
        ENDDO

! DEPENDS ON: calc_zw
        CALL CALC_ZW(NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX                    &
     &   ,BEXP,SATHH,SMCL,SMCLZW,SMCLSAT,SMCLSATZW,V_SAT,ZW)

!-----------------------------------------------------------------------
! Dont allow negative base flows:
!-----------------------------------------------------------------------
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          QBASE(I)=0.0
          DO N=1,NSHYD+1
             QBASE_L(I,N)=MAX(QBASE_L(I,N),0.0)
             QBASE(I)=QBASE(I)+QBASE_L(I,N)
          ENDDO
        ENDDO

      ENDIF
!-----------------------------------------------------------------------
! Output slow runoff (drainage) diagnostic.
!-----------------------------------------------------------------------
      IF (STF_SLOW_RUNOFF) THEN
        DO I=1,NPNTS
          SLOW_RUNOFF(I)=0.0
        ENDDO
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
! Ensure correct field is output as deep runoff: dependant on L_TOP
          IF(L_TOP)THEN
            SLOW_RUNOFF(I)=QBASE(I)
            DRAIN(I)=W_FLUX(I,NSHYD)
          ELSE
            SLOW_RUNOFF(I)=W_FLUX(I,NSHYD)
          ENDIF

        ENDDO
      ENDIF

!-----------------------------------------------------------------------
! Update surface runoff diagnostic.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        SURF_ROFF(I)=SURF_ROFF(I)+(FW(I)-W_FLUX(I,0))
      ENDDO
      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('SOILHYD ',104)
      ENDIF
      RETURN
      END SUBROUTINE SOIL_HYD
