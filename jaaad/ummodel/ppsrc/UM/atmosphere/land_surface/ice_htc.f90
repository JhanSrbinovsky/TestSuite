
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE ICE_HTC------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE ICE_HTC (                                              &
     & NPNTS,NSHYD,LICE_PTS,LICE_INDEX,DZ                               &
     &,SURF_HT_FLUX,TIMESTEP                                            &
     &,TSOIL                                                            &
! LOGICAL LTIMER
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE
!
! Description:
!     Updates deep soil temperatures for ice. No external subroutines
!     are called.
!                                                     (Cox, 10/95)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1       6/96    New deck.   Peter Cox
!LL   4.5   18/06/98  Changed Timer calls to indicate non-barrier
!LL                                                   P.Burton
!  6.0        06/02/04    Added compiler directive            A. Malcolm
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
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

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & LICE_PTS                                                         &
                            ! IN Number of land ice points.
     &,NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil moisture levels.

      REAL                                                              &
     & TIMESTEP             ! IN Model timestep (s).


!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & LICE_INDEX(NPNTS)    ! IN Array of ice points.

      REAL                                                              &
     & DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,SURF_HT_FLUX(NPNTS)  ! IN Net downward surface heat flux (W/m2).
!
      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(INOUT) :
      REAL                                                              &
     & TSOIL(NPNTS,NSHYD)   ! INOUT Sub-surface temperatures (K).

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.

! Local arrays:
      REAL                                                              &
     & H_FLUX(NPNTS,0:NSHYD)! WORK The fluxes of heat between layers
!                           !      (W/m2).

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('ICEHTC  ',103)
      ENDIF

!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries
!--------------------------------------------------------------------
      DO N=1,NSHYD-1
        DO J=1,LICE_PTS
          I=LICE_INDEX(J)
          H_FLUX(I,N)=-SNOW_HCON*2.0*(TSOIL(I,N+1)-TSOIL(I,N))          &
     &                             /(DZ(N+1)+DZ(N))
        ENDDO
      ENDDO

!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,LICE_PTS
        I=LICE_INDEX(J)
        H_FLUX(I,NSHYD)=0.0
        H_FLUX(I,0)=SURF_HT_FLUX(I)
      ENDDO


!--------------------------------------------------------------------
! Update the sub-surface temperatures
!--------------------------------------------------------------------
      DO N=1,NSHYD
!CDIR NOVECTOR
! CDIR$ IVDEP here would force vectorization but changes results!
        DO J=1,LICE_PTS
          I=LICE_INDEX(J)

          TSOIL(I,N)=TSOIL(I,N)                                         &
     &     +1.0/(SNOW_HCAP*DZ(N))*(H_FLUX(I,N-1)                        &
     &     -H_FLUX(I,N))*TIMESTEP

        ENDDO
      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('ICEHTC  ',104)
      ENDIF

      RETURN
      END SUBROUTINE ICE_HTC
