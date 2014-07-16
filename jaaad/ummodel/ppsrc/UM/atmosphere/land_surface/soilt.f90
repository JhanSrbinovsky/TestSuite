
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE SOILT--------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SOILT(NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX                  &
     &  ,DZ,TSOIL,TSOIL_D )

      IMPLICIT NONE
!
! Description:
!     Diagnoses the soil temperature in a layer at the surface
!
!
! Documentation :
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5    03/02/03   New Deck         Nic Gedney
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!
      INTEGER                                                           &
     & NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,NSHYD                                                            &
                            ! IN Number of soil moisture levels.
     &,SOIL_PTS                                                         &
                            ! IN Number of soil points.
     &,SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL                                                              &
     & DZ(NSHYD)                                                        &
                            ! IN Thicknesses of the soil layers (m).
     &,TSOIL(NPNTS,NSHYD)   ! IN Unfrozen soil moisture content of

      REAL                                                              &
     & TSOIL_D(NPNTS)           ! OUT Soil moisture (kg/m2).


      REAL                                                              &
     & Z1,Z2                                                            &
                            ! WORK Depth of the top and bottom of the
!                           !      soil layers (m).
     &,ZST                  ! WORK Depth of layer for soil moisture
!                           !      diagnostic (m).
      PARAMETER ( ZST = 1. )

      INTEGER                                                           &
     & I,J,N

      Z2 = 0.

      DO I=1,NPNTS
        TSOIL_D(I) = 0.0
      ENDDO

      DO N=1,NSHYD
        Z1 = Z2
        Z2 = Z2 + DZ(N)
        IF ( Z2 <  ZST ) THEN
          DO J=1,SOIL_PTS
            I=SOIL_INDEX(J)
            TSOIL_D(I) = TSOIL_D(I) + DZ(N) * TSOIL(I,N)
          ENDDO
        ELSEIF ( Z2 >= ZST .AND. Z1 <  ZST ) THEN
          DO J=1,SOIL_PTS
            I = SOIL_INDEX(J)
            TSOIL_D(I) = TSOIL_D(I) + ( ZST - Z1 ) * TSOIL(I,N)
          ENDDO
        ENDIF

      ENDDO

      RETURN
      END SUBROUTINE SOILT
