
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE CH4_WETL-----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE CH4_WETL(NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX               &
     &  ,TSOIL_D,CS,F_WETL,FCH4_WETL)

      IMPLICIT NONE
!
! Description:
!     Calculates methane emissions from wetland area.
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
     &,SOIL_INDEX(NPNTS)   ! IN Array of soil points.

      REAL                                                              &
     & TSOIL_D(NPNTS)                                                   &
                           ! IN Diagnosted soil temp to 1metre.
     &,CS(NPNTS)                                                        &
                           ! IN Soil carbon (kg C/m2).
     &,F_WETL(NPNTS)       ! IN Wetland fraction

      REAL                                                              &
     & FCH4_WETL(NPNTS)    ! OUT Scaled methane flux (10^-9 kg C/m2/s)

      REAL                                                              &
     & Q10T_CH4                                                         &
                           ! Q10 value at T
     &,CONST_TDEP          ! T and Q10(0) dependent function

! C_CH4 start
! 5.5 17/02/03    Required for large-scale hydrology L_TOP code.
!                 Used in calculation of methane flux from wetlands.
!
      REAL,PARAMETER :: T0_CH4 = 273.15               ! T0 value
      REAL,PARAMETER :: Q10_CH4 = 3.7                 ! Q10 value
      REAL,PARAMETER :: CONST_CH4 = 7.41E-12          ! Scale factor
!
! C_CH4 end

      INTEGER                                                           &
     & I,J

      CONST_TDEP=T0_CH4*LOG(Q10_CH4)
      DO J=1,SOIL_PTS
         I=SOIL_INDEX(J)
         IF(TSOIL_D(I) >  273.15.AND.F_WETL(I) >  0.0)THEN
           Q10T_CH4=EXP(CONST_TDEP/TSOIL_D(I))
           FCH4_WETL(I)=1.E9*CONST_CH4*CS(I)*                           &
     &      F_WETL(I)*Q10T_CH4**(0.1*(TSOIL_D(I)-T0_CH4))
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE CH4_WETL

