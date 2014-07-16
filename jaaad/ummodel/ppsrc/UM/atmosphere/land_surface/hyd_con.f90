
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE HYD_CON------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,THETAK,K       &
     &,                   DK_DTHK                                       &
! LOGICAL LTIMER
     &,LTIMER                                                           &
     &)

      IMPLICIT NONE
!
! Description:
!     Calculates the hydraulic conductivity
!
!
! Documentation : UM Documentation Paper 25
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.2   15/11/00   New Deck         M. Best
!  6.1  08/12/03  Add !CDIR NODEP to force vectorisation. R Barnes
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!
! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NPNTS                                                            &
                        ! IN points in grid
     &,SOIL_PTS         ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)! IN Array of soil points.

      REAL                                                              &
     & B(NPNTS)                                                         &
                        ! IN Exponent in conductivity and soil water
!                       !    suction fits.
     &,KS(NPNTS)                                                        &
                        ! IN The saturated hydraulic conductivity (kg/m2
     &,THETAK(NPNTS)    ! IN Fractional saturation.
!
      LOGICAL LTIMER    ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & K(NPNTS)                                                         &
                        ! OUT The hydraulic conductivity (kg/m2/s).
     &,DK_DTHK(NPNTS)   ! OUT The rate of change of K with THETAK
!                       !     (kg/m2/s).

!---------------------------------------------------------------------
!  External routines called :-

      EXTERNAL TIMER

!-----------------------------------------------------------------------

! Local scalars:
      INTEGER                                                           &
     & I,J              ! WORK Loop counter.

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HYDCON  ',103)
      ENDIF

!CDIR NODEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)

        DK_DTHK(I)=0.0
        IF (THETAK(I) >= 0.0.AND.THETAK(I) <  1.0) THEN
          K(I)=KS(I)*THETAK(I)**(2*B(I)+3)
          DK_DTHK(I)=(2*B(I)+3)*KS(I)*(THETAK(I)**(2*B(I)+2))
        ELSEIF (THETAK(I) <  0.0) THEN
          K(I)=0.0
        ELSE
          K(I)=KS(I)
        ENDIF

      ENDDO

      IF (LTIMER) THEN
! DEPENDS ON: timer
        CALL TIMER('HYDCON  ',104)
      ENDIF

      RETURN
      END SUBROUTINE HYD_CON
