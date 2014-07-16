#if defined(A25_1A)
      MODULE VELOC_MOD
      IMPLICIT NONE

      INTERFACE VELOC
        MODULE PROCEDURE VELOC
        MODULE PROCEDURE VELOC_V
      END INTERFACE ! VELOC

      CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      SUBROUTINE VELOC(pos,u,v,w,su,sv,sw,v1,hdiff,vdiff)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Return interpolated wind vector
!-
!-   Inputs  : POS,U,V,W,S
!-            (Longm,Latm,Longm_half,Latm_half,Dlongm,Dlatm are used
!-             from module ND_STOCH_MOD1)
!-   Outputs : 3-D interpolated wind vector: V1
!-   Controls:
!-
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.4  11/04/94   Created. W.J. Collins
!  4.1  10/07/96   No S.D.s - use diffusion coeffs. W.J. Collins
!  4.1  06/08/96   Parameters in INCLUDE. W.J. Collins
!  4.1  07/08/96   No longer have surface eta level for u,v,w.
!                  Moved I2 (I4) calculation to INTERP. Do
!                  interpolation over poles explicitly. W.J.Collins
!  4.1  30/09/96   u,v output in degrees/second. C.E.Johnson
!  4.4  08/10/97   Converted to F90 and parallelised. W.J.Collins
!  6.1  07/09/04   Added '&' for F90 compatability as needed.
!                  M.G. Sanderson
!  6.2  31/01/06   Placed names of interface statement end blocks
!                  behind comments for portability. T. Edwards
!
!-
!VVV  V5.0  VELOC 31/V/01
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE INTERP_MOD
      IMPLICIT NONE
! ----------------------------------------------------------------------

      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: u
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: v
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: w
      REAL,                           INTENT(IN) :: su  ! Random numbers
      REAL,                           INTENT(IN) :: sv  ! Random numbers
      REAL,                           INTENT(IN) :: sw  ! Random numbers
      REAL, DIMENSION(4),             INTENT(IN) :: pos ! Cell position
      REAL, INTENT(IN)                             :: hdiff
      REAL, INTENT(IN)                             :: vdiff

! Interpolated wind vector
      REAL, DIMENSION(3),              INTENT(OUT) :: v1

      REAL, PARAMETER :: deg_per_m=360.0/(2.0*pi*earth_radius)
      REAL            :: u7
      REAL            :: v7
      REAL            :: w7
      REAL            :: x
      REAL            :: y
      REAL            :: lambda  ! longitude angle
      REAL            :: theta   ! latitude angle

!kk      CHARACTER(LEN=72) :: cmessage

! DEPENDS ON: interp_s
      u7 = INTERP_S(pos,u,.FALSE.,.TRUE.,.TRUE.)
! DEPENDS ON: interp_s
      v7 = INTERP_S(pos,v,.FALSE.,.TRUE.,.FALSE.)
! DEPENDS ON: interp_s
      w7 = INTERP_S(pos,w,.TRUE.,.FALSE.)

      IF (pos(2)<=dlatm) THEN ! South pole
! convert to x and y coords, with point lying along y-axis
        x = 0.0
        y = pos(2)
! advect in this coordinate frame
        x = x + u7*stochem_advection_step*deg_per_m
        y = y + v7*stochem_advection_step*deg_per_m
! convert back to lat and long
        lambda = pos(1) + ATAN2(x,y)/pi_over_180
        theta = SQRT(x**2+y**2)
! work out velocities
        v1(1) = (lambda-pos(1))/stochem_advection_step
        v1(2) = (theta-pos(2))/stochem_advection_step
        v1(3) = w7
      ELSE IF (pos(2) >= 180.0-dlatm) THEN ! North pole
! convert to x and y coords, with point lying along y-axis
        x = 0.0
        y = 180.0 - pos(2)
! advect in this coordinate frame
        x = x - u7*stochem_advection_step*deg_per_m
        y = y - v7*stochem_advection_step*deg_per_m
! convert back to lat and long
        lambda = pos(1)-ATAN2(x,y)/pi_over_180
        theta = 180.0-SQRT(x**2+y**2)
! work out velocities
        v1(1) = (lambda-pos(1))/stochem_advection_step
        v1(2) = (theta-pos(2))/stochem_advection_step
        v1(3) = w7
      ELSE
! Add in random wind component
!       write(6,*) 'V1: ', u7,v7,w7
!       write(6,*) 'S: ', s,vdiff
!       write(6,*) 'POS: ', pos
!       write(6,*)
!     & 'Stochem_advection_step,HDIFF,VDIFF,pi_over_180,DPMY: ',
!     &  Stochem_advection_step,hdiff,vdiff,pi_over_180,dpmy

! Convert horizontal wind components to degrees per second.
! Debug: no random component
!        v1(1) = u7*deg_per_m/COS((pos(2)-90.0)*pi_over_180)
!        v1(2) = v7*deg_per_m
!        v1(3) = w7
        v1(1) = (u7 + su * SQRT(2.0*hdiff/stochem_advection_step))      &
     &    * deg_per_m / COS((pos(2)-90.0)*pi_over_180)
        v1(2) = (v7 + sv * SQRT(2.0*hdiff/stochem_advection_step))      &
     &    *deg_per_m
        v1(3) = w7 + sw * SQRT(2.0*vdiff/stochem_advection_step)
      END IF

      END SUBROUTINE VELOC

!kk   Vector Version of VELOC

      SUBROUTINE VELOC_V(nfill,todo,pos,u,v,w,su,sv,sw,v1,hdiff,vdiff)
!----------------------------------------------------------------------
!-
!-   Purpose and Methods : Return interpolated wind vector
!-
!-   Inputs  : POS,U,V,W,S
!-            (Longm,Latm,Longm_half,Latm_half,Dlongm,Dlatm are used
!-             from module ND_STOCH_MOD1)
!-   Outputs : 3-D interpolated wind vector: V1
!-   Controls:
!
! Current Code Owner: Bill Collins
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.5  21/01/04   Created. K. Ketelsen
!  6.1  07/09/04   Added '&' at end of lines for F90 compatability.
!                  M.G. Sanderson
!-
!VVV  V6.1  VELOC_V 07/09/04
! ----------------------------------------------------------------------
      USE IN_STOCHEM_GRD
      USE IN_STOCHEM_INTF
      USE INTERP_MOD              !kk
      IMPLICIT NONE
! ----------------------------------------------------------------------

      INTEGER, INTENT(IN)                :: nfill
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: u
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: v
      REAL, DIMENSION(nlonpe,nlatpe,0:nmetlev), INTENT(IN) :: w
      LOGICAL, DIMENSION(nfill), INTENT(IN) :: todo
      REAL, DIMENSION(4,nfill),  INTENT(IN) :: pos ! Cell position
      REAL, DIMENSION(nfill), INTENT(IN) :: su  ! Random numbers
      REAL, DIMENSION(nfill), INTENT(IN) :: sv  ! Random numbers
      REAL, DIMENSION(nfill), INTENT(IN) :: sw  ! Random numbers
      REAL, DIMENSION(nfill), INTENT(IN) :: hdiff
      REAL, DIMENSION(nfill), INTENT(IN) :: vdiff

! Interpolated wind vector
      REAL, DIMENSION(3,nfill), INTENT(OUT) :: v1

      INTEGER         :: k
      REAL, PARAMETER :: deg_per_m = 360.0 / (2.0*pi*earth_radius)
      REAL, DIMENSION(nfill) :: u7
      REAL, DIMENSION(nfill) :: v7
      REAL, DIMENSION(nfill) :: w7
      REAL            :: x
      REAL            :: y
      REAL            :: lambda  ! longitude angle
      REAL            :: theta   ! latitude angle
!
! Interpolate wind velocity vectors to parcel positions
!
      CALL INTERP_V(nfill,todo,pos,u,u7,.FALSE.,.TRUE.,.TRUE.)
      CALL INTERP_V(nfill,todo,pos,v,v7,.FALSE.,.TRUE.,.FALSE.)
      CALL INTERP_V(nfill,todo,pos,w,w7,.TRUE.,.FALSE.)

      DO k=1,nfill
        IF (todo(k)) THEN
          IF (pos(2,k) <= dlatm) THEN ! South pole
! convert to x and y coords, with point lying along y-axis
            x=0.0
            y=pos(2,k)
! advect in this coordinate frame
            x=x+u7(k)*stochem_advection_step*deg_per_m
            y=y+v7(k)*stochem_advection_step*deg_per_m
! convert back to lat and long
            lambda=pos(1,k)+ATAN2(x,y)/pi_over_180
            theta=SQRT(x**2+y**2)
! work out velocities
            v1(1,k)=(lambda-pos(1,k))/stochem_advection_step
            v1(2,k)=(theta-pos(2,k))/stochem_advection_step
            v1(3,k)=w7(k)
          ELSE IF(pos(2,k)>=180.-dlatm) THEN ! North pole
! convert to x and y coords, with point lying along y-axis
            x=0.0
            y=180.0-pos(2,k)
! advect in this coordinate frame
            x=x-u7(k)*stochem_advection_step*deg_per_m
            y=y-v7(k)*stochem_advection_step*deg_per_m
! convert back to lat and long
            lambda=pos(1,k)-ATAN2(x,y)/pi_over_180
            theta=180.0-SQRT(x**2+y**2)
! work out velocities
            v1(1,k)=(lambda-pos(1,k))/stochem_advection_step
            v1(2,k)=(theta-pos(2,k))/stochem_advection_step
            v1(3,k)=w7(k)
          ELSE

! Convert horizontal wind components to degrees per second.
! Debug: no random component

            v1(1,k)=(u7(k)+su(k)*SQRT(2.0*hdiff(k)/                     &
     &        stochem_advection_step))*                                 &
     &        deg_per_m/COS((pos(2,k)-90.0)*pi_over_180)
            v1(2,k)=(v7(k)+sv(k)*SQRT(2.0*hdiff(k)/                     &
     &        stochem_advection_step))*deg_per_m
            v1(3,k)=w7(k)+sw(k)*SQRT(2.0*vdiff(k)/                      &
     &        stochem_advection_step)
          END IF
        END IF
      END DO

      END SUBROUTINE VELOC_V

      END MODULE VELOC_MOD
#endif
