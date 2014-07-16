#if defined(FLDCALC)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate Max Wind value and height


!=======================================================================


! ====================================================================



! ====================================================================

SUBROUTINE MaxWindSplineV(NumCols, NumRows, NumLevs,           &
                          UFields, VFields, PFields, MaxWLev,  &
                          Uinc, Vinc, Pinc )

! Description:
!
! Method:
!
! Owner: Dave Robinson
!
! History:
! Version Date     Comment
! ------- ----     -------
! 6.2     19/12/05 New subroutine. Vectorised version of existing
!                  MaxWindSpline subroutine. J-C Rioual (NEC)
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6


USE IO_Mod,  ONLY:     &
  PP_Field_type

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumCols, NumRows, NumLevs

TYPE(PP_Field_type), INTENT(IN) :: UFields(NumLevs) ! U-wind on B-grid
TYPE(PP_Field_type), INTENT(IN) :: VFields(NumLevs) ! V-wind on B-grid
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! P on wind B-grid

INTEGER, INTENT(IN) :: MaxWLev(PFields(1)%Hdr%NumCols,  &
                               PFields(1)%Hdr%NumRows)

INTEGER, PARAMETER :: ninc=16           ! number of increments used
REAL, INTENT(OUT) :: Uinc(NumCols,NumRows,2*ninc)  ! u at increment pts
REAL, INTENT(OUT) :: Vinc(NumCols,NumRows,2*ninc)  ! v at increment pts
REAL, INTENT(OUT) :: Pinc(NumCols,NumRows,2*ninc)  ! p at increment pts

! Local constants

INTEGER, PARAMETER :: kmax=5            ! total levels needed for spline
INTEGER, PARAMETER :: khalf=(kmax+1)/2  ! number of levels each side

! Local variables:
INTEGER :: i, j, k, kint
REAL, DIMENSION(NumCols,NumRows)         :: P_Lwr, P_Mid, P_Upr
REAL, DIMENSION(NumCols,NumRows,NumLevs) :: Arr_3D

! Variables involving spline interpolation
INTEGER :: spl(2*ninc)
REAL, DIMENSION(NumCols,NumRows,kmax) :: xsp, ysp, bsp, csp, dsp
REAL :: dx(NumCols,NumRows,2*ninc)


! Use levels above and below to create splines of U and V
! Use the splines to evaluate U and V at intervals between
! level-1 and level+1      OLD : level-1/2 and level+1/2

! Since MaxWLev may be 0 on places where U or V is missing
! we have to define an array which has save values there


! Store PFields into a conventional 3D array for vectorization
DO k=1,NumLevs
  Arr_3D(:,:,k) = PFields(k) % RData(:,:)
ENDDO

! P_Lwr : P on Rho Level below maximum
! P_Mid : P on Rho Level at    maximum
! P_Upr : P on Rho Level above maximum

DO j=1,NumRows
  DO i=1,NumCols
    if (maxwlev(i,j) /= 0) then
      P_Lwr(i,j) = Arr_3D(i,j,MaxWLev(i,j)-1)
      P_Mid(i,j) = Arr_3D(i,j,MaxWLev(i,j)  )
      P_Upr(i,j) = Arr_3D(i,j,MaxWLev(i,j)+1)
    end if
  END DO
END DO

DO k = 1, ninc                               ! Calc p increments
   DO j=1,NumRows
      DO i=1,NumCols
         if (maxwlev(i,j) /= 0) then
            Pinc(i,j,k)      = P_Lwr(i,j) + &
    &                   (P_Mid(i,j)-P_Lwr(i,j))*FLOAT(k)/FLOAT(ninc)
            Pinc(i,j,k+ninc) = P_Mid(i,j) + &
    &                   (P_Upr(i,j)-P_Mid(i,j))*FLOAT(k)/FLOAT(ninc)
         end if
      END DO
   END DO
END DO

!-- Pressure is x -- going down through levels so it is increasing
DO k = 1,kmax                      ! P on surrounding levels
  DO j=1,NumRows
    DO i=1,NumCols
      if (maxwlev(i,j) /= 0) then
         xsp(i,j,k) = Arr_3D(i,j,MaxWLev(i,j)+khalf-k)
      end if
    END DO
  END DO
END DO


spl (1     :  ninc) = khalf
spl (1+ninc:2*ninc) = khalf-1
!DO k=1,2*ninc
!  kint = spl(k)
!  dx(:,:,k) = Pinc(:,:,k) - xsp(:,:,kint)
!ENDDO

do k=1,2*ninc
   kint = spl(k)
   do j=1,NumRows
      do i=1,NumCols
         if (maxwlev(i,j) /= 0) then
            dx(i,j,k) = Pinc(i,j,k) - xsp(i,j,kint)
         end if
      end do
   end do
end do

!-- Store UFields into a conventional 3D array for vectorization
DO k=1,NumLevs
  Arr_3D(:,:,k) = UFields(k) % RData(:,:)
ENDDO

!-- Get U values --
DO k = 1,kmax                      ! U on surrounding levels
  DO j=1,NumRows
    DO i=1,NumCols
      if (maxwlev(i,j) /= 0) then
         ysp(i,j,k) = Arr_3D(i,j,MaxWLev(i,j)+khalf-k)
      end if
    END DO
  END DO
END DO

! DEPENDS ON: splinesetupv
CALL SplineSetupV(NumCols, NumRows, kmax, xsp, ysp, bsp, csp, &
&                  dsp, MaxWLev )
!DO k=1,2*ninc
!  kint = spl(k)
!  Uinc(:,:,k) = ysp(:,:,kint) + dx(:,:,k)*( bsp(:,:,kint) + &
!      dx(:,:,k)*( csp(:,:,kint) + dx(:,:,k)*dsp(:,:,kint)))
!ENDDO

do k=1,2*ninc
   kint = spl(k)
   do j=1,NumRows
      do i=1,NumCols
         if (maxwlev(i,j) /= 0) then
              Uinc(i,j,k) = ysp(i,j,kint) + dx(i,j,k)*( bsp(i,j,kint) + &
             &      dx(i,j,k)*(csp(i,j,kint) + dx(i,j,k)*dsp(i,j,kint)))
         end if
      end do
   end do
end do

!-- Store VFields into a conventional 3D array for vectorization
DO k=1,NumLevs
  Arr_3D(:,:,k) = VFields(k) % RData(:,:)
ENDDO

!-- Get V values --
DO k = 1,kmax                      ! V on surrounding levels
  DO j=1,NumRows
    DO i=1,NumCols
      if (maxwlev(i,j) /= 0) then
         ysp(i,j,k) = Arr_3D(i,j,MaxWLev(i,j)+khalf-k)
      end if
    END DO
  END DO
END DO

! DEPENDS ON: splinesetupv
CALL SplineSetupV(NumCols, NumRows, kmax, xsp, ysp, bsp, csp, &
&                 dsp, MaxWLev )

DO k=1,2*ninc
  kint = spl(k)
  DO j=1,NumRows
    DO i=1,NumCols
      if (maxwlev(i,j) /= 0) then
         Vinc(i,j,k) = ysp(i,j,kint) + dx(i,j,k)*( bsp(i,j,kint) + &
         dx(i,j,k)*( csp(i,j,kint) + dx(i,j,k)*dsp(i,j,kint)))
      end if
    END DO
  END DO
ENDDO

END SUBROUTINE MaxWindSplineV


!=======================================================================


!=======================================================================


!=======================================================================

#endif
