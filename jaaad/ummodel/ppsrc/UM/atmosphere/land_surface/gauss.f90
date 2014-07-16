
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE GAUSS--------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE GAUSS (NLEVS,NPNTS,SOIL_PTS,SOIL_INDEX,A,B,C,D         &
     &,                 XMIN,XMAX,X)

      IMPLICIT NONE
!
! Description:
!     Solves a tridiagnonal matrix equation of the form:
!
!             A(n) X(n-1) + B(n) X(n) + C(n) X(n+1) = D(n)
!
!     by Gausian elimination, assuming boundary conditions:
!
!             A(1) = 0.0    at the top
!             C(N) = 0.0    at the bottom.
!                                                          (Cox, 2/99)
!
!
! Documentation :
!
! Current Code Owner : Peter Cox
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  5.2    15/11/00   New Deck         M. Best
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!

! CONSTANTS:
      CHARACTER (LEN=17), PARAMETER :: SubName = "GAUSS (hydrology)"   

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER                                                           &
     & NLEVS                                                            &
                            ! IN Number of levels.
     &,NPNTS                                                            &
                            ! IN Number of gridpoints.
     &,SOIL_PTS             ! IN Number of soil points.


!   Array arguments with intent(IN) :
      INTEGER                                                           &
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL                                                              &
     & A(NPNTS,NLEVS)                                                   &
                            ! IN Matrix elements corresponding
!                           !    to the coefficients of X(n-1).
     &,B(NPNTS,NLEVS)                                                   &
                            ! IN Matrix elements corresponding
!                           !    to the coefficients of X(n).
     &,C(NPNTS,NLEVS)                                                   &
                            ! IN Matrix elements corresponding
!                           !    to the coefficients of X(n+1).
     &,D(NPNTS,NLEVS)                                                   &
                            ! IN Matrix elements corresponding
!                           !    to the RHS of the equation.
     &,XMIN(NPNTS,NLEVS)                                                &
                            ! IN Minimum permitted value of X.
     &,XMAX(NPNTS,NLEVS)    ! IN Maximum permitted value of X.

!   Array arguments with intent(OUT) :
      REAL                                                              &
     & X(NPNTS,NLEVS)       ! OUT Solution.

! Local scalars:
      INTEGER                                                           &
     & I,J,N                ! WORK Loop counters.
     
      INTEGER :: ErrCount, ErrCode   !error reporting
      CHARACTER (LEN=80) :: ErrMsg

! Local arrays:
      REAL                                                              &
     & ADASH(NPNTS,NLEVS),BDASH(NPNTS,NLEVS),CDASH(NPNTS,NLEVS)         &
     &,DDASH(NPNTS,NLEVS)   ! WORK Transformed matrix elements

!-----------------------------------------------------------------------
! By default set the implicit increment to the explicit increment
! (for when denominators vanish).
!-----------------------------------------------------------------------
      DO N=1,NLEVS
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          X(I,N)=D(I,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Upward Sweep: eliminate "C" elements by replacing nth equation with:
!                  B'(n+1)*Eq(n)-C(n)*Eq'(n+1)
! where "'" denotes a previously tranformed equation. The resulting
! equations take the form:
!                A'(n) X(n-1) + B'(n) X(n) = D'(n)
! (NB. The bottom boundary condition implies that the NLEV equation does
!  not need transforming.)
!-----------------------------------------------------------------------
!CDIR NODEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        ADASH(I,NLEVS)=A(I,NLEVS)
        BDASH(I,NLEVS)=B(I,NLEVS)
        DDASH(I,NLEVS)=D(I,NLEVS)
      ENDDO

      DO N=NLEVS-1,1,-1
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          ADASH(I,N)=BDASH(I,N+1)*A(I,N)
          BDASH(I,N)=BDASH(I,N+1)*B(I,N)-C(I,N)*ADASH(I,N+1)
          DDASH(I,N)=BDASH(I,N+1)*D(I,N)-C(I,N)*DDASH(I,N+1)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Top boundary condition: A(1) = 0.0 , allows X(1) to be diagnosed
!-----------------------------------------------------------------------
      ErrCount=0
!CDIR NODEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        IF (BDASH(I,1) /= 0.0) THEN
          X(I,1)=DDASH(I,1)/BDASH(I,1)
        ELSE
          ErrCount = ErrCount + 1   ! occurence count
        ENDIF
        X(I,1)=MAX(X(I,1),XMIN(I,1))
        X(I,1)=MIN(X(I,1),XMAX(I,1))
      ENDDO
      
      !output any error(s)
      IF (ErrCount > 0) THEN
        ErrCode = -1   ! probs with first layer
        WRITE (ErrMsg,*) ErrCount, ' divide by zero warnings at layer 1'
! DEPENDS ON: ereport
        CALL ereport( SubName, ErrCode, ErrMsg )
      END IF

!-----------------------------------------------------------------------
! Downward Sweep: calculate X(n) from X(n-1):
!                X(n) = (D'(n) - A'(n) X(n-1)) / B'(n)
!-----------------------------------------------------------------------
      DO N=2,NLEVS
        ErrCount=0
!CDIR NODEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          IF (BDASH(I,N) /= 0.0) THEN
            X(I,N)=(DDASH(I,N)-ADASH(I,N)*X(I,N-1))/BDASH(I,N)
          ELSE
            ErrCount = ErrCount + 1
          ENDIF
          X(I,N)=MAX(X(I,N),XMIN(I,N))
          X(I,N)=MIN(X(I,N),XMAX(I,N))
        ENDDO
        
        !output any error(s)
        IF (ErrCount > 0) THEN
          ErrCode = -2   ! probs with other layers
          WRITE (ErrMsg,*) ErrCount, ' divide by zero warnings at layer', N
! DEPENDS ON: ereport
          CALL ereport( SubName, ErrCode, ErrMsg )
        END IF

      ENDDO

      RETURN
      END SUBROUTINE GAUSS
