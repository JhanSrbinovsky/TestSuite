#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
!
! (c) [University of Leeds] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Solve the coagulation/nucleation equation.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Current Code Owner:       Graham Mann
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_SOLVECOAGNUCL_V(NV,MASK,A,B,C,ND,DTZ,DELN)
!---------------------------------------------------------------
!
! Purpose
! -------
! Solve the coagulation/nucleation equation which is of the form
! dN/dt=A*N^2 + B*N + C
!
! Note that the A term refers to intra-modal coagulation
!           the B term refers to inter-modal coagulation
!           the C term refers to new particle formation
!
! Rearrange this to be in the form of an indefinite integral
!
! int_{N_0}^{N} dx/X = int_{t_0}^{t_0+deltat} dt = deltat
!
! where X= A*x^2 + B*x + C
!
! Use analytical solutions to solve this indefinite integral for
! different cases below following pg 28 of Bronshtein & Semendyayev,
! "Handbook of Mathematics" (1997, 3rd edition) Equation 40.
!
! D = -discriminant = 4*A*C - B^2
!
! 1) [LOGIC1A is true] : A /= 0, and D < 0
!
! int dx/X = {2/sqrt(D)} * arctan{ (2*A*x+B)/sqrt(D) }
!
! 2) [LOGIC1B is true] : A /= 0, and D > 0
!
! int dx/X = {1/sqrt(-D)}*log{ (2*A*x+B-sqrt(-D)) / (2*A*x+B+sqrt(-D) }
!
! 3) [LOGIC1CA is true] : A /= 0, D = 0 and B/=0 (causes error message)
!
! In this case, B^2 = 4AC (B=/0) and there is no solution given.
!
! 4) [LOGIC1CB is true] : A /= 0, D = 0 and B=0
!
! In this case, B^2 = 4AC, and B=0, so C must also = 0
! in which case we have dN/dt = A*N^2 which has solution
!
!       N=1/(1/N0-3*A*deltat)
!
! 5) [LOGIC2A is true]  : A == 0 and B /= 0
!
!    In this case dN/dt = B*N + C which has the solution
!
!       N={(B*N0+C)*exp(B*deltat)-C}/B
!
!     n.b. if no nucleation (C=0) then N=N0*exp{B*deltat}
!
! 6) [LOGIC2B is true]  : A == 0 and B == 0
!
!    In this case dN/dt = C and N=N0+C*deltat.
!
! Parameters
! ----------
! None
!
! Inputs:
! ------
! A      : constant in coag/nucl equation
! B      : constant in coag/nucl equation
! C      : constant in coag/nucl equation
! ND     : initial aerosol particle no. conc. (ptcls/cc)
! DTZ    : nucl/coag/cond time step (s)
!
! Outputs:
! -------
! DELN   : change in ND due to nucl & coag (ptcls/cc)
!
!
! Local variables:
! ---------------
! D      : 4*A*C - B^2
! SQD    : SQRT(D) if D>0 or SQRT(-D) if D<0
! NDNEW  : Number concentration after combined coag-nucl.
! LOGIC1,LOGIC2,LOGIC3,etc. : Logicals for various cases (see above)
! IERR   : stores which boxes had no solution (should not occur).
!
!--------------------------------------------------------------------
      IMPLICIT NONE
!
!     Subroutine interface
      INTEGER :: NV
      REAL    :: A(NV)
      REAL    :: B(NV)
      REAL    :: C(NV)
      REAL    :: ND(NV)
      REAL    :: DELN(NV)
      REAL    :: DTZ
      LOGICAL :: MASK(NV)
!
!     Local variables
      INTEGER :: I
      INTEGER :: IERR(NV)
      LOGICAL :: LOGIC1(NV)
      LOGICAL :: LOGIC2(NV)
      LOGICAL :: LOGIC3(NV)
      LOGICAL :: LOGIC1A(NV)
      LOGICAL :: LOGIC1B(NV)
      LOGICAL :: LOGIC1C(NV)
      LOGICAL :: LOGIC1CA(NV)
      LOGICAL :: LOGIC1CB(NV)
      LOGICAL :: LOGIC2A(NV)
      LOGICAL :: LOGIC2B(NV)
      REAL    :: NDNEW(NV)
      REAL    :: D(NV)
      REAL    :: SQD(NV)
      REAL    :: TERM1(NV)
      REAL    :: TERM2(NV)
      REAL    :: TERM3(NV)
      REAL    :: TERM4(NV)

      IERR(:) = 0

      NDNEW(:)=0.
      D(:)=4.0*A(:)*C(:)-B(:)*B(:)

      LOGIC1(:)=MASK(:).AND.(A(:) /= 0.0)
      LOGIC2(:)=MASK(:).AND.(A(:) == 0.0)
      LOGIC3(:)=B(:) == 0.0

      LOGIC1A(:)=LOGIC1(:).AND.(D(:) < 0.0)
      LOGIC1B(:)=LOGIC1(:).AND.(D(:) > 0.0)
      LOGIC1C(:)=LOGIC1(:).AND.(D(:) == 0.0)

      LOGIC1CA(:)=LOGIC1C(:).AND.(.NOT.LOGIC3(:))
      LOGIC1CB(:)=LOGIC1C(:).AND.(LOGIC3(:))

      LOGIC2A(:)=LOGIC2(:).AND.(.NOT.LOGIC3(:))
      LOGIC2B(:)=LOGIC2(:).AND.LOGIC3(:)

! Below is for case where A /= 0 & D < 0
      WHERE(LOGIC1A(:))
       SQD(:)=SQRT(-D(:))
       TERM1(:)=2.0*A(:)*ND(:)+B(:)+SQD(:)
       TERM2(:)=2.0*A(:)*ND(:)+B(:)-SQD(:)
       TERM3(:)=TERM1(:)/TERM2(:)
       TERM4(:)=1.0-EXP(SQD(:)*DTZ)/TERM3(:)
       NDNEW(:)=(2.0*SQD(:)/TERM4(:)-B(:)-SQD(:))/2.0/A(:)
      ENDWHERE
!
! Below is for case where A /= 0 & D > 0
      WHERE(LOGIC1B(:))
       SQD(:)=SQRT(D(:))
       TERM1(:)=(2.0*A(:)*ND(:)+B(:))/SQD(:)
       TERM2(:)=ATAN(TERM1(:))+SQD(:)*DTZ/2.0
       NDNEW(:)=(SQD(:)*TAN(TERM2(:))-B(:))/2.0/A(:)
      ENDWHERE
!
! Below is for case where A /= 0, D = 0 and B/=0 (error)
      WHERE(LOGIC1CA(:)) IERR(:)=1 ! ! A /= 0, D=0 & B=0
!
! Below is for case where A /= 0, D = 0 and B=0
      WHERE(LOGIC1CB(:)) NDNEW(:)=1.0/(1.0/ND(:)-3.0*A(:)*DTZ)

! Below is for case where A == 0, B /= 0 (==> D<0 )
      WHERE(LOGIC2A(:))
       TERM1(:)=EXP(B(:)*DTZ)
       TERM2(:)=B(:)*ND(:)+C(:)
       NDNEW(:)=(TERM1(:)*TERM2(:)-C(:))/B(:)
      ENDWHERE

! Below is for case where A == 0, B == 0 (==> D=0 )
      WHERE(LOGIC2B(:)) NDNEW(:)=ND(:)+C(:)*DTZ

! Calculate increment to ND due to combined coag & nucl.
      DELN(:)=NDNEW(:)-ND(:)

      IF(SUM(IERR(:)) > 0) THEN
       DO I=1,NV
        IF(IERR(I) == 0) CYCLE
        D=4.0*A(I)*C(I)-B(I)*B(I)
        write(6,*) 'D=0, A not = 0, B not = 0, C not = 0'
        write(6,*) 'so B^2=4AC>0',A(I),B(I),C(I),D
! Note that if C=0 & D=0 then B=0
       ENDDO
       pause
      ENDIF

      RETURN
      END SUBROUTINE UKCA_SOLVECOAGNUCL_V
#endif
