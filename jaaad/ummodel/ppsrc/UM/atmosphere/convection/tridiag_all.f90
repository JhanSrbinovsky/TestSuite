
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  solves tridiagonal matrix -  convection scheme
!
      SUBROUTINE TRIDIAG_ALL(N,nvec,nvec_len,A,B,C,R,U)
!
! Purpose: Solves the equations A.X = Y,  where A is a tridiagnol matrix
!
!          for several matrices A at once.
!
!
!   Called by shconv_grad_h
!
! Current owners of code: _Convection code owner
!
! History:
! Version     Date     Comment
! -------     ----     -------
!
!   6.2    22/02/05   New verion of code to work on many columns
!                      R A Stratton
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

      IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!
      integer, intent(in) ::                                            &
     &   N                                                              &
                    ! maximum size of vectors X and Y
     &,  nvec                                                           &
                   ! Number of vectors/matrices to solve
     &,  nvec_len(nvec)  ! length of each vector.

      real, intent(in) ::                                               &
     &    A(nvec,N)                                                     &
                       ! Components of tridiagonal matrix
     &,   B(nvec,N)                                                     &
                       ! Components of tridiagonal matrix
     &,   C(nvec,N)                                                     &
                       ! Components of tridiagonal matrix
     &,   R(nvec,N)    ! vector Y, R.H.S of linear equation

      real, intent(out) ::                                              &
     &    U(nvec,N)    ! solution vectors


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
!
      integer ::                                                        &
     &  i,j            ! loop counter

      real ::                                                           &
     &  gam(nvec,n)                                                     &
                      ! work array
     &, bet(nvec)

!
! External routines called:  NONE
!

!-----------------------------------------------------------------------

      Do i=1,nvec
        bet(i) = B(i,1)
        U(i,1) = R(i,1)/bet(i)
      End do

      Do J=2,N
        Do i=1, nvec
          If (j <= nvec_len(i)) then
            gam(i,j) = C(i,j-1)/bet(i)
            bet(i)   = B(i,j) - A(i,j)*gam(i,j)
!for info !    if (bet(i) == 0.0) stop   ! was in original code
            u(i,j)   = (R(i,j) - A(i,j)*U(i,j-1))/bet(i)
          End if
        End do
      End do

      Do J=N-1,1,-1
        Do i=1, nvec
          If (j <= (nvec_len(i)-1)) then
            U(i,j) = U(i,j) - Gam(i,j+1)*U(i,j+1)
          Endif
        End do
      End do

!-----------------------------------------------------------------------

      RETURN
      END SUBROUTINE TRIDIAG_ALL
