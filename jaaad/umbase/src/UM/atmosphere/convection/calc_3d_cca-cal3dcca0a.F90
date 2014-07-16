#if defined(A05_0A)
!
! SUBROUTINE CALC_3D_CCA------------------------------------------------------
!
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! ******************************COPYRIGHT*************************************
!

subroutine calc_3d_cca(np_field, npnts, nlev, nbl, cld_base, cld_top,         &
           p_lyr_bnds, frz_lev, cca_2d, cca_3d, z_theta, z_rho, l_q_interact, &
           l_use_sh_mask, l_pc2_diag_sh_pts  )

!-----------------------------------------------------------------------------
! Description:
!   Dummy stub of CALC_3D_CCA which is called if Convection Scheme is disabled
!
!   Current code owner:  R.A. Stratton
!
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.
!
!-----------------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------------
!   Scalar arguments with intent(in):
!-----------------------------------------------------------------------------
  INTEGER, intent(in) :: npnts         ! Number of points
  INTEGER, intent(in) :: np_field      ! Full data length
  INTEGER, intent(in) :: nlev          ! Number of levels
  INTEGER, intent(in) :: nbl           ! Number of Boundary layer levels

  LOGICAL, intent(in) :: l_q_interact  ! .TRUE. : PC2 cloud scheme is in use

  LOGICAL, intent(in) :: l_pc2_diag_sh_pts(npnts) 
                                       ! .TRUE. : If replacing pc2 cloud with
                                       !          diagnostic shallow Cloud

  LOGICAL, intent(in) :: l_use_sh_mask ! Logical to indicate whether or not
                                       ! to use l_pc2_diag_sh_pts as a test
                                       ! criteria


!-----------------------------------------------------------------------------
!   Array  arguments with intent(in):
!-----------------------------------------------------------------------------

  REAL,    intent(in) :: z_theta    (np_field,   nlev) ! z (th layer centres)
  REAL,    intent(in) :: z_rho      (np_field,   nlev) ! z (rh level  bounds)
  REAL,    intent(in) :: p_lyr_bnds (np_field, 0:nlev)
                                                 ! Pressure on layer
                                                 ! boundaries (rho levels-1)

  REAL,    intent(in) :: cca_2d   (npnts) ! 2D convective cloud amount
  INTEGER, intent(in) :: cld_top  (npnts) ! Conv. cloud top  (theta level)
  INTEGER, intent(in) :: cld_base (npnts) ! Conv. cloud base (theta level)
  INTEGER, intent(in) :: frz_lev  (npnts) ! Freezing level   (theta level)



!-----------------------------------------------------------------------------
!   Array  arguments with intent(out):
!-----------------------------------------------------------------------------

  REAL,    intent(in) :: cca_3d(np_field, nlev) ! Convective cloud amount
                                                ! (Theta levels)





! Local Variables
!-----------------

  CHARACTER(Len=52)   :: Message
  CHARACTER(Len=11)   :: RoutineName
  INTEGER             :: ErrorStat            ! Return code:
                                              !   0 = Normal exit
                                              ! +ve = Fatal Error
                                              ! -ve = Warning

!-----------------------------------------------------------------------------
! Code Statements
!-----------------------------------------------------------------------------
        
  ErrorStat   = 1
  RoutineName = 'CALC_3D_CCA'
  Message     = 'Convection Scheme Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: CALC_3D_CCA called but is unavailable.'
  WRITE (6,*) '  Sections 5: Convection Scheme is required '


! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)

  RETURN

!-----------------------------------------------------------------------------
END SUBROUTINE calc_3d_cca
#endif
