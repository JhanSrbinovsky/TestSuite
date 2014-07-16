#if defined(A05_0A) 
!
! SUBROUTINE SCNSCV2------------------------------------------------------------
!
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
!

SUBROUTINE scnscv2(row_length, rows, off_x, off_y, halo_i, halo_j, model_levels &
         , wet_model_levels, r_rho_levels, r_theta_levels, timestep             &
         , substep_number, rho, q, qcl, qcf, tracer, ccldbase, ccldtop          &
         , rainrate, snowrate, l_scav_below_cloud, k_rain, k_snow, accu_scav_tr &
         )

!-------------------------------------------------------------------------------
! Description:
!   Dummy stub of SCNSCV2 which is called if Convection Scheme is disabled.
!
!   Current code owner:  R.A. Stratton
!
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.
!
! Documentation:  UMDP 20
!
!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Arguments
! ---------
  INTEGER, INTENT(IN) :: row_length               ! Number of point in row_length
  INTEGER, INTENT(IN) :: rows                     ! Number of rows
  INTEGER, INTENT(IN) :: off_x                    ! EW size of std. halo
  INTEGER, INTENT(IN) :: off_y                    ! NS size of std. halo
  INTEGER, INTENT(IN) :: halo_i                   ! EW extended halo
  INTEGER, INTENT(IN) :: halo_j                   ! NS extended halo
  INTEGER, INTENT(IN) :: model_levels             
  INTEGER, INTENT(IN) :: wet_model_levels 
  INTEGER, INTENT(IN) :: substep_number   

  INTEGER, INTENT(IN) :: ccldbase (row_length, rows)   ! Convective cloud base
  INTEGER, INTENT(IN) :: ccldtop  (row_length, rows)   ! Convective cloud top

  REAL, INTENT(IN) :: k_rain                    ! Scavenging rate coeff for rain
  REAL, INTENT(IN) :: k_snow                    ! Scavenging rate coeff for snow
  REAL, INTENT(IN) :: timestep                         ! Timestep in secs
  REAL, INTENT(IN) :: rainrate      (row_length, rows) ! Conv. srf
                                                       ! rainrate(kg/m2)
  REAL, INTENT(IN) :: snowrate      (row_length, rows) ! Conv. srf
                                                       ! snowrate(kg/m2)
  REAL, INTENT(IN) :: accu_scav_tr  (row_length, rows) ! Column total of
                                                       ! scavenged tracer

  REAL, INTENT(IN) :: q  (row_length, rows, wet_model_levels) ! Water vap(kg/kg)
  REAL, INTENT(IN) :: qcl(row_length, rows, wet_model_levels) ! Water liq(kg/kg)
  REAL, INTENT(IN) :: qcf(row_length, rows, wet_model_levels) ! Water ice(kg/kg)


  REAL, INTENT(IN) :: r_rho_levels   (1-halo_i : row_length+halo_i,             &
                                      1-halo_j : rows      +halo_j,             &
                                                 model_levels)
                                               ! Radius to Rho levels

  REAL, INTENT(IN) :: r_theta_levels (1-halo_i : row_length+halo_i,             &
                                      1-halo_j : rows      +halo_j,             &
                                             0 : model_levels)
                                               ! Radius to Theta levels

  REAL, INTENT(IN) :: rho            (1-off_x  : row_length+off_x,              &
                                      1-off_y  : rows      +off_y,              &
                                                 model_levels)
                                               ! Density*r*r

  REAL, INTENT(IN) :: tracer         (1-off_x  : row_length+off_x,              &
                                      1-off_y  : rows      +off_y,              &
                                                 model_levels)

  LOGICAL, INTENT(IN) :: l_scav_below_cloud    ! Control for scavenging levels


! Local Variables
!-----------------

  CHARACTER(Len=52)   :: Message
  CHARACTER(Len=7 )   :: RoutineName
  INTEGER             :: ErrorStat        ! Return code:
                                          !   0 = Normal exit
                                          ! +ve = Fatal Error
                                          ! -ve = Warning

!-------------------------------------------------------------------------------
! Code Statements
!-------------------------------------------------------------------------------
        
  ErrorStat   = 1
  RoutineName = 'SCNSCV2'
  Message     = 'Convection Scheme Routines unavailable - see output.'

  WRITE (6,*) '**ERROR**: SCNSCV2 called but is unavailable.    '
  WRITE (6,*) '  Sections 5: Convection Scheme is required '


! DEPENDS ON: ereport
  CALL ereport(RoutineName, ErrorStat, Message)
        
  RETURN

!-------------------------------------------------------------------------------
END SUBROUTINE SCNSCV2
#endif
