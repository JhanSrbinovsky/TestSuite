#if defined(A34_1A) || defined(A34_1C) || defined(A34_1G)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to define fields from D1 needed for UKCA
!  these are set in UKCA_SETD1DEFS
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Current Code Owner:       Colin Johnson/Olaf Morgenstern
!                           Fiona O'Connor
!
!  Code Description:
!    Language:  FORTRAN 90 (formatted)
!
! ######################################################################
!
      MODULE UKCA_D1_DEFS
      IMPLICIT NONE

!     No of Prognostics and Diagnostics required

      INTEGER, SAVE :: Nukca_D1items      ! Size of UkcaD1codes array
      INTEGER, SAVE :: n_all_tracers=150  ! max number of tracers
      INTEGER, SAVE :: n_use_tracers      ! no of tracers used
      INTEGER, SAVE :: n_use_emissions    ! no of emissions used
      INTEGER, SAVE :: n_3d_emissions     ! no of 3-D emissions used
      INTEGER, SAVE :: n_in_progs         ! No of prognostics required
      INTEGER, SAVE :: n_in_diags0        ! No of diagnostics (sect 0)
      INTEGER, SAVE :: n_in_diags1        ! No of diagnostics (sect 1)
      INTEGER, SAVE :: n_in_diags3        ! No of diagnostics (sect 3)
      INTEGER, SAVE :: n_in_diags4        ! No of diagnostics (sect 4)
      INTEGER, SAVE :: n_in_diags5        ! No of diagnostics (sect 5)
      INTEGER, SAVE :: n_in_diags8        ! No of diagnostics (sect 8)
      INTEGER, SAVE :: n_in_diags15       ! No of diagnostics (sect 15)
      INTEGER, SAVE :: n_in_diags30       ! No of diagnostics (sect 30)
      INTEGER, SAVE :: n_in_diags33       ! No of diagnostics (sect 33)
      INTEGER, SAVE :: n_out_diags        ! Not used
      INTEGER, SAVE :: n_emiss_first      ! Position of emissions in section
      INTEGER, SAVE :: n_emiss_last       ! Position of emissions in section
      INTEGER, SAVE :: idiag_first        ! Position of first chem diag
      INTEGER, SAVE :: idiag_last         ! Position of last chem diag
      INTEGER, SAVE :: iflux_first        ! item for 1st flux diag
      INTEGER, SAVE :: iflux_last         ! item for last flux diag
      INTEGER, SAVE :: istrat_first       ! item for 1st strat flux diag
      INTEGER, SAVE :: istrat_last        ! item for last strat flux diag
      INTEGER, SAVE :: UKCA_sect=34       ! stash section for UKCA
      INTEGER, SAVE :: UKCA_ems_sect=0    ! stash section for UKCA ems
      INTEGER, SAVE :: MODE_diag_sect=33  ! stash section for MODE diags

      TYPE CODE
       INTEGER :: section     ! section code
       INTEGER :: item        ! item code
       INTEGER :: n_levels    ! number of levels
       INTEGER :: address     ! address in D1
       INTEGER :: length      ! length of field
       INTEGER :: halo_type   ! field halo type
       INTEGER :: grid_type   ! grid type
       INTEGER :: field_type  ! field grid type
       INTEGER :: len_dim1    ! length of array dim1
       INTEGER :: len_dim2    ! length of array dim2
       INTEGER :: len_dim3    ! length of array dim3
       LOGICAL :: prognostic  ! prognostic t/f
       LOGICAL :: required    ! t/f
      ENDTYPE CODE

!     Number of tracers, emissions and diagnostics, set according
!     to UMUI choices in UKCA_SETD1DEFS

      INTEGER, SAVE   :: n_chem_tracers     ! No. tracers for chemistry
      INTEGER, SAVE   :: n_aero_tracers     ! No. tracers for chemistry
      INTEGER, SAVE   :: n_chem_emissions   ! "   emissions        "
      INTEGER, SAVE   :: n_chem_diags       ! "   diagnostics      "
      INTEGER, SAVE   :: n_BE_fluxdiags     ! "   BE flux diagnostics "
      INTEGER, SAVE   :: nmax_strat_fluxdiags  ! Max no strat flux diags
      INTEGER, SAVE   :: n_strat_fluxdiags  ! No. strat flux diags
      INTEGER, SAVE   :: n_MODE_tracers     ! No. tracers for MODE
      INTEGER, SAVE   :: n_MODE_emissions   ! "   emissions   "
      INTEGER, SAVE   :: n_MODE_diags       ! "   diagnostics "
      INTEGER, SAVE   :: n_dust_tracers     ! No. tracers for dust
      INTEGER, SAVE   :: n_dust_emissions   ! "   emissions   "
      INTEGER, SAVE   :: n_dust_diags       ! "   diagnostics "
      INTEGER, SAVE   :: n_RnPb_tracers     ! No. tracers for Rn-Pb
      INTEGER, SAVE   :: n_RnPb_emissions   ! "   emissions   "
      INTEGER, SAVE   :: n_RnPb_diags       ! "   diagnostics "
      INTEGER, SAVE   :: nr_therm           ! "   thermal reactions
      INTEGER, SAVE   :: nr_phot            ! "   photolytic reactions

!     list of prognostic/diagnostics from/to D1

      TYPE(CODE), DIMENSION(:), ALLOCATABLE, SAVE :: UkcaD1Codes

!     Names for tracers which have surface emissions

      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_chem_spec
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_MODE_spec
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_dust_spec
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_RnPb_spec

!     Names for defined tracers

      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: nm_spec

!     Index arrays for tracers and emissions

      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE           :: tr_index
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE           :: em_index

!     UKCA Logicals - will eventually be replaced with umui/namelist

      LOGICAL, SAVE :: L_ukca_BEflux    = .false.  ! T for BE chem fluxes
      LOGICAL, SAVE :: L_ukca_stratflux = .false.  ! T for stratospheric fluxes
      LOGICAL, SAVE :: L_ukca_Tbias     = .false.  ! T for temp bias correction
      LOGICAL, SAVE :: L_ukca_Qbias     = .false.  ! T for sp.hum. bias correction

      END MODULE UKCA_D1_DEFS

#endif
