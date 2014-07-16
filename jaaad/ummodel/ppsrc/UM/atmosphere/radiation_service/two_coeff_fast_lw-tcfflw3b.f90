

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate coefficients in the two-stream equations.
!
! Method:
!
! The basic two-stream coefficients in the differential equations
! are calculated. These are then used to determine the
! transmission and reflection coefficients. Coefficients for
! determining the solar or infra-red source terms are calculated.
!
! Current owner of code: J. M. Edwards
!
! History:
! Version  Date      Comment
! -------  ----      -------
! 5.3      04-10-01  Original Code
!                    (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      Subroutine two_coeff_fast_lw(ierr                                 &
     &   , n_profile, i_layer_first, i_layer_last                       &
     &   , l_ir_source_quad, tau                                        &
     &   , trans, source_coeff                                          &
     &   , npd_profile, npd_layer                                       &
     &   )
!
!
!
      Implicit None
!
!
!     Sizes of dummy arrays.
      Integer, Intent(IN) :: npd_profile
!                              Maximum number of profiles
      Integer, Intent(IN) :: npd_layer
!                              Maximum number of layers
!
!     Include header files.
! DIMFIX3A defines internal dimensions tied to algorithms for
! two-stream radiation code, mostly for clouds

      ! number of components of clouds
      INTEGER,PARAMETER:: NPD_CLOUD_COMPONENT=4

      ! number of permitted types of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_TYPE=4

      ! number of permitted representations of clouds.
      INTEGER,PARAMETER:: NPD_CLOUD_REPRESENTATION=4

      ! number of overlap coefficients for clouds
      INTEGER,PARAMETER:: NPD_OVERLAP_COEFF=18

      ! number of coefficients for two-stream sources
      INTEGER,PARAMETER:: NPD_SOURCE_COEFF=2

      ! number of regions in a layer
      INTEGER,PARAMETER:: NPD_REGION=3

! DIMFIX3A end
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET ERROR FLAGS IN THE RADIATION CODE.
!
      INTEGER                                                           &
     &     I_NORMAL                                                     &
!             ERROR FREE CONDITION
     &   , I_ERR_FATAL                                                  &
!             FATAL ERROR: IMMEDIATE RETURN
     &   , I_ABORT_CALCULATION                                          &
!             CALCULATION ABORTED
     &   , I_MISSING_DATA                                               &
!             MISSING DATA ERROR: CONDITIONAL
     &   , I_ERR_IO                                                     &
!             I/O ERROR
     &   , I_ERR_RANGE                                                  &
!             INTERPOLATION RANGE ERROR
     &   , I_ERR_EXIST
!             EXISTENCE ERROR
!
      PARAMETER(                                                        &
     &     I_NORMAL=0                                                   &
     &   , I_ERR_FATAL=1                                                &
     &   , I_ABORT_CALCULATION=2                                        &
     &   , I_MISSING_DATA=3                                             &
     &   , I_ERR_IO=4                                                   &
     &   , I_ERR_RANGE=5                                                &
     &   , I_ERR_EXIST=6                                                &
     &   )
!
!     ------------------------------------------------------------------
! SCFPT3A defines pointers to source coefficients in two-stream
! radiation code.

      ! pointer to source coeficient for upward solar beam
      INTEGER,PARAMETER::IP_SCF_SOLAR_UP=1

      ! pointer to source coeficient for downward solar beam
      INTEGER,PARAMETER:: IP_SCF_SOLAR_DOWN=2

      ! pointer to source coeficient for 1st difference of planckian
      INTEGER,PARAMETER:: IP_SCF_IR_1D=1

      ! pointer to source coeficient for 2nd difference of planckian
      INTEGER,PARAMETER:: IP_SCF_IR_2D=2

! SCFPT3A end
!
!
!
!     Dummy arguments.
      Integer, Intent(OUT) :: ierr
!                               Error flag
!
      Integer, Intent(IN)  :: n_profile
!                               Number of profiles used
      Integer, Intent(IN)  :: i_layer_first
!                               First layer to process
      Integer, Intent(IN)  :: i_layer_last
!                               Last layer to process
      Logical, Intent(IN)  :: l_ir_source_quad
!                               Use a quadratic source function
!
!     Optical properties of layer:
      Real, Intent(IN) :: tau(npd_profile, npd_layer)
!                               Optical depth
!
!
!     Coefficients in the two-stream equations:
      Real, Intent(OUT) :: trans(npd_profile, npd_layer)
!                            Diffuse transmission coefficient
      Real, Intent(OUT) :: source_coeff(npd_profile, npd_layer          &
     &                                 , npd_source_coeff)
!                            Source coefficients in two-stream equations
!
!
!     Local variables.
      Integer :: i    ! Loop variable
      Integer :: l    ! Loop variable
!
!
!
      Do i=i_layer_first, i_layer_last
        Do l=1, n_profile
          trans(l, i)=exp(-1.66e+00*tau(l, i))
        Enddo
      Enddo
!
      Do i=i_layer_first, i_layer_last
        Do l=1, n_profile
          source_coeff(l, i, ip_scf_ir_1d)                              &
     &      =(1.0e+00-trans(l, i)+sqrt(epsilon(trans)))                 &
     &      /(1.66e+00*tau(l, i)+sqrt(epsilon(tau)))
        Enddo
      Enddo
!
      If (l_ir_source_quad) then
        Do i=i_layer_first, i_layer_last
          Do l=1, n_profile
            source_coeff(l, i, ip_scf_ir_2d)                            &
     &        =-(1.0e+00+trans(l, i)                                    &
     &        -2.0e+00*source_coeff(l, i, ip_scf_ir_1d))                &
     &        /(1.66e+00*tau(l, i)+sqrt(epsilon(tau)))
          Enddo
        Enddo
      Endif
!
!
!
      Return
      END SUBROUTINE two_coeff_fast_lw
