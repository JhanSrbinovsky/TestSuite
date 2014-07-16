#if defined(A70_1B) || defined(A70_1C)
#if defined(A02_3A) || defined(A02_3C)
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
#include "dimfix3a.h"
#include "error3a.h"
#include "scfpt3a.h"
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
#if defined(VECTLIB)
      Integer :: n_in ! Number of elements for vector function
#endif
!
!
!
#if defined(VECTLIB)
      Do i=i_layer_first, i_layer_last
        Do l=1, n_profile
          trans(l, i)=-1.66e+00*tau(l, i)
        Enddo
        Do l=n_profile+1, npd_profile
          trans(l, i)=0.0e+00
        Enddo
      Enddo
      n_in=npd_profile*(i_layer_last-i_layer_first+1)
! DEPENDS ON: exp_v
      call exp_v(n_in                                                   &
     &   , trans(1, i_layer_first), trans(1, i_layer_first))
#else
      Do i=i_layer_first, i_layer_last
        Do l=1, n_profile
          trans(l, i)=exp(-1.66e+00*tau(l, i))
        Enddo
      Enddo
#endif
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
#endif
#endif
