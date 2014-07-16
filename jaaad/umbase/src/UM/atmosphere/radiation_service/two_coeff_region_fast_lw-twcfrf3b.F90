#if defined(A70_1B) || defined(A70_1C)
#if defined(A02_3A) || defined(A02_3C)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate two-stream coefficients in the regions.
!
! Method:
!
! The coefficients for each region are determined and averaged.
!
! Current owner of code: J. M. Edwards
!
! History:
! Version  Date      Comment
! -------  ----      -------
! 5.3      04-10-01  Original Code
!                    (J. M. Edwards)
!  6.0  21/08/03  NEC SX-6 optimisation - add vectorisation
!                 !CDIR NODEP directives.  R Barnes & J-C Rioual.
!
!
! Description of code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      Subroutine two_coeff_region_fast_lw(ierr                          &
     &   , n_profile, n_layer, n_cloud_top                              &
     &   , i_2stream, l_ir_source_quad, n_source_coeff                  &
     &   , n_cloud_type, frac_cloud                                     &
     &   , n_region, i_region_cloud, frac_region                        &
     &   , tau_free, tau_cloud                                          &
     &   , isolir                                                       &
     &   , trans, reflect, source_coeff                                 &
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
#include "spcrg3a.h"
#include "error3a.h"
#include "cldreg3a.h"
#include "stdio3a.h"
!
!
!
!     Dummy arguments.
      Integer, Intent(OUT) :: ierr
!                               Error flag
      Integer, Intent(IN)  :: n_profile
!                               Number of profiles
      Integer, Intent(IN)  :: n_layer
!                               Number of layers
      Integer, Intent(IN)  :: n_cloud_top
!                               Topmost cloudy layer
      Integer, Intent(IN)  :: isolir
!                               Spectral region
      Integer, Intent(IN)  :: n_cloud_type
!                               Number of types of clouds
      Integer, Intent(IN)  :: i_2stream
!                               Two stream scheme
      Integer, Intent(IN)  :: n_source_coeff
!                               Number of source coefficients
!
      Integer, Intent(IN)  :: n_region
!                               Number of cloudy regions
      Integer, Intent(IN)  :: i_region_cloud(npd_cloud_type)
!                               Regions in which types of clouds fall
!
      Logical, Intent(IN)  :: l_ir_source_quad
!                               Use a quadratic source in the infra-red
!
!     Optical properties of the layer:
      Real, Intent(IN) :: frac_cloud(npd_profile, npd_layer             &
     &                             , npd_cloud_type)
!                           Fractions of different types of clouds
      Real, Intent(IN) :: frac_region(npd_profile, npd_layer            &
     &                              , npd_region)
!                           Fractions of total cloud occupied
!                           by each region
      Real, Intent(IN) :: tau_free(npd_profile, npd_layer)
!                           Clear-sky optical depth
      Real, Intent(IN) :: tau_cloud(npd_profile, npd_layer              &
     &                            , npd_cloud_type)
!                           Optical depth
!
!     Coefficients in the two-stream equations:
      Real, Intent(OUT) :: trans(npd_profile, npd_layer, npd_region)
!                            Diffuse transmission coefficient
      Real, Intent(OUT) :: reflect(npd_profile, npd_layer, npd_region)
!                            Diffuse reflection coefficient
      Real, Intent(OUT) :: source_coeff(npd_profile, npd_layer          &
     &                                , npd_source_coeff, npd_region)
!                            Source coefficients in two-stream equations
!
!     Local variables.
      integer :: i        !             loop variable
      integer :: j        !             loop variable
      integer :: k        !             loop variable
      integer :: l        !             loop variable
      integer :: i_region !             loop variable over regions
!
!     Coefficients in the two-stream equations:
      Real :: trans_temp(npd_profile, npd_layer)
!               Temporary diffuse transmission coefficient
      Real :: source_coeff_temp(npd_profile, npd_layer                  &
     &                        , npd_source_coeff)
!               Temporary source coefficients in two-stream equations
!
!     Variables for gathering:
      Integer :: n_list
!                  Number of points in list
      Integer :: l_list(npd_profile)
!                  List of collected points
      Integer :: ll
!                  Loop variable
      Real :: tau_gathered(npd_profile, npd_layer)
!               Gathered optical depth
      Real :: tmp_inv(npd_profile)
!               Temporary work array
!
!     Subroutines called:
      External                                                          &
     &     two_coeff_fast_lw
!
!     Cray directives for the whole routine:
!     points are not repeated in the indexing array, so it is safe
!     to vectorize over indirectly addressed arrays.
!FPP$ NODEPCHK R
!
!
!
!     This routine should not be used outside the IR.
      If (isolir /= ip_infra_red) then
         Write(iu_err, '(/a)')                                          &
     &     '*** Erroneous use of fast non-scattering code.'
         ierr=i_err_fatal
         Return
      Endif
!
!     Determine the optical properties of the clear-sky regions of
!     the layers.
!
! DEPENDS ON: two_coeff_fast_lw
      Call two_coeff_fast_lw(ierr                                       &
     &   , n_profile, 1, n_layer                                        &
     &   , l_ir_source_quad, tau_free                                   &
     &   , trans(1, 1, ip_region_clear)                                 &
     &   , source_coeff(1, 1, 1, ip_region_clear)                       &
     &   , npd_profile, npd_layer                                       &
     &   )
      If (ierr /= i_normal) Return
      Do i=1, n_layer
        Do l=1, n_profile
          reflect(l, i, ip_region_clear)=0.0
        Enddo
      Enddo
!
!
!     Now deal with clouds.
!
!     Initialize the full arrays for cloudy regions.
!
      Do i_region=1, n_region
        If (i_region /= ip_region_clear) Then
          Do i=n_cloud_top, n_layer
            Do l=1, n_profile
              trans(l, i, i_region)=0.0
              reflect(l, i, i_region)=0.0
            Enddo
          Enddo
          Do j=1, n_source_coeff
            Do i=n_cloud_top, n_layer
              Do l=1, n_profile
                source_coeff(l, i, j, i_region)=0.0
              Enddo
            Enddo
          Enddo
!
        Endif
!
      Enddo
!
!
!
!     Consider each type of cloud in turn, checking which region it
!     contrubutes to and form weighted sums of cloud properties.
!
      Do k=1, n_cloud_type
!
!
!       Set the region in which clouds of this type are included.
        i_region=i_region_cloud(k)
!
        Do i=n_cloud_top, n_layer
!
!         Form a list of points where cloud of this type exists
!         on this row for gathering.
          n_list=0
          Do l=1, n_profile
             if (frac_cloud(l, i, k) >  0.0e+00) then
                n_list=n_list+1
                l_list(n_list)=l
             Endif
          Enddo
!
!
          If (n_list /= 0) Then
!
!           Gather the optical properties. Though we consider only
!           one layer at a time the lower routines will operate on
!           arrays with vertical structure, so the gathered arrays
!           are two-dimensional.
!
            Do l=1, n_list
              tau_gathered(l, i)                                        &
     &          =tau_cloud(l_list(l), i, k)
            Enddo
!
!
! DEPENDS ON: two_coeff_fast_lw
            Call two_coeff_fast_lw(ierr                                 &
     &        , n_list, i, i                                            &
     &        , l_ir_source_quad, tau_gathered                          &
     &        , trans_temp                                              &
     &        , source_coeff_temp                                       &
     &        , npd_profile, npd_layer                                  &
     &        )
            if (ierr /= i_normal) return
!
!
!CDIR NODEP
            Do l=1, n_list
              ll=l_list(l)
              trans(ll, i, i_region)=trans(ll, i, i_region)             &
     &           +frac_cloud(ll, i, k)*trans_temp(l, i)
            Enddo
            Do j=1, n_source_coeff
!CDIR NODEP
              Do l=1, n_list
                ll=l_list(l)
                source_coeff(ll, i, j, i_region)                        &
     &            =source_coeff(ll, i, j, i_region)                     &
     &            +frac_cloud(ll, i, k)                                 &
     &            *source_coeff_temp(l, i, j)
              Enddo
            Enddo
!
          Endif
!
        Enddo
      Enddo
!
!
!     Finally, scale the weighted sums by the cloud fractions.
      Do i_region=1, n_region
        If (i_region /= ip_region_clear) Then
          Do i=n_cloud_top, n_layer
!
!           Gather points within this region.
            n_list=0
            Do l=1,n_profile
              If (frac_region(l, i, i_region) > 0.0) Then
                n_list=n_list+1
                l_list(n_list)=l
              Endif
            Enddo
!CDIR NODEP
            Do l=1, n_list
              ll=l_list(l)
              tmp_inv(l)=1.0e+00/frac_region(ll, i, i_region)
              trans(ll, i, i_region)=trans(ll, i, i_region)             &
     &          *tmp_inv(l)
            Enddo
            Do j=1, n_source_coeff
!CDIR NODEP
              Do l=1, n_list
                ll=l_list(l)
                source_coeff(ll, i, j, i_region)                        &
     &            =source_coeff(ll, i, j, i_region)                     &
     &            *tmp_inv(l)
              Enddo
            Enddo
          Enddo
        Endif
      Enddo
!
!
!
      Return
      END SUBROUTINE two_coeff_region_fast_lw
#endif
#endif
