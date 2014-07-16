#if defined(A70_1Z)
#if defined(A01_3Z) ||  defined(A02_3Z)
! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************
!
!+ subroutine to set the mixing ratios of gases.
!
! Purpose:
!   the full array of mass mixing ratios of gases is filled.
!
! Method:
!   the arrays of supplied mixing ratios are inverted and fed
!   into the array to pass to the radiation code. for well-mixed
!   gases the constant mixing ratios are fed into this array.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- --------------------------------------------------------------------
      subroutine r2_set_gas_mix_ratio(ierr, i_call                      &
     &   , n_profile, nlevs, n_layer, nwet, nozone                      &
     &   , i_gather, l_extra_top                                        &
     &   , n_absorb, type_absorb                                        &
     &   , l_n2o, l_ch4, l_cfc11, l_cfc12, l_o2                         &
     &   , l_cfc113, l_hcfc22, l_hfc125, l_hfc134a                      &
     &   , h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio                   &
     &   , c11_mix_ratio, c12_mix_ratio, o2_mix_ratio                   &
     &   , c113_mix_ratio, hcfc22_mix_ratio                             &
     &   , hfc125_mix_ratio, hfc134a_mix_ratio                          &
     &   , gas_mix_ratio                                                &
     &   , co2_dim1, co2_dim2, co2_3d, l_co2_3d                         &
     &   , nd_field, nd_profile, nd_layer, nd_species                   &
     &   )
!
!
!     comdecks included
#include "gas_list_pcf3z.h"
#include "def_std_io_icf3z.h"
#include "error_pcf3z.h"
!
!     dummy arguments.
!
      integer                                                           &
                !, intent(out)
     &     ierr
!             error flag
!
!     sizes of arrays:
      integer                                                           &
                !, intent(in)
     &     nd_field                                                     &
!             size of array from um
     &   , nd_profile                                                   &
!             size of array
     &   , nd_layer                                                     &
!             size of array
     &   , nd_species
!             size of array
!
!     sizes used:
!       5.1             04-04-00                tolerances replaced by
!                                               f90 intrinsics.
!                                               (J.-C. Thelen)
      integer                                                           &
                !, intent(in)
     &     n_profile                                                    &
!             number of profiles
     &   , nlevs                                                        &
!             number of layers in the main model
     &   , n_layer                                                      &
!             number of radiative layers
     &   , nwet                                                         &
!             number of wet levels
     &   , nozone
!             number of ozone levels
!
      integer :: i_call
!         Number of call

!     gathering array:
      integer                                                           &
                !, intent(in)
     &     i_gather(nd_field)
!             list of points to be gathered
      logical                                                           &
                !, intent(in)
     &     l_extra_top
!             flag to use an extra top layer in radiative calculations
!
!     types of gases:
      integer                                                           &
                !, intent(in)
     &     n_absorb                                                     &
!             number of absorbers
     &   , type_absorb(nd_species)
!             types of absorbers
!
!     flags for minor gases:
      logical                                                           &
                !,intent(in)
     &     l_n2o                                                        &
!             flag for nitrous oxide
     &   , l_ch4                                                        &
!             flag for methane
     &   , l_cfc11                                                      &
!             flag for cfc11
     &   , l_cfc12                                                      &
!             flag for cfc12
     &   , l_o2                                                         &
!             flag for o2
     &   , l_cfc113                                                     &
!             flag for cfc113
     &   , l_hcfc22                                                     &
!             flag for hcfc22
     &   , l_hfc125                                                     &
!             flag for hfc125
     &   , l_hfc134a
!             flag for hfc134a
!
!     mixing ratios supplied:
      integer  co2_dim1, co2_dim2   ! dimensions of co2_3d field
      logical  l_co2_3d    !  controls use of 3d co2 field
      real                                                              &
                !, intent(in)
     &     h2o(nd_field, nwet)                                          &
!             mass mixing ratio of water vapour
     &   , co2                                                          &
!             mass mixing ratio of carbon dioxide
     &   , co2_3d(co2_dim1, co2_dim2)                                   &
!             3d mass mixing ratio of co2 (full field)
     &   , o3(nd_field, nozone)                                         &
!             mass mixing ratio of ozone
     &   , n2o_mix_ratio                                                &
!             mass mixing ratio of nitrous oxide
     &   , ch4_mix_ratio                                                &
!             mass mixing ratio of methane
     &   , c11_mix_ratio                                                &
!             mass mixing ratio of cfc11
     &   , c12_mix_ratio                                                &
!             mass mixing ratio of cfc12
     &   , o2_mix_ratio                                                 &
!             mass mixing ratio of o2
     &   , c113_mix_ratio                                               &
!             mass mixing ratio of cfc113
     &   , hcfc22_mix_ratio                                             &
!             mass mixing ratio of hcfc22
     &   , hfc125_mix_ratio                                             &
!             mass mixing ratio of hfc125
     &   , hfc134a_mix_ratio
!             mass mixing ratio of hfc134a
!
!     array of assigned mxing ratios:
      real                                                              &
                !, intent(out)
     &     gas_mix_ratio(nd_profile, nd_layer, nd_species)
!             mixing ratios
!
!     local variables.
!
!     pointers to gases:
      integer                                                           &
     &     iump_h2o                                                     &
!             pointer to water vapour
     &   , iump_co2                                                     &
!             pointer to carbon dioxide
     &   , iump_o3                                                      &
!             pointer to ozone
     &   , iump_n2o                                                     &
!             pointer to nitous oxide
     &   , iump_ch4                                                     &
!             pointer to methane
     &   , iump_cfc11                                                   &
!             pointer to cfc11
     &   , iump_cfc12                                                   &
!             pointer to cfc12
     &   , iump_o2                                                      &
!             pointer to o2
     &   , iump_cfc113                                                  &
!             pointer to cfc113
     &   , iump_hcfc22                                                  &
!             pointer to hcfc22
     &   , iump_hfc125                                                  &
!             pointer to hfc125
     &   , iump_hfc134a
!             pointer to hfc134a
      integer                                                           &
     &     i                                                            &
!             loop variable
     &   , l                                                            &
!             loop variable
     &   , lg                                                           &
!             corresponding ungathered index
     &   , i_top_copy
!             topmost layer where properties are set by copying the
!             input fields.
!
!
!
!
!     match the indexing numbers of gaseous species in the spectral
!     file with actual types of gases known to the um.
!
!     set all pointers to 0 initially to flag missing gases.
      iump_h2o=0
      iump_co2=0
      iump_o3=0
      iump_n2o=0
      iump_ch4=0
      iump_cfc11=0
      iump_cfc12=0
      iump_o2=0
      iump_cfc113=0
      iump_hcfc22=0
      iump_hfc125=0
      iump_hfc134a=0
!
!
      do i=1, n_absorb
!
         if (type_absorb(i) == ip_h2o) then
            iump_h2o=i
         else if (type_absorb(i) == ip_co2) then
            iump_co2=i
         else if (type_absorb(i) == ip_o3) then
            iump_o3=i
         else if (type_absorb(i) == ip_n2o) then
            iump_n2o=i
         else if (type_absorb(i) == ip_ch4) then
            iump_ch4=i
         else if (type_absorb(i) == ip_cfc11) then
            iump_cfc11=i
         else if (type_absorb(i) == ip_cfc12) then
            iump_cfc12=i
         else if (type_absorb(i) == ip_o2) then
            iump_o2=i
         else if (type_absorb(i) == ip_cfc113) then
            iump_cfc113=i
         else if (type_absorb(i) == ip_hcfc22) then
            iump_hcfc22=i
         else if (type_absorb(i) == ip_hfc125) then
            iump_hfc125=i
         else if (type_absorb(i) == ip_hfc134a) then
            iump_hfc134a=i
         endif
!
      enddo
!
!
      if (l_extra_top) then
!       the second radiative layer will be the first to have properties
!       set by copying input fields.
        i_top_copy=2
      else
!       the first radiative layer will be the first to have properties
!       set by copying input fields.
        i_top_copy=1
      endif
!
!
!     assign mixing ratios of the gases to the main arrays.
!
!     water vapour:
!
      if (iump_h2o >  0) then
!        no water exists above the wet levels.
         do i=1, n_layer-nwet
            do l=1, n_profile
               gas_mix_ratio(l, i, iump_h2o)=0.0e+00
            enddo
         enddo
         do i=n_layer-nwet+1, n_layer
            do l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_h2o)=h2o(lg, n_layer-i+1)
            enddo
         enddo
      else if (i_call == 1) then
         write(iu_err, '(/a)')                                          &
     &      '*** error: water vapour is not in the spectral file.'
         ierr=i_err_fatal
         return
      endif
!
!     carbon dioxide:
!
      if (iump_co2 >  0) then
         do i=1, n_layer
           if (l_co2_3d) then
             do l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_co2)=co2_3d(lg, n_layer-i+1)
             enddo
           else
             do l=1, n_profile
               gas_mix_ratio(l, i, iump_co2)=co2
             enddo
           endif
         enddo
      else if (i_call == 1) then
         write(iu_err, '(/a)')                                          &
     &      '*** error: carbon dioxide is not in the spectral file.'
         ierr=i_err_fatal
         return
      endif
!
!     ozone:
!
      if (iump_o3 >  0) then
!        the input field of ozone is supplied on nozone levels.
!        these values apply to the upper layers used by the full um.
!        if nozone is smaller than nlevs, the mixing ratio on the
!        bottom level supplied is copied to lower levels. if an
!        extra top level is used its mixing ratio is set by copying
!        the value for the top non-radiative level.
         if (l_extra_top) then
           do l=1, n_profile
             lg=i_gather(l)
             gas_mix_ratio(l, 1, iump_o3)=o3(lg, nozone)
           enddo
         endif
         do i=i_top_copy, nozone+i_top_copy-1
            do l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_o3)=o3(lg, nozone+i_top_copy-i)
            enddo
         enddo
         do i=nozone+i_top_copy, n_layer
            do l=1, n_profile
               lg=i_gather(l)
               gas_mix_ratio(l, i, iump_o3)=o3(lg, 1)
            enddo
         enddo
!      else
!         write(iu_err, '(/a)')
!     &      '*** error: ozone is not in the spectral file.'
!         ierr=i_err_fatal
!         return
      endif
!
!
!
!     other trace gases:
!
!     these gases are not always included in the calculation.
!     testing is therefore more intricate.
!
      if (iump_n2o >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_n2o) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_n2o)=n2o_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_n2o)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_n2o) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: nitrous oxide is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_ch4 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_ch4) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_ch4)=ch4_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_ch4)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_ch4) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: methane is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_cfc11 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_cfc11) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc11)=c11_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc11)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_cfc11) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: cfc11 is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_cfc12 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_cfc12) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc12)=c12_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc12)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_cfc12) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: cfc12 is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_o2 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_o2) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_o2)=o2_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_o2)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_o2) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: o2 is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_cfc113 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_cfc113) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc113)=c113_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_cfc113)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_cfc113) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: cfc113 is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_hcfc22 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_hcfc22) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_hcfc22)=hcfc22_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_hcfc22)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_hcfc22) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: hcfc22 is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_hfc125 >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_hfc125) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc125)=hfc125_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc125)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_hfc125) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: hfc125 is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
      if (iump_hfc134a >  0) then
!        the gas is in the spectral file. if it has been selected
!        from the ui its mixing ratio must be set. if it is in the
!        file but not selected the mixing ratio is set to 0.
         if (l_hfc134a) then
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc134a)=hfc134a_mix_ratio
               enddo
            enddo
         else
            do i=1, n_layer
               do l=1, n_profile
                  gas_mix_ratio(l, i, iump_hfc134a)=0.0e+00
               enddo
            enddo
         endif
      else if (i_call == 1) then
!        the gas is not in the spectral file. an error results if
!        it was to be included in the calculation.
         if (l_hfc134a) then
            write(iu_err, '(/a)')                                       &
     &         '*** error: hfc134a is not in the spectral file.'
            ierr=i_err_fatal
            return
         endif
      endif
!
!
!
      return
      END SUBROUTINE r2_set_gas_mix_ratio
!+ subroutine to set thermodynamic properties
!
! Purpose:
!   pressures, temperatures at the centres and edges of layers
!   and the masses in layers are set.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to assign properties of clouds.
!
! Purpose:
!   the fractions of different types of clouds and their microphysical
!   preoperties are set.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set the parametrization schemes for clouds.
!
! Purpose:
!   the parametrization schemes for each component within a cloud
!   are set.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set fields of aerosols.
!
! Purpose:
!   the mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set fields of climatological aerosols in hadcm3.
!
! Purpose:
!   this routine sets the mixing ratios of climatological aerosols.
!   a separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   the climatoogy used here is the one devised for hadcm3.
!
! Method:
!   straightforward.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
!
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to calculate the total cloud cover.
!
! Purpose:
!   the total cloud cover at all grid-points is determined.
!
! Method:
!   a separate calculation is made for each different assumption about
!   the overlap.
!
! Current owner of code: J.-C. Thelen
!
! History:
!
!     Version            Date                   Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to implement the mrf umist parametrization.
!
! Purpose:
!   effective radii are calculated in accordance with this
!   parametrization.
!
! Method:
!   the number density of ccn is found from the concentration
!   of aerosols, if available. this yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. effective radii are calculated from the number of
!   droplets and the lwc. limits are applied to these values. in
!   deep convective clouds fixed values are assumed.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
!+ subroutine to set the actual process options for the radiation code.
!
! Purpose:
!   to set a consistent set of process options for the radiation.
!
! Method:
!   the global options for the spectral region are compared with the
!   contents of the spectral file. the global options should be set
!   to reflect the capabilities of the code enabled in the model.
!
! Current owner of code: J.-C. Thelen
!
! History:
!       Version         Date                    Comment
!
!       6.2             13-02-06                Original Code included
!                                               into the UM build
!                                               (J.-C. Thelen)
! description of code:
!   fortran 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
#endif
#endif
