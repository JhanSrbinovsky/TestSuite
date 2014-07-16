! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine for setting internal dimensions and arrays associated 
!  with the aerosol climatology for NWP.
!
! Purpose:
!   Translate model switches from CNTLATM to an internal array of
!   switches. Compute the number of aerosol species and components
!   that are requested.
!
! Method:
!   Straightforward.
!
! Current owner of code: N. Bellouin
!
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
subroutine set_arcl_dimensions(                                         &
                               ! CNTLATM model switches
                               L_USE_ARCLBIOM,                          &
                               L_USE_ARCLBLCK,                          &
                               L_USE_ARCLSSLT,                          &
                               L_USE_ARCLSULP,                          &
                               L_USE_ARCLDUST,                          &
                               L_USE_ARCLOCFF,                          &
                               L_USE_ARCLDLTA,                          &
                               ! Number of aerosol species used
                               n_arcl_species,                          &
                               ! Corresponding number of components
                               n_arcl_compnts,                          &
                               ! Internal model switches
                               L_USE_ARCL                               &
                              )

  implicit none

#include "arcl_dim.h"
  
  !
  ! Arguments with intent(in)
  !
  
  !
  ! CNTLATM model switches
  !
  logical, intent(in) ::                                                &
                            ! Biomass-burning aerosol climatology
                         L_USE_ARCLBIOM,                                &
                            ! Black-carbon aerosol climatology
                         L_USE_ARCLBLCK,                                &
                            ! Sea-salt aerosol climatology
                         L_USE_ARCLSSLT,                                &
                            ! Sulphate aerosol climatology
                         L_USE_ARCLSULP,                                &
                            ! Mineral dust aerosol climatology
                         L_USE_ARCLDUST,                                &
                            ! Fossil-fuel organic carbon aerosol clim.
                         L_USE_ARCLOCFF,                                &
                            ! Delta aerosol climatology
                         L_USE_ARCLDLTA

  !
  ! Arguments with intent(out)
  !
  
  !
  ! Number of species used
  !
  integer, intent(out) :: n_arcl_species
  
  !
  ! Corresponding number of components
  !
  integer, intent(out) :: n_arcl_compnts
  
  !
  ! Internal model switches
  !
  logical, intent(out), dimension(NPD_ARCL_SPECIES) :: L_USE_ARCL

  !
  ! Local variables
  !
  
  !
  ! IDs of aerosol species and components in the climatology for NWP
  !
#include "arcl_ids.h"

  !
  ! Number of components for each aerosol species
  !
  integer, dimension(NPD_ARCL_SPECIES), parameter ::                    &
        i_arcl_compnts_per_species = (/                                 &
            ! IP_ARCL_SULP: Sulphate (accumulation, Aitken, dissolved)
            3,                                                          &
            ! IP_ARCL_DUST: Six size bins                          
            6,                                                          &
            ! IP_ARCL_SSLT: Sea-salt (film and jet)
            2,                                                          &
            ! IP_ARCL_BLCK: Black-carbon (fresh and aged)
            2,                                                          &
            ! IP_ARCL_BIOM: Biomass (fresh, aged, in-cloud)
            3,                                                          &
            ! IP_ARCL_OCFF: Fossil-fuel org. carb. (fresh, aged, in-cld)
            3,                                                          &
            ! IP_ARCL_DLTA: Delta aerosol
            1                                                           &
        /)

  
  n_arcl_species = 0
  n_arcl_compnts = 0
  L_USE_ARCL = .false.
  
  !
  ! For each CNTLATM model switch that is set, set the corresponding
  ! internal model switch, increase the number of requested species
  ! by 1, and increase the number of components by the corresponding
  ! amount.
  !
  if (L_USE_ARCLBIOM) then
    
    L_USE_ARCL(IP_ARCL_BIOM) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_BIOM)
    
  end if

  if (L_USE_ARCLBLCK) then
    
    L_USE_ARCL(IP_ARCL_BLCK) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_BLCK)
    
  end if

  if (L_USE_ARCLSSLT) then
    
    L_USE_ARCL(IP_ARCL_SSLT) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_SSLT)
    
  end if

  if (L_USE_ARCLSULP) then
    
    L_USE_ARCL(IP_ARCL_SULP) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_SULP)
    
  end if
  
  if (L_USE_ARCLDUST) then
    
    L_USE_ARCL(IP_ARCL_DUST) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_DUST)
    
  end if

  if (L_USE_ARCLOCFF) then
    
    L_USE_ARCL(IP_ARCL_OCFF) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_OCFF)
    
  end if
  
  if (L_USE_ARCLDLTA) then
  
    L_USE_ARCL(IP_ARCL_DLTA) = .true.
    n_arcl_species = n_arcl_species + 1
    n_arcl_compnts = n_arcl_compnts +                                   &
                     i_arcl_compnts_per_species(IP_ARCL_DLTA)
  
  end if

end subroutine set_arcl_dimensions
