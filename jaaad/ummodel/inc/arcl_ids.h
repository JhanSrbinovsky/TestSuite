! arcl_ids.h
!
! Constants used by the aerosol climatology for NWP
!
      !
      ! Available species:
      !
      !   1. Sulphate
      !   2. Mineral dust
      !   3. Sea-salt
      !   4. Black-carbon (also named soot)
      !   5. Biomass-burning
      !   6. Fossil-fuel organic carbon
      !   7. Delta aerosol
      !
      ! Note: When adding a species, increase NPD_ARCL_SPECIES in
      !       arcl_dim.h. NPD_ARCL_SPECIES must be equal to the
      !       largest IP_ARCL_????

      integer, parameter :: IP_ARCL_SULP = 1
      integer, parameter :: IP_ARCL_DUST = 2
      integer, parameter :: IP_ARCL_SSLT = 3
      integer, parameter :: IP_ARCL_BLCK = 4
      integer, parameter :: IP_ARCL_BIOM = 5
      integer, parameter :: IP_ARCL_OCFF = 6
      integer, parameter :: IP_ARCL_DLTA = 7
      
      !
      ! List of components.
      !
      !   The number of components depends on the species:
      !
      !   1. Sulphate: 3 (accumulation, Aitken, dissolved)
      !   2. Mineral dust: 6 size bins
      !   3. Sea-salt: 2 (film and jet)
      !   4. Black-carbon: 2 (fresh and aged)
      !   5. Biomass-burning: 3 (fresh, aged, in-cloud)
      !   6. Fossil-fuel organic carbon: 3 (fresh, aged, in-cloud)
      !   7. Delta aerosol: 1
      !
      ! Note: when adding a component, increase NPD_ARCL_COMPNTS
      !       in arcl_dim.h. NPD_ARCL_COMPNTS must be equal to the
      !       largest IP_ARCL_????_??.
      !       Array i_arcl_compnts_per_species in set_arcl_dimensions()
      !       will also have to be updated.
      !
      
      integer, parameter :: IP_ARCL_SULP_AC = 1
      integer, parameter :: IP_ARCL_SULP_AK = 2
      integer, parameter :: IP_ARCL_SULP_DI = 3
      integer, parameter :: IP_ARCL_DUST_B1 = 4
      integer, parameter :: IP_ARCL_DUST_B2 = 5
      integer, parameter :: IP_ARCL_DUST_B3 = 6
      integer, parameter :: IP_ARCL_DUST_B4 = 7
      integer, parameter :: IP_ARCL_DUST_B5 = 8
      integer, parameter :: IP_ARCL_DUST_B6 = 9
      integer, parameter :: IP_ARCL_SSLT_FI = 10
      integer, parameter :: IP_ARCL_SSLT_JT = 11
      integer, parameter :: IP_ARCL_BLCK_FR = 12
      integer, parameter :: IP_ARCL_BLCK_AG = 13
      integer, parameter :: IP_ARCL_BIOM_FR = 14
      integer, parameter :: IP_ARCL_BIOM_AG = 15
      integer, parameter :: IP_ARCL_BIOM_IC = 16
      integer, parameter :: IP_ARCL_OCFF_FR = 17
      integer, parameter :: IP_ARCL_OCFF_AG = 18
      integer, parameter :: IP_ARCL_OCFF_IC = 19
      integer, parameter :: IP_ARCL_DLTA_DL = 20

! end of arcl_ids.h
