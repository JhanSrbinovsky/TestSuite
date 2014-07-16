

!+ ---------------------------------------------------------------------
!  Sructure containing new boundary layer diagnostics.
!  This permits easier addition of new boundary layer
!  diagnostics without additional passing of arguments
!  though the boundary layer tree.
!  It also does not require the addition
!  of extra subroutine arguments when adding a new diagnostic.
!- ----------------------------------------------------------------------

      MODULE bl_diags_mod
      TYPE Strnewbldiag

! Need to create a flag and a pointer

        Logical ::  L_dscbase
        Logical ::  L_cldbase
        Logical ::  L_weparm
        Logical ::  L_weparm_dsc
        Logical ::  L_oblen
        Logical ::  L_ustar
        Logical ::  L_wbsurf
        Logical ::  L_gradrich
        Logical ::  L_wstar
        Logical ::  L_dbdz
        Logical ::  L_dvdzm
        Logical ::  L_rhokm
        Logical ::  L_rhokh
        Logical ::  L_tke
        Logical ::  L_ostressx
        Logical ::  L_ostressy
!
        Real, Pointer :: dscbase(:, :)
!                    360      Height of decoupled layer base
        Real, Pointer :: cldbase(:, :)
!                    361      Height of stratocumulus cloud base
        Real, Pointer :: weparm(:, :)
!                    362      Entrainment rate for SML
        Real, Pointer :: weparm_dsc(:, :)
!                    363      Entrainment rate for DSC
        Real, Pointer :: oblen(:, :)
!                    464      Surface Obukhov length
        Real, Pointer :: ustar(:, :)
!                    465      Friction velocity
        Real, Pointer :: wbsurf(:, :)
!                    467      Surface buoyancy flux
        Real, Pointer :: gradrich(:, :, :)
!                    468      Gradient Richardson number
        Real, Pointer :: wstar(:, :)
!                    466      Convective velocity scale
        Real, Pointer :: dbdz(:, :, :)
!                    469      Vertical buoyancy gradient
        Real, Pointer :: dvdzm(:, :, :)
!                    470      Modulus of wind shear
        Real, Pointer :: rhokm(:, :, :)
!                    471      BL Momentum diffusivity
        Real, Pointer :: rhokh(:, :, :)
!                    472      BL Thermal diffusivity
        Real, Pointer :: tke(:, :, :)
!                    473      Turbulent kinetic energy
        Real, Pointer :: ostressx(:, :, :)
!                    474      Orographic stress (x-component)
        Real, Pointer :: ostressy(:, :, :)
!                    475      Orographic stress (y-component)

      END TYPE Strnewbldiag
! ----------------------------------------------------------------------
      END MODULE bl_diags_mod
