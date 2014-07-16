!     ------------------------------------------------------------------
!     Module to set representations of clouds.
!
!     Representations:
      INTEGER, PARAMETER :: IP_CLOUD_HOMOGEN    = 1
!           All components are mixed homogeneously
      INTEGER, PARAMETER :: IP_CLOUD_ICE_WATER  = 2
!           Ice and water clouds are treated separately
      INTEGER, PARAMETER :: IP_CLOUD_CONV_STRAT = 3
!           Clouds are divided into homogeneously mixed
!           stratiform and convective parts
      INTEGER, PARAMETER :: IP_CLOUD_CSIW       = 4
!           Clouds divided into ice and water phases and
!           into stratiform and convective components.
#if defined(XICE)
      INTEGER, PARAMETER :: IP_CLOUD_CSIW_CRYS  = 5
!           Clouds divided into ice and water phases and
!           into stratiform and convective components, with
!           the ice clouds containing a mixture of crystal
!           shapes
#endif
!     ------------------------------------------------------------------
