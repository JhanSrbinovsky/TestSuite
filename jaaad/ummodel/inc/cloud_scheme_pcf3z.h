!     ------------------------------------------------------------------
!     Module to define reference numbers for cloud schemes.
!
      INTEGER, Parameter :: IP_cloud_mix_max         = 2
!           Maximum/random overlap in a mixed column
      INTEGER, Parameter :: IP_cloud_mix_random      = 4
!           Random overlap in a mixed column
      INTEGER, Parameter :: IP_cloud_column_max      = 3
!           Maximum overlap in a column model
      INTEGER, Parameter :: IP_cloud_clear           = 5
!           Clear column
      INTEGER, Parameter :: IP_cloud_triple          = 6
!           Mixed column with split between convective and layer cloud
      INTEGER, Parameter :: IP_cloud_part_corr       = 7
!           Coupled overlap with partial correlation of cloud
      INTEGER, Parameter :: IP_cloud_part_corr_cnv   = 8
!           Coupled overlap with partial correlation of cloud
!           with a separate treatment of convective cloud
#if defined(EXPCLOP)
      INTEGER, Parameter :: IP_cloud_ovlptest        = 9
!           Version of cloud scheme for testing
#endif
!
!     ------------------------------------------------------------------
