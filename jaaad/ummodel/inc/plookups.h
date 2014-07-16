!----------------------------------------------------------------------
! comdeck: PLOOKUPS
! Purpose: declares parameters used to declare lookup tables, Int_Head
!          and Real_Head. *CALLed by CLOOKUPS
! History:
! version  date         change
! 4.5      21/09/98     New code
! 6.1      30/07/04     Increase max number of input fluxes. A. Hines.
! Author:  M. J. Bell
!----------------------------------------------------------------------

! Parameters

! lengths of tables
      INTEGER,PARAMETER:: Len_FixHd   = 256 ! length of fixed headers

      ! length of first dimension of lookups
      INTEGER,PARAMETER:: Len1_Lookup = 64

      ! length of integer part of each lookup
      INTEGER,PARAMETER:: Len_IntHd   =  45

      ! length of real part of each lookup
      INTEGER,PARAMETER:: Len_RealHd  = 19

      ! grid types
      INTEGER,PARAMETER:: ITGrid = 0 ! tracer grid
      INTEGER,PARAMETER:: IUGrid = 1 !  velocity grid

! values used to construct the max. number of lookup tables
      integer Max_Num_FC_times   ! max. number of f/c times to store
      integer Max_Num_Clim_times ! max. number of climate times (12)
      integer Max_Num_In_Flux    ! max. number of input fluxes at
                                 ! one validity time

      parameter ( Max_Num_FC_times = 30   ) ! from (T-30 to T-24)
                                            ! to (T+144 to T+150)
      parameter ( Max_Num_Clim_times = 12 ) ! 12 months in year
      parameter ( Max_Num_In_Flux    = 25 ) !

      ! max. number of lookup tables

      ! for preferred NWP file
      INTEGER,PARAMETER:: Len2_LookupPreferred = 4096

      ! for previous NWP file
      INTEGER,PARAMETER:: Len2_LookupPrevious  = 4096

      ! for climate file
      INTEGER,PARAMETER:: Len2_LookupClimate =                          &
     &  Max_Num_Clim_times * Max_Num_In_Flux
! PLOOKUPS end
