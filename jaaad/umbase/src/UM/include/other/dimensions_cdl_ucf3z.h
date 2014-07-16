!     ------------------------------------------------------------------
!     Module setting maximum dimensions for CDL-files.
!
      INTEGER                                                           &
     &     NPD_CDL_DIMEN                                                &
!             Maximum number of CDL dimensions
     &   , NPD_CDL_DIMEN_SIZE                                           &
!             Maximum size of a dimension
     &   , NPD_CDL_DATA                                                 &
!             Maximum size of CDL data
     &   , NPD_CDL_VAR
!             Maximum number of CDL variables in a file
!
      PARAMETER(                                                        &
     &     NPD_CDL_DIMEN=5                                              &
     &   , NPD_CDL_DIMEN_SIZE=502                                       &
     &   , NPD_CDL_DATA=10000                                           &
     &   , NPD_CDL_VAR=4                                                &
     &   )
!
!     ------------------------------------------------------------------
