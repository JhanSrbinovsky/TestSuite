! AERCMP3A start
!     ------------------------------------------------------------------
!     MODULE TO SET INDICES OF AEROSOL COMPONENTS.
!   4.5   Aug 1998     Set indices for two soot aerosol species.
!                                                 Luke Robinson
!
!
!   5.2   Dec 2000     Set indices for two sea-salt aerosol modes.
!                                                 Andy Jones
!
!   5.4   May 2002     Correct value of NPD_AEROSOL_COMPONENT
!                      and add warning comment.
!                                                 Andy Jones
!   5.5   Feb 2003     Set indices for two modes of biomass
!                      aerosol.                   P Davison
!   5.5   Feb 2003     Set indices for mineral dust aerosol.
!                                        S Woodward
!
      ! maximum number of aerosol components
      ! N.B: this must be at least as large as the
      !      largest value in the list below
      INTEGER,PARAMETER:: NPD_AEROSOL_COMPONENT=28

      INTEGER,PARAMETER:: IP_WATER_SOLUBLE=1
      INTEGER,PARAMETER:: IP_DUST_LIKE=2
      INTEGER,PARAMETER:: IP_OCEANIC=3
      INTEGER,PARAMETER:: IP_SOOT=4
      INTEGER,PARAMETER:: IP_ASH=5
      INTEGER,PARAMETER:: IP_SULPHURIC=6
      INTEGER,PARAMETER:: IP_ACCUM_SULPHATE=10
      INTEGER,PARAMETER:: IP_AITKEN_SULPHATE=11
      INTEGER,PARAMETER:: IP_FRESH_SOOT=12
      INTEGER,PARAMETER:: IP_AGED_SOOT=13
      INTEGER,PARAMETER:: IP_SEASALT_FILM=15
      INTEGER,PARAMETER:: IP_SEASALT_JET=16
      INTEGER,PARAMETER:: IP_DUST_1=17
      INTEGER,PARAMETER:: IP_DUST_2=18
      INTEGER,PARAMETER:: IP_DUST_3=19
      INTEGER,PARAMETER:: IP_DUST_4=20
      INTEGER,PARAMETER:: IP_DUST_5=21
      INTEGER,PARAMETER:: IP_DUST_6=22
      INTEGER,PARAMETER:: IP_BIOMASS_1=23
      INTEGER,PARAMETER:: IP_BIOMASS_2=24
      INTEGER,PARAMETER:: IP_BIOGENIC=25
      INTEGER,PARAMETER:: IP_OCFF_FRESH=26
      INTEGER,PARAMETER:: IP_OCFF_AGED=27
      INTEGER,PARAMETER:: IP_DELTA=28
! AERCMP3A end
