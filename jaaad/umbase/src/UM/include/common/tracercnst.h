! ----------------------- Comdeck: TRACERCNST  -----------------------
! Description: COMDECK containing a common block of constants used in
!              the routines to model the air-sea flux of chemical
!              tracers
!
!
! Author : S.Liddicoat
!
! History:
! Version  Date      Comment.
!  5.5  26/02/03  New code. S.Liddicoat
!
! Constants used in the calculations of air-sea flux of chemical
! tracers
!
      INTEGER, PARAMETER :: MASTF  = 10  ! Nmbr of air-sea fluxes

! Pointers to air-sea (and virtual) fluxes.
      INTEGER :: ASTF_ALK
      INTEGER :: ASTF_TCO2
      INTEGER :: ASTF_CEX
      INTEGER :: ASTF_C14
      INTEGER :: ASTF_BOMC14
      INTEGER :: ASTF_CFC11
      INTEGER :: ASTF_CFC12
      INTEGER :: ASTF_HE3
      INTEGER :: ASTF_HE4
      INTEGER :: ASTF_O2
      INTEGER :: ASTF_DMS

!  Reference values for the virtual fluxes of each tracer
      REAL    ::                                                        &
     &   he3_ref                                                        &
     & , he4_ref                                                        &
     & , c14_ref                                                        &
     & , co2_ref                                                        &
     & , alk_ref                                                        &
     & , cfc11_ref                                                      &
     & , cfc12_ref                                                      &
     & , yr_ofst_co2                                                    &
     & , yr_ofst_c14                                                    &
     & , yr_ofst_cfc11                                                  &
     & , yr_ofst_cfc12

      COMMON /TRACERCNST/                                               &
     & ASTF_ALK, ASTF_TCO2, ASTF_CEX, ASTF_C14, ASTF_BOMC14, ASTF_O2,   &
     & ASTF_CFC11, ASTF_CFC12, ASTF_HE3, ASTF_HE4, ASTF_DMS,            &
     & he3_ref, he4_ref, c14_ref, co2_ref, alk_ref, cfc11_ref,          &
     & cfc12_ref,                                                       &
     & yr_ofst_co2, yr_ofst_c14, yr_ofst_cfc11, yr_ofst_cfc12
