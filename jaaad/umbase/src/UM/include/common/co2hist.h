! ----------------------- Comdeck: CO2HIST  ----------------------------
! Description: COMDECK containing a common block for the time evolution
!              of atmospheric CO2.
!
! Author : S.Spall based on code of C.D.Jones
!
! History:
! Version  Date      Comment.
!  5.3  13/09/01  New code. S.Spall
!  5.5  28/02/03  Allow similar time evolution of atmospheric
!                 carbon 14.                      S.Liddicoat

      INTEGER, PARAMETER :: MSCNRC = 2501
      INTEGER, PARAMETER :: MSCNLT = 19

      INTEGER :: nrecs_co2
      INTEGER :: nlats_co2

      INTEGER :: nrecs_c14
      INTEGER :: nlats_c14

      INTEGER :: nrecs_cfc11
      INTEGER :: nlats_cfc11

      INTEGER :: nrecs_cfc12
      INTEGER :: nlats_cfc12

      REAL :: atmscn_n_co2(MSCNRC,0:MSCNLT)
      REAL :: atmscn_s_co2(MSCNRC,0:MSCNLT)
      REAL :: atmscn_t_co2(MSCNRC)

      REAL :: atmscn_n_c14(MSCNRC,0:MSCNLT)
      REAL :: atmscn_s_c14(MSCNRC,0:MSCNLT)
      REAL :: atmscn_t_c14(MSCNRC)

      REAL :: atmscn_n_cfc11(MSCNRC,0:MSCNLT)
      REAL :: atmscn_s_cfc11(MSCNRC,0:MSCNLT)
      REAL :: atmscn_t_cfc11(MSCNRC)

      REAL :: atmscn_n_cfc12(MSCNRC,0:MSCNLT)
      REAL :: atmscn_s_cfc12(MSCNRC,0:MSCNLT)
      REAL :: atmscn_t_cfc12(MSCNRC)

      CHARACTER*(120) file_cscen
      CHARACTER*(120) file_c14scen
      CHARACTER*(120) file_cfc11scen
      CHARACTER*(120) file_cfc12scen

      COMMON /CO2HIST/                                                  &
     &  nrecs_co2, nlats_co2                                            &
     &, atmscn_n_co2, atmscn_s_co2, atmscn_t_co2                        &
     &, nrecs_c14, nlats_c14                                            &
     &, atmscn_n_c14, atmscn_s_c14, atmscn_t_c14                        &
     &, nrecs_cfc11, nlats_cfc11                                        &
     &, atmscn_n_cfc11, atmscn_s_cfc11, atmscn_t_cfc11                  &
     &, nrecs_cfc12, nlats_cfc12                                        &
     &, atmscn_n_cfc12, atmscn_s_cfc12, atmscn_t_cfc12

      COMMON /FILECSN/ file_cscen, file_c14scen                         &
     &, file_cfc11scen, file_cfc12scen
