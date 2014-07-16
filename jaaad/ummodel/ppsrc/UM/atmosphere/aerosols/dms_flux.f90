
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      SUBROUTINE DMS_FLUX(                                              &
!
! Arguments IN
     &             row_length                                           &
     &,            rows                                                 &
     &,            u_1                                                  &
     &,            v_1                                                  &
     &,            height                                               &
     &,            Tstar                                                &
     &,            land_fract                                           &
     &,            DMS_conc                                             &
     &,            L_Liss_Merlivat                                      &
     &,            L_Wanninkhof                                         &
     &,            L_Nightingale                                        &
! Arguments OUT
     &,            f_DMS                                                &
     &            )
!---------------------------------------------------------------------
! Purpose: To calculate the flux of DMS (as kg m-2 s-1 of sulphur)
!          from the ocean surface as a function of its concentration
!          in seawater and of windspeed. The sea-air exchange can
!          be determined according to one of three commonly-used
!          parametrization schemes, those of Liss & Merlivat (1986),
!          Wanninkhof (1992) or Nightingale et al. (2000). The routine
!          is called by Aero_Ctl.
!
! Method:  10m windspeed is calculated assuming a neutral profile and a
!          roughness length appropriate for water. The Schmidt number
!          for DMS is then calculated as in Saltzman et al. (1993), and
!          used with the windspeed to determine the mass transfer (or
!          "piston") velocity according to the desired parametrization.
!          This is then used to determine the sea-air mass flux of DMS
!          as a function of sea-water DMS concentration. High surface
!          temperatures (caused by the land portion of a gridbox when
!          coastal tiling is not active) cause negative Sc values which
!          would give a floating-point error in the k_DMS calculation,
!          so the Tstar values are capped. This shouldn't be a problem
!          when coastal tiling is on as then the Tstar values passed in
!          are those for sea only.
!
! Current owner of code: A. Jones
!
! History:
!
! Version     Date     Comment
! -------     ----     -------
!   6.1     07/04/04   Original code.
!                      (A. Jones)
!   6.2     20/10/05   DEF for a17_2b added.
!                      (A. Jones)
!
! Code Description:
!   FORTRAN77 + common extensions
!
!---------------------------------------------------------------------
!
!


      USE auscom_cpl_data_mod,                                          &
      &   Only : auscom_salinity, access_tfs, ocn_sss


      IMPLICIT NONE
!
!
!*L------------------COMDECK C_O_DG_C-----------------------------------
! ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
! TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
! TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS

      Real, Parameter :: ZeroDegC = 273.15
      Real, Parameter :: TFS      = 271.35
      Real, Parameter :: TM       = 273.15

!*----------------------------------------------------------------------
!
!
! Arguments with intent IN:
!
      INTEGER                                                           &
     &     row_length                                                   &
     &,    rows

      REAL                                                              &
     &     u_1(row_length, rows)                                        &
                                         ! Level 1 u-wind (ms-1)
     &,    v_1(row_length, rows)                                        &
                                         ! Level 1 v-wind (ms-1)
     &,    height(row_length, rows)                                     &
                                         ! Level 1 centre height (m)
     &,    Tstar(row_length, rows)                                      &
                                         ! Surface temperature (K)
     &,    land_fract(row_length, rows)                                 &
                                         ! Fraction of land in gridbox
     &,    DMS_conc(row_length, rows)    ! Concentration of DMS in
                                         !         seawater (nmol l-1)
      LOGICAL                                                           &
     &     L_Liss_Merlivat                                              &
                                         ! Switches to determine which
     &,    L_Wanninkhof                                                 &
                                         !  scheme to use to calculate
     &,    L_Nightingale                 !  mass transfer velocity
!
! Arguments with intent OUT:
!
      REAL                                                              &
     &     f_DMS(row_length, rows)       ! Sea-air flux of DMS
!                                        !          (kg[S] m-2 s-1)
!
! Local variables:
!
!
      INTEGER                                                           &
     &     i                                                            &
                                         ! Loop counters
     &,    j

      REAL                                                              &
     &     Sc(row_length, rows)                                         &
                                         ! Schmidt number
     &,    wind_10m(row_length, rows)                                   &
                                         ! 10m windspeed (ms-1)
     &,    k_DMS(row_length, rows)       ! Piston velocity of DMS
                                         !                  (cm h-1)
      REAL                                                              &
     &     T_C                                                          &
                          ! Surface temperature in degrees Celsius
     &,    T_max                                                        &
                          ! Max T to avoid breaking the Sc fit (C)
     &,    k_600                                                        &
                          ! Piston velocities for gases with Schmidt
     &,    k_660                                                        &
                          !      numbers of 600 & 660 resp. (cm h-1)
     &,    n                                                            &
                          ! Schmidt number exponent
     &,    z0_sea         ! Roughness length over sea (m)

      PARAMETER(                                                        &
     &     z0_sea=2.5e-04                                               &
     &,    T_max=47.0                                                   &
     &   )

      REAL ltfs

!
! Calculate 10m windspeed and the Schmidt number (Sc):
!
      ltfs = access_tfs

      Do j = 1, rows
        Do i = 1, row_length
          wind_10m(i, j) = (sqrt((u_1(i, j)**2)+(v_1(i, j)**2)))        &
     &                *(alog(10.0/z0_sea))/(alog(height(i, j)/z0_sea))
          if (ocn_sss) then
              ltfs = ZeroDegC - 0.054 * auscom_salinity(I,J)
          end if
          T_C = min((max(Tstar(i, j), lTFS) - ZERODEGC), T_max)
          Sc(i, j) = 2674.0 - (147.12*T_C) + (3.726*T_C**2)             &
     &                                                - (0.038*T_C**3)
        End Do
      End Do
!
! Determine the mass transfer (or "piston") velocity (k_DMS) over sea
! according to the specified parametrization scheme:
!
      If (L_Liss_Merlivat) Then
        Do j = 1, rows
          Do i = 1, row_length
            If (wind_10m(i, j)  <=  3.6) Then
              k_600 = 0.17 * wind_10m(i, j)
              n = -2.0/3.0
            End If
            If (wind_10m(i, j)  >   3.6 .AND.                           &
     &          wind_10m(i, j)  <=  13.0) Then
              k_600 = (2.85*wind_10m(i, j)) - 9.65
              n = -0.5
            End If
            If (wind_10m(i, j)  >   13.0) Then
              k_600 = (5.9*wind_10m(i, j)) - 49.3
              n = -0.5
            End If
            If (land_fract(i, j)  <   1.0) Then
              k_DMS(i, j) = k_600 * (Sc(i, j)/600.0)**n
            Else
              k_DMS(i, j) = 0.0
            End If
          End Do
        End Do
      End If

      If (L_Wanninkhof) Then
        Do j = 1, rows
          Do i = 1, row_length
            k_660 = 0.31 * wind_10m(i, j)**2
            n = -0.5
            If (land_fract(i, j)  <   1.0) Then
              k_DMS(i, j) = k_660 * (Sc(i, j)/660.0)**n
            Else
              k_DMS(i, j) = 0.0
            End If
          End Do
        End Do
      End If

      If (L_Nightingale) Then
        Do j = 1, rows
          Do i = 1, row_length
            k_600 = (0.222*wind_10m(i, j)**2) + (0.333*wind_10m(i, j))
            n = -0.5
            If (land_fract(i, j)  <   1.0) Then
              k_DMS(i, j) = k_600 * (Sc(i, j)/600.0)**n
            Else
              k_DMS(i, j) = 0.0
            End If
          End Do
        End Do
      End If
!
! Finally, calculate the sea-air flux of DMS as a function of k_DMS
! and dissolved DMS concentration. The former requires a conversion
! from cm hour-1 to ms-1, and the latter from nanomoles per litre to
! kg[S] m-3, to return the flux in kg[S] m-2 sec-1.
!
      Do j = 1, rows
        Do i = 1, row_length
          f_DMS(i, j) = (k_DMS(i, j) / 3.6e5)                           &
     &                * (DMS_conc(i, j) * 32.0e-9)
        End Do
      End Do
!
!-----------------------------------------------------------------------
!
      Return
      END SUBROUTINE DMS_FLUX
!
